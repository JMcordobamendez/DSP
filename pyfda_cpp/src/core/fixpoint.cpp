#include "fixpoint.hpp"

#include <algorithm>
#include <cmath>

namespace pyfda {

const char *quant_key(Quant q) {
    switch (q) {
    case Quant::Floor: return "floor";
    case Quant::Round: return "round";
    case Quant::Fix: return "fix";
    case Quant::Ceil: return "ceil";
    case Quant::None: return "none";
    }
    return "";
}

const char *ovfl_key(Ovfl o) {
    switch (o) {
    case Ovfl::Wrap: return "wrap";
    case Ovfl::Sat: return "sat";
    case Ovfl::None: return "none";
    }
    return "";
}

Quant quant_from_key(const std::string &k) {
    for (Quant q : {Quant::Floor, Quant::Round, Quant::Fix, Quant::Ceil, Quant::None})
        if (k == quant_key(q)) return q;
    throw DesignError("Unknown quantization '" + k + "'");
}

Ovfl ovfl_from_key(const std::string &k) {
    for (Ovfl o : {Ovfl::Wrap, Ovfl::Sat, Ovfl::None})
        if (k == ovfl_key(o)) return o;
    throw DesignError("Unknown overflow behaviour '" + k + "'");
}

double QFormat::lsb() const { return std::ldexp(1.0, -WF); }
double QFormat::min() const { return -std::ldexp(1.0, WI); }
double QFormat::max() const { return std::ldexp(1.0, WI) - lsb(); }

double Quantizer::fixp(double y) {
    // (2) quantization in integer scale
    y = y * std::ldexp(1.0, m_q.WF);
    double yq = y;
    switch (m_q.quant) {
    case Quant::Floor: yq = std::floor(y); break;
    case Quant::Round: yq = std::nearbyint(y); break;  // numpy.round: half to even
    case Quant::Fix: yq = std::trunc(y); break;
    case Quant::Ceil: yq = std::ceil(y); break;
    case Quant::None: break;
    }
    // (3) overflow w.r.t. MSB = 2^(W - 2)
    if (m_q.ovfl != Ovfl::None) {
        const double MSB = std::ldexp(1.0, m_q.WI + m_q.WF - 1);
        const double MAX = 2 * MSB - 1, MIN = -2 * MSB;
        const bool neg = yq < MIN, pos = yq > MAX;
        if (neg || pos) {
            ++m_n_over;
            if (m_q.ovfl == Ovfl::Sat) yq = pos ? MAX : MIN;
            else {
                const double sgn = yq > 0 ? 1.0 : (yq < 0 ? -1.0 : 0.0);
                yq = yq - 4. * MSB * std::trunc((sgn * 2 * MSB + yq) / (4 * MSB));
            }
        }
    }
    // (4) back to fractional scale
    return yq / std::ldexp(1.0, m_q.WF);
}

Vec Quantizer::fixp(const Vec &y) {
    Vec o(y.size());
    for (size_t i = 0; i < y.size(); ++i) o[i] = fixp(y[i]);
    return o;
}

long long to_int(double v_q, const QFormat &q) { return std::llround(v_q * std::ldexp(1.0, q.WF)); }

std::string to_base(long long v, int W, int base) {
    if (W < 1 || W > 64) throw DesignError("Word length must be 1 ... 64 bits.");
    const unsigned long long mask = W == 64 ? ~0ULL : ((1ULL << W) - 1);
    unsigned long long u = static_cast<unsigned long long>(v) & mask;
    const int bits = base == 2 ? 1 : base == 8 ? 3 : 4;
    const int digits = (W + bits - 1) / bits;
    static const char *hex = "0123456789ABCDEF";
    std::string s(size_t(digits), '0');
    for (int i = digits - 1; i >= 0; --i) {
        s[size_t(i)] = hex[u & ((1ULL << bits) - 1)];
        u >>= bits;
    }
    return s;
}

std::string to_csd(long long v, int W) {
    // non-adjacent form, most significant digit first, padded to W digits
    std::string s;
    long long n = v;
    while (n != 0) {
        if (n & 1) {
            const long long d = 2 - (((n % 4) + 4) % 4);  // +1 or -1
            s += d > 0 ? '+' : '-';
            n -= d;
        } else {
            s += '0';
        }
        n /= 2;
    }
    while (int(s.size()) < W) s += '0';
    std::reverse(s.begin(), s.end());
    return s;
}

int coeff_wi(const Vec &c) {
    double m = 0;
    for (double v : c) m = std::max(m, std::fabs(v));
    if (!(m > 0) || !std::isfinite(m)) return 0;
    // smallest WI with max |c| < 2^WI
    return std::max(0, int(std::floor(std::log2(m))) + 1);
}

void update_auto_formats(FxSpec &s, const Ba &ba, const Sos &sos, bool fir) {
    if (fir) {
        if (s.coeff_auto) s.qcb.WI = coeff_wi(ba.b);
        if (s.acc_auto) {
            double sum = 0;
            for (double v : ba.b) sum += std::fabs(v);
            const int guard = sum > 1 ? int(std::ceil(std::log2(sum))) : 0;
            s.qacc.WF = s.qi.WF + s.qcb.WF;
            s.qacc.WI = s.qi.WI + s.qcb.WI + guard;
        }
        return;
    }
    if (s.coeff_auto) {
        Vec b, a;
        for (const auto &sec : sos) {
            b.insert(b.end(), {sec[0], sec[1], sec[2]});
            a.insert(a.end(), {sec[4], sec[5]});
        }
        s.qcb.WI = coeff_wi(b);
        s.qca.WI = coeff_wi(a);
    }
    if (s.acc_auto) {
        // section inputs are the filter input (first section) or the output format
        const int wf_in = std::max(s.qi.WF, s.qo.WF), wi_in = std::max(s.qi.WI, s.qo.WI);
        s.qacc.WF = std::max(wf_in + s.qcb.WF, s.qo.WF + s.qca.WF);
        s.qacc.WI = std::max(wi_in + s.qcb.WI, s.qo.WI + s.qca.WI) + 3;  // 5 products + sign growth
    }
}

Vec quant_coeffs(const Vec &c, const QFormat &q, bool recursive, long long *n_over) {
    Quantizer Q(q);
    Vec o(c.size());
    for (size_t i = 0; i < c.size(); ++i) o[i] = (recursive && i == 0) ? 1.0 : Q.fixp(c[i]);
    if (n_over) *n_over = Q.overflows();
    return o;
}

FxResult fx_filter_fir(const Vec &b, const FxSpec &spec, const Vec &x) {
    FxResult r;
    Quantizer Qi(spec.qi), Qacc(spec.qacc), Qo(spec.qo);
    r.b_q = quant_coeffs(b, spec.qcb, false, &r.n_over_coeff);
    r.x_q = Qi.fixp(x);
    r.y.resize(x.size());
    const size_t L = r.b_q.size();
    for (size_t k = 0; k < x.size(); ++k) {
        double acc = 0;
        for (size_t i = 0; i < L && i <= k; ++i) acc += Qacc.fixp(r.x_q[k - i] * r.b_q[i]);
        r.y[k] = Qo.fixp(Qacc.fixp(acc));
    }
    r.n_over_i = Qi.overflows();
    r.n_over_acc = Qacc.overflows();
    r.n_over_o = Qo.overflows();
    return r;
}

FxResult fx_filter_sos(const Sos &sos, const FxSpec &spec, const Vec &x) {
    FxResult r;
    Quantizer Qi(spec.qi), Qacc(spec.qacc), Qo(spec.qo), Qb(spec.qcb), Qa(spec.qca);
    r.x_q = Qi.fixp(x);
    Vec s = r.x_q;
    for (const auto &sec : sos) {
        if (sec[3] != 1.0) throw DesignError("Second-order sections must be normalized to a0 = 1.");
        const double b0 = Qb.fixp(sec[0]), b1 = Qb.fixp(sec[1]), b2 = Qb.fixp(sec[2]);
        const double a1 = Qa.fixp(sec[4]), a2 = Qa.fixp(sec[5]);
        r.b_q.insert(r.b_q.end(), {b0, b1, b2});
        r.a_q.insert(r.a_q.end(), {1.0, a1, a2});
        double x1 = 0, x2 = 0, y1 = 0, y2 = 0;
        for (double &v : s) {
            const double acc_b = Qacc.fixp(b0 * v) + Qacc.fixp(b1 * x1) + Qacc.fixp(b2 * x2);
            const double acc_a = Qacc.fixp(a1 * y1) + Qacc.fixp(a2 * y2);
            const double y = Qo.fixp(Qacc.fixp(acc_b - acc_a));
            x2 = x1;
            x1 = v;
            y2 = y1;
            y1 = y;
            v = y;
        }
    }
    r.y = s;
    r.n_over_coeff = Qb.overflows() + Qa.overflows();
    r.n_over_i = Qi.overflows();
    r.n_over_acc = Qacc.overflows();
    r.n_over_o = Qo.overflows();
    return r;
}

}  // namespace pyfda
