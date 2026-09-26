#include "filter_info.hpp"

#include <limits>

#include "conversions.hpp"

#include <algorithm>
#include <cmath>

namespace pyfda {

namespace {
// |H| in dB on a dense grid between the normalized frequencies F0 and F1 (F = f / f_S)
Vec mag_db(const FilterDesign &d, double F0, double F1, int n = 2000) {
    Vec w(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) w[size_t(i)] = 2 * PI * (F0 + (F1 - F0) * i / (n - 1));
    const CVec H = d.sos.empty() ? freqz(d.ba, w) : freqz(d.sos, w);
    Vec m(H.size());
    for (size_t i = 0; i < H.size(); ++i) m[i] = 20 * std::log10(std::max(std::abs(H[i]), 1e-15));
    return m;
}
}  // namespace

FilterInfo filter_info(const FilterDesign &d) {
    FilterInfo fi;
    const FilterSpec &s = d.spec;
    fi.fir = d.fir;
    fi.n_b = d.ba.b.size();
    fi.n_a = d.ba.a.size();
    fi.n_sos = d.sos.size();
    fi.order = int(std::max(fi.n_b, fi.n_a)) - 1;
    for (const cplx &p : d.zpk.p) fi.max_pole_radius = std::max(fi.max_pole_radius, std::abs(p));
    fi.stable = fi.max_pole_radius < 1.0;
    for (const cplx &z : d.zpk.z)
        if (std::abs(z) > 1.0 + 1e-9) fi.min_phase = false;
    if (d.fir && fi.n_b > 1) {
        const Vec &b = d.ba.b;
        double scale = 0;
        for (double v : b) scale = std::max(scale, std::fabs(v));
        bool sym = true, anti = true;
        for (size_t i = 0; i < b.size(); ++i) {
            sym = sym && std::fabs(b[i] - b[b.size() - 1 - i]) <= 1e-9 * scale;
            anti = anti && std::fabs(b[i] + b[b.size() - 1 - i]) <= 1e-9 * scale;
        }
        fi.linear_phase = sym || anti;
    }
    const CVec H0 = d.sos.empty() ? freqz(d.ba, {0.0, PI}) : freqz(d.sos, {0.0, PI});
    fi.gain_dc = std::abs(H0[0]);
    fi.gain_ny = std::abs(H0[1]);
    const Vec all = mag_db(d, 0, 0.5, 4096);
    fi.h_max_db = *std::max_element(all.begin(), all.end());

    // like pyfda's Info tab, the band edge specs are also checked for manual order designs
    // (there they are only a reference); not for manual filters and MA band pass / stop
    if (is_manual(s.method)) return fi;
    const bool ma = s.method == DesignMethod::MovingAverage;
    if (ma && (s.rt == RespType::BP || s.rt == RespType::BS)) return fi;
    const double fs = s.f_s;
    auto add = [&](const char *name, double f0, double f1, bool pass) {
        if (!(f1 > f0)) return;
        const Vec m = mag_db(d, f0 / fs, f1 / fs);
        const auto mm = std::minmax_element(m.begin(), m.end());
        BandCheck b;
        b.name = name;
        b.f0 = f0;
        b.f1 = f1;
        b.pass = pass;
        b.spec_db = pass ? s.A_PB : s.A_SB;
        b.achieved_db = pass ? *mm.second - *mm.first : -*mm.second;
        // small tolerance: the minimum order formulas are estimates, the grid is finite
        b.ok = pass ? b.achieved_db <= b.spec_db * 1.001 + 1e-9 : b.achieved_db >= b.spec_db * 0.999 - 1e-9;
        fi.bands.push_back(b);
    };
    const double ny = fs / 2;
    if (ma) {  // only the stop band is specified
        if (s.rt == RespType::LP) add("SB", s.f_sb, ny, false);
        else add("SB", 0, s.f_sb, false);
        return fi;
    }
    switch (s.rt) {
    case RespType::LP:
        add("PB", 0, s.f_pb, true);
        add("SB", s.f_sb, ny, false);
        break;
    case RespType::HP:
        add("SB", 0, s.f_sb, false);
        add("PB", s.f_pb, ny, true);
        break;
    case RespType::BP:
        add("SB", 0, s.f_sb, false);
        add("PB", s.f_pb, s.f_pb2, true);
        add("SB2", s.f_sb2, ny, false);
        break;
    case RespType::BS:
        add("PB", 0, s.f_pb, true);
        add("SB", s.f_sb, s.f_sb2, false);
        add("PB2", s.f_pb2, ny, true);
        break;
    }
    return fi;
}

double h_mag_z(const Ba &ba, cplx z) {
    if (z == cplx(0, 0)) {  // limit z -> 0: ratio of the highest order terms
        auto last = [](const Vec &c) {
            size_t n = c.size();
            while (n > 0 && c[n - 1] == 0) --n;
            return n;  // index of the highest nonzero coefficient + 1
        };
        const size_t nb = last(ba.b), na = last(ba.a);
        if (nb == 0) return 0.0;
        if (na == 0 || nb > na) return std::numeric_limits<double>::infinity();
        return nb < na ? 0.0 : std::fabs(ba.b[nb - 1] / ba.a[na - 1]);
    }
    // Horner scheme in w = 1 / z
    const cplx w = cplx(1.0, 0.0) / z;
    auto poly = [&](const Vec &c) {
        cplx acc = 0;
        for (size_t i = c.size(); i-- > 0;) acc = acc * w + c[i];
        return acc;
    };
    const cplx num = poly(ba.b), den = poly(ba.a);
    if (den == cplx(0, 0)) return std::numeric_limits<double>::infinity();
    const double m = std::abs(num / den);
    return std::isnan(m) ? std::numeric_limits<double>::infinity() : m;
}

}  // namespace pyfda
