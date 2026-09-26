#include "filter_design.hpp"

#include "conversions.hpp"

#include <cmath>
#include <sstream>

namespace pyfda {

const char *resp_type_name(RespType rt) {
    switch (rt) {
    case RespType::LP: return "Lowpass";
    case RespType::HP: return "Highpass";
    case RespType::BP: return "Bandpass";
    case RespType::BS: return "Bandstop";
    }
    return "";
}

const char *method_name(DesignMethod m) {
    switch (m) {
    case DesignMethod::Butter: return "Butterworth";
    case DesignMethod::Cheby1: return "Chebyshev 1";
    case DesignMethod::Cheby2: return "Chebyshev 2";
    case DesignMethod::Ellip: return "Elliptic";
    case DesignMethod::Bessel: return "Bessel";
    case DesignMethod::Firwin: return "Windowed FIR";
    case DesignMethod::Equiripple: return "Equiripple";
    }
    return "";
}

bool is_fir(DesignMethod m) { return m == DesignMethod::Firwin || m == DesignMethod::Equiripple; }

double fir_a_pb_lin(double A) {
    const double g = std::pow(10.0, A / 20.0);
    return (g - 1) / (g + 1);
}

double sb_lin(double A) { return std::pow(10.0, -A / 20.0); }

namespace {

BType btype(RespType rt) {
    switch (rt) {
    case RespType::LP: return BType::Lowpass;
    case RespType::HP: return BType::Highpass;
    case RespType::BP: return BType::Bandpass;
    case RespType::BS: return BType::Bandstop;
    }
    return BType::Lowpass;
}

void check_freqs(const FilterSpec &s) {
    if (!(s.f_s > 0)) throw DesignError("The sampling frequency f_S must be > 0.");
    auto chk = [&](double f, const char *name) {
        if (!(f > 0 && f < s.f_s / 2))
            throw DesignError(std::string(name) + " must be between 0 and f_S / 2.");
    };
    const bool two = s.rt == RespType::BP || s.rt == RespType::BS;
    if (s.fo == OrderMode::Min) {
        chk(s.f_pb, "F_PB");
        chk(s.f_sb, "F_SB");
        if (two) {
            chk(s.f_pb2, "F_PB2");
            chk(s.f_sb2, "F_SB2");
        }
        switch (s.rt) {
        case RespType::LP:
            if (!(s.f_pb < s.f_sb)) throw DesignError("Lowpass: F_PB < F_SB is required.");
            break;
        case RespType::HP:
            if (!(s.f_sb < s.f_pb)) throw DesignError("Highpass: F_SB < F_PB is required.");
            break;
        case RespType::BP:
            if (!(s.f_sb < s.f_pb && s.f_pb < s.f_pb2 && s.f_pb2 < s.f_sb2))
                throw DesignError("Bandpass: F_SB < F_PB < F_PB2 < F_SB2 is required.");
            break;
        case RespType::BS:
            if (!(s.f_pb < s.f_sb && s.f_sb < s.f_sb2 && s.f_sb2 < s.f_pb2))
                throw DesignError("Bandstop: F_PB < F_SB < F_SB2 < F_PB2 is required.");
            break;
        }
        if (!(s.A_PB > 0 && s.A_SB > s.A_PB)) throw DesignError("0 < A_PB < A_SB is required.");
    } else {
        chk(s.f_c, "F_C");
        if (two) {
            chk(s.f_c2, "F_C2");
            if (!(s.f_c < s.f_c2)) throw DesignError("F_C < F_C2 is required.");
        }
        if (s.N < 1) throw DesignError("The filter order must be >= 1.");
    }
}

FilterDesign design_iir(FilterSpec s) {
    const IirType ft = s.method == DesignMethod::Butter   ? IirType::Butter
                       : s.method == DesignMethod::Cheby1 ? IirType::Cheby1
                       : s.method == DesignMethod::Cheby2 ? IirType::Cheby2
                       : s.method == DesignMethod::Ellip  ? IirType::Ellip
                                                          : IirType::Bessel;
    const bool two = s.rt == RespType::BP || s.rt == RespType::BS;
    // normalize to f_Ny
    const double k = 2.0 / s.f_s;
    int N;
    Vec Wn;
    if (s.fo == OrderMode::Min) {
        Vec wp{s.f_pb * k}, ws{s.f_sb * k};
        if (two) {
            wp.push_back(s.f_pb2 * k);
            ws.push_back(s.f_sb2 * k);
        }
        OrdResult o;
        switch (ft) {
        case IirType::Butter:
        case IirType::Bessel: o = buttord(wp, ws, s.A_PB, s.A_SB); break;
        case IirType::Cheby1: o = cheb1ord(wp, ws, s.A_PB, s.A_SB); break;
        case IirType::Cheby2: o = cheb2ord(wp, ws, s.A_PB, s.A_SB); break;
        case IirType::Ellip: o = ellipord(wp, ws, s.A_PB, s.A_SB); break;
        }
        N = o.N;
        Wn = o.Wn;
        s.N = two ? 2 * N : N;
        s.f_c = Wn[0] / k;
        if (two) s.f_c2 = Wn[1] / k;
    } else {
        N = two ? s.N / 2 : s.N;
        Wn = {s.f_c * k};
        if (two) Wn.push_back(s.f_c2 * k);
    }
    if (N < 1) throw DesignError("The calculated filter order is < 1, check the specifications.");
    if (N > 50) throw DesignError("Filter order N = " + std::to_string(two ? 2 * N : N) +
                                  " is too high for a reasonable IIR design.");
    FilterDesign d;
    d.fir = false;
    d.zpk = iirfilter(N, Wn, s.A_PB, s.A_SB, btype(s.rt), ft);
    d.ba = zpk2tf(d.zpk);
    d.sos = zpk2sos(d.zpk);
    d.spec = s;
    return d;
}

FilterDesign design_fir(FilterSpec s) {
    const double f_pb = s.f_pb / s.f_s, f_sb = s.f_sb / s.f_s, f_pb2 = s.f_pb2 / s.f_s,
                 f_sb2 = s.f_sb2 / s.f_s;  // normalized to f_S
    const double a_pb = fir_a_pb_lin(s.A_PB), a_sb = sb_lin(s.A_SB);
    const bool min = s.fo == OrderMode::Min;
    int numtaps = s.N + 1;
    Vec h;

    if (s.method == DesignMethod::Firwin) {
        if (min) {
            if (s.window == WindowType::Kaiser) {
                const double delta_f = (s.rt == RespType::LP || s.rt == RespType::HP)
                                           ? std::fabs(f_sb - f_pb) * 2
                                           : std::min(std::fabs(f_pb - f_sb), std::fabs(f_sb2 - f_pb2)) * 2;
                double beta;
                kaiserord(s.A_SB, delta_f, numtaps, beta);
                s.win_par = beta;
            } else {
                RemezOrd o;
                switch (s.rt) {
                case RespType::LP: o = remezord({f_pb, f_sb}, {1, 0}, {a_pb, a_sb}, 1, s.order_alg); break;
                case RespType::HP: o = remezord({f_sb, f_pb}, {0, 1}, {a_sb, a_pb}, 1, s.order_alg); break;
                case RespType::BP:
                    o = remezord({f_sb, f_pb, f_pb2, f_sb2}, {0, 1, 0}, {a_sb, a_pb, a_sb}, 1, s.order_alg);
                    break;
                case RespType::BS:
                    o = remezord({f_pb, f_sb, f_sb2, f_pb2}, {1, 0, 1}, {a_pb, a_sb, a_pb}, 1, s.order_alg);
                    break;
                }
                numtaps = o.numtaps;
            }
            // cutoff in the middle of the transition band(s)
            s.f_c = (s.f_pb + s.f_sb) / 2;
            s.f_c2 = (s.f_pb2 + s.f_sb2) / 2;
        }
        // HP and BS need a passband at f_S / 2, i.e. an odd number of taps (even order)
        if ((s.rt == RespType::HP || s.rt == RespType::BS) && numtaps % 2 == 0) ++numtaps;
        if (numtaps < 2) numtaps = 2;
        if (numtaps > 5001) throw DesignError("Filter order N = " + std::to_string(numtaps - 1) + " is too high.");
        const Vec win = get_window(s.window, numtaps, s.win_par);
        const double fc = s.f_c / s.f_s * 2, fc2 = s.f_c2 / s.f_s * 2;  // normalized to f_Ny
        switch (s.rt) {
        case RespType::LP: h = firwin(numtaps, {fc}, win, true); break;
        case RespType::HP: h = firwin(numtaps, {fc}, win, false); break;
        case RespType::BP: h = firwin(numtaps, {fc, fc2}, win, false); break;
        case RespType::BS: h = firwin(numtaps, {fc, fc2}, win, true); break;
        }
    } else {  // equiripple
        Vec bands, desired, weight;
        RemezType type = RemezType::Bandpass;
        if (min) {
            RemezOrd o;
            switch (s.rt) {
            case RespType::LP: o = remezord({f_pb, f_sb}, {1, 0}, {a_pb, a_sb}, 1, s.order_alg); break;
            case RespType::HP: o = remezord({f_sb, f_pb}, {0, 1}, {a_sb, a_pb}, 1, s.order_alg); break;
            case RespType::BP:
                o = remezord({f_sb, f_pb, f_pb2, f_sb2}, {0, 1, 0}, {a_sb, a_pb, a_sb}, 1, s.order_alg);
                break;
            case RespType::BS:
                o = remezord({f_pb, f_sb, f_sb2, f_pb2}, {1, 0, 1}, {a_pb, a_sb, a_pb}, 1, s.order_alg);
                break;
            }
            numtaps = o.numtaps;
            bands = o.bands;
            desired = o.desired;
            weight = o.weight;
            // store pass / stop band weights for a subsequent manual design
            if (s.rt == RespType::LP || s.rt == RespType::BS) {
                s.W_PB = weight[0];
                s.W_SB = weight[1];
            } else {
                s.W_SB = weight[0];
                s.W_PB = weight[1];
            }
            if ((s.rt == RespType::HP || s.rt == RespType::BS) && numtaps % 2 == 0) ++numtaps;
        } else {
            switch (s.rt) {
            case RespType::LP:
                bands = {0, f_pb, f_sb, 0.5};
                desired = {1, 0};
                weight = {s.W_PB, s.W_SB};
                break;
            case RespType::HP:
                bands = {0, f_sb, f_pb, 0.5};
                desired = {0, 1};
                weight = {s.W_SB, s.W_PB};
                // odd order (even number of taps): antisymmetric (type IV) filter
                if (numtaps % 2 == 0) type = RemezType::Hilbert;
                break;
            case RespType::BP:
                bands = {0, f_sb, f_pb, f_pb2, f_sb2, 0.5};
                desired = {0, 1, 0};
                weight = {s.W_SB, s.W_PB, s.W_SB};
                break;
            case RespType::BS:
                bands = {0, f_pb, f_sb, f_sb2, f_pb2, 0.5};
                desired = {1, 0, 1};
                weight = {s.W_PB, s.W_SB, s.W_PB};
                if (numtaps % 2 == 0) ++numtaps;
                break;
            }
            // band edges have to be specified for equiripple filters also in manual mode
            FilterSpec chk = s;
            chk.fo = OrderMode::Min;
            check_freqs(chk);
        }
        if (numtaps < 3) numtaps = 3;
        if (numtaps > 2001) throw DesignError("Filter order N = " + std::to_string(numtaps - 1) + " is too high.");
        h = remez(numtaps, bands, desired, weight, type, 1.0, 25, s.grid_density);
    }
    s.N = numtaps - 1;
    FilterDesign d;
    d.fir = true;
    d.ba = {h, {1.0}};
    d.zpk = tf2zpk(d.ba);
    d.spec = s;
    return d;
}

}  // namespace

FilterDesign design_filter(const FilterSpec &spec) {
    check_freqs(spec);
    FilterDesign d = is_fir(spec.method) ? design_fir(spec) : design_iir(spec);
    std::ostringstream ss;
    ss << method_name(spec.method) << " " << resp_type_name(spec.rt) << ", N = " << d.spec.N;
    d.info = ss.str();
    return d;
}

}  // namespace pyfda
