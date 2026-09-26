// Command line interface to the pyfda C++ core, used by tests/verify_scipy.py to
// compare the results with scipy. Reads one command per line from stdin and writes
// one JSON object per line to stdout.
//
// Commands (tokens separated by blanks, lists separated by commas):
//   iir <butter|cheby1|cheby2|ellip|bessel> <low|high|bandpass|bandstop> N rp rs Wn
//   ord <butter|cheby1|cheby2|ellip> wp ws gpass gstop
//   window <name> M par
//   firwin numtaps cutoff <0|1 pass_zero> <window> par
//   kaiserord ripple width
//   remezord <ichige|kaiser|herrmann> freqs amps rips
//   remez numtaps bands desired weight <bandpass|hilbert|differentiator> grid_density
//   filter <lfilter|filtfilt> b a x          (x as comma separated list)
//   sosfilter <sosfilt|sosfiltfilt> sos x    (sos as 6*n comma separated values)
//   freqz b a w | gd b a w | gdsos sos w
//   roots coeffs
//   spectrum x
//   csv <file>
//   design key=value ...                     (see parse_spec)
//   tojson unit key=value ...                filter file (JSON text) of a design
//   fromjson <file>                          design from a filter file
//   export <csv|matlab|c|python> key=value ...  exported coefficients as text
#include "../src/core/conversions.hpp"
#include "../src/core/data_io.hpp"
#include "../src/core/filter_design.hpp"
#include "../src/core/filter_io.hpp"
#include "../src/core/filtering.hpp"
#include "../src/core/fir_design.hpp"
#include "../src/core/iir_design.hpp"
#include "../src/core/poly.hpp"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

using namespace pyfda;

namespace {

std::string num(double v) {
    if (std::isnan(v)) return "NaN";
    if (std::isinf(v)) return v > 0 ? "Infinity" : "-Infinity";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

std::string arr(const Vec &v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + num(v[i]);
    return s + "]";
}

std::string carr(const CVec &v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? "," : "") + std::string("[") + num(v[i].real()) + "," + num(v[i].imag()) + "]";
    return s + "]";
}

std::string sosarr(const Sos &sos) {
    std::string s = "[";
    for (size_t i = 0; i < sos.size(); ++i) s += (i ? "," : "") + arr(Vec(sos[i].begin(), sos[i].end()));
    return s + "]";
}

std::string str(const std::string &s) {
    std::string o = "\"";
    for (char c : s) {
        if (c == '\n') o += "\\n";
        else if (c == '\t') o += "\\t";
        else if (c == '\r') o += "\\r";
        else {
            if (c == '"' || c == '\\') o += '\\';
            o += c;
        }
    }
    return o + "\"";
}

Vec list(const std::string &s) {
    Vec v;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!tok.empty()) v.push_back(std::stod(tok));
    return v;
}

Sos to_sos(const Vec &v) {
    Sos sos(v.size() / 6);
    for (size_t i = 0; i < sos.size(); ++i)
        for (int k = 0; k < 6; ++k) sos[i][k] = v[6 * i + k];
    return sos;
}

WindowType win_type(const std::string &n) {
    for (const auto &w : window_list()) {
        std::string name = w.name;
        for (char &c : name) c = char(std::tolower(c));
        if (name == n) return w.type;
    }
    throw DesignError("unknown window " + n);
}

FilterSpec parse_spec(std::istringstream &in) {
    FilterSpec s;
    std::string kv;
    while (in >> kv) {
        const auto p = kv.find('=');
        const std::string k = kv.substr(0, p), v = kv.substr(p + 1);
        if (k == "rt") s.rt = v == "LP" ? RespType::LP : v == "HP" ? RespType::HP : v == "BP" ? RespType::BP : RespType::BS;
        else if (k == "method") {
            static const std::map<std::string, DesignMethod> m = {
                {"butter", DesignMethod::Butter}, {"cheby1", DesignMethod::Cheby1}, {"cheby2", DesignMethod::Cheby2},
                {"ellip", DesignMethod::Ellip}, {"bessel", DesignMethod::Bessel}, {"firwin", DesignMethod::Firwin},
                {"equiripple", DesignMethod::Equiripple}};
            s.method = m.at(v);
        } else if (k == "fo") s.fo = v == "min" ? OrderMode::Min : OrderMode::Manual;
        else if (k == "N") s.N = std::stoi(v);
        else if (k == "f_s") s.f_s = std::stod(v);
        else if (k == "f_pb") s.f_pb = std::stod(v);
        else if (k == "f_pb2") s.f_pb2 = std::stod(v);
        else if (k == "f_sb") s.f_sb = std::stod(v);
        else if (k == "f_sb2") s.f_sb2 = std::stod(v);
        else if (k == "f_c") s.f_c = std::stod(v);
        else if (k == "f_c2") s.f_c2 = std::stod(v);
        else if (k == "A_PB") s.A_PB = std::stod(v);
        else if (k == "A_SB") s.A_SB = std::stod(v);
        else if (k == "W_PB") s.W_PB = std::stod(v);
        else if (k == "W_SB") s.W_SB = std::stod(v);
        else if (k == "window") s.window = win_type(v);
        else if (k == "win_par") s.win_par = std::stod(v);
        else if (k == "alg") s.order_alg = v == "kaiser" ? RemezAlg::Kaiser : v == "herrmann" ? RemezAlg::Herrmann : RemezAlg::Ichige;
        else throw DesignError("unknown key " + k);
    }
    return s;
}

std::string run(const std::string &line) {
    std::istringstream in(line);
    std::string cmd;
    in >> cmd;
    if (cmd == "iir") {
        std::string ft, bt;
        int N;
        double rp, rs;
        std::string wn;
        in >> ft >> bt >> N >> rp >> rs >> wn;
        static const std::map<std::string, IirType> fm = {{"butter", IirType::Butter}, {"cheby1", IirType::Cheby1},
                                                          {"cheby2", IirType::Cheby2}, {"ellip", IirType::Ellip},
                                                          {"bessel", IirType::Bessel}};
        static const std::map<std::string, BType> bm = {{"low", BType::Lowpass}, {"high", BType::Highpass},
                                                        {"bandpass", BType::Bandpass}, {"bandstop", BType::Bandstop}};
        const Zpk zpk = iirfilter(N, list(wn), rp, rs, bm.at(bt), fm.at(ft));
        const Ba ba = zpk2tf(zpk);
        return "{\"z\":" + carr(zpk.z) + ",\"p\":" + carr(zpk.p) + ",\"k\":" + num(zpk.k) + ",\"b\":" + arr(ba.b) +
               ",\"a\":" + arr(ba.a) + ",\"sos\":" + sosarr(zpk2sos(zpk)) + "}";
    }
    if (cmd == "ord") {
        std::string ft, wp, ws;
        double gp, gs;
        in >> ft >> wp >> ws >> gp >> gs;
        OrdResult o;
        if (ft == "butter") o = buttord(list(wp), list(ws), gp, gs);
        else if (ft == "cheby1") o = cheb1ord(list(wp), list(ws), gp, gs);
        else if (ft == "cheby2") o = cheb2ord(list(wp), list(ws), gp, gs);
        else o = ellipord(list(wp), list(ws), gp, gs);
        return "{\"N\":" + std::to_string(o.N) + ",\"Wn\":" + arr(o.Wn) + "}";
    }
    if (cmd == "window") {
        std::string n;
        int M;
        double par;
        in >> n >> M >> par;
        return "{\"w\":" + arr(get_window(win_type(n), M, par)) + "}";
    }
    if (cmd == "firwin") {
        int numtaps, pz;
        std::string cutoff, w;
        double par;
        in >> numtaps >> cutoff >> pz >> w >> par;
        return "{\"h\":" + arr(firwin(numtaps, list(cutoff), get_window(win_type(w), numtaps, par), pz != 0)) + "}";
    }
    if (cmd == "kaiserord") {
        double r, w, beta;
        int n;
        in >> r >> w;
        kaiserord(r, w, n, beta);
        return "{\"N\":" + std::to_string(n) + ",\"beta\":" + num(beta) + "}";
    }
    if (cmd == "remezord") {
        std::string alg, f, a, r;
        in >> alg >> f >> a >> r;
        const RemezAlg ra = alg == "kaiser" ? RemezAlg::Kaiser : alg == "herrmann" ? RemezAlg::Herrmann : RemezAlg::Ichige;
        const RemezOrd o = remezord(list(f), list(a), list(r), 1.0, ra);
        return "{\"N\":" + std::to_string(o.numtaps) + ",\"bands\":" + arr(o.bands) + ",\"weight\":" + arr(o.weight) + "}";
    }
    if (cmd == "remez") {
        int numtaps, grid;
        std::string b, d, w, t;
        in >> numtaps >> b >> d >> w >> t >> grid;
        const RemezType rt = t == "hilbert" ? RemezType::Hilbert : t == "differentiator" ? RemezType::Differentiator : RemezType::Bandpass;
        return "{\"h\":" + arr(remez(numtaps, list(b), list(d), list(w), rt, 1.0, 25, grid)) + "}";
    }
    if (cmd == "filter") {
        std::string m, b, a, x;
        in >> m >> b >> a >> x;
        const Vec y = m == "filtfilt" ? filtfilt(list(b), list(a), list(x)) : lfilter(list(b), list(a), list(x));
        return "{\"y\":" + arr(y) + "}";
    }
    if (cmd == "sosfilter") {
        std::string m, s, x;
        in >> m >> s >> x;
        const Sos sos = to_sos(list(s));
        const Vec y = m == "sosfiltfilt" ? sosfiltfilt(sos, list(x)) : sosfilt(sos, list(x));
        return "{\"y\":" + arr(y) + "}";
    }
    if (cmd == "freqz" || cmd == "gd") {
        std::string b, a, w;
        in >> b >> a >> w;
        if (cmd == "gd") return "{\"gd\":" + arr(group_delay(Ba{list(b), list(a)}, list(w))) + "}";
        return "{\"H\":" + carr(freqz(Ba{list(b), list(a)}, list(w))) + "}";
    }
    if (cmd == "gdsos") {
        std::string s, w;
        in >> s >> w;
        return "{\"gd\":" + arr(group_delay(to_sos(list(s)), list(w))) + "}";
    }
    if (cmd == "roots") {
        std::string c;
        in >> c;
        return "{\"r\":" + carr(roots(list(c))) + "}";
    }
    if (cmd == "spectrum") {
        std::string x;
        in >> x;
        return "{\"A\":" + arr(amplitude_spectrum(list(x))) + "}";
    }
    if (cmd == "csv") {
        std::string f;
        std::getline(in >> std::ws, f);
        const DataTable t = load_data_file(f);
        std::string names = "[";
        for (size_t i = 0; i < t.names.size(); ++i) names += (i ? "," : "") + str(t.names[i]);
        names += "]";
        return "{\"rows\":" + std::to_string(t.n_rows) + ",\"cols\":" + std::to_string(t.n_cols) +
               ",\"names\":" + names + ",\"values\":" + arr(t.values) + ",\"fs\":" + num(t.fs) + "}";
    }
    if (cmd == "design") {
        const FilterDesign d = design_filter(parse_spec(in));
        return "{\"N\":" + std::to_string(d.spec.N) + ",\"f_c\":" + num(d.spec.f_c) + ",\"f_c2\":" + num(d.spec.f_c2) +
               ",\"W_PB\":" + num(d.spec.W_PB) + ",\"W_SB\":" + num(d.spec.W_SB) + ",\"win_par\":" + num(d.spec.win_par) +
               ",\"b\":" + arr(d.ba.b) + ",\"a\":" + arr(d.ba.a) + ",\"sos\":" + sosarr(d.sos) + "}";
    }
    if (cmd == "tojson") {
        std::string unit;
        in >> unit;
        return "{\"json\":" + str(filter_to_json(design_filter(parse_spec(in)), unit)) + "}";
    }
    if (cmd == "fromjson") {
        std::string f;
        std::getline(in >> std::ws, f);
        const FilterFile ff = load_filter(f);
        const FilterDesign d = design_filter(ff.spec);
        return "{\"unit\":" + str(ff.unit) + ",\"rt\":" + str(resp_type_key(ff.spec.rt)) + ",\"method\":" +
               str(method_key(ff.spec.method)) + ",\"N\":" + std::to_string(d.spec.N) + ",\"b\":" + arr(d.ba.b) +
               ",\"a\":" + arr(d.ba.a) + ",\"sos\":" + sosarr(d.sos) + ",\"file_b\":" + arr(ff.ba.b) +
               ",\"file_a\":" + arr(ff.ba.a) + ",\"file_sos\":" + sosarr(ff.sos) + "}";
    }
    if (cmd == "export") {
        std::string f;
        in >> f;
        static const std::map<std::string, CoeffFormat> fm = {{"csv", CoeffFormat::Csv}, {"matlab", CoeffFormat::Matlab},
                                                             {"c", CoeffFormat::CHeader}, {"python", CoeffFormat::Python}};
        return "{\"text\":" + str(export_coeffs(design_filter(parse_spec(in)), fm.at(f), "lp filter")) + "}";
    }
    throw DesignError("unknown command " + cmd);
}

}  // namespace

int main() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty()) continue;
        try {
            std::cout << run(line) << "\n";
        } catch (const std::exception &e) {
            std::cout << "{\"error\":" << str(e.what()) << "}\n";
        }
        std::cout.flush();
    }
    return 0;
}
