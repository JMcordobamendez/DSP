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
//   stim n_end key=value ...                 stimulus (see parse_stim)
//   wfft <window> par x                      windowed FFT (pyfda transient analysis)
//   fixp WI WF quant ovfl x                  quantize, returns values and overflows
//   fxbase v W                               v in bin, oct, hex and csd
//   fxfir b x key=WI,WF,quant,ovfl ...       fixpoint FIR (keys qi qcb qacc qo)
//   fxsos sos x key=WI,WF,quant,ovfl ...     fixpoint SOS cascade (keys qi qcb qca qacc qo)
//   fxauto <fir|sos> coeffs key=... ...      automatic formats
//   coe radix b WI,WF,quant,ovfl             Xilinx COE file
//   vhdl <fir|sos> coeffs key=... ...        VHDL entity (name fir_filter / iir_filter)
#include "../src/core/filter_info.hpp"
#include "../src/core/conversions.hpp"
#include "../src/core/data_io.hpp"
#include "../src/core/filter_design.hpp"
#include "../src/core/filter_io.hpp"
#include "../src/core/filtering.hpp"
#include "../src/core/fixpoint.hpp"
#include "../src/core/hdl_export.hpp"
#include "../src/core/fir_design.hpp"
#include "../src/core/iir_design.hpp"
#include "../src/core/poly.hpp"
#include "../src/core/stimulus.hpp"

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
        if (name == n || name.substr(0, name.find(' ')) == n) return w.type;
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
                {"equiripple", DesignMethod::Equiripple}, {"ma", DesignMethod::MovingAverage},
                {"delay", DesignMethod::Delay}, {"manual", DesignMethod::ManualIIR}};
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
        else if (k == "stages") s.ma_stages = std::stoi(v);
        else if (k == "norm") s.ma_norm = v == "1";
        else if (k == "b") s.manual_ba.b = list(v);
        else if (k == "a") s.manual_ba.a = list(v);
        else if (k == "z" || k == "p") {  // re:im,re:im,...
            CVec c;
            std::stringstream ss(v);
            std::string tok;
            while (std::getline(ss, tok, ',')) {
                const auto col = tok.find(':');
                c.emplace_back(std::stod(tok.substr(0, col)), col == std::string::npos ? 0.0 : std::stod(tok.substr(col + 1)));
            }
            (k == "z" ? s.manual_zpk.z : s.manual_zpk.p) = c;
            s.manual_from_zpk = true;
        } else if (k == "k") s.manual_zpk.k = std::stod(v);
        else if (k == "alg") s.order_alg = v == "kaiser" ? RemezAlg::Kaiser : v == "herrmann" ? RemezAlg::Herrmann : RemezAlg::Ichige;
        else throw DesignError("unknown key " + k);
    }
    return s;
}

StimParams parse_stim(std::istringstream &in) {
    StimParams p;
    std::string kv;
    while (in >> kv) {
        const auto pos = kv.find('=');
        const std::string k = kv.substr(0, pos), v = kv.substr(pos + 1);
        const std::map<std::string, double *> nums = {{"a1", &p.a1},   {"a2", &p.a2},   {"f1", &p.f1},   {"f2", &p.f2},
                                                      {"phi1", &p.phi1}, {"phi2", &p.phi2}, {"t1", &p.t1}, {"t2", &p.t2},
                                                      {"tw", &p.tw},   {"bw1", &p.bw1}, {"bw2", &p.bw2}, {"duty", &p.duty},
                                                      {"noi", &p.noi}, {"dc", &p.dc}, {"a1_im", &p.a1_im},
                                                      {"a2_im", &p.a2_im}, {"noi_im", &p.noi_im}, {"dc_im", &p.dc_im}};
        if (auto it = nums.find(k); it != nums.end()) *it->second = std::stod(v);
        else if (k == "stim") {
            bool found = false;
            for (const auto &i : stim_list())
                if (v == i.key) {
                    p.stim = i.stim;
                    found = true;
                }
            if (!found) throw DesignError("unknown stimulus " + v);
        } else if (k == "n1") p.n1 = std::stoi(v);
        else if (k == "mls_b") p.mls_b = std::stoi(v);
        else if (k == "bl") p.bl = v == "1";
        else if (k == "chirp") {
            static const std::map<std::string, ChirpType> cm = {{"linear", ChirpType::Linear}, {"quadratic", ChirpType::Quadratic},
                                                               {"logarithmic", ChirpType::Logarithmic}, {"hyperbolic", ChirpType::Hyperbolic}};
            p.chirp = cm.at(v);
        } else if (k == "noise") {
            static const std::map<std::string, Noise> nm = {{"none", Noise::None}, {"gauss", Noise::Gauss}, {"uniform", Noise::Uniform},
                                                           {"randint", Noise::RandInt}, {"mls", Noise::MLS}, {"brownian", Noise::Brownian}};
            p.noise = nm.at(v);
        } else if (k == "fs") p.f_s = std::stod(v);
        else if (k == "formula") {  // the rest of the line
            std::string rest;
            std::getline(in, rest);
            p.formula = v + rest;
            p.stim = Stim::Formula;
            break;
        } else throw DesignError("unknown key " + k);
    }
    return p;
}

QFormat parse_q(const std::string &v) {
    std::stringstream ss(v);
    std::string wi, wf, q, o;
    std::getline(ss, wi, ',');
    std::getline(ss, wf, ',');
    std::getline(ss, q, ',');
    std::getline(ss, o, ',');
    return QFormat{std::stoi(wi), std::stoi(wf), quant_from_key(q), ovfl_from_key(o)};
}

FxSpec parse_fx(std::istringstream &in) {
    FxSpec s;
    s.acc_auto = s.coeff_auto = false;
    std::string kv;
    while (in >> kv) {
        const auto p = kv.find('=');
        const std::string k = kv.substr(0, p), v = kv.substr(p + 1);
        if (k == "qi") s.qi = parse_q(v);
        else if (k == "qcb") s.qcb = parse_q(v);
        else if (k == "qca") s.qca = parse_q(v);
        else if (k == "qacc") s.qacc = parse_q(v);
        else if (k == "qo") s.qo = parse_q(v);
        else if (k == "auto") s.acc_auto = s.coeff_auto = v == "1";
        else throw DesignError("unknown key " + k);
    }
    return s;
}

std::string qstr(const QFormat &q) { return "[" + std::to_string(q.WI) + "," + std::to_string(q.WF) + "]"; }

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
               ",\"names\":" + names + ",\"values\":" + arr(t.values) + ",\"imag\":" + arr(t.imag) +
               ",\"fs\":" + num(t.fs) + "}";
    }
    if (cmd == "design") {
        const FilterDesign d = design_filter(parse_spec(in));
        return "{\"N\":" + std::to_string(d.spec.N) + ",\"f_c\":" + num(d.spec.f_c) + ",\"f_c2\":" + num(d.spec.f_c2) +
               ",\"W_PB\":" + num(d.spec.W_PB) + ",\"W_SB\":" + num(d.spec.W_SB) + ",\"win_par\":" + num(d.spec.win_par) +
               ",\"b\":" + arr(d.ba.b) + ",\"a\":" + arr(d.ba.a) + ",\"sos\":" + sosarr(d.sos) +
               ",\"z\":" + carr(d.zpk.z) + ",\"p\":" + carr(d.zpk.p) + ",\"k\":" + num(d.zpk.k) +
               ",\"fir\":" + (d.fir ? "true" : "false") + ",\"method\":" + str(method_key(d.spec.method)) +
               ",\"info\":" + str(d.info) + "}";
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
               ",\"file_a\":" + arr(ff.ba.a) + ",\"file_sos\":" + sosarr(ff.sos) + ",\"has_fx\":" +
               (ff.has_fx ? "true" : "false") + ",\"fx_sim\":" + (ff.fx_sim ? "true" : "false") + ",\"qacc\":[" +
               std::to_string(ff.fx.qacc.WI) + "," + std::to_string(ff.fx.qacc.WF) + "," + str(quant_key(ff.fx.qacc.quant)) +
               "," + str(ovfl_key(ff.fx.qacc.ovfl)) + "]}";
    }
    if (cmd == "export") {
        std::string f;
        in >> f;
        static const std::map<std::string, CoeffFormat> fm = {{"csv", CoeffFormat::Csv}, {"matlab", CoeffFormat::Matlab},
                                                             {"c", CoeffFormat::CHeader}, {"python", CoeffFormat::Python}};
        return "{\"text\":" + str(export_coeffs(design_filter(parse_spec(in)), fm.at(f), "lp filter")) + "}";
    }
    if (cmd == "stim") {
        int n;
        in >> n;
        const StimParams p = parse_stim(in);
        bool cmplx = false;
        const CVec x = calc_stimulus_c(p, n, &cmplx);
        Vec xr, xi;
        for (const cplx &v : x) {
            xr.push_back(v.real());
            xi.push_back(v.imag());
        }
        return "{\"x\":" + arr(xr) + ",\"xi\":" + arr(xi) + ",\"complex\":" + (cmplx ? "true" : "false") +
               ",\"title\":" + str(stim_title(p)) + ",\"scale\":" + num(impulse_scale(p)) + "}";
    }
    if (cmd == "wprops") {  // wprops <window> <N> <par>
        std::string w;
        int n;
        double par;
        in >> w >> n >> par;
        const WindowProps p = window_props(fft_window(win_type(w), n, par));
        return "{\"cgain\":" + num(p.cgain) + ",\"nenbw\":" + num(p.nenbw) + ",\"scallop\":" + num(p.scallop_db) +
               ",\"bw3\":" + num(p.bw3_bins) + ",\"bw6\":" + num(p.bw6_bins) + ",\"sidelobe\":" + num(p.sidelobe_db) + "}";
    }
    if (cmd == "info") {
        const FilterDesign d = design_filter(parse_spec(in));
        const FilterInfo fi = filter_info(d);
        std::string bands = "[";
        for (size_t i = 0; i < fi.bands.size(); ++i) {
            const BandCheck &b = fi.bands[i];
            bands += std::string(i ? "," : "") + "{\"name\":" + str(b.name) + ",\"f0\":" + num(b.f0) + ",\"f1\":" +
                     num(b.f1) + ",\"pass\":" + (b.pass ? "true" : "false") + ",\"spec\":" + num(b.spec_db) +
                     ",\"achieved\":" + num(b.achieved_db) + ",\"ok\":" + (b.ok ? "true" : "false") + "}";
        }
        return "{\"order\":" + std::to_string(fi.order) + ",\"stable\":" + (fi.stable ? "true" : "false") +
               ",\"min_phase\":" + (fi.min_phase ? "true" : "false") + ",\"linear_phase\":" +
               (fi.linear_phase ? "true" : "false") + ",\"rmax\":" + num(fi.max_pole_radius) + ",\"b\":" + arr(d.ba.b) +
               ",\"a\":" + arr(d.ba.a) + ",\"bands\":" + bands + "]}";
    }
    if (cmd == "amp") {  // amp <to|from> <value> <dB|V|W> <fir 0/1> <pb 0/1>
        std::string dir, u;
        double v;
        int fir, pb;
        in >> dir >> v >> u >> fir >> pb;
        const AmpUnit au = u == "V" ? AmpUnit::V : u == "W" ? AmpUnit::W : AmpUnit::dB;
        return "{\"v\":" + num(dir == "to" ? amp_to_db(v, au, fir, pb) : amp_from_db(v, au, fir, pb)) + "}";
    }
    if (cmd == "spgr") {  // spgr <psd|magnitude|angle> <density 0/1> <fs> <window> <par> <nperseg> <noverlap> <x> [<x imag>]
        std::string mode, w, x, xi;
        int dens, nper, novl;
        double fs, par;
        in >> mode >> dens >> fs >> w >> par >> nper >> novl >> x >> xi;
        const SpgrMode m = mode == "psd" ? SpgrMode::PSD : mode == "magnitude" ? SpgrMode::Magnitude : SpgrMode::Angle;
        const Vec win = fft_window(win_type(w), nper, par);
        Spectrogram r;
        if (xi.empty()) {
            r = spectrogram(list(x), fs, win, novl, m, dens == 1);
        } else {
            const Vec re = list(x), im = list(xi);
            CVec xc(re.size());
            for (size_t i = 0; i < re.size(); ++i) xc[i] = cplx(re[i], im[i]);
            r = spectrogram(xc, fs, win, novl, m, dens == 1);
        }
        std::string s = "[";
        for (size_t i = 0; i < r.s.size(); ++i) s += (i ? "," : "") + arr(r.s[i]);
        return "{\"f\":" + arr(r.f) + ",\"t\":" + arr(r.t) + ",\"s\":" + s + "]}";
    }
    if (cmd == "wfft") {  // wfft <window> <par> <x> [<x imag>]
        std::string w, x, xi;
        double par;
        in >> w >> par >> x >> xi;
        const Vec xv = list(x);
        const Vec win = fft_window(win_type(w), int(xv.size()), par);
        if (!xi.empty()) {
            const Vec im = list(xi);
            CVec xc(xv.size());
            for (size_t i = 0; i < xv.size(); ++i) xc[i] = cplx(xv[i], im[i]);
            return "{\"X\":" + carr(windowed_fft(xc, win)) + ",\"shifted\":" + carr(fftshift(windowed_fft(xc, win))) + "}";
        }
        return "{\"win\":" + arr(win) + ",\"X\":" + carr(windowed_fft(xv, win)) + ",\"ssb\":" +
               carr(ssb_spectrum(windowed_fft(xv, win))) + ",\"cgain\":" + num(window_cgain(win)) + ",\"nenbw\":" +
               num(window_nenbw(win)) + "}";
    }
    if (cmd == "fixp") {
        std::string q, o, x;
        int wi, wf;
        in >> wi >> wf >> q >> o >> x;
        Quantizer Q(QFormat{wi, wf, quant_from_key(q), ovfl_from_key(o)});
        const Vec y = Q.fixp(list(x));
        return "{\"y\":" + arr(y) + ",\"n_over\":" + std::to_string(Q.overflows()) + "}";
    }
    if (cmd == "fxbase") {
        long long v;
        int W;
        in >> v >> W;
        return "{\"bin\":" + str(to_base(v, W, 2)) + ",\"oct\":" + str(to_base(v, W, 8)) + ",\"hex\":" +
               str(to_base(v, W, 16)) + ",\"csd\":" + str(to_csd(v, W)) + "}";
    }
    if (cmd == "fxfir" || cmd == "fxsos") {
        std::string c, x;
        in >> c >> x;
        const FxSpec s = parse_fx(in);
        const FxResult r = cmd == "fxfir" ? fx_filter_fir(list(c), s, list(x)) : fx_filter_sos(to_sos(list(c)), s, list(x));
        return "{\"y\":" + arr(r.y) + ",\"x_q\":" + arr(r.x_q) + ",\"b_q\":" + arr(r.b_q) + ",\"a_q\":" +
               arr(r.a_q) + ",\"ov\":[" + std::to_string(r.n_over_i) + "," + std::to_string(r.n_over_acc) + "," +
               std::to_string(r.n_over_o) + "," + std::to_string(r.n_over_coeff) + "]}";
    }
    if (cmd == "coe") {
        int radix;
        std::string b, q;
        in >> radix >> b >> q;
        return "{\"text\":" + str(export_coe(list(b), parse_q(q), radix)) + "}";
    }
    if (cmd == "vhdl") {
        std::string t, c;
        in >> t >> c;
        const FxSpec s = parse_fx(in);
        return "{\"text\":" + str(t == "fir" ? export_vhdl_fir(list(c), s) : export_vhdl_sos(to_sos(list(c)), s)) + "}";
    }
    if (cmd == "verilog") {
        std::string t, c;
        in >> t >> c;
        const FxSpec s = parse_fx(in);
        return "{\"text\":" + str(t == "fir" ? export_verilog_fir(list(c), s) : export_verilog_sos(to_sos(list(c)), s)) + "}";
    }
    if (cmd == "tb") {  // tb vhdl|verilog fir|sos <coeffs> <x or -> <fx spec>
        std::string lang, t, c, xs;
        in >> lang >> t >> c >> xs;
        const FxSpec s = parse_fx(in);
        const bool fir = t == "fir";
        const Vec b = fir ? list(c) : Vec();
        const Sos sos = fir ? Sos() : to_sos(list(c));
        const Vec x = xs == "-" ? hdl_test_stimulus(s, fir ? b.size() : 3 * sos.size()) : list(xs);
        const std::string name = fir ? "fir_filter" : "iir_filter";
        return "{\"text\":" + str(lang == "vhdl" ? export_vhdl_testbench(b, sos, s, x, name)
                                                  : export_verilog_testbench(b, sos, s, x, name)) + "}";
    }
    if (cmd == "fxauto") {
        std::string t, c;
        in >> t >> c;
        FxSpec s = parse_fx(in);
        s.acc_auto = s.coeff_auto = true;
        const bool fir = t == "fir";
        update_auto_formats(s, Ba{list(c), {1.0}}, fir ? Sos() : to_sos(list(c)), fir);
        return "{\"qcb\":" + qstr(s.qcb) + ",\"qca\":" + qstr(s.qca) + ",\"qacc\":" + qstr(s.qacc) + "}";
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
