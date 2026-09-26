#include "filter_io.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <memory>
#include <sstream>
#include <vector>

namespace pyfda {

// ---------------------------------------------------------------------------
// names

const char *resp_type_key(RespType rt) {
    switch (rt) {
    case RespType::LP: return "LP";
    case RespType::HP: return "HP";
    case RespType::BP: return "BP";
    case RespType::BS: return "BS";
    }
    return "";
}

const char *method_key(DesignMethod m) {
    switch (m) {
    case DesignMethod::Butter: return "butter";
    case DesignMethod::Cheby1: return "cheby1";
    case DesignMethod::Cheby2: return "cheby2";
    case DesignMethod::Ellip: return "ellip";
    case DesignMethod::Bessel: return "bessel";
    case DesignMethod::Firwin: return "firwin";
    case DesignMethod::Equiripple: return "equiripple";
    case DesignMethod::MovingAverage: return "ma";
    case DesignMethod::Delay: return "delay";
    case DesignMethod::ManualFIR: return "manual_fir";
    case DesignMethod::ManualIIR: return "manual_iir";
    }
    return "";
}

const char *order_alg_key(RemezAlg a) {
    switch (a) {
    case RemezAlg::Ichige: return "ichige";
    case RemezAlg::Kaiser: return "kaiser";
    case RemezAlg::Herrmann: return "herrmann";
    }
    return "";
}

RespType resp_type_from_key(const std::string &k) {
    for (RespType rt : {RespType::LP, RespType::HP, RespType::BP, RespType::BS})
        if (k == resp_type_key(rt)) return rt;
    throw DesignError("Unknown response type '" + k + "'");
}

DesignMethod method_from_key(const std::string &k) {
    for (DesignMethod m : {DesignMethod::Butter, DesignMethod::Cheby1, DesignMethod::Cheby2, DesignMethod::Ellip,
                           DesignMethod::Bessel, DesignMethod::Firwin, DesignMethod::Equiripple,
                           DesignMethod::MovingAverage, DesignMethod::Delay, DesignMethod::ManualFIR,
                           DesignMethod::ManualIIR})
        if (k == method_key(m)) return m;
    throw DesignError("Unknown design method '" + k + "'");
}

RemezAlg order_alg_from_key(const std::string &k) {
    for (RemezAlg a : {RemezAlg::Ichige, RemezAlg::Kaiser, RemezAlg::Herrmann})
        if (k == order_alg_key(a)) return a;
    throw DesignError("Unknown order estimation algorithm '" + k + "'");
}

namespace {
std::string lower(std::string s) {
    for (char &c : s) c = char(std::tolower(static_cast<unsigned char>(c)));
    return s;
}
}  // namespace

WindowType window_from_name(const std::string &n) {
    for (const auto &w : window_list())
        if (lower(w.name) == lower(n)) return w.type;
    throw DesignError("Unknown window '" + n + "'");
}

const char *window_name(WindowType t) {
    for (const auto &w : window_list())
        if (w.type == t) return w.name;
    return "";
}

// ---------------------------------------------------------------------------
// minimal JSON reader / writer

namespace {

std::string num(double v) {
    if (!std::isfinite(v)) return "null";
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

std::string quote(const std::string &s) {
    std::string o = "\"";
    for (char c : s) {
        switch (c) {
        case '"': o += "\\\""; break;
        case '\\': o += "\\\\"; break;
        case '\n': o += "\\n"; break;
        case '\t': o += "\\t"; break;
        default:
            if (static_cast<unsigned char>(c) < 0x20) {
                char buf[8];
                std::snprintf(buf, sizeof buf, "\\u%04x", c);
                o += buf;
            } else {
                o += c;
            }
        }
    }
    return o + "\"";
}

std::string vec(const Vec &v) {
    std::string s = "[";
    for (size_t i = 0; i < v.size(); ++i) s += (i ? ", " : "") + num(v[i]);
    return s + "]";
}

struct Json {
    enum Type { Null, Bool, Number, String, Array, Object } type = Null;
    double number = 0;
    bool boolean = false;
    std::string string;
    std::vector<Json> array;
    std::vector<std::pair<std::string, Json>> object;

    const Json *get(const std::string &key) const {
        for (const auto &kv : object)
            if (kv.first == key) return &kv.second;
        return nullptr;
    }
};

class Parser {
public:
    explicit Parser(const std::string &t) : m_t(t) {}

    Json parse() {
        Json v = value();
        ws();
        if (m_i != m_t.size()) fail("unexpected characters after the end");
        return v;
    }

private:
    [[noreturn]] void fail(const std::string &msg) const {
        // line number for the error message
        int line = 1;
        for (size_t i = 0; i < m_i && i < m_t.size(); ++i) line += m_t[i] == '\n';
        throw DesignError("Invalid filter file (line " + std::to_string(line) + "): " + msg);
    }
    void ws() {
        while (m_i < m_t.size() && std::isspace(static_cast<unsigned char>(m_t[m_i]))) ++m_i;
    }
    bool eat(char c) {
        ws();
        if (m_i < m_t.size() && m_t[m_i] == c) {
            ++m_i;
            return true;
        }
        return false;
    }
    void expect(char c) {
        if (!eat(c)) fail(std::string("expected '") + c + "'");
    }
    bool word(const char *w) {
        const std::string s(w);
        if (m_t.compare(m_i, s.size(), s) == 0) {
            m_i += s.size();
            return true;
        }
        return false;
    }

    Json value() {
        if (++m_depth > 64) fail("nested too deeply");
        ws();
        if (m_i >= m_t.size()) fail("unexpected end of file");
        Json v;
        const char c = m_t[m_i];
        if (c == '{') {
            ++m_i;
            v.type = Json::Object;
            if (!eat('}')) {
                do {
                    ws();
                    if (m_i >= m_t.size() || m_t[m_i] != '"') fail("expected a key");
                    std::string k = str();
                    expect(':');
                    v.object.emplace_back(std::move(k), value());
                } while (eat(','));
                expect('}');
            }
        } else if (c == '[') {
            ++m_i;
            v.type = Json::Array;
            if (!eat(']')) {
                do v.array.push_back(value());
                while (eat(','));
                expect(']');
            }
        } else if (c == '"') {
            v.type = Json::String;
            v.string = str();
        } else if (word("true")) {
            v.type = Json::Bool;
            v.boolean = true;
        } else if (word("false")) {
            v.type = Json::Bool;
        } else if (word("null")) {
            v.type = Json::Null;
        } else if (word("NaN")) {
            v.type = Json::Number;
            v.number = NAN;
        } else {
            const char *start = m_t.c_str() + m_i;
            char *end = nullptr;
            v.number = std::strtod(start, &end);
            if (end == start) fail("invalid value");
            v.type = Json::Number;
            m_i += size_t(end - start);
        }
        --m_depth;
        return v;
    }

    std::string str() {
        ++m_i;  // opening quote
        std::string s;
        while (m_i < m_t.size() && m_t[m_i] != '"') {
            char c = m_t[m_i++];
            if (c == '\\' && m_i < m_t.size()) {
                c = m_t[m_i++];
                switch (c) {
                case 'n': s += '\n'; break;
                case 't': s += '\t'; break;
                case 'r': s += '\r'; break;
                case 'b': s += '\b'; break;
                case 'f': s += '\f'; break;
                case 'u': {
                    if (m_i + 4 > m_t.size()) fail("invalid escape");
                    const unsigned cp = unsigned(std::stoul(m_t.substr(m_i, 4), nullptr, 16));
                    m_i += 4;
                    // UTF-8 (surrogate pairs are not needed for filter files)
                    if (cp < 0x80) s += char(cp);
                    else if (cp < 0x800) {
                        s += char(0xC0 | (cp >> 6));
                        s += char(0x80 | (cp & 0x3F));
                    } else {
                        s += char(0xE0 | (cp >> 12));
                        s += char(0x80 | ((cp >> 6) & 0x3F));
                        s += char(0x80 | (cp & 0x3F));
                    }
                    break;
                }
                default: s += c;
                }
            } else {
                s += c;
            }
        }
        if (m_i >= m_t.size()) fail("unterminated string");
        ++m_i;
        return s;
    }

    const std::string &m_t;
    size_t m_i = 0;
    int m_depth = 0;
};

double as_num(const Json &v, const std::string &key) {
    if (v.type == Json::Null) return NAN;
    if (v.type != Json::Number) throw DesignError("Invalid filter file: '" + key + "' must be a number");
    return v.number;
}

std::string as_str(const Json &v, const std::string &key) {
    if (v.type != Json::String) throw DesignError("Invalid filter file: '" + key + "' must be a string");
    return v.string;
}

Vec as_vec(const Json &v, const std::string &key) {
    if (v.type != Json::Array) throw DesignError("Invalid filter file: '" + key + "' must be a list");
    Vec out;
    for (const Json &e : v.array) out.push_back(as_num(e, key));
    return out;
}

const char *FORMAT_ID = "pyfda_cpp filter";
const int FORMAT_VERSION = 1;

}  // namespace

std::string filter_to_json(const FilterDesign &d, const std::string &unit, const FxSpec *fx, bool fx_sim) {
    const FilterSpec &s = d.spec;
    std::ostringstream o;
    o << "{\n  \"format\": " << quote(FORMAT_ID) << ",\n  \"version\": " << FORMAT_VERSION << ",\n";
    o << "  \"info\": " << quote(d.info) << ",\n";
    o << "  \"spec\": {\n";
    o << "    \"rt\": " << quote(resp_type_key(s.rt)) << ",\n";
    o << "    \"method\": " << quote(method_key(s.method)) << ",\n";
    o << "    \"order\": " << quote(s.fo == OrderMode::Min ? "min" : "manual") << ",\n";
    o << "    \"N\": " << s.N << ",\n";
    o << "    \"unit\": " << quote(unit) << ",\n";
    o << "    \"f_S\": " << num(s.f_s) << ",\n";
    o << "    \"F_PB\": " << num(s.f_pb) << ", \"F_PB2\": " << num(s.f_pb2) << ",\n";
    o << "    \"F_SB\": " << num(s.f_sb) << ", \"F_SB2\": " << num(s.f_sb2) << ",\n";
    o << "    \"F_C\": " << num(s.f_c) << ", \"F_C2\": " << num(s.f_c2) << ",\n";
    o << "    \"A_PB\": " << num(s.A_PB) << ", \"A_SB\": " << num(s.A_SB) << ",\n";
    o << "    \"W_PB\": " << num(s.W_PB) << ", \"W_SB\": " << num(s.W_SB) << ",\n";
    o << "    \"window\": " << quote(window_name(s.window)) << ", \"win_par\": " << num(s.win_par) << ",\n";
    o << "    \"order_alg\": " << quote(order_alg_key(s.order_alg)) << ",\n";
    o << "    \"grid_density\": " << s.grid_density << ",\n";
    o << "    \"ma_stages\": " << s.ma_stages << ", \"ma_norm\": " << (s.ma_norm ? "true" : "false") << "\n";
    o << "  },\n";
    if (fx) {
        auto q = [](const QFormat &f) {
            return "{\"WI\": " + std::to_string(f.WI) + ", \"WF\": " + std::to_string(f.WF) + ", \"quant\": " +
                   quote(quant_key(f.quant)) + ", \"ovfl\": " + quote(ovfl_key(f.ovfl)) + "}";
        };
        o << "  \"fixpoint\": {\n    \"simulate\": " << (fx_sim ? "true" : "false") << ",\n";
        o << "    \"coeff_auto\": " << (fx->coeff_auto ? "true" : "false") << ", \"acc_auto\": "
          << (fx->acc_auto ? "true" : "false") << ",\n";
        o << "    \"QI\": " << q(fx->qi) << ",\n    \"QCB\": " << q(fx->qcb) << ",\n    \"QCA\": " << q(fx->qca)
          << ",\n    \"QACC\": " << q(fx->qacc) << ",\n    \"QO\": " << q(fx->qo) << "\n  },\n";
    }
    if (is_manual(s.method) && s.manual_from_zpk) {
        // poles / zeros entered by hand are the reference, store them exactly
        auto cv = [](const CVec &v) {
            std::string r = "[";
            for (size_t i = 0; i < v.size(); ++i)
                r += (i ? ", " : "") + std::string("[") + num(v[i].real()) + ", " + num(v[i].imag()) + "]";
            return r + "]";
        };
        o << "  \"zpk\": {\"z\": " << cv(s.manual_zpk.z) << ",\n          \"p\": " << cv(s.manual_zpk.p)
          << ",\n          \"k\": " << num(s.manual_zpk.k) << "},\n";
    }
    o << "  \"b\": " << vec(d.ba.b) << ",\n";
    o << "  \"a\": " << vec(d.ba.a) << ",\n";
    o << "  \"sos\": [";
    for (size_t i = 0; i < d.sos.size(); ++i)
        o << (i ? ",\n          " : "") << vec(Vec(d.sos[i].begin(), d.sos[i].end()));
    o << "]\n}\n";
    return o.str();
}

FilterFile filter_from_json(const std::string &text) {
    Json root = Parser(text).parse();
    if (root.type != Json::Object) throw DesignError("Invalid filter file: expected a JSON object");
    const Json *fmt = root.get("format");
    if (!fmt || fmt->type != Json::String || fmt->string != FORMAT_ID)
        throw DesignError("Not a pyfda_cpp filter file (missing \"format\": \"pyfda_cpp filter\")");
    if (const Json *v = root.get("version"); v && as_num(*v, "version") > FORMAT_VERSION)
        throw DesignError("The filter file was written by a newer version of pyfda_cpp");
    const Json *spec = root.get("spec");
    if (!spec || spec->type != Json::Object) throw DesignError("Invalid filter file: missing \"spec\"");

    FilterFile f;
    FilterSpec &s = f.spec;
    const std::map<std::string, double *> nums = {
        {"f_S", &s.f_s},   {"F_PB", &s.f_pb}, {"F_PB2", &s.f_pb2}, {"F_SB", &s.f_sb}, {"F_SB2", &s.f_sb2},
        {"F_C", &s.f_c},   {"F_C2", &s.f_c2}, {"A_PB", &s.A_PB},   {"A_SB", &s.A_SB}, {"W_PB", &s.W_PB},
        {"W_SB", &s.W_SB}, {"win_par", &s.win_par}};
    for (const auto &kv : spec->object) {
        const std::string &k = kv.first;
        const Json &v = kv.second;
        if (auto it = nums.find(k); it != nums.end()) {
            *it->second = as_num(v, k);
            if (!std::isfinite(*it->second)) throw DesignError("Invalid filter file: '" + k + "' is not finite");
        } else if (k == "rt") s.rt = resp_type_from_key(as_str(v, k));
        else if (k == "method") s.method = method_from_key(as_str(v, k));
        else if (k == "order") {
            const std::string o = as_str(v, k);
            if (o != "min" && o != "manual") throw DesignError("Invalid filter file: order must be 'min' or 'manual'");
            s.fo = o == "min" ? OrderMode::Min : OrderMode::Manual;
        } else if (k == "N") s.N = int(as_num(v, k));
        else if (k == "grid_density") s.grid_density = int(as_num(v, k));
        else if (k == "ma_stages") s.ma_stages = int(as_num(v, k));
        else if (k == "ma_norm") {
            if (v.type != Json::Bool) throw DesignError("Invalid filter file: 'ma_norm' must be true or false");
            s.ma_norm = v.boolean;
        }
        else if (k == "unit") f.unit = as_str(v, k);
        else if (k == "window") s.window = window_from_name(as_str(v, k));
        else if (k == "order_alg") s.order_alg = order_alg_from_key(as_str(v, k));
        // unknown keys are ignored for forward compatibility
    }
    if (f.unit != "f_S" && f.unit != "f_Ny" && f.unit != "mHz" && f.unit != "Hz" && f.unit != "kHz" &&
        f.unit != "MHz" && f.unit != "GHz")
        throw DesignError("Invalid filter file: unknown unit '" + f.unit + "'");
    if (!(s.f_s > 0)) throw DesignError("Invalid filter file: f_S must be > 0");
    if (s.N < (is_manual(s.method) ? 0 : 1)) throw DesignError("Invalid filter file: N must be >= 1");
    if (s.ma_stages < 1 || s.ma_stages > 100) throw DesignError("Invalid filter file: ma_stages out of range");
    if (const Json *fx = root.get("fixpoint"); fx && fx->type == Json::Object) {
        f.has_fx = true;
        auto flag = [&](const char *k, bool &target) {
            if (const Json *v = fx->get(k)) {
                if (v->type != Json::Bool) throw DesignError(std::string("Invalid filter file: '") + k + "' must be true or false");
                target = v->boolean;
            }
        };
        flag("simulate", f.fx_sim);
        flag("coeff_auto", f.fx.coeff_auto);
        flag("acc_auto", f.fx.acc_auto);
        const std::map<std::string, QFormat *> qs = {
            {"QI", &f.fx.qi}, {"QCB", &f.fx.qcb}, {"QCA", &f.fx.qca}, {"QACC", &f.fx.qacc}, {"QO", &f.fx.qo}};
        for (const auto &kv : qs) {
            const Json *q = fx->get(kv.first);
            if (!q) continue;
            if (q->type != Json::Object) throw DesignError("Invalid filter file: '" + kv.first + "' must be an object");
            if (const Json *v = q->get("WI")) kv.second->WI = int(as_num(*v, "WI"));
            if (const Json *v = q->get("WF")) kv.second->WF = int(as_num(*v, "WF"));
            if (const Json *v = q->get("quant")) kv.second->quant = quant_from_key(as_str(*v, "quant"));
            if (const Json *v = q->get("ovfl")) kv.second->ovfl = ovfl_from_key(as_str(*v, "ovfl"));
            if (kv.second->WI < 0 || kv.second->WI > 32 || kv.second->WF < 0 || kv.second->WF > 60)
                throw DesignError("Invalid filter file: word lengths of '" + kv.first + "' out of range");
        }
    }
    if (const Json *b = root.get("b")) f.ba.b = as_vec(*b, "b");
    if (const Json *a = root.get("a")) f.ba.a = as_vec(*a, "a");
    if (const Json *sos = root.get("sos"); sos && sos->type == Json::Array) {
        for (const Json &row : sos->array) {
            const Vec r = as_vec(row, "sos");
            if (r.size() != 6) throw DesignError("Invalid filter file: every sos row needs 6 values");
            std::array<double, 6> sec;
            std::copy(r.begin(), r.end(), sec.begin());
            f.sos.push_back(sec);
        }
    }
    if (s.method == DesignMethod::ManualFIR || s.method == DesignMethod::ManualIIR) {
        if (const Json *zpk = root.get("zpk"); zpk && zpk->type == Json::Object) {
            auto cv = [](const Json *j, const char *key) {
                CVec out;
                if (!j) return out;
                if (j->type != Json::Array) throw DesignError(std::string("Invalid filter file: '") + key + "' must be an array");
                for (const Json &e : j->array) {
                    const Vec r = as_vec(e, key);
                    if (r.size() != 2 || !std::isfinite(r[0]) || !std::isfinite(r[1]))
                        throw DesignError(std::string("Invalid filter file: '") + key + "' needs [re, im] pairs");
                    out.emplace_back(r[0], r[1]);
                }
                return out;
            };
            s.manual_zpk.z = cv(zpk->get("z"), "z");
            s.manual_zpk.p = cv(zpk->get("p"), "p");
            s.manual_zpk.k = zpk->get("k") ? as_num(*zpk->get("k"), "k") : 1.0;
            s.manual_from_zpk = true;
        } else {
            if (f.ba.b.empty() || f.ba.a.empty()) throw DesignError("Invalid filter file: manual filter without coefficients");
            for (double v : f.ba.b)
                if (!std::isfinite(v)) throw DesignError("Invalid filter file: 'b' is not finite");
            for (double v : f.ba.a)
                if (!std::isfinite(v)) throw DesignError("Invalid filter file: 'a' is not finite");
            s.manual_ba = f.ba;
            s.manual_from_zpk = false;
        }
    }
    return f;
}

namespace {
std::string read_file(const std::string &fn) {
    std::ifstream in(fn, std::ios::binary);
    if (!in) throw DesignError("Couldn't open '" + fn + "'");
    std::ostringstream ss;
    ss << in.rdbuf();
    std::string t = ss.str();
    if (t.size() >= 3 && t.compare(0, 3, "\xEF\xBB\xBF") == 0) t.erase(0, 3);  // UTF-8 BOM
    return t;
}

void write_file(const std::string &fn, const std::string &text) {
    std::ofstream out(fn, std::ios::binary);
    if (!out) throw DesignError("Couldn't write '" + fn + "'");
    out << text;
    if (!out) throw DesignError("Couldn't write '" + fn + "'");
}
}  // namespace

void save_filter(const std::string &file_name, const FilterDesign &d, const std::string &unit, const FxSpec *fx,
                 bool fx_sim) {
    write_file(file_name, filter_to_json(d, unit, fx, fx_sim));
}

FilterFile load_filter(const std::string &file_name) { return filter_from_json(read_file(file_name)); }

// ---------------------------------------------------------------------------
// coefficient export

const char *coeff_format_filter(CoeffFormat f) {
    switch (f) {
    case CoeffFormat::Csv: return "CSV (*.csv)";
    case CoeffFormat::Matlab: return "MATLAB / Octave (*.m)";
    case CoeffFormat::CHeader: return "C header (*.h)";
    case CoeffFormat::Python: return "Python / NumPy (*.py)";
    }
    return "";
}

const char *coeff_format_suffix(CoeffFormat f) {
    switch (f) {
    case CoeffFormat::Csv: return "csv";
    case CoeffFormat::Matlab: return "m";
    case CoeffFormat::CHeader: return "h";
    case CoeffFormat::Python: return "py";
    }
    return "";
}

namespace {

std::string identifier(const std::string &name) {
    std::string id;
    for (char c : name) id += std::isalnum(static_cast<unsigned char>(c)) ? c : '_';
    if (id.empty() || std::isdigit(static_cast<unsigned char>(id[0]))) id = "f_" + id;
    return id;
}

std::string upper(std::string s) {
    for (char &c : s) c = char(std::toupper(static_cast<unsigned char>(c)));
    return s;
}

// values separated by sep, line break (with prefix) after every `per_line` values
std::string join(const Vec &v, const std::string &sep, const std::string &brk, int per_line = 4) {
    std::string s;
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) s += (i % size_t(per_line) == 0) ? sep + brk : sep + " ";
        s += num(v[i]);
    }
    return s;
}

std::string description(const FilterDesign &d) {
    const FilterSpec &s = d.spec;
    std::ostringstream o;
    o << d.info << ", f_S = " << num(s.f_s);
    return o.str();
}

}  // namespace

std::string export_coeffs(const FilterDesign &d, CoeffFormat f, const std::string &name) {
    const Vec &b = d.ba.b, &a = d.ba.a;
    const std::string desc = description(d);
    const std::string id = identifier(name);
    std::ostringstream o;
    switch (f) {
    case CoeffFormat::Csv: {
        o << "b,a\n";
        const size_t n = std::max(b.size(), a.size());
        for (size_t i = 0; i < n; ++i)
            o << (i < b.size() ? num(b[i]) : "") << "," << (i < a.size() ? num(a[i]) : "") << "\n";
        break;
    }
    case CoeffFormat::Matlab: {
        o << "% " << desc << "\n% exported by pyfda_cpp; y = filter(b, a, x)";
        if (!d.sos.empty()) o << " or y = sosfilt(sos, x)";
        o << "\n";
        o << "b = [" << join(b, ",", "...\n     ") << "];\n";
        o << "a = [" << join(a, ",", "...\n     ") << "];\n";
        if (!d.sos.empty()) {
            o << "sos = [";
            for (size_t i = 0; i < d.sos.size(); ++i)
                o << (i ? ";\n       " : "") << join(Vec(d.sos[i].begin(), d.sos[i].end()), ",", "", 6);
            o << "];\n";
        }
        break;
    }
    case CoeffFormat::CHeader: {
        const std::string guard = upper(id) + "_H";
        o << "/* " << desc << "\n * exported by pyfda_cpp */\n";
        o << "#ifndef " << guard << "\n#define " << guard << "\n\n";
        o << "#define " << upper(id) << "_NB " << b.size() << "\n";
        o << "#define " << upper(id) << "_NA " << a.size() << "\n";
        if (!d.sos.empty()) o << "#define " << upper(id) << "_NSOS " << d.sos.size() << "\n";
        o << "\nstatic const double " << id << "_b[" << upper(id) << "_NB] = {\n    " << join(b, ",", "\n   ") << "};\n";
        o << "static const double " << id << "_a[" << upper(id) << "_NA] = {\n    " << join(a, ",", "\n   ") << "};\n";
        if (!d.sos.empty()) {
            o << "/* second-order sections: b0, b1, b2, a0, a1, a2 */\n";
            o << "static const double " << id << "_sos[" << upper(id) << "_NSOS][6] = {\n";
            for (size_t i = 0; i < d.sos.size(); ++i)
                o << "    {" << join(Vec(d.sos[i].begin(), d.sos[i].end()), ",", "", 6) << "}"
                  << (i + 1 < d.sos.size() ? "," : "") << "\n";
            o << "};\n";
        }
        o << "\n#endif /* " << guard << " */\n";
        break;
    }
    case CoeffFormat::Python: {
        o << "# " << desc << "\n# exported by pyfda_cpp; y = scipy.signal.lfilter(b, a, x)";
        if (!d.sos.empty()) o << " or sosfilt(sos, x)";
        o << "\nimport numpy as np\n\n";
        o << "fs = " << num(d.spec.f_s) << "\n";
        o << "b = np.array([" << join(b, ",", "\n              ") << "])\n";
        o << "a = np.array([" << join(a, ",", "\n              ") << "])\n";
        if (!d.sos.empty()) {
            o << "sos = np.array([";
            for (size_t i = 0; i < d.sos.size(); ++i)
                o << (i ? ",\n                " : "") << "[" << join(Vec(d.sos[i].begin(), d.sos[i].end()), ",", "", 6) << "]";
            o << "])\n";
        }
        break;
    }
    }
    return o.str();
}

void export_coeffs(const std::string &file_name, const FilterDesign &d, CoeffFormat f) {
    // name of the arrays in C headers from the file name
    std::string base = file_name;
    const auto slash = base.find_last_of("/\\");
    if (slash != std::string::npos) base = base.substr(slash + 1);
    const auto dot = base.find_last_of('.');
    if (dot != std::string::npos && dot > 0) base = base.substr(0, dot);
    write_file(file_name, export_coeffs(d, f, base));
}

}  // namespace pyfda
