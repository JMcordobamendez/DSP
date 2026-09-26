#include "expr.hpp"

#include <cctype>
#include <cmath>
#include <functional>

namespace pyfda {

struct ExprNode {
    enum Kind { Num, Var, Unary, Binary, Call, Imag } kind = Num;  // Imag: value * j
    double value = 0;
    int var = 0;
    std::string op;  // operator or function name
    std::vector<std::unique_ptr<ExprNode>> args;
};

namespace {

using Node = std::unique_ptr<ExprNode>;

Node make(ExprNode::Kind k, const std::string &op = {}) {
    auto n = std::make_unique<ExprNode>();
    n->kind = k;
    n->op = op;
    return n;
}

const std::map<std::string, int> FUNCS = {
    {"sin", 1},    {"cos", 1},     {"tan", 1},     {"arcsin", 1}, {"arccos", 1}, {"arctan", 1},
    {"arctan2", 2}, {"sinh", 1},   {"cosh", 1},    {"tanh", 1},   {"arcsinh", 1}, {"arccosh", 1},
    {"arctanh", 1}, {"exp", 1},    {"expm1", 1},   {"log", 1},    {"log10", 1},  {"log1p", 1},
    {"sqrt", 1},   {"abs", 1},     {"sign", 1},    {"floor", 1},  {"ceil", 1},   {"round", 1},
    {"where", 3},  {"minimum", 2}, {"maximum", 2}, {"real", 1},   {"imag", 1},   {"conj", 1},
    {"complex", 2}};

class Parser {
public:
    Parser(const std::string &t, const std::vector<std::string> &vars) : s(t), vars(vars) {}
    bool complex = false;  // imaginary literal or complex() found

    Node parse() {
        Node n = parse_or();
        skip();
        if (i < s.size()) error("unexpected '" + std::string(1, s[i]) + "'");
        return n;
    }

private:
    [[noreturn]] void error(const std::string &msg) const {
        throw DesignError("Formula, position " + std::to_string(i + 1) + ": " + msg);
    }
    void skip() {
        while (i < s.size() && std::isspace(static_cast<unsigned char>(s[i]))) ++i;
    }
    bool accept(const std::string &tok) {
        skip();
        if (s.compare(i, tok.size(), tok) == 0) {
            i += tok.size();
            return true;
        }
        return false;
    }
    Node binary(const std::string &op, Node a, Node b) {
        Node n = make(ExprNode::Binary, op);
        n->args.push_back(std::move(a));
        n->args.push_back(std::move(b));
        return n;
    }
    Node parse_or() {
        Node a = parse_and();
        while (accept("|")) a = binary("|", std::move(a), parse_and());
        return a;
    }
    Node parse_and() {
        Node a = parse_cmp();
        while (accept("&")) a = binary("&", std::move(a), parse_cmp());
        return a;
    }
    Node parse_cmp() {
        Node a = parse_add();
        for (;;) {
            std::string op;
            for (const char *o : {"<=", ">=", "==", "!=", "<", ">"})
                if (accept(o)) {
                    op = o;
                    break;
                }
            if (op.empty()) return a;
            a = binary(op, std::move(a), parse_add());
        }
    }
    Node parse_add() {
        Node a = parse_mul();
        for (;;) {
            if (accept("+")) a = binary("+", std::move(a), parse_mul());
            else if (accept("-")) a = binary("-", std::move(a), parse_mul());
            else return a;
        }
    }
    Node parse_mul() {
        Node a = parse_unary();
        for (;;) {
            skip();
            if (s.compare(i, 2, "**") == 0) return a;  // handled in parse_pow
            if (accept("*")) a = binary("*", std::move(a), parse_unary());
            else if (accept("/")) a = binary("/", std::move(a), parse_unary());
            else if (accept("%")) a = binary("%", std::move(a), parse_unary());
            else return a;
        }
    }
    Node parse_unary() {
        for (const char *o : {"-", "+", "~"})
            if (accept(o)) {
                Node n = make(ExprNode::Unary, o);
                n->args.push_back(parse_unary());
                return n;
            }
        return parse_pow();
    }
    Node parse_pow() {
        Node a = parse_atom();
        if (accept("**")) return binary("**", std::move(a), parse_unary());  // right associative, -2**2 = -4
        return a;
    }
    Node parse_atom() {
        skip();
        if (i >= s.size()) error("unexpected end of the formula");
        const char c = s[i];
        if (accept("(")) {
            Node n = parse_or();
            if (!accept(")")) error("')' expected");
            return n;
        }
        if (std::isdigit(static_cast<unsigned char>(c)) || c == '.') {
            size_t j = i;
            while (j < s.size() && (std::isdigit(static_cast<unsigned char>(s[j])) || s[j] == '.')) ++j;
            if (j < s.size() && (s[j] == 'e' || s[j] == 'E')) {
                size_t k = j + 1;
                if (k < s.size() && (s[k] == '+' || s[k] == '-')) ++k;
                if (k < s.size() && std::isdigit(static_cast<unsigned char>(s[k]))) {
                    j = k;
                    while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j]))) ++j;
                }
            }
            std::string num = s.substr(i, j - i);
            if (num == ".") error("invalid number");
            const bool imag = j < s.size() && (s[j] == 'j' || s[j] == 'J') &&
                              !(j + 1 < s.size() && (std::isalnum(static_cast<unsigned char>(s[j + 1])) || s[j + 1] == '_'));
            Node n = make(imag ? ExprNode::Imag : ExprNode::Num);
            complex = complex || imag;
            try {
                size_t used = 0;
                n->value = std::stod(num, &used);
                if (used != num.size()) error("invalid number '" + num + "'");
            } catch (const std::logic_error &) {
                error("invalid number '" + num + "'");
            }
            i = imag ? j + 1 : j;
            return n;
        }
        if (std::isalpha(static_cast<unsigned char>(c)) || c == '_') {
            size_t j = i;
            while (j < s.size() && (std::isalnum(static_cast<unsigned char>(s[j])) || s[j] == '_')) ++j;
            const std::string name = s.substr(i, j - i);
            const size_t start = i;
            i = j;
            if (accept("(")) {
                auto f = FUNCS.find(name);
                if (f == FUNCS.end()) {
                    i = start;
                    error("unknown function '" + name + "'");
                }
                Node n = make(ExprNode::Call, name);
                complex = complex || name == "complex";
                if (!accept(")")) {
                    do n->args.push_back(parse_or());
                    while (accept(","));
                    if (!accept(")")) error("')' expected");
                }
                if (int(n->args.size()) != f->second)
                    error(name + "() needs " + std::to_string(f->second) + " argument" + (f->second > 1 ? "s" : ""));
                return n;
            }
            for (size_t k = 0; k < vars.size(); ++k)
                if (vars[k] == name) {
                    Node n = make(ExprNode::Var);
                    n->var = int(k);
                    return n;
                }
            if (name == "j" || name == "J") {  // imaginary unit like in pyfda
                complex = true;
                Node n = make(ExprNode::Imag);
                n->value = 1;
                return n;
            }
            i = start;
            error("unknown name '" + name + "'");
        }
        error("unexpected '" + std::string(1, c) + "'");
    }

    const std::string &s;
    const std::vector<std::string> &vars;
    size_t i = 0;
};

double eval_node(const ExprNode &n, const std::vector<double> &v) {
    switch (n.kind) {
    case ExprNode::Num: return n.value;
    case ExprNode::Imag: return std::nan("");  // not used, complex expressions are evaluated by eval_c
    case ExprNode::Var: return v[size_t(n.var)];
    case ExprNode::Unary: {
        const double a = eval_node(*n.args[0], v);
        if (n.op == "-") return -a;
        if (n.op == "~") return a != 0 ? 0.0 : 1.0;
        return a;
    }
    case ExprNode::Binary: {
        const double a = eval_node(*n.args[0], v), b = eval_node(*n.args[1], v);
        const std::string &o = n.op;
        if (o == "+") return a + b;
        if (o == "-") return a - b;
        if (o == "*") return a * b;
        if (o == "/") return a / b;
        if (o == "%") return b == 0 ? std::nan("") : a - std::floor(a / b) * b;  // numpy / Python modulo
        if (o == "**") return std::pow(a, b);
        if (o == "<") return a < b;
        if (o == "<=") return a <= b;
        if (o == ">") return a > b;
        if (o == ">=") return a >= b;
        if (o == "==") return a == b;
        if (o == "!=") return a != b;
        if (o == "&") return (a != 0) && (b != 0);
        if (o == "|") return (a != 0) || (b != 0);
        return 0;
    }
    case ExprNode::Call: {
        const std::string &f = n.op;
        if (f == "where") return eval_node(*n.args[0], v) != 0 ? eval_node(*n.args[1], v) : eval_node(*n.args[2], v);
        const double a = eval_node(*n.args[0], v);
        if (f == "real") return a;
        if (f == "imag") return 0.0;
        if (f == "conj") return a;
        if (n.args.size() == 2) {
            const double b = eval_node(*n.args[1], v);
            if (f == "arctan2") return std::atan2(a, b);
            if (f == "minimum") return std::fmin(a, b);
            return std::fmax(a, b);
        }
        static const std::map<std::string, std::function<double(double)>> fn = {
            {"sin", [](double x) { return std::sin(x); }},     {"cos", [](double x) { return std::cos(x); }},
            {"tan", [](double x) { return std::tan(x); }},     {"arcsin", [](double x) { return std::asin(x); }},
            {"arccos", [](double x) { return std::acos(x); }}, {"arctan", [](double x) { return std::atan(x); }},
            {"sinh", [](double x) { return std::sinh(x); }},   {"cosh", [](double x) { return std::cosh(x); }},
            {"tanh", [](double x) { return std::tanh(x); }},   {"arcsinh", [](double x) { return std::asinh(x); }},
            {"arccosh", [](double x) { return std::acosh(x); }}, {"arctanh", [](double x) { return std::atanh(x); }},
            {"exp", [](double x) { return std::exp(x); }},     {"expm1", [](double x) { return std::expm1(x); }},
            {"log", [](double x) { return std::log(x); }},     {"log10", [](double x) { return std::log10(x); }},
            {"log1p", [](double x) { return std::log1p(x); }}, {"sqrt", [](double x) { return std::sqrt(x); }},
            {"abs", [](double x) { return std::fabs(x); }},
            {"sign", [](double x) { return x > 0 ? 1.0 : x < 0 ? -1.0 : x; }},
            {"floor", [](double x) { return std::floor(x); }}, {"ceil", [](double x) { return std::ceil(x); }},
            {"round", [](double x) { return std::nearbyint(x); }}};
        return fn.at(f)(a);
    }
    }
    return 0;
}

// complex power like numpy: repeated squaring for small integer exponents
cplx cpow(cplx a, cplx b) {
    if (b.imag() == 0 && b.real() == std::nearbyint(b.real()) && std::fabs(b.real()) <= 100) {
        long long e = (long long)(std::fabs(b.real()));
        cplx r(1, 0), x = a;
        while (e) {
            if (e & 1) r *= x;
            x *= x;
            e >>= 1;
        }
        return b.real() < 0 ? cplx(1, 0) / r : r;
    }
    if (a == cplx(0, 0)) return b.real() > 0 ? cplx(0, 0) : cplx(std::nan(""), std::nan(""));
    return std::pow(a, b);
}

cplx eval_c(const ExprNode &n, const std::vector<cplx> &v) {
    switch (n.kind) {
    case ExprNode::Num: return n.value;
    case ExprNode::Imag: return cplx(0, n.value);
    case ExprNode::Var: return v[size_t(n.var)];
    case ExprNode::Unary: {
        const cplx a = eval_c(*n.args[0], v);
        if (n.op == "-") return -a;
        if (n.op == "~") return a.real() != 0 ? 0.0 : 1.0;
        return a;
    }
    case ExprNode::Binary: {
        const cplx a = eval_c(*n.args[0], v), b = eval_c(*n.args[1], v);
        const double ar = a.real(), br = b.real();
        const std::string &o = n.op;
        if (o == "+") return a + b;
        if (o == "-") return a - b;
        if (o == "*") return a * b;
        if (o == "/") return a / b;
        if (o == "%") return br == 0 ? std::nan("") : ar - std::floor(ar / br) * br;
        if (o == "**") return cpow(a, b);
        if (o == "<") return ar < br;
        if (o == "<=") return ar <= br;
        if (o == ">") return ar > br;
        if (o == ">=") return ar >= br;
        if (o == "==") return a == b;
        if (o == "!=") return a != b;
        if (o == "&") return (ar != 0) && (br != 0);
        if (o == "|") return (ar != 0) || (br != 0);
        return 0.0;
    }
    case ExprNode::Call: {
        const std::string &f = n.op;
        if (f == "where") return eval_c(*n.args[0], v).real() != 0 ? eval_c(*n.args[1], v) : eval_c(*n.args[2], v);
        const cplx a = eval_c(*n.args[0], v);
        if (n.args.size() == 2) {
            const cplx b = eval_c(*n.args[1], v);
            if (f == "complex") return cplx(a.real(), b.real());
            if (f == "arctan2") return std::atan2(a.real(), b.real());
            if (f == "minimum") return std::fmin(a.real(), b.real());
            return std::fmax(a.real(), b.real());
        }
        static const std::map<std::string, std::function<cplx(cplx)>> fn = {
            {"sin", [](cplx x) { return std::sin(x); }},     {"cos", [](cplx x) { return std::cos(x); }},
            {"tan", [](cplx x) { return std::tan(x); }},     {"arcsin", [](cplx x) { return std::asin(x); }},
            {"arccos", [](cplx x) { return std::acos(x); }}, {"arctan", [](cplx x) { return std::atan(x); }},
            {"sinh", [](cplx x) { return std::sinh(x); }},   {"cosh", [](cplx x) { return std::cosh(x); }},
            {"tanh", [](cplx x) { return std::tanh(x); }},   {"arcsinh", [](cplx x) { return std::asinh(x); }},
            {"arccosh", [](cplx x) { return std::acosh(x); }}, {"arctanh", [](cplx x) { return std::atanh(x); }},
            {"exp", [](cplx x) { return std::exp(x); }},
            {"expm1", [](cplx x) { return x.imag() == 0 ? cplx(std::expm1(x.real())) : std::exp(x) - 1.0; }},
            {"log", [](cplx x) { return std::log(x); }},     {"log10", [](cplx x) { return std::log10(x); }},
            {"log1p", [](cplx x) { return x.imag() == 0 && x.real() > -1 ? cplx(std::log1p(x.real())) : std::log(1.0 + x); }},
            {"sqrt", [](cplx x) { return std::sqrt(x); }},
            {"abs", [](cplx x) { return cplx(std::abs(x)); }},
            {"sign", [](cplx x) { return cplx(x.real() > 0 ? 1.0 : x.real() < 0 ? -1.0 : x.real()); }},
            {"floor", [](cplx x) { return cplx(std::floor(x.real())); }},
            {"ceil", [](cplx x) { return cplx(std::ceil(x.real())); }},
            {"round", [](cplx x) { return cplx(std::nearbyint(x.real()), std::nearbyint(x.imag())); }},
            {"real", [](cplx x) { return cplx(x.real()); }}, {"imag", [](cplx x) { return cplx(x.imag()); }},
            {"conj", [](cplx x) { return std::conj(x); }}};
        return fn.at(f)(a);
    }
    }
    return 0.0;
}

}  // namespace

Expr::Expr(const std::string &text, const std::vector<std::string> &vars) {
    std::string t = text;
    // like pyfda's safe_numexpr_eval: ',' as decimal separator (function arguments
    // are separated by ', ' or ',' followed by a non-digit)
    for (size_t k = 0; k + 1 < t.size(); ++k)
        if (t[k] == ',' && k > 0 && std::isdigit(static_cast<unsigned char>(t[k - 1])) &&
            std::isdigit(static_cast<unsigned char>(t[k + 1])) && t.find('(') == std::string::npos)
            t[k] = '.';
    bool blank = true;
    for (char c : t) blank = blank && std::isspace(static_cast<unsigned char>(c));
    if (blank) throw DesignError("Formula: enter an expression, e.g. A1 * sin(2 * pi * f1 * n)");
    Parser parser(t, vars);
    m_root = parser.parse();
    m_complex = parser.complex;
}

Expr::~Expr() = default;
Expr::Expr(Expr &&) noexcept = default;
Expr &Expr::operator=(Expr &&) noexcept = default;

double Expr::eval(const std::vector<double> &values) const { return eval_node(*m_root, values); }

cplx Expr::evalc(const std::vector<cplx> &values) const { return eval_c(*m_root, values); }

Vec eval_formula(const std::string &text, int n_end, double f_s, const std::map<std::string, double> &consts) {
    std::vector<std::string> names{"n", "t"};
    std::vector<double> vals{0, 0};
    for (const auto &kv : consts) {
        names.push_back(kv.first);
        vals.push_back(kv.second);
    }
    const Expr e(text, names);
    Vec x(size_t(std::max(n_end, 0)));
    for (int n = 0; n < n_end; ++n) {
        vals[0] = n;
        vals[1] = n / f_s;
        x[size_t(n)] = e.eval(vals);
        if (!std::isfinite(x[size_t(n)]))
            throw DesignError("Formula: the result is not finite for n = " + std::to_string(n));
    }
    return x;
}

CVec eval_formula_c(const std::string &text, int n_end, double f_s, const std::map<std::string, cplx> &consts,
                    bool *cmplx) {
    std::vector<std::string> names{"n", "t"};
    std::vector<cplx> vals{0.0, 0.0};
    bool c = false;
    for (const auto &kv : consts) {
        names.push_back(kv.first);
        vals.push_back(kv.second);
        c = c || kv.second.imag() != 0;
    }
    const Expr e(text, names);
    c = c || e.isComplex();
    if (cmplx) *cmplx = c;
    CVec x(size_t(std::max(n_end, 0)));
    if (!c) {  // real evaluation, identical to eval_formula
        std::vector<double> rv;
        for (const cplx &v : vals) rv.push_back(v.real());
        for (int n = 0; n < n_end; ++n) {
            rv[0] = n;
            rv[1] = n / f_s;
            x[size_t(n)] = e.eval(rv);
            if (!std::isfinite(x[size_t(n)].real()))
                throw DesignError("Formula: the result is not finite for n = " + std::to_string(n));
        }
        return x;
    }
    for (int n = 0; n < n_end; ++n) {
        vals[0] = double(n);
        vals[1] = n / f_s;
        x[size_t(n)] = e.evalc(vals);
        if (!std::isfinite(x[size_t(n)].real()) || !std::isfinite(x[size_t(n)].imag()))
            throw DesignError("Formula: the result is not finite for n = " + std::to_string(n));
    }
    return x;
}

}  // namespace pyfda
