// Small expression evaluator for formula stimuli, a real valued subset of the
// numexpr syntax used by pyfda (plot_tran_stim.py, "Formula"):
//   numbers (also with ',' as decimal separator), variables, + - * / % **,
//   comparisons < <= > >= == != (1 / 0), & | ~, parentheses and the functions
//   sin cos tan arcsin arccos arctan arctan2 sinh cosh tanh arcsinh arccosh
//   arctanh exp expm1 log log10 log1p sqrt abs sign floor ceil round where
//   minimum maximum real imag conj complex
// Like numexpr, an expression is evaluated with complex numbers when it contains
// an imaginary literal (2j, j) or complex(), or when a constant is complex;
// comparisons, logic, %, floor, ceil, sign, arctan2, minimum and maximum then use
// the real parts. Real expressions are evaluated with real numbers (sqrt(-1) = nan).
#pragma once

#include "types.hpp"

#include <map>
#include <memory>
#include <string>

namespace pyfda {

struct ExprNode;

class Expr {
public:
    /// Parse `text`, throws DesignError with the position of syntax errors and
    /// unknown names; `vars` are the names that can be used
    Expr(const std::string &text, const std::vector<std::string> &vars);
    ~Expr();
    Expr(Expr &&) noexcept;
    Expr &operator=(Expr &&) noexcept;
    /// Evaluate with the variable values in the order of `vars`
    double eval(const std::vector<double> &values) const;
    /// Evaluate with complex numbers
    cplx evalc(const std::vector<cplx> &values) const;
    /// True when the expression contains imaginary literals or complex()
    bool isComplex() const { return m_complex; }

private:
    std::unique_ptr<ExprNode> m_root;
    bool m_complex = false;
};

/// Evaluate `text` for n = 0 ... n_end - 1 with the variables n, t = n / f_S and
/// the constants in `consts` (e.g. A1, f1, pi)
Vec eval_formula(const std::string &text, int n_end, double f_s, const std::map<std::string, double> &consts);
/// Same with complex constants; evaluated with complex numbers when the formula or a
/// constant is complex (then `*cmplx` is set to true), else like eval_formula
CVec eval_formula_c(const std::string &text, int n_end, double f_s, const std::map<std::string, cplx> &consts,
                    bool *cmplx = nullptr);

}  // namespace pyfda
