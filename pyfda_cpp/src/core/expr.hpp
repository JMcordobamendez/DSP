// Small expression evaluator for formula stimuli, a real valued subset of the
// numexpr syntax used by pyfda (plot_tran_stim.py, "Formula"):
//   numbers (also with ',' as decimal separator), variables, + - * / % **,
//   comparisons < <= > >= == != (1 / 0), & | ~, parentheses and the functions
//   sin cos tan arcsin arccos arctan arctan2 sinh cosh tanh arcsinh arccosh
//   arctanh exp expm1 log log10 log1p sqrt abs sign floor ceil round where
//   minimum maximum
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

private:
    std::unique_ptr<ExprNode> m_root;
};

/// Evaluate `text` for n = 0 ... n_end - 1 with the variables n, t = n / f_S and
/// the constants in `consts` (e.g. A1, f1, pi)
Vec eval_formula(const std::string &text, int n_end, double f_s, const std::map<std::string, double> &consts);

}  // namespace pyfda
