// Fixpoint quantization and simulation, ported from pyfda
// (libs/pyfda_fix_lib.py: class Fixed, quant_coeffs; fixpoint_widgets:
// fir_df_pyfixp.py, iir_df1_pyfixp.py).
//
// Numbers are signed fixpoint numbers in the format WI.WF: WI integer bits,
// WF fractional bits and a sign bit, i.e. W = WI + WF + 1 bits in total, range
// -2^WI ... 2^WI - 2^-WF. Values are handled in pyfda's fractional format
// ('qfrac'): the quantized values are returned as doubles (exact up to 53 bits).
#pragma once

#include "types.hpp"

#include <cstdint>
#include <string>

namespace pyfda {

enum class Quant { Floor, Round, Fix, Ceil, None };  // Round = numpy.round (half to even)
enum class Ovfl { Wrap, Sat, None };

const char *quant_key(Quant q);
const char *ovfl_key(Ovfl o);
Quant quant_from_key(const std::string &k);
Ovfl ovfl_from_key(const std::string &k);

struct QFormat {
    int WI = 0, WF = 15;
    Quant quant = Quant::Round;
    Ovfl ovfl = Ovfl::Sat;
    int W() const { return WI + WF + 1; }
    double lsb() const;
    double min() const;  // -2^WI
    double max() const;  // 2^WI - LSB
};

/// Quantizer with overflow counter (pyfda Fixed.fixp with in_frmt = out_frmt = 'qfrac')
class Quantizer {
public:
    explicit Quantizer(const QFormat &q = QFormat()) : m_q(q) {}
    const QFormat &format() const { return m_q; }
    void setFormat(const QFormat &q) { m_q = q; }
    double fixp(double y);
    Vec fixp(const Vec &y);
    long long overflows() const { return m_n_over; }
    void reset() { m_n_over = 0; }

private:
    QFormat m_q;
    long long m_n_over = 0;
};

/// Integer representation of a quantized value (value * 2^WF)
long long to_int(double v_q, const QFormat &q);
/// Two's complement string of an integer with W bits in base 2, 8 (octal digits
/// of the W bit pattern), 16 or CSD (canonical signed digit, "+0-" digits)
std::string to_base(long long v, int W, int base);
std::string to_csd(long long v, int W);

/// Settings of a fixpoint filter
struct FxSpec {
    QFormat qi{0, 15, Quant::Round, Ovfl::Sat};    // input
    QFormat qcb{0, 15, Quant::Round, Ovfl::Sat};   // coefficients b
    QFormat qca{1, 14, Quant::Round, Ovfl::Sat};   // coefficients a (IIR)
    QFormat qacc{2, 30, Quant::Floor, Ovfl::Wrap}; // accumulator
    QFormat qo{0, 15, Quant::Round, Ovfl::Sat};    // output (and between IIR sections)
    bool acc_auto = true;    // accumulator word length from input and coefficient formats
    bool coeff_auto = true;  // integer bits of the coefficients from their maximum
};

/// Integer bits needed for the largest magnitude of the coefficients: smallest WI >= 0
/// with max |c| < 2^WI (pyfda uses ceil(log2(max |c|)) which saturates c = 1.0)
int coeff_wi(const Vec &c);
/// Update the automatic parts of `spec` for the design (pyfda update_accu_settings):
/// FIR: WF_acc = WF_i + WF_b, WI_acc = WI_i + WI_b + ceil(log2(sum |b|))
/// IIR (per section): WF_acc = max(WF_i + WF_b, WF_o + WF_a), WI_acc = max(WI_i + WI_b, WI_o + WI_a) + guard bits
void update_auto_formats(FxSpec &spec, const Ba &ba, const Sos &sos, bool fir);

/// Quantize coefficients; the leading a0 = 1 of recursive coefficients is kept
Vec quant_coeffs(const Vec &c, const QFormat &q, bool recursive = false, long long *n_over = nullptr);

struct FxResult {
    Vec x_q;   // quantized input
    Vec y;     // fixpoint output
    long long n_over_i = 0, n_over_acc = 0, n_over_o = 0, n_over_coeff = 0;
    Vec b_q, a_q;  // quantized FIR coefficients / flattened SOS coefficients
};

/// FIR direct form (pyfda fir_df_pyfixp): partial products and their sum are
/// quantized with the accumulator format, the sum is requantized to the output format.
/// Unlike pyfda (which multiplies b[0] with the oldest sample, which only matters for
/// asymmetric b) this is the regular convolution y[n] = sum b[k] x[n-k].
FxResult fx_filter_fir(const Vec &b, const FxSpec &spec, const Vec &x);
/// IIR as a cascade of second-order sections in direct form 1 (pyfda's iir_df1_pyfixp
/// works on the whole transfer function, which is unstable for higher orders):
/// acc = Q_acc(sum Q_acc(b_k x[n-k]) - sum Q_acc(a_k y[n-k])), y = Q_o(acc);
/// the output of every section is the input of the next one
FxResult fx_filter_sos(const Sos &sos, const FxSpec &spec, const Vec &x);

}  // namespace pyfda
