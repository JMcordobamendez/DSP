// Export of fixpoint filters for hardware: Xilinx COE coefficient files and
// synthesizable VHDL (numeric_std) that is bit exact with fx_filter_fir() /
// fx_filter_sos() (checked with GHDL by tests/verify_vhdl.py).
//
// pyfda generates Verilog with amaranth (formerly nMigen); here a small VHDL
// generator is used instead so that no Python is needed.
#pragma once

#include "fixpoint.hpp"

#include <string>

namespace pyfda {

/// Xilinx COE file (FIR Compiler) with the quantized coefficients, radix 10 or 16
std::string export_coe(const Vec &b, const QFormat &qcb, int radix = 16);

/// VHDL entity `name` with ports clk, rst (synchronous, active high), en (sample
/// enable), x (signed, input format) and y (signed, output format, registered:
/// y is valid one clock after x). Throws DesignError for settings that can't be
/// implemented (quantization 'none').
std::string export_vhdl_fir(const Vec &b, const FxSpec &spec, const std::string &name = "fir_filter");
std::string export_vhdl_sos(const Sos &sos, const FxSpec &spec, const std::string &name = "iir_filter");

/// Valid VHDL identifier from a file name ("my filter-1" -> "my_filter_1")
std::string vhdl_identifier(const std::string &name);

}  // namespace pyfda
