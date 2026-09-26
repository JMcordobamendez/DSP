// Export of fixpoint filters for hardware: Xilinx COE coefficient files and
// synthesizable VHDL (numeric_std) that is bit exact with fx_filter_fir() /
// fx_filter_sos() (checked with GHDL by tests/verify_vhdl.py).
//
// pyfda generates Verilog with amaranth (formerly nMigen); here small VHDL and
// Verilog generators are used instead so that no Python is needed. The Verilog
// modules are checked with Icarus Verilog (tests/verify_vhdl.py).
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

/// Verilog-2001 module with the same ports and timing as the VHDL entities
std::string export_verilog_fir(const Vec &b, const FxSpec &spec, const std::string &name = "fir_filter");
std::string export_verilog_sos(const Sos &sos, const FxSpec &spec, const std::string &name = "iir_filter");

/// Self-checking testbench for the entity / module `name` (VHDL-2008 or Verilog):
/// applies the stimulus x (quantized to the input format), compares y with the
/// fixpoint model sample by sample and reports "PASSED" or the number of errors.
/// For FIR filters pass the coefficients b and an empty sos.
std::string export_vhdl_testbench(const Vec &b, const Sos &sos, const FxSpec &spec, const Vec &x,
                                  const std::string &name);
std::string export_verilog_testbench(const Vec &b, const Sos &sos, const FxSpec &spec, const Vec &x,
                                     const std::string &name);
/// Default test signal for the testbenches: impulse, positive and negative step and
/// pseudo random noise, scaled to the input format
Vec hdl_test_stimulus(const FxSpec &spec, size_t n_taps);

/// Valid VHDL identifier from a file name ("my filter-1" -> "my_filter_1")
std::string vhdl_identifier(const std::string &name);

}  // namespace pyfda
