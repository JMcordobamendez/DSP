// High level filter design from specifications, modelled on pyfda's filter
// widgets (pyfda/filter_widgets/*.py): response type, design method, minimum
// or manual order, frequency and amplitude specs.
#pragma once

#include "fir_design.hpp"
#include "iir_design.hpp"
#include "types.hpp"

#include <string>

namespace pyfda {

enum class RespType { LP, HP, BP, BS };
enum class DesignMethod { Butter, Cheby1, Cheby2, Ellip, Bessel, Firwin, Equiripple };
enum class OrderMode { Min, Manual };

const char *resp_type_name(RespType rt);
const char *method_name(DesignMethod m);
bool is_fir(DesignMethod m);

struct FilterSpec {
    RespType rt = RespType::LP;
    DesignMethod method = DesignMethod::Ellip;
    OrderMode fo = OrderMode::Min;
    int N = 10;          // filter order
    double f_s = 1.0;    // sampling frequency, all frequencies below use the same unit
    // pass / stop band edges and corner frequencies (second value for BP / BS)
    double f_pb = 0.1, f_pb2 = 0.35, f_sb = 0.15, f_sb2 = 0.3;
    double f_c = 0.125, f_c2 = 0.325;
    // maximum pass band ripple and minimum stop band attenuation in dB
    double A_PB = 1.0, A_SB = 60.0;
    // weights for manual order equiripple designs
    double W_PB = 1.0, W_SB = 1.0;
    // window for the window method
    WindowType window = WindowType::Kaiser;
    double win_par = 10.0;
    RemezAlg order_alg = RemezAlg::Ichige;
    int grid_density = 16;
};

struct FilterDesign {
    FilterSpec spec;   // specs used, with order / corner frequencies / weights calculated
    bool fir = false;
    Ba ba;
    Sos sos;           // empty for FIR filters
    Zpk zpk;
    std::string info;  // short description
};

/// Design a filter; throws DesignError with a user readable message on failure
FilterDesign design_filter(const FilterSpec &spec);

/// Linear ripples from dB for FIR filters (pyfda.libs.special_functions.unit2lin)
double fir_a_pb_lin(double A_PB_dB);
double sb_lin(double A_SB_dB);

}  // namespace pyfda
