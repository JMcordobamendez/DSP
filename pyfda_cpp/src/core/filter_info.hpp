// Properties of a designed filter and a comparison of the achieved magnitude
// response with the specifications (pyfda's "Info" tab, input_widgets/input_info.py)
#pragma once

#include "filter_design.hpp"

#include <string>
#include <vector>

namespace pyfda {

struct BandCheck {
    std::string name;      // "PB", "SB", "PB2", "SB2"
    double f0, f1;         // band edges in the unit of f_S
    bool pass;             // pass band (ripple) or stop band (attenuation)
    double spec_db;        // A_PB (max. ripple) or A_SB (min. attenuation)
    double achieved_db;    // ripple max - min in the pass band, -max |H| / dB in the stop band
    bool ok;
};

struct FilterInfo {
    int order = 0;
    size_t n_b = 0, n_a = 0, n_sos = 0;
    bool fir = false;
    bool stable = true;          // all poles inside the unit circle
    double max_pole_radius = 0;
    bool min_phase = true;       // all zeros inside or on the unit circle
    bool linear_phase = false;   // FIR with (anti)symmetric coefficients
    double gain_dc = 0, gain_ny = 0;  // |H| at f = 0 and f_S / 2
    double h_max_db = 0;         // maximum of |H| in dB
    std::vector<BandCheck> bands;  // empty when the design has no band edge specs
};

FilterInfo filter_info(const FilterDesign &d);

}  // namespace pyfda
