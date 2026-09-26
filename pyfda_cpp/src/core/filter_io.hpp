// Saving and loading filter designs (JSON) and exporting the coefficients in
// formats for other tools (MATLAB / Octave, C, Python / NumPy, CSV).
//
// pyfda saves filters as pickled Python dicts (*.pkl / *.npz), which can't be read
// without Python; pyfda_cpp uses a JSON file with the specifications instead. The
// coefficients are stored as well for reference and for other tools.
#pragma once

#include "filter_design.hpp"

#include <string>

namespace pyfda {

/// Short names used in files and in the CLI ("LP", "ellip", "kaiser", "ichige")
const char *resp_type_key(RespType rt);
const char *method_key(DesignMethod m);
const char *order_alg_key(RemezAlg a);
RespType resp_type_from_key(const std::string &k);
DesignMethod method_from_key(const std::string &k);
RemezAlg order_alg_from_key(const std::string &k);
/// Window from its display name, case insensitive ("Kaiser", "blackmanharris" ...)
WindowType window_from_name(const std::string &n);
const char *window_name(WindowType w);

/// Contents of a filter file
struct FilterFile {
    FilterSpec spec;
    std::string unit = "f_S";  // frequency unit of the specs: "f_S" (normalized), "Hz", "kHz", "MHz"
    Ba ba;                     // stored coefficients (only for reference)
    Sos sos;
};

/// JSON text of a design (the specs used and the resulting coefficients)
std::string filter_to_json(const FilterDesign &d, const std::string &unit = "f_S");
/// Parse a filter file, throws DesignError with a readable message on errors.
/// Missing keys keep their default values.
FilterFile filter_from_json(const std::string &text);

void save_filter(const std::string &file_name, const FilterDesign &d, const std::string &unit = "f_S");
FilterFile load_filter(const std::string &file_name);

enum class CoeffFormat { Csv, Matlab, CHeader, Python };
/// File name filter for dialogs, e.g. "MATLAB / Octave (*.m)"
const char *coeff_format_filter(CoeffFormat f);
const char *coeff_format_suffix(CoeffFormat f);
/// Coefficients of the design as text in the given format. Values are written
/// with 17 significant digits, i.e. they read back exactly.
std::string export_coeffs(const FilterDesign &d, CoeffFormat f, const std::string &name = "filter");
void export_coeffs(const std::string &file_name, const FilterDesign &d, CoeffFormat f);

}  // namespace pyfda
