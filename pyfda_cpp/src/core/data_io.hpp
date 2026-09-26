// Import of measured data (csv / txt / wav / npy) for the "Data Filt" tab and
// CSV export of the results. Port of pyfda/plot_widgets/plot_data_filt.py
#pragma once

#include "types.hpp"

#include <string>

namespace pyfda {

struct DataTable {
    size_t n_rows = 0;
    size_t n_cols = 0;
    Vec values;                      // row-major, NaN for non-numeric cells
    std::vector<std::string> names;  // column names (UTF-8)
    double fs = 0.0;                 // sampling rate stored in the file (wav), 0 if unknown
    std::vector<std::string> warnings;

    double at(size_t row, size_t col) const { return values[row * n_cols + col]; }
    Vec column(size_t col) const;
};

/// Convert a string to a number, NaN for empty or non-numeric strings.
/// A decimal comma is accepted when there is no dot.
double str2num(const std::string &s);
bool is_num(const std::string &s);

/// Tolerant reader for csv / txt files: detects encoding, delimiter, header
/// and skips comment (#, %, //) and metadata lines. Throws DesignError on failure.
DataTable read_text_table(const std::string &file_name);
/// Same, reading from the file contents
DataTable parse_text_table(const std::string &raw_bytes);

DataTable read_wav(const std::string &file_name);
DataTable read_npy(const std::string &file_name);

/// Load a file, the type is derived from the extension
DataTable load_data_file(const std::string &file_name);

/// True when column `col` looks like a time axis (monotonic increasing, no NaN)
bool is_time_column(const DataTable &t, size_t col);

/// Write time axis, original and filtered data as CSV
void write_csv(const std::string &file_name, const std::vector<std::string> &header,
               const std::vector<const Vec *> &columns);

}  // namespace pyfda
