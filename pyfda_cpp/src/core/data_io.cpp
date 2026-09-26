#include "data_io.hpp"

#include <algorithm>
#include <cctype>
#include <clocale>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <limits>
#include <locale>
#include <map>
#include <sstream>

namespace pyfda {

namespace {
const double NaN = std::numeric_limits<double>::quiet_NaN();

std::string read_file(const std::string &file_name) {
    std::ifstream f(file_name, std::ios::binary);
    if (!f) throw DesignError("Couldn't open '" + file_name + "'.");
    std::ostringstream ss;
    ss << f.rdbuf();
    return ss.str();
}

std::string strip(const std::string &s) {
    size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

bool valid_utf8(const std::string &s) {
    size_t i = 0;
    while (i < s.size()) {
        const unsigned char c = static_cast<unsigned char>(s[i]);
        int n = 0;
        if (c < 0x80) n = 0;
        else if ((c >> 5) == 0x6) n = 1;
        else if ((c >> 4) == 0xE) n = 2;
        else if ((c >> 3) == 0x1E) n = 3;
        else return false;
        for (int k = 1; k <= n; ++k) {
            if (i + k >= s.size()) return false;
            if ((static_cast<unsigned char>(s[i + k]) >> 6) != 0x2) return false;
        }
        i += n + 1;
    }
    return true;
}

void append_utf8(std::string &out, uint32_t cp) {
    if (cp < 0x80) out += char(cp);
    else if (cp < 0x800) {
        out += char(0xC0 | (cp >> 6));
        out += char(0x80 | (cp & 0x3F));
    } else if (cp < 0x10000) {
        out += char(0xE0 | (cp >> 12));
        out += char(0x80 | ((cp >> 6) & 0x3F));
        out += char(0x80 | (cp & 0x3F));
    } else {
        out += char(0xF0 | (cp >> 18));
        out += char(0x80 | ((cp >> 12) & 0x3F));
        out += char(0x80 | ((cp >> 6) & 0x3F));
        out += char(0x80 | (cp & 0x3F));
    }
}

// Decode to UTF-8: UTF-16 (with BOM), UTF-8 (BOM stripped), else cp1252 / latin-1
std::string decode(const std::string &raw) {
    auto u = [&](size_t i) { return static_cast<unsigned char>(raw[i]); };
    if (raw.size() >= 2 && ((u(0) == 0xFF && u(1) == 0xFE) || (u(0) == 0xFE && u(1) == 0xFF))) {
        const bool le = u(0) == 0xFF;
        std::string out;
        for (size_t i = 2; i + 1 < raw.size(); i += 2) {
            uint32_t cp = le ? (u(i) | (u(i + 1) << 8)) : ((u(i) << 8) | u(i + 1));
            if (cp >= 0xD800 && cp < 0xDC00 && i + 3 < raw.size()) {
                const uint32_t lo = le ? (u(i + 2) | (u(i + 3) << 8)) : ((u(i + 2) << 8) | u(i + 3));
                cp = 0x10000 + ((cp - 0xD800) << 10) + (lo - 0xDC00);
                i += 2;
            }
            append_utf8(out, cp);
        }
        return out;
    }
    if (raw.size() >= 3 && u(0) == 0xEF && u(1) == 0xBB && u(2) == 0xBF) return raw.substr(3);
    if (valid_utf8(raw)) return raw;
    static const uint16_t cp1252[32] = {
        0x20AC, 0x81, 0x201A, 0x0192, 0x201E, 0x2026, 0x2020, 0x2021, 0x02C6, 0x2030, 0x0160,
        0x2039, 0x0152, 0x8D, 0x017D, 0x8F, 0x90, 0x2018, 0x2019, 0x201C, 0x201D, 0x2022,
        0x2013, 0x2014, 0x02DC, 0x2122, 0x0161, 0x203A, 0x0153, 0x9D, 0x017E, 0x0178};
    std::string out;
    for (unsigned char c : raw) {
        if (c < 0x80) out += char(c);
        else if (c < 0xA0) append_utf8(out, cp1252[c - 0x80]);
        else append_utf8(out, c);
    }
    return out;
}

std::vector<std::string> split_lines(const std::string &text) {
    std::vector<std::string> lines;
    std::string cur;
    for (size_t i = 0; i < text.size(); ++i) {
        const char c = text[i];
        if (c == '\n' || c == '\r') {
            lines.push_back(cur);
            cur.clear();
            if (c == '\r' && i + 1 < text.size() && text[i + 1] == '\n') ++i;
        } else {
            cur += c;
        }
    }
    if (!cur.empty()) lines.push_back(cur);
    return lines;
}

bool starts_with_comment(const std::string &line) {
    size_t i = 0;
    while (i < line.size() && std::isspace(static_cast<unsigned char>(line[i]))) ++i;
    const std::string s = line.substr(i);
    return s.rfind("#", 0) == 0 || s.rfind("%", 0) == 0 || s.rfind("//", 0) == 0;
}

std::vector<std::string> split_ws(const std::string &line) {
    std::vector<std::string> out;
    std::istringstream ss(line);
    std::string tok;
    while (ss >> tok) out.push_back(tok);
    return out;
}

// csv.reader with the given delimiter (double quotes as quote char), fields stripped
std::vector<std::string> split_csv(const std::string &line, char delim) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;
    for (size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];
        if (in_quotes) {
            if (c == '"') {
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    cur += '"';
                    ++i;
                } else {
                    in_quotes = false;
                }
            } else {
                cur += c;
            }
        } else if (c == '"' && strip(cur).empty()) {
            in_quotes = true;
            cur.clear();
        } else if (c == delim) {
            out.push_back(strip(cur));
            cur.clear();
        } else {
            cur += c;
        }
    }
    out.push_back(strip(cur));
    return out;
}

template <typename T>
T most_frequent(const std::vector<T> &v) {
    std::map<T, int> cnt;
    for (const T &x : v) ++cnt[x];
    T best = v.front();
    int best_n = -1;
    for (const auto &kv : cnt)
        if (kv.second > best_n) {
            best_n = kv.second;
            best = kv.first;
        }
    return best;
}
}  // namespace

Vec DataTable::column(size_t col) const {
    Vec c(n_rows);
    for (size_t r = 0; r < n_rows; ++r) c[r] = at(r, col);
    return c;
}

double str2num(const std::string &s_in) {
    std::string s = strip(s_in);
    if (s.empty()) return NaN;
    if (s.find(',') != std::string::npos && s.find('.') == std::string::npos)
        std::replace(s.begin(), s.end(), ',', '.');
    // strtod is locale dependent, parse in the "C" locale
    std::istringstream ss(s);
    ss.imbue(std::locale::classic());
    double v;
    ss >> v;
    if (ss.fail()) {
        // handle inf / nan spelled like in Python
        std::string l = s;
        std::transform(l.begin(), l.end(), l.begin(), [](unsigned char c) { return char(std::tolower(c)); });
        if (l == "inf" || l == "+inf" || l == "infinity") return std::numeric_limits<double>::infinity();
        if (l == "-inf" || l == "-infinity") return -std::numeric_limits<double>::infinity();
        return NaN;
    }
    ss >> std::ws;
    if (!ss.eof()) return NaN;  // trailing garbage
    return v;
}

bool is_num(const std::string &s) { return !std::isnan(str2num(s)); }

DataTable parse_text_table(const std::string &raw) {
    const std::string text = decode(raw);
    std::vector<std::string> lines;
    for (const std::string &ln : split_lines(text))
        if (!strip(ln).empty() && !starts_with_comment(ln)) lines.push_back(ln);
    if (lines.empty()) throw DesignError("The file doesn't contain any data.");

    const size_t n_tail = std::min<size_t>(lines.size(), 50);
    const std::vector<std::string> tail(lines.end() - long(n_tail), lines.end());

    // a single column with decimal commas ("0,5") and a header without comma
    bool single_col_dec_comma = !is_num(lines[0]) && lines[0].find(',') == std::string::npos;
    for (const std::string &ln : tail)
        if (std::count(ln.begin(), ln.end(), ',') != 1 || !is_num(ln)) single_col_dec_comma = false;

    char delim = 0;
    for (char d : {'\t', ';', ',', '|'}) {
        if (d == ',' && single_col_dec_comma) continue;
        std::vector<long> counts;
        for (const std::string &ln : tail) counts.push_back(std::count(ln.begin(), ln.end(), d));
        const long mode = most_frequent(counts);
        const long n_mode = std::count(counts.begin(), counts.end(), mode);
        if (mode > 0 && double(n_mode) >= 0.8 * double(counts.size())) {
            delim = d;
            break;
        }
    }

    std::vector<std::vector<std::string>> rows;
    for (const std::string &ln : lines) rows.push_back(delim ? split_csv(ln, delim) : split_ws(ln));

    // drop a trailing empty field caused by a delimiter at the end of each line
    {
        size_t n_trailing = 0;
        for (size_t i = rows.size() - n_tail; i < rows.size(); ++i)
            if (!rows[i].empty() && rows[i].back().empty()) ++n_trailing;
        if (double(n_trailing) >= 0.8 * double(n_tail))
            for (auto &r : rows)
                if (!r.empty() && r.back().empty()) r.pop_back();
    }

    auto has_num = [](const std::vector<std::string> &r) {
        return std::any_of(r.begin(), r.end(), [](const std::string &c) { return is_num(c); });
    };
    std::vector<size_t> lens;
    for (size_t i = rows.size() - n_tail; i < rows.size(); ++i)
        if (has_num(rows[i])) lens.push_back(rows[i].size());
    if (lens.empty()) throw DesignError("The file doesn't contain any numeric data.");
    const size_t n_cols = most_frequent(lens);

    size_t i_data = 0;
    while (i_data < rows.size() && !(rows[i_data].size() == n_cols && has_num(rows[i_data]))) ++i_data;

    DataTable t;
    t.n_cols = n_cols;
    for (size_t c = 0; c < n_cols; ++c) t.names.push_back("Col " + std::to_string(c + 1));
    for (size_t i = 0; i < i_data; ++i) {
        const auto &r = rows[i];
        if (r.size() == n_cols && !has_num(r)) {
            for (size_t c = 0; c < n_cols; ++c)
                if (!r[c].empty()) t.names[c] = r[c];
            break;  // use the first text line (names), not a following units line
        }
    }
    size_t n_skipped = 0;
    for (size_t i = i_data; i < rows.size(); ++i) {
        if (rows[i].size() != n_cols) {
            ++n_skipped;
            continue;
        }
        for (const std::string &cell : rows[i]) t.values.push_back(str2num(cell));
        ++t.n_rows;
    }
    if (n_skipped)
        t.warnings.push_back("Skipped " + std::to_string(n_skipped) +
                             " line(s) with a number of fields different from " + std::to_string(n_cols) + ".");
    if (std::all_of(t.values.begin(), t.values.end(), [](double v) { return std::isnan(v); }))
        throw DesignError("The data doesn't contain any numeric values.");
    return t;
}

DataTable read_text_table(const std::string &file_name) { return parse_text_table(read_file(file_name)); }

namespace {
uint32_t le32(const std::string &s, size_t i) {
    return uint32_t(uint8_t(s[i])) | (uint32_t(uint8_t(s[i + 1])) << 8) | (uint32_t(uint8_t(s[i + 2])) << 16) |
           (uint32_t(uint8_t(s[i + 3])) << 24);
}
uint16_t le16(const std::string &s, size_t i) { return uint16_t(uint8_t(s[i]) | (uint8_t(s[i + 1]) << 8)); }
}  // namespace

DataTable read_wav(const std::string &file_name) {
    const std::string s = read_file(file_name);
    if (s.size() < 12 || s.compare(0, 4, "RIFF") != 0 || s.compare(8, 4, "WAVE") != 0)
        throw DesignError("Not a RIFF / WAVE file.");
    size_t pos = 12;
    uint16_t fmt = 0, n_ch = 0, bits = 0;
    uint32_t fs = 0;
    size_t data_pos = 0, data_len = 0;
    while (pos + 8 <= s.size()) {
        const std::string id = s.substr(pos, 4);
        size_t len = le32(s, pos + 4);
        const size_t body = pos + 8;
        if (id == "fmt ") {
            if (body + 16 > s.size()) throw DesignError("Malformed wav header.");
            fmt = le16(s, body);
            n_ch = le16(s, body + 2);
            fs = le32(s, body + 4);
            bits = le16(s, body + 14);
            if (fmt == 0xFFFE && len >= 26) fmt = le16(s, body + 24);  // WAVE_FORMAT_EXTENSIBLE
        } else if (id == "data") {
            data_pos = body;
            data_len = std::min(len, s.size() - body);
            break;
        }
        pos = body + len + (len & 1);
    }
    if (!data_pos || !n_ch) throw DesignError("No audio data found in wav file.");
    if (!((fmt == 1 && (bits == 8 || bits == 16 || bits == 24 || bits == 32)) ||
          (fmt == 3 && (bits == 32 || bits == 64))))
        throw DesignError("Unsupported wav format (only PCM 8/16/24/32 bit and float 32/64 bit).");
    const size_t bps = bits / 8;
    const size_t n_frames = data_len / (bps * n_ch);
    DataTable t;
    t.n_rows = n_frames;
    t.n_cols = n_ch;
    t.fs = fs;
    t.values.resize(n_frames * n_ch);
    for (size_t i = 0; i < n_frames * n_ch; ++i) {
        const size_t p = data_pos + i * bps;
        double v = 0.0;
        if (fmt == 3) {
            if (bits == 32) {
                float f;
                std::memcpy(&f, s.data() + p, 4);
                v = f;
            } else {
                std::memcpy(&v, s.data() + p, 8);
            }
        } else if (bits == 8) {
            v = (double(uint8_t(s[p])) - 128.0) / 128.0;
        } else if (bits == 16) {
            v = double(int16_t(le16(s, p))) / 32768.0;
        } else if (bits == 24) {
            int32_t x = int32_t(uint8_t(s[p]) | (uint8_t(s[p + 1]) << 8) | (uint8_t(s[p + 2]) << 16));
            if (x & 0x800000) x -= 0x1000000;
            v = double(x) / 8388608.0;
        } else {
            v = double(int32_t(le32(s, p))) / 2147483648.0;
        }
        t.values[i] = v;
    }
    for (size_t c = 0; c < n_ch; ++c) t.names.push_back("Ch " + std::to_string(c + 1));
    return t;
}

DataTable read_npy(const std::string &file_name) {
    const std::string s = read_file(file_name);
    if (s.size() < 10 || s.compare(0, 6, "\x93NUMPY") != 0) throw DesignError("Not a .npy file.");
    const int major = uint8_t(s[6]);
    size_t hlen, hpos;
    if (major == 1) {
        hlen = le16(s, 8);
        hpos = 10;
    } else {
        hlen = le32(s, 8);
        hpos = 12;
    }
    const std::string header = s.substr(hpos, hlen);
    auto field = [&](const std::string &key) {
        const size_t k = header.find("'" + key + "'");
        if (k == std::string::npos) throw DesignError("Malformed .npy header.");
        return header.substr(header.find(':', k) + 1);
    };
    std::string descr = field("descr");
    descr = descr.substr(descr.find('\'') + 1);
    descr = descr.substr(0, descr.find('\''));
    const std::string fortran = strip(field("fortran_order")).substr(0, 4);
    std::string shape = field("shape");
    shape = shape.substr(shape.find('(') + 1, shape.find(')') - shape.find('(') - 1);
    std::vector<size_t> dims;
    std::istringstream ss(shape);
    std::string tok;
    while (std::getline(ss, tok, ','))
        if (!strip(tok).empty()) dims.push_back(std::stoul(strip(tok)));
    if (dims.empty() || dims.size() > 2) throw DesignError("Only 1D and 2D arrays are supported.");
    const char order = descr[0];
    const char kind = descr[1];
    const size_t size = std::stoul(descr.substr(2));
    if (order == '>') throw DesignError("Big endian .npy files are not supported.");
    const size_t n_rows = dims[0], n_cols = dims.size() == 2 ? dims[1] : 1;
    const size_t data_pos = hpos + hlen;
    if (data_pos + n_rows * n_cols * size > s.size()) throw DesignError("Truncated .npy file.");
    DataTable t;
    t.n_rows = n_rows;
    t.n_cols = n_cols;
    t.values.resize(n_rows * n_cols);
    for (size_t i = 0; i < n_rows * n_cols; ++i) {
        const char *p = s.data() + data_pos + i * size;
        double v;
        if (kind == 'f' && size == 8) std::memcpy(&v, p, 8);
        else if (kind == 'f' && size == 4) { float f; std::memcpy(&f, p, 4); v = f; }
        else if (kind == 'i' && size == 2) { int16_t x; std::memcpy(&x, p, 2); v = x; }
        else if (kind == 'i' && size == 4) { int32_t x; std::memcpy(&x, p, 4); v = x; }
        else if (kind == 'i' && size == 8) { int64_t x; std::memcpy(&x, p, 8); v = double(x); }
        else if (kind == 'u' && size == 1) { v = uint8_t(*p); }
        else throw DesignError("Unsupported .npy data type '" + descr + "'.");
        // store row-major
        const size_t r = fortran == "True" ? i % n_rows : i / n_cols;
        const size_t c = fortran == "True" ? i / n_rows : i % n_cols;
        t.values[r * n_cols + c] = v;
    }
    for (size_t c = 0; c < n_cols; ++c) t.names.push_back("Col " + std::to_string(c + 1));
    return t;
}

DataTable load_data_file(const std::string &file_name) {
    std::string ext = file_name.substr(file_name.find_last_of('.') + 1);
    std::transform(ext.begin(), ext.end(), ext.begin(), [](unsigned char c) { return char(std::tolower(c)); });
    if (ext == "wav") return read_wav(file_name);
    if (ext == "npy") return read_npy(file_name);
    return read_text_table(file_name);
}

bool is_time_column(const DataTable &t, size_t col) {
    if (t.n_cols < 2 || t.n_rows < 2) return false;
    for (size_t r = 0; r < t.n_rows; ++r) {
        if (std::isnan(t.at(r, col))) return false;
        if (r > 0 && t.at(r, col) < t.at(r - 1, col)) return false;
    }
    return t.at(t.n_rows - 1, col) > t.at(0, col);
}

void write_csv(const std::string &file_name, const std::vector<std::string> &header,
               const std::vector<const Vec *> &columns) {
    std::ofstream f(file_name, std::ios::binary);
    if (!f) throw DesignError("Couldn't write '" + file_name + "'.");
    f.imbue(std::locale::classic());
    f.precision(17);
    for (size_t i = 0; i < header.size(); ++i) f << (i ? "," : "") << header[i];
    f << "\r\n";
    const size_t n = columns.empty() ? 0 : columns[0]->size();
    for (size_t r = 0; r < n; ++r) {
        for (size_t c = 0; c < columns.size(); ++c) f << (c ? "," : "") << (*columns[c])[r];
        f << "\r\n";
    }
    if (!f) throw DesignError("Couldn't write '" + file_name + "'.");
}

}  // namespace pyfda
