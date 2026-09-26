#include "filtering.hpp"

#include <algorithm>
#include <cmath>
#include <string>

namespace pyfda {

Vec lfilter(const Vec &b_in, const Vec &a_in, const Vec &x, const Vec &zi, Vec *zf) {
    Vec a(a_in), b(b_in);
    while (a.size() > 1 && a[0] == 0.0) a.erase(a.begin());
    if (a.empty() || a[0] == 0.0) throw DesignError("Filter denominator must not be zero.");
    if (b.empty()) b = {0.0};
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, 0.0);
    b.resize(n, 0.0);
    const double a0 = a[0];
    for (double &v : a) v /= a0;
    for (double &v : b) v /= a0;
    Vec z(n - 1, 0.0);
    if (!zi.empty()) {
        if (zi.size() != n - 1) throw DesignError("Initial state has the wrong length.");
        z = zi;
    }
    Vec y(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        const double xi = x[i];
        const double yi = (n > 1 ? z[0] : 0.0) + b[0] * xi;
        for (size_t k = 1; k + 1 < n; ++k) z[k - 1] = z[k] + b[k] * xi - a[k] * yi;
        if (n > 1) z[n - 2] = b[n - 1] * xi - a[n - 1] * yi;
        y[i] = yi;
    }
    if (zf) *zf = z;
    return y;
}

Vec lfilter_zi(const Vec &b_in, const Vec &a_in) {
    Vec a(a_in), b(b_in);
    while (a.size() > 1 && a[0] == 0.0) a.erase(a.begin());
    if (a.empty() || a[0] == 0.0) throw DesignError("There must be at least one nonzero `a` coefficient.");
    if (a[0] != 1.0) {
        const double a0 = a[0];
        for (double &v : a) v /= a0;
        for (double &v : b) v /= a0;
    }
    const size_t n = std::max(a.size(), b.size());
    a.resize(n, 0.0);
    b.resize(n, 0.0);
    if (n == 1) return {};
    // (I - companion(a).T) zi = b[1:] - a[1:] b[0]; row i reads
    // a[i+1] z0 + z[i] - z[i+1] = B[i], solved as z[i] = u[i] + v[i] z0
    const size_t m = n - 1;
    Vec B(m), u(m + 1, 0.0), v(m + 1, 0.0);
    for (size_t i = 0; i < m; ++i) B[i] = b[i + 1] - a[i + 1] * b[0];
    for (size_t i = m; i-- > 0;) {
        u[i] = B[i] + u[i + 1];
        v[i] = -a[i + 1] + v[i + 1];
    }
    const double z0 = u[0] / (1.0 - v[0]);
    Vec zi(m);
    for (size_t i = 0; i < m; ++i) zi[i] = u[i] + v[i] * z0;
    return zi;
}

Vec sosfilt(const Sos &sos, const Vec &x, const std::vector<std::array<double, 2>> &zi) {
    Vec y(x);
    for (size_t s = 0; s < sos.size(); ++s) {
        const auto &c = sos[s];
        const double a0 = c[3];
        if (a0 == 0.0) throw DesignError("Second-order section with a0 = 0.");
        const double b0 = c[0] / a0, b1 = c[1] / a0, b2 = c[2] / a0, a1 = c[4] / a0, a2 = c[5] / a0;
        double z0 = zi.empty() ? 0.0 : zi[s][0];
        double z1 = zi.empty() ? 0.0 : zi[s][1];
        for (double &v : y) {
            const double xi = v;
            const double yi = b0 * xi + z0;
            z0 = b1 * xi - a1 * yi + z1;
            z1 = b2 * xi - a2 * yi;
            v = yi;
        }
    }
    return y;
}

std::vector<std::array<double, 2>> sosfilt_zi(const Sos &sos) {
    std::vector<std::array<double, 2>> zi(sos.size());
    double scale = 1.0;
    for (size_t s = 0; s < sos.size(); ++s) {
        const Vec b{sos[s][0], sos[s][1], sos[s][2]};
        const Vec a{sos[s][3], sos[s][4], sos[s][5]};
        const Vec z = lfilter_zi(b, a);
        zi[s] = {scale * z[0], scale * z[1]};
        scale *= (b[0] + b[1] + b[2]) / (a[0] + a[1] + a[2]);
    }
    return zi;
}

namespace {
Vec odd_ext(const Vec &x, size_t n) {
    if (n < 1) return x;
    if (n > x.size() - 1)
        throw DesignError("The data must contain more than " + std::to_string(n) +
                          " samples for zero-phase filtering with this filter.");
    Vec ext;
    ext.reserve(x.size() + 2 * n);
    for (size_t i = n; i >= 1; --i) ext.push_back(2 * x[0] - x[i]);
    ext.insert(ext.end(), x.begin(), x.end());
    const size_t L = x.size();
    for (size_t i = 0; i < n; ++i) ext.push_back(2 * x[L - 1] - x[L - 2 - i]);
    return ext;
}
}  // namespace

Vec filtfilt(const Vec &b, const Vec &a, const Vec &x) {
    const size_t edge = 3 * std::max(a.size(), b.size());
    if (x.size() <= edge)
        throw DesignError("The data must contain more than " + std::to_string(edge) +
                          " samples for zero-phase filtering with this filter.");
    const Vec ext = odd_ext(x, edge);
    const Vec zi = lfilter_zi(b, a);
    Vec z(zi);
    for (double &v : z) v *= ext.front();
    Vec y = lfilter(b, a, ext, z);
    std::reverse(y.begin(), y.end());
    z = zi;
    for (double &v : z) v *= y.front();
    y = lfilter(b, a, y, z);
    std::reverse(y.begin(), y.end());
    return Vec(y.begin() + long(edge), y.end() - long(edge));
}

Vec sosfiltfilt(const Sos &sos, const Vec &x) {
    size_t ntaps = 2 * sos.size() + 1;
    size_t nb0 = 0, na0 = 0;
    for (const auto &s : sos) {
        nb0 += s[2] == 0.0;
        na0 += s[5] == 0.0;
    }
    ntaps -= std::min(nb0, na0);
    const size_t edge = 3 * ntaps;
    if (x.size() <= edge)
        throw DesignError("The data must contain more than " + std::to_string(edge) +
                          " samples for zero-phase filtering with this filter.");
    const Vec ext = odd_ext(x, edge);
    const auto zi = sosfilt_zi(sos);
    auto scaled = [&](double f) {
        auto z = zi;
        for (auto &r : z) r = {r[0] * f, r[1] * f};
        return z;
    };
    Vec y = sosfilt(sos, ext, scaled(ext.front()));
    std::reverse(y.begin(), y.end());
    y = sosfilt(sos, y, scaled(y.front()));
    std::reverse(y.begin(), y.end());
    return Vec(y.begin() + long(edge), y.end() - long(edge));
}

// ---------------------------------------------------------------------------
namespace {
void fft_radix2(CVec &a, bool inverse) {
    const size_t n = a.size();
    for (size_t i = 1, j = 0; i < n; ++i) {
        size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) std::swap(a[i], a[j]);
    }
    for (size_t len = 2; len <= n; len <<= 1) {
        const double ang = 2 * PI / double(len) * (inverse ? 1 : -1);
        for (size_t i = 0; i < n; i += len) {
            for (size_t j = 0; j < len / 2; ++j) {
                const cplx w = std::polar(1.0, ang * double(j));
                const cplx u = a[i + j], v = a[i + j + len / 2] * w;
                a[i + j] = u + v;
                a[i + j + len / 2] = u - v;
            }
        }
    }
    if (inverse)
        for (cplx &v : a) v /= double(n);
}
}  // namespace

CVec fft(const CVec &x) {
    const size_t n = x.size();
    if (n == 0) return {};
    if ((n & (n - 1)) == 0) {
        CVec a(x);
        fft_radix2(a, false);
        return a;
    }
    // Bluestein's algorithm
    size_t m = 1;
    while (m < 2 * n - 1) m <<= 1;
    CVec w(n);
    for (size_t k = 0; k < n; ++k) {
        const double kk = double((k * k) % (2 * n));
        w[k] = std::polar(1.0, -PI * kk / double(n));
    }
    CVec a(m, 0.0), b(m, 0.0);
    for (size_t k = 0; k < n; ++k) a[k] = x[k] * w[k];
    b[0] = std::conj(w[0]);
    for (size_t k = 1; k < n; ++k) b[k] = b[m - k] = std::conj(w[k]);
    fft_radix2(a, false);
    fft_radix2(b, false);
    for (size_t i = 0; i < m; ++i) a[i] *= b[i];
    fft_radix2(a, true);
    CVec X(n);
    for (size_t k = 0; k < n; ++k) X[k] = a[k] * w[k];
    return X;
}

Vec amplitude_spectrum(const Vec &x) {
    const size_t N = x.size();
    if (N == 0) return {};
    CVec cx(x.begin(), x.end());
    const CVec X = fft(cx);
    Vec A(N / 2 + 1);
    for (size_t k = 0; k < A.size(); ++k) A[k] = std::abs(X[k]) / double(N) * (k > 0 ? 2.0 : 1.0);
    return A;
}

}  // namespace pyfda
