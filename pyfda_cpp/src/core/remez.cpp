// Parks-McClellan (Remez exchange) FIR design, ported from scipy's
// _sigtoolsmodule.c (SciPy license, BSD), which is based on remez.c by
// Egil Kvaleberg, converted from the original FORTRAN program by
// J. H. McClellan, T. W. Parks and L. R. Rabiner.
#include "fir_design.hpp"

#include <cmath>
#include <string>

namespace pyfda {

namespace {

constexpr double TWOPI = 2 * PI;

double lagrange_interp(int k, int n, int m, const double *x) {
    double retval = 1.0;
    const double q = x[k];
    for (int l = 1; l <= m; ++l)
        for (int j = l; j <= n; j += m)
            if (j != k) retval *= 2.0 * (q - x[j]);
    return 1.0 / retval;
}

double freq_eval(int k, int n, const double *grid, const double *x, const double *y, const double *ad) {
    double d = 0.0, p = 0.0;
    const double xf = std::cos(TWOPI * grid[k]);
    for (int j = 1; j <= n; ++j) {
        const double c = ad[j] / (xf - x[j]);
        d += c;
        p += c * y[j];
    }
    return p / d;
}

// returns 0 on success, -1 on failure to converge
int remez_core(double *dev, double des[], double grid[], double edge[], double wt[], int ngrid,
               int nbands, int iext[], double alpha[], int nfcns, int itrmax, double *work,
               int dimsize, int *niter_out) {
    int k, k1, kkk, kn, knz, klow, kup, nz, nzz, nm1;
    int cn;
    int j, jchnge, jet, jm1, jp1;
    int l, luck = 0, nu, nut, nut1 = 0, niter;
    double ynz = 0.0, comp = 0.0, devl, gtemp, fsh, y1 = 0.0, err, dtemp, delf, dnum, dden;
    double aa = 0.0, bb = 0.0, ft, xe, xt;
    double *a, *p, *q, *ad, *x, *y;

    a = work;
    p = a + dimsize + 1;
    q = p + dimsize + 1;
    ad = q + dimsize + 1;
    x = ad + dimsize + 1;
    y = x + dimsize + 1;
    devl = -1.0;
    nz = nfcns + 1;
    nzz = nfcns + 2;
    niter = 0;

    do {
    L100:
        iext[nzz] = ngrid + 1;
        ++niter;
        if (niter > itrmax) break;

        for (j = 1; j <= nz; ++j) x[j] = std::cos(grid[iext[j]] * TWOPI);
        jet = (nfcns - 1) / 15 + 1;
        for (j = 1; j <= nz; ++j) ad[j] = lagrange_interp(j, nz, jet, x);

        dnum = 0.0;
        dden = 0.0;
        k = 1;
        for (j = 1; j <= nz; ++j) {
            l = iext[j];
            dnum += ad[j] * des[l];
            dden += double(k) * ad[j] / wt[l];
            k = -k;
        }
        *dev = dnum / dden;

        nu = 1;
        if ((*dev) > 0.0) nu = -1;
        (*dev) = -double(nu) * (*dev);
        k = nu;
        for (j = 1; j <= nz; ++j) {
            l = iext[j];
            y[j] = des[l] + double(k) * (*dev) / wt[l];
            k = -k;
        }
        if ((*dev) <= devl) {
            *niter_out = niter;
            return -1;
        }
        devl = (*dev);
        jchnge = 0;
        k1 = iext[1];
        knz = iext[nz];
        klow = 0;
        nut = -nu;
        j = 1;

    L200:
        if (j == nzz) ynz = comp;
        if (j >= nzz) goto L300;
        kup = iext[j + 1];
        l = iext[j] + 1;
        nut = -nut;
        if (j == 2) y1 = comp;
        comp = (*dev);
        if (l >= kup) goto L220;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L220;
        comp = double(nut) * err;
    L210:
        if (++l >= kup) goto L215;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L215;
        comp = double(nut) * err;
        goto L210;

    L215:
        iext[j++] = l - 1;
        klow = l - 1;
        ++jchnge;
        goto L200;

    L220:
        --l;
    L225:
        if (--l <= klow) goto L250;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) > 0.0) goto L230;
        if (jchnge <= 0) goto L225;
        goto L260;

    L230:
        comp = double(nut) * err;
    L235:
        if (--l <= klow) goto L240;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L240;
        comp = double(nut) * err;
        goto L235;
    L240:
        klow = iext[j];
        iext[j] = l + 1;
        ++j;
        ++jchnge;
        goto L200;

    L250:
        l = iext[j] + 1;
        if (jchnge > 0) goto L215;

    L255:
        if (++l >= kup) goto L260;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L255;
        comp = double(nut) * err;
        goto L210;
    L260:
        klow = iext[j++];
        goto L200;

    L300:
        if (j > nzz) goto L320;
        if (k1 > iext[1]) k1 = iext[1];
        if (knz < iext[nz]) knz = iext[nz];
        nut1 = nut;
        nut = -nu;
        l = 0;
        kup = k1;
        comp = ynz * (1.00001);
        luck = 1;
    L310:
        if (++l >= kup) goto L315;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L310;
        comp = double(nut) * err;
        j = nzz;
        goto L210;

    L315:
        luck = 6;
        goto L325;

    L320:
        if (luck > 9) goto L350;
        if (comp > y1) y1 = comp;
        k1 = iext[nzz];
    L325:
        l = ngrid + 1;
        klow = knz;
        nut = -nut1;
        comp = y1 * (1.00001);
    L330:
        if (--l <= klow) goto L340;
        err = (freq_eval(l, nz, grid, x, y, ad) - des[l]) * wt[l];
        if ((double(nut) * err - comp) <= 0.0) goto L330;
        j = nzz;
        comp = double(nut) * err;
        luck = luck + 10;
        goto L235;
    L340:
        if (luck == 6) goto L370;
        for (j = 1; j <= nfcns; ++j) iext[nzz - j] = iext[nz - j];
        iext[1] = k1;
        goto L100;
    L350:
        kn = iext[nzz];
        for (j = 1; j <= nfcns; ++j) iext[j] = iext[j + 1];
        iext[nz] = kn;
        goto L100;
    L370:;
    } while (jchnge > 0);

    // Calculation of the coefficients of the best approximation using the inverse DFT
    nm1 = nfcns - 1;
    fsh = 1.0e-06;
    gtemp = grid[1];
    x[nzz] = -2.0;
    cn = 2 * nfcns - 1;
    delf = 1.0 / cn;
    l = 1;
    kkk = 0;

    if (edge[1] == 0.0 && edge[2 * nbands] == 0.5) kkk = 1;
    if (nfcns <= 3) kkk = 1;
    if (kkk != 1) {
        dtemp = std::cos(TWOPI * grid[1]);
        dnum = std::cos(TWOPI * grid[ngrid]);
        aa = 2.0 / (dtemp - dnum);
        bb = -(dtemp + dnum) / (dtemp - dnum);
    }

    for (j = 1; j <= nfcns; ++j) {
        ft = (j - 1) * delf;
        xt = std::cos(TWOPI * ft);
        if (kkk != 1) {
            xt = (xt - bb) / aa;
            ft = std::acos(xt) / TWOPI;
        }
    L410:
        xe = x[l];
        if (xt > xe) goto L420;
        if ((xe - xt) < fsh) goto L415;
        ++l;
        goto L410;
    L415:
        a[j] = y[l];
        goto L425;
    L420:
        if ((xt - xe) < fsh) goto L415;
        grid[1] = ft;
        a[j] = freq_eval(1, nz, grid, x, y, ad);
    L425:
        if (l > 1) l = l - 1;
    }

    grid[1] = gtemp;
    dden = TWOPI / cn;
    for (j = 1; j <= nfcns; ++j) {
        dtemp = 0.0;
        dnum = (j - 1) * dden;
        if (nm1 >= 1)
            for (k = 1; k <= nm1; ++k) dtemp += a[k + 1] * std::cos(dnum * k);
        alpha[j] = 2.0 * dtemp + a[1];
    }

    for (j = 2; j <= nfcns; ++j) alpha[j] *= 2.0 / cn;
    alpha[1] /= cn;

    if (kkk != 1) {
        p[1] = 2.0 * alpha[nfcns] * bb + alpha[nm1];
        p[2] = 2.0 * aa * alpha[nfcns];
        q[1] = alpha[nfcns - 2] - alpha[nfcns];
        for (j = 2; j <= nm1; ++j) {
            if (j >= nm1) {
                aa *= 0.5;
                bb *= 0.5;
            }
            p[j + 1] = 0.0;
            for (k = 1; k <= j; ++k) {
                a[k] = p[k];
                p[k] = 2.0 * bb * a[k];
            }
            p[2] += a[1] * 2.0 * aa;
            jm1 = j - 1;
            for (k = 1; k <= jm1; ++k) p[k] += q[k] + aa * a[k + 1];
            jp1 = j + 1;
            for (k = 3; k <= jp1; ++k) p[k] += aa * a[k - 1];

            if (j != nm1) {
                for (k = 1; k <= j; ++k) q[k] = -a[k];
                q[1] += alpha[nfcns - 1 - j];
            }
        }
        for (j = 1; j <= nfcns; ++j) alpha[j] = p[j];
    }

    if (nfcns <= 3) alpha[nfcns + 1] = alpha[nfcns + 2] = 0.0;
    return 0;
}

double eff(double freq, const double *fx, int lband, int jtype) {
    if (jtype != 2) return fx[lband];
    return fx[lband] * freq;
}

double wate(double freq, const double *fx, const double *wtx, int lband, int jtype) {
    if (jtype != 2) return wtx[lband];
    if (fx[lband] >= 0.0001) return wtx[lband] / freq;
    return wtx[lband];
}

}  // namespace

Vec remez(int numtaps, const Vec &bands_in, const Vec &desired, const Vec &weight_in, RemezType type,
          double fs, int maxiter, int grid_density) {
    const int numbands = int(desired.size());
    Vec weight = weight_in.empty() ? Vec(numbands, 1.0) : weight_in;
    if (numtaps < 2) throw DesignError("The number of taps must be greater than 1.");
    if (int(bands_in.size()) != 2 * numbands || int(weight.size()) != numbands)
        throw DesignError("The inputs desired and weight must have same length, bands twice this length.");
    Vec bands(bands_in);
    double oldvalue = 0;
    for (double &b : bands) {
        if (b < oldvalue) throw DesignError("Bands must be monotonic starting at zero.");
        if (b * 2 > fs) throw DesignError("Band edges should be less than 1/2 the sampling frequency.");
        oldvalue = b;
        b = b / fs;
    }
    const int jtype = int(type);
    Vec h2(numtaps, 0.0);

    const int lgrid = grid_density;
    const int dimsize = int(std::ceil(numtaps / 2.0 + 2));
    const int wrksize = grid_density * dimsize;
    const int nfilt = numtaps;
    const int nbands = numbands;
    // 1-based views
    double *edge = bands.data() - 1;
    double *h = h2.data() - 1;
    const double *fx = desired.data() - 1;
    const double *wtx = weight.data() - 1;

    // a few elements of slack: the original C code relies on the arrays being contiguous
    Vec des_v(wrksize + 1), grid_v(wrksize + 1), wt_v(wrksize + 1), alpha_v(dimsize + 3),
        work_v((dimsize + 1) * 6 + 8);
    std::vector<int> iext_v(dimsize + 3);
    double *des = des_v.data(), *grid = grid_v.data(), *wt = wt_v.data(), *alpha = alpha_v.data(),
           *work = work_v.data();
    int *iext = iext_v.data();

    int neg = 1;
    if (jtype == 1) neg = 0;
    const int nodd = nfilt % 2;
    int nfcns = nfilt / 2;
    if (nodd == 1 && neg == 0) nfcns = nfcns + 1;

    grid[1] = edge[1];
    double delf = lgrid * nfcns;
    delf = 0.5 / delf;
    if (neg != 0 && edge[1] < delf) grid[1] = delf;
    int j = 1, l = 1, lband = 1;
    double temp, fup, change;

    for (;;) {
        fup = edge[l + 1];
        do {
            temp = grid[j];
            des[j] = eff(temp, fx, lband, jtype);
            wt[j] = wate(temp, fx, wtx, lband, jtype);
            if (++j > wrksize) throw DesignError("remez: too many grid points, reduce grid density.");
            grid[j] = temp + delf;
        } while (grid[j] <= fup);
        grid[j - 1] = fup;
        des[j - 1] = eff(fup, fx, lband, jtype);
        wt[j - 1] = wate(fup, fx, wtx, lband, jtype);
        ++lband;
        l += 2;
        if (lband > nbands) break;
        grid[j] = edge[l];
    }

    int ngrid = j - 1;
    if (neg == nodd && grid[ngrid] > (0.5 - delf)) --ngrid;
    if (ngrid < nfcns + 1)
        throw DesignError("Band edges are too close together to build the dense frequency grid. "
                          "Widen the bands, reduce the order or increase the grid density.");

    if (neg <= 0) {
        if (nodd != 1)
            for (j = 1; j <= ngrid; ++j) {
                change = std::cos(PI * grid[j]);
                des[j] = des[j] / change;
                wt[j] = wt[j] * change;
            }
    } else {
        if (nodd != 1) {
            for (j = 1; j <= ngrid; ++j) {
                change = std::sin(PI * grid[j]);
                des[j] = des[j] / change;
                wt[j] = wt[j] * change;
            }
        } else {
            for (j = 1; j <= ngrid; ++j) {
                change = std::sin(TWOPI * grid[j]);
                des[j] = des[j] / change;
                wt[j] = wt[j] * change;
            }
        }
    }

    temp = double(ngrid - 1) / double(nfcns);
    for (j = 1; j <= nfcns; ++j) iext[j] = int((j - 1) * temp) + 1;
    iext[nfcns + 1] = ngrid;
    const int nm1 = nfcns - 1;
    const int nz = nfcns + 1;

    double dev = 0.0;
    int niter = -1;
    if (remez_core(&dev, des, grid, edge, wt, ngrid, numbands, iext, alpha, nfcns, maxiter, work,
                   dimsize, &niter) < 0)
        throw DesignError("remez: failure to converge at iteration " + std::to_string(niter) +
                          ", try reducing the transition band width.");

    if (neg <= 0) {
        if (nodd != 0) {
            for (j = 1; j <= nm1; ++j) h[j] = 0.5 * alpha[nz - j];
            h[nfcns] = alpha[1];
        } else {
            h[1] = 0.25 * alpha[nfcns];
            for (j = 2; j <= nm1; ++j) h[j] = 0.25 * (alpha[nz - j] + alpha[nfcns + 2 - j]);
            h[nfcns] = 0.5 * alpha[1] + 0.25 * alpha[2];
        }
    } else {
        if (nodd != 0) {
            h[1] = 0.25 * alpha[nfcns];
            h[2] = 0.25 * alpha[nm1];
            for (j = 3; j <= nm1; ++j) h[j] = 0.25 * (alpha[nz - j] - alpha[nfcns + 3 - j]);
            h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[3];
            h[nz] = 0.0;
        } else {
            h[1] = 0.25 * alpha[nfcns];
            for (j = 2; j <= nm1; ++j) h[j] = 0.25 * (alpha[nz - j] - alpha[nfcns + 2 - j]);
            h[nfcns] = 0.5 * alpha[1] - 0.25 * alpha[2];
        }
    }
    for (j = 1; j <= nfcns; ++j) {
        const int k = nfilt + 1 - j;
        h[k] = (neg == 0) ? h[j] : -h[j];
    }
    if (neg == 1 && nodd == 1) h[nz] = 0.0;
    return h2;
}

}  // namespace pyfda
