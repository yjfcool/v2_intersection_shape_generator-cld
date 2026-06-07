// ─────────────────────────────────────────────────────────────────────────────
//  L-BFGS with strong-Wolfe line search
// ─────────────────────────────────────────────────────────────────────────────
#include "lbfgs_solver.h"
#include <cmath>
#include <deque>
#include <cassert>

// ── Cubic interpolation for interval [a, b] ──────────────────────────────────
// Finds minimiser of cubic through (a, fa, ga) and (b, fb, gb) in [a,b].
static double cubicMinimiser(double a, double fa, double ga, double b, double fb, double gb) {
    double d1 = ga + gb - 3.0 * (fb - fa) / (b - a);
    double d2sq = d1 * d1 - ga * gb;
    if (d2sq < 0) return 0.5 * (a + b); // fallback: bisect
    double d2 = std::sqrt(d2sq);
    double t = b - (b - a) * (gb + d2 - d1) / (gb - ga + 2.0 * d2);
    return std::max(a + 0.1 * (b - a), std::min(b - 0.1 * (b - a), t)); // clamp in (a,b)
}

// ── Strong-Wolfe line search ──────────────────────────────────────────────────
double LBFGSSolver::lineSearch(
    CostFn& fn, const VecXd& x, const VecXd& dir,
    double f0, const VecXd& g0, VecXd& xn, double& fn_val, VecXd& gn, int& evals) {
    const double c1 = cfg_.wolfe_c1;
    const double c2 = cfg_.wolfe_c2;
    double derphi0 = g0.dot(dir);

    // Non-descent direction guard — revert to steepest descent step
    if (derphi0 >= 0) {
        xn = x - 1e-4 * g0;
        fn_val = fn(xn, gn);
        ++evals;
        return 1e-4;
    }

    double alpha_prev = 0, f_prev = f0;
    double alpha = 1.0; // try unit step first (Newton-like)
    const double alpha_max = 50.0;
    VecXd g_cur(g0.size());
    double f_cur = f0; // track last evaluated f

    auto phi = [&](double a, VecXd& gout) -> double {
        xn = x + a * dir;
        double v = fn(xn, gout);
        ++evals;
        f_cur = v; // keep f_cur in sync
        return v;
    };

    for (int iter = 0; iter < cfg_.max_ls_iter; ++iter) {
        double f_cur = phi(alpha, g_cur);

        bool armijo_fail = (f_cur > f0 + c1 * alpha * derphi0);
        bool non_descent = (iter > 0 && f_cur >= f_prev);

        if (armijo_fail || non_descent) {
            // Zoom into [alpha_prev, alpha]
            double a_lo = alpha_prev, f_lo = f_prev;
            double a_hi = alpha;
            VecXd g_lo(g0.size()), g_hi(g_cur);
            g_lo = (iter == 0) ? g0 : g_cur; // approximation; good enough

            for (int z = 0; z < 10; ++z) {
                double dlo = g_lo.dot(dir);
                double dhi = g_hi.dot(dir);
                double a_j = cubicMinimiser(a_lo, f_lo, dlo, a_hi, fn_val, dhi);
                double f_j = phi(a_j, g_cur);

                if (f_j > f0 + c1 * a_j * derphi0 || f_j >= f_lo) {
                    a_hi = a_j;
                    g_hi = g_cur;
                    fn_val = f_j;
                } else {
                    double dp = g_cur.dot(dir);
                    if (std::abs(dp) <= -c2 * derphi0) {
                        gn = g_cur;
                        fn_val = f_j;
                        return a_j;
                    }
                    f_lo = f_j;
                    g_lo = g_cur;
                    if (dp * (a_hi - a_lo) >= 0)
                        a_hi = a_lo;
                    a_lo = a_j;
                }
            }
            // Return best found so far
            gn = g_cur;
            fn_val = f_lo;
            return a_lo;
        }

        double dp = g_cur.dot(dir);
        if (std::abs(dp) <= -c2 * derphi0) {
            // Strong Wolfe satisfied
            gn = g_cur;
            fn_val = f_cur;
            return alpha;
        }
        if (dp >= 0) {
            // Zoom into [alpha, alpha_prev]
            double a_lo = alpha, f_lo = f_cur;
            VecXd g_lo = g_cur;
            for (int z = 0; z < 8; ++z) {
                double dlo = g_lo.dot(dir);
                double a_j = 0.5 * (a_lo + alpha_prev);
                double f_j = phi(a_j, g_cur);
                if (f_j > f0 + c1 * a_j * derphi0 || f_j >= f_lo) {
                    alpha_prev = a_j;
                } else {
                    double dpj = g_cur.dot(dir);
                    if (std::abs(dpj) <= -c2 * derphi0) {
                        gn = g_cur;
                        fn_val = f_j;
                        return a_j;
                    }
                    f_lo = f_j;
                    g_lo = g_cur;
                    a_lo = a_j;
                }
            }
            gn = g_lo;
            fn_val = f_lo;
            return a_lo;
        }

        alpha_prev = alpha;
        f_prev = f_cur;
        // Expand: try larger step (cubic extrapolation)
        double new_alpha = std::min(alpha * 2.5, alpha_max);
        alpha = new_alpha;
        if (alpha < 1e-15)
            break;
    }
    // Fallback: return whatever we have
    gn = g_cur;
    fn_val = f_cur;
    return alpha;
}

// ── L-BFGS two-loop recursion ──────────────────────────────────────────────
SolveResult LBFGSSolver::solve(CostFn fn, const VecXd& x0) {
    const int n = (int)x0.size();
    VecXd x = x0, g(n);
    double f = fn(x, g);
    int evals = 1;

    std::deque<VecXd> sv, yv;
    std::deque<double> rho;

    SolveResult res;
    res.converged = false;

    for (int iter = 0; iter < cfg_.max_iter; ++iter) {
        res.iterations = iter;

        if (g.norm() < cfg_.grad_tol) {
            res.converged = true;
            break;
        }

        // ── Two-loop L-BFGS direction ─────────────────────────────────────
        const int m = (int)sv.size();
        VecXd q = g;
        std::vector<double> alpha_arr(m);
        for (int i = m - 1; i >= 0; --i) {
            alpha_arr[i] = rho[i] * sv[i].dot(q);
            q -= alpha_arr[i] * yv[i];
        }

        // Initial Hessian scaling: γ = sᵀy / yᵀy
        VecXd r = q;
        if (m > 0) {
            double sy = sv.back().dot(yv.back());
            double yy = yv.back().dot(yv.back());
            double gamma = (yy > 1e-20) ? sy / yy : 1.0;
            r *= std::max(1e-10, std::min(gamma, 1e3)); // clamp Hessian scaling
        }
        for (int i = 0; i < m; ++i) {
            double beta = rho[i] * yv[i].dot(r);
            r += sv[i] * (alpha_arr[i] - beta);
        }

        VecXd dir = -r;

        // ── Line search ───────────────────────────────────────────────────
        VecXd xn(n), gn(n);
        double fn_val = f;
        double alpha = lineSearch(fn, x, dir, f, g, xn, fn_val, gn, evals);
        (void)alpha;

        // ── Update history ────────────────────────────────────────────────
        VecXd s = xn - x;
        VecXd y = gn - g;
        double sy = s.dot(y);

        if (sy > 1e-20) {
            // skip update if curvature condition not met
            if ((int)sv.size() == cfg_.history_size) {
                sv.pop_front();
                yv.pop_front();
                rho.pop_front();
            }
            sv.push_back(s);
            yv.push_back(y);
            rho.push_back(1.0 / sy);
        }

        double f_diff = std::abs(fn_val - f);
        x = xn, g = gn, f = fn_val;
        if (f_diff < cfg_.func_tol) {
            res.converged = true;
            break;
        }
    }

    res.x = x;
    res.final_cost = f;
    res.fn_evals = evals;
    return res;
}

SolveResult LBFGSSolver::solveWarm(CostFn fn, const VecXd& x_warm) {
    return solve(fn, x_warm); // history naturally resets; warm start just from x
}
