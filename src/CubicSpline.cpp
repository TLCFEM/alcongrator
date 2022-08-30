#include "CubicSpline.h"

using namespace arma;

CubicSpline::CubicSpline(const vec& T, const vec& M, const bool N)
    : natural(N), time(T), magnitude(M) {
    dt = diff(time);
    dy = diff(magnitude);

    const auto np = time.n_elem;
    const auto n = np - 1llu;
    const auto nm = n - 1llu;

    sp_mat system(np, np);

    system.diag().fill(2.);

    vec b(np, fill::none);

    if(natural) {
        b(0) = b(n) = 0.;

        system.at(0, 1) = system.at(n, nm) = 0.;
    }
    else {
        b(0) = dy(0) / dt(0) / dt(0);
        b(n) = -dy(nm) / dt(nm) / dt(nm);

        system.at(0llu, 1llu) = 1.;
        system.at(n, nm) = 1.;
    }

    for(auto I = 1llu; I < n; ++I) {
        const auto J = I - 1llu, K = I + 1llu;
        const auto denom = time(K) - time(J);
        const auto mu = dt(J) / denom;
        system.at(I, J) = mu;
        system.at(I, I) = 2.;
        system.at(I, K) = 1. - mu;
        b(I) = (dy(I) / dt(I) - dy(J) / dt(J)) / denom;
    }

    m = spsolve(system, b);
}

double CubicSpline::evaluate(const double step_time) const {
    if(step_time <= time.front()) return magnitude.front();

    if(step_time >= time.back()) return magnitude.back();

    const auto I = std::distance(time.cbegin(), std::lower_bound(time.cbegin(), time.cend(), step_time));
    const auto J = I - 1;

    double y = (m(J) * pow(time(I) - step_time, 3.) + m(I) * pow(step_time - time(J), 3.)) / dt(J);
    y += (magnitude(J) / dt(J) - m(J) * dt(J)) * (time(I) - step_time);
    y += (magnitude(I) / dt(J) - m(I) * dt(J)) * (step_time - time(J));

    return y;
}
