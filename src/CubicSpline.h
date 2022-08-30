#ifndef CUBICSPLINE_H
#define CUBICSPLINE_H

#include <armadillo>

class CubicSpline {
    const bool natural;
    arma::vec time;      // time
    arma::vec magnitude; // magnitude
    arma::vec dt, dy, m;

public:
    CubicSpline(const arma::vec&, const arma::vec&, bool = false);

    [[nodiscard]] double evaluate(double) const;
};

#endif // CUBICSPLINE_H
