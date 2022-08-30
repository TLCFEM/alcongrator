#ifndef INTEGRATIONSCHEME_H
#define INTEGRATIONSCHEME_H

#include <armadillo>

class IntegrationScheme {
protected:
    const double DT;

public:
    explicit IntegrationScheme(double);

    virtual void update_from_acceleration(arma::vec&, arma::vec&, arma::vec&) = 0;
    virtual void update_from_velocity(arma::vec&, arma::vec&, arma::vec&) = 0;
    virtual void update_from_displacement(arma::vec&, arma::vec&, arma::vec&) = 0;
};

class Newmark : public IntegrationScheme {
    const double beta;  /**< parameter = .25 */
    const double gamma; /**< parameter = .5 */

    double C0 = 0., C2 = 0., C3 = 0., C4 = 0., C5 = 0.; /**< parameters */
public:
    Newmark(double, double, double);

    void update_from_acceleration(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_velocity(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_displacement(arma::vec&, arma::vec&, arma::vec&) override;
};

class BatheTwoStep : public IntegrationScheme {
    double C0 = 0., C1 = 0., C2 = 0., C3 = 0., C4 = 0., C6 = 0.;

public:
    explicit BatheTwoStep(double);

    void update_from_acceleration(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_velocity(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_displacement(arma::vec&, arma::vec&, arma::vec&) override;
};

class GeneralizedAlpha : public IntegrationScheme {
    const double gamma;
    const double beta;

    const double F9;

    double F7 = 0., F8 = 0., F10 = 0., F11 = 0.;

public:
    GeneralizedAlpha(double, double);

    void update_from_acceleration(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_velocity(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_displacement(arma::vec&, arma::vec&, arma::vec&) override;
};

class GSSSS : public IntegrationScheme {
protected:
    double L1 = 0., L2 = 0., L3 = 0., L4 = 0., L5 = 0.;

    // ReSharper disable once CppMemberFunctionMayBeStatic
    template<typename T> void generate_constants(double, double, double) { throw std::invalid_argument("need a proper scheme"); }

public:
    explicit GSSSS(double);

    void update_from_acceleration(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_velocity(arma::vec&, arma::vec&, arma::vec&) override;
    void update_from_displacement(arma::vec&, arma::vec&, arma::vec&) override;
};

class GSSSSU0 final : public GSSSS {
public:
    GSSSSU0(double, arma::vec&&);
};

class GSSSSV0 final : public GSSSS {
public:
    GSSSSV0(double, arma::vec&&);
};

#endif // INTEGRATIONSCHEME_H
