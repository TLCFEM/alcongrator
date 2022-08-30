#include "IntegrationScheme.h"

IntegrationScheme::IntegrationScheme(const double step_size)
    : DT(step_size) {
}

Newmark::Newmark(const double SS, const double A, const double B)
    : IntegrationScheme(SS), beta(A), gamma(B) {
    C5 = DT;
    C2 = 1. / beta / C5;
    C0 = C2 / C5;
    C4 = .5 / beta;
    C3 = C5 * gamma;
}

void Newmark::update_from_acceleration(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = next_a - current_a;
        const auto incre_v = C5 * current_a + C3 * incre_a;
        const auto incre_u = (incre_a + C2 * current_v + C4 * current_a) / C0;

        next_v = current_v + incre_v;
        next_u = current_u + incre_u;
    }
}

void Newmark::update_from_velocity(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = (next_v - current_v - C5 * current_a) / C3;
        const auto incre_u = (incre_a + C2 * current_v + C4 * current_a) / C0;

        next_a = current_a + incre_a;
        next_u = current_u + incre_u;
    }
}

void Newmark::update_from_displacement(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_u = next_u - current_u;
        const auto incre_a = C0 * incre_u - C2 * current_v - C4 * current_a;
        const auto incre_v = C5 * current_a + C3 * incre_a;

        next_v = current_v + incre_v;
        next_a = current_a + incre_a;
    }
}

BatheTwoStep::BatheTwoStep(const double SS)
    : IntegrationScheme(SS) {
    C0 = DT;
    C1 = .5 / C0;
    C2 = 3. * C1;
    C3 = 4. * C1;
    C4 = 2. * C3;
    C6 = C4 / C0;
}

void BatheTwoStep::update_from_acceleration(arma::vec& u, arma::vec& v, arma::vec& a) {
    enum class FLAG {
        TRAP,
        EULER
    };

    FLAG step_flag = FLAG::TRAP;

    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        if(FLAG::TRAP == step_flag) {
            const auto incre_a = next_a - current_a;
            const auto incre_u = incre_a / C6 + C0 * current_v + 2. / C6 * current_a;
            const auto incre_v = C3 / C6 * (incre_a + C4 * current_v + 2. * current_a) - 2. * current_v;

            next_v = current_v + incre_v;
            next_u = current_u + incre_u;

            step_flag = FLAG::EULER;
        }
        else {
            const auto pre_u = u(J - 1);
            const auto pre_v = v(J - 1);

            next_v = (next_a - C1 * pre_v + C3 * current_v) / C2;
            next_u = current_u + (next_v + C1 * (current_u - pre_u)) / C2;

            step_flag = FLAG::TRAP;
        }
    }
}

void BatheTwoStep::update_from_velocity(arma::vec& u, arma::vec& v, arma::vec& a) {
    enum class FLAG {
        TRAP,
        EULER
    };

    FLAG step_flag = FLAG::TRAP;

    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        if(FLAG::TRAP == step_flag) {
            const auto incre_a = C6 / C3 * next_v + (C6 / C3 - C4) * current_v - 2. * current_a;
            const auto incre_u = incre_a / C6 + C0 * current_v + 2. / C6 * current_a;

            next_a = current_a + incre_a;
            next_u = current_u + incre_u;

            step_flag = FLAG::EULER;
        }
        else {
            const auto pre_u = u(J - 1);
            const auto pre_v = v(J - 1);

            next_a = C2 * next_v + C1 * pre_v - C3 * current_v;
            next_u = current_u + (next_v + C1 * (current_u - pre_u)) / C2;

            step_flag = FLAG::TRAP;
        }
    }
}

void BatheTwoStep::update_from_displacement(arma::vec& u, arma::vec& v, arma::vec& a) {
    enum class FLAG {
        TRAP,
        EULER
    };

    FLAG step_flag = FLAG::TRAP;

    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        if(FLAG::TRAP == step_flag) {
            const auto incre_u = next_u - current_u;
            const auto incre_a = C6 * incre_u - C0 * C6 * current_v - 2. * current_a;
            const auto incre_v = C3 / C6 * (incre_a + C4 * current_v + 2. * current_a) - 2. * current_v;

            next_v = current_v + incre_v;
            next_a = current_a + incre_a;

            step_flag = FLAG::EULER;
        }
        else {
            const auto pre_u = u(J - 1);
            const auto pre_v = v(J - 1);

            next_v = C2 * (next_u - current_u) + C1 * (pre_u - current_u);
            next_a = C2 * next_v + C1 * pre_v - C3 * current_v;

            step_flag = FLAG::TRAP;
        }
    }
}

GeneralizedAlpha::GeneralizedAlpha(const double SS, const double R)
    : IntegrationScheme(SS)
    , gamma(.5 - (R - 1.) / (R + 1.))
    , beta(pow(R + 1., -2.))
    , F9(-.5 / beta) {
    F10 = DT;
    F11 = F10 * gamma;
    F8 = -1. / beta / F10;
    F7 = -F8 / F10;
}

void GeneralizedAlpha::update_from_acceleration(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = next_a - current_a;
        const auto incre_v = F10 * current_a + F11 * incre_a;
        const auto incre_u = incre_a / F7 + F10 * current_v - F9 / F7 * current_a;

        next_v = current_v + incre_v;
        next_u = current_u + incre_u;
    }
}

void GeneralizedAlpha::update_from_velocity(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = (next_v - current_v - F10 * current_a) / F11;
        const auto incre_u = incre_a / F7 + F10 * current_v - F9 / F7 * current_a;

        next_a = current_a + incre_a;
        next_u = current_u + incre_u;
    }
}

void GeneralizedAlpha::update_from_displacement(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_u = next_u - current_u;
        const auto incre_a = F7 * incre_u + F8 * current_v + F9 * current_a;
        const auto incre_v = F10 * current_a + F11 * incre_a;

        next_v = current_v + incre_v;
        next_a = current_a + incre_a;
    }
}

GSSSS::GSSSS(const double SS)
    : IntegrationScheme(SS) {
}

void GSSSS::update_from_acceleration(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = next_a - current_a;
        const auto incre_v = DT * (L4 * current_a + L5 * incre_a);
        const auto incre_u = DT * DT * (L1 / DT * current_v + L2 * current_a + L3 * incre_a);

        next_v = current_v + incre_v;
        next_u = current_u + incre_u;
    }
}

void GSSSS::update_from_velocity(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_a = ((next_v - current_v) / DT - L4 * current_a) / L5;
        const auto incre_u = DT * DT * (L1 / DT * current_v + L2 * current_a + L3 * incre_a);

        next_a = current_a + incre_a;
        next_u = current_u + incre_u;
    }
}

void GSSSS::update_from_displacement(arma::vec& u, arma::vec& v, arma::vec& a) {
    for(auto I = 1llu, J = 0llu; I < a.n_elem; ++I, ++J) {
        const auto current_u = u(J);
        const auto current_v = v(J);
        const auto current_a = a(J);
        auto& next_u = u(I);
        auto& next_v = v(I);
        auto& next_a = a(I);

        const auto incre_u = next_u - current_u;
        const auto incre_a = (incre_u / DT / DT - L1 / DT * current_v - L2 * current_a) / L3;
        const auto incre_v = DT * (L4 * current_a + L5 * incre_a);

        next_v = current_v + incre_v;
        next_a = current_a + incre_a;
    }
}

template<> void GSSSS::generate_constants<GSSSSU0>(double, const double R1, const double R2) {
    L1 = 1.;
    L2 = .5;
    L3 = 1. / (1. + R1) / (1. + R2);
    L4 = 1.;
    L5 = .5 * (3. + R1 + R2 - R1 * R2) * L3;
}

template<> void GSSSS::generate_constants<GSSSSV0>(const double R3, double, double) {
    L1 = 1.;
    L2 = .5;
    L3 = .5 / (1. + R3);
    L4 = 1.;
    L5 = 2. * L3;
}

GSSSSU0::GSSSSU0(const double SS, arma::vec&& R)
    : GSSSS(SS) {
    R = arma::sort(R.clamp(0., 1.));
    generate_constants<GSSSSU0>(R(0), R(1), R(2));
}

GSSSSV0::GSSSSV0(const double SS, arma::vec&& R)
    : GSSSS(SS) {
    R = arma::sort(R.clamp(0., 1.));
    generate_constants<GSSSSV0>(R(0), R(1), R(2));
}
