/*******************************************************************************
 * Copyright (C) 2017-2022 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/

#ifndef RESAMPLING_H
#define RESAMPLING_H

#include <armadillo>

enum class WindowType {
    Hamming,
    Hann,
    Blackman,
    BlackmanNuttall,
    BlackmanHarris,
    FlatTop
};

arma::uword gcd(arma::uword, arma::uword);

arma::vec hamming(arma::uword);
arma::vec hann(arma::uword);
arma::vec blackman(arma::uword);
arma::vec blackman_nuttall(arma::uword);
arma::vec blackman_harris(arma::uword);
arma::vec flat_top(arma::uword);

arma::vec fir_low_pass(arma::uword, double, arma::vec (*)(arma::uword));
arma::vec fir_high_pass(arma::uword, double, arma::vec (*)(arma::uword));
arma::vec fir_band_pass(arma::uword, double, double, arma::vec (*)(arma::uword));
arma::vec fir_band_stop(arma::uword, double, double, arma::vec (*)(arma::uword));

template<WindowType T> arma::vec fir_low_pass(arma::uword, double) { throw std::invalid_argument("not supported"); }

template<> arma::vec fir_low_pass<WindowType::Hamming>(arma::uword, double);

template<> arma::vec fir_low_pass<WindowType::Hann>(arma::uword, double);

template<> arma::vec fir_low_pass<WindowType::Blackman>(arma::uword, double);

template<> arma::vec fir_low_pass<WindowType::BlackmanNuttall>(arma::uword, double);

template<> arma::vec fir_low_pass<WindowType::BlackmanHarris>(arma::uword, double);

template<> arma::vec fir_low_pass<WindowType::FlatTop>(arma::uword, double);

template<WindowType T> arma::vec fir_high_pass(arma::uword, double) { throw std::invalid_argument("not supported"); }

template<> arma::vec fir_high_pass<WindowType::Hamming>(arma::uword, double);

template<> arma::vec fir_high_pass<WindowType::Hann>(arma::uword, double);

template<> arma::vec fir_high_pass<WindowType::Blackman>(arma::uword, double);

template<> arma::vec fir_high_pass<WindowType::BlackmanNuttall>(arma::uword, double);

template<> arma::vec fir_high_pass<WindowType::BlackmanHarris>(arma::uword, double);

template<> arma::vec fir_high_pass<WindowType::FlatTop>(arma::uword, double);

template<WindowType T> arma::vec fir_band_pass(arma::uword, double, double) { throw std::invalid_argument("not supported"); }

template<> arma::vec fir_band_pass<WindowType::Hamming>(arma::uword, double, double);

template<> arma::vec fir_band_pass<WindowType::Hann>(arma::uword, double, double);

template<> arma::vec fir_band_pass<WindowType::Blackman>(arma::uword, double, double);

template<> arma::vec fir_band_pass<WindowType::BlackmanNuttall>(arma::uword, double, double);

template<> arma::vec fir_band_pass<WindowType::BlackmanHarris>(arma::uword, double, double);

template<> arma::vec fir_band_pass<WindowType::FlatTop>(arma::uword, double, double);

template<WindowType T> arma::vec fir_band_stop(arma::uword, double, double) { throw std::invalid_argument("not supported"); }

template<> arma::vec fir_band_stop<WindowType::Hamming>(arma::uword, double, double);

template<> arma::vec fir_band_stop<WindowType::Hann>(arma::uword, double, double);

template<> arma::vec fir_band_stop<WindowType::Blackman>(arma::uword, double, double);

template<> arma::vec fir_band_stop<WindowType::BlackmanNuttall>(arma::uword, double, double);

template<> arma::vec fir_band_stop<WindowType::BlackmanHarris>(arma::uword, double, double);

template<> arma::vec fir_band_stop<WindowType::FlatTop>(arma::uword, double, double);

template<WindowType T> arma::vec upsampling(const arma::vec& in, const arma::uword up_rate) {
    arma::vec out(up_rate * in.n_elem, arma::fill::none);

    const auto coef = fir_low_pass<T>(8llu * up_rate, 1. / static_cast<double>(up_rate));

    const auto buffer_size = coef.n_elem;

    arma::vec buffer(buffer_size, arma::fill::zeros);
    auto index = 0llu;

    auto feed = [&](const double next) {
        buffer(index) = next;

        arma::uword p = 0;
        double result = 0.;
        for(auto m = index; m < buffer_size; ++m) result += coef(p++) * buffer(m);
        for(auto m = 0llu; m < index; ++m) result += coef(p++) * buffer(m);

        index = (index == 0 ? buffer_size : index) - 1;

        return result;
    };

    for(auto n = 0llu; n < out.n_elem; ++n) out(n) = static_cast<double>(up_rate) * feed(0 == n % up_rate ? in(n / up_rate) : 0.);

    return out;
}

arma::vec apply_filter(const arma::vec&, const arma::vec&);

#endif
