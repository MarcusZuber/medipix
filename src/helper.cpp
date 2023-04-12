/*
 * Simple Medipix simulation
 * Copyright (C) 2023  Marcus Zuber
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
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <random>
#include "helper.h"
#include "Medipix.h"
#include <ctime>
#include <iostream>


void exposure(const std::shared_ptr<Medipix> &medipix, float energy, float exposure_time, float flux_density,
              const std::function<bool(float, float)> &photon_interacting) {
    //flux density in photons per second per square mm

    float total_area =
            float(medipix->get_num_pixels_x() * medipix->get_num_pixels_y()) * (medipix->get_pixel_pitch() * 0.001f) *
            (medipix->get_pixel_pitch() * 0.001f);
    float duration = exposure_time * float(1E6);
    unsigned int number_of_photons = int(flux_density * total_area * exposure_time);
    unsigned int seed = time(nullptr);
    std::default_random_engine generator_x(seed);
    std::default_random_engine generator_y(seed + 1);
    std::default_random_engine generator_t(seed + 2);
    std::uniform_real_distribution<float> distribution_x(medipix->get_min_x(), medipix->get_max_x());
    std::uniform_real_distribution<float> distribution_y(medipix->get_min_y(), medipix->get_max_y());
    std::uniform_real_distribution<float> distribution_t(0, duration);
    if (medipix->get_timed()) {
        for (unsigned int i = 0; i < number_of_photons; ++i) {
            float x = distribution_x(generator_x);
            float y = distribution_y(generator_y);
            float t = distribution_t(generator_t);
            if (photon_interacting(x, y)) {
                medipix->add_photon(energy, x, y, 3, t);
            }
        }
    } else {
        #pragma omp parallel for default(none) shared(generator_x, generator_y, generator_t, medipix, energy, distribution_x, distribution_y, distribution_t, number_of_photons, photon_interacting)
        for (unsigned int i = 0; i < number_of_photons; ++i) {

            float x = distribution_x(generator_x);
            float y = distribution_y(generator_y);
            float t = distribution_t(generator_t);
            if (photon_interacting(x, y)) {
                medipix->add_photon(energy, x, y, 3, t);
            }
        }
    }
}

[[maybe_unused]] void
homogeneous_exposure(const std::shared_ptr<Medipix> &medipix, float energy, float exposure_time, float flux_density) {
    exposure(medipix, energy, exposure_time, flux_density, [](float x, float y) { return true; });
}

[[maybe_unused]] void
edge_exposure(const std::shared_ptr<Medipix> &medipix, float energy, float m, float c, float exposure_time,
              float flux_density) {
    std::function<bool(float, float)> g = [m, c](float _x, float _y) { return edge(_x, _y, m, c); };
    exposure(medipix, energy, exposure_time, flux_density, g);
}

[[maybe_unused]] void
frequency_exposure(const std::shared_ptr<Medipix> &medipix, float energy, float exposure_time, float period,
                   float phase, float n_x, float n_y, float flux_density) {
    float r = std::sqrt(n_x * n_x + n_y * n_y);
    n_x /= r;
    n_y /= r;
    std::function<bool(float, float)> g = [period, phase, n_x, n_y](float _x, float _y) {
        return frequency(_x, _y, period, phase, n_x, n_y);
    };
    exposure(medipix, energy, exposure_time, flux_density, g);
}

bool edge(float x, float y, float m, float c) {
    return y > m * x + c;
}

bool frequency(float x, float y, float period, float phase, float n_x, float n_y) {
    float a1 = n_x * x + n_y * y;
    float probability = sinf(2.f * M_PIf * a1 / period + phase);
    std::uniform_real_distribution<float> distribution(0, 1);
    unsigned int seed = time(nullptr);
    std::default_random_engine freq_generator(seed);
    return (distribution(freq_generator) < probability);
}
