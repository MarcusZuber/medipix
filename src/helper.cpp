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
#include <omp.h>
#include <chrono>


void exposure(const std::shared_ptr<Medipix> &medipix, float energy, double exposure_time, double flux_density,
              const std::function<bool(float, float)> &photon_interacting) {
    //flux density in photons per second per square mm

    float total_area =
            float(medipix->get_num_pixels_x() * medipix->get_num_pixels_y()) * (medipix->get_pixel_pitch() * 0.001f) *
            (medipix->get_pixel_pitch() * 0.001f);
    double duration = exposure_time * double(1E6);
    unsigned int number_of_photons = int(flux_density * total_area * exposure_time);

    unsigned int seed = time(nullptr);
    #pragma omp parallel for default(none) shared(medipix, energy, number_of_photons, photon_interacting, duration, seed)
    for (unsigned int i = 0; i < number_of_photons; ++i) {
        std::random_device randomDevice;
        std::mt19937 twisterEngine(randomDevice());
        twisterEngine.seed(seed + i);
        std::uniform_real_distribution<float> distribution_x(medipix->get_min_x(), medipix->get_max_x());
        std::uniform_real_distribution<float> distribution_y(medipix->get_min_y(), medipix->get_max_y());
        std::uniform_real_distribution<float> distribution_t(0.f, float(duration));
        float x = distribution_x(twisterEngine);
        float y = distribution_y(twisterEngine);
        float t = distribution_t(twisterEngine);

        if (photon_interacting(x, y)) {
            medipix->add_photon(energy, x, y, 3, t);
        }
    }

}

[[maybe_unused]] void
homogeneous_exposure(const std::shared_ptr<Medipix> &medipix, float energy, double exposure_time, double flux_density) {
    exposure(medipix, energy, exposure_time, flux_density, [](float x, float y) { return true; });
}

[[maybe_unused]] void
edge_exposure(const std::shared_ptr<Medipix> &medipix, float energy, float m, float c, double exposure_time,
              double flux_density) {
    std::function<bool(float, float)> g = [m, c](float _x, float _y) { return edge(_x, _y, m, c); };
    exposure(medipix, energy, exposure_time, flux_density, g);
}

[[maybe_unused]] void
frequency_exposure(const std::shared_ptr<Medipix> &medipix, float energy, double exposure_time, float period,
                   float phase, float n_x, float n_y, double flux_density) {
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

float get_random_free_path(float mu_norm, float rho) {
    float mu = mu_norm * rho;

    // From beer-lambert:
    // PDF(x) = exp(-mu * x) / mu
    // CTF(x) =  (-1/mu) * (exp(-mu * x) - 1)
    // CTF^-1(y) = -ln(1 - y * mu) / mu.

    // With the inverse transfer method we sample now the interaction depth.

    static std::random_device rd;
    static std::mt19937 gen(rd());
    float max_range = std::min(1.0f, 1.0f/mu);
    static std::uniform_real_distribution<float> dis(0.0, max_range);
    float u = dis(gen);
    float x = -std::log(1 - u * mu) / mu;
    //while (x <= 0 || x!=x) {
    //    u = dis(gen);
    //    x = -std::log(1 - u * mu) / mu;
    //}
    return x * 1E4f; // in um
}


bool frequency(float x, float y, float period, float phase, float n_x, float n_y) {
    float a1 = n_x * x + n_y * y;
    float probability = sinf(2.f * M_PIf * a1 / period + phase);
    std::uniform_real_distribution<float> distribution(0, 1);
    unsigned int seed = time(nullptr);
    std::default_random_engine freq_generator(seed);
    return (distribution(freq_generator) < probability);
}
