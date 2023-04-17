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

#ifndef MEDIPIX_HELPER_H
#define MEDIPIX_HELPER_H

#include <memory>
#include <functional>

class Medipix;

/**
 * @brief Simulates a homogeneous exposure of the Medipix
 * @param medipix
 * @param energy in keV
 * @param exposure_time in s
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted in photons / (s mm^2)
 */
[[maybe_unused]] void homogeneous_exposure(const std::shared_ptr<Medipix>& medipix, float energy, double exposure_time, double flux_density);

/**
 *
 * @param medipix
 * @param energy energy in keV
 * @param m slope of the edge
 * @param c y-intercept of the edge in um
 * @param exposure_time in s
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted in photons / (s mm^2)
 */
[[maybe_unused]] void edge_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float m, float c, double exposure_time, double flux_density);

/**
 *
 * @param medipix
 * @param energy in keV
 * @param exposure_time in s
 * @param period in um
 * @param phase in um
 * @param n_x Vector of the normal of the sin-wave in x-direction
 * @param n_y Vector of the normal of the sin-wave in y-direction
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted in photons / (s mm^2)
 */
[[maybe_unused]] void frequency_exposure(const std::shared_ptr<Medipix>& medipix, float energy, double exposure_time, float period, float phase, float n_x, float n_y, double flux_density);

/**
 *
 * @param medipix
 * @param energy in keV
 * @param exposure_time in s
 * @param flux_density in photons / (s mm^2)
 * @param photon_interacting Function that returns true if the photon is interacting with the detector.
 *  The function takes the x and y position of the photon in um as arguments.
 */
void exposure(const std::shared_ptr<Medipix>& medipix, float energy, double exposure_time, double flux_density, const std::function<bool (float, float)>& photon_interacting);


/**
 * Returns true if (x, y) is right of the line y = m * x + c
 * @param x
 * @param y
 * @param m
 * @param c
 * @return
 */
bool edge(float x, float y, float m, float c);

/**
 * Makes a random guess if the photon is interacting with the detector depending on the "attenuation" of the sin-wave.
 *
 * @param x
 * @param y
 * @param period
 * @param phase
 * @param n_x Normalized vector of the normal of the sin-wave in x-direction
 * @param n_y Normalized vector of the normal of the sin-wave in y-direction
 * @return
 */
bool frequency(float x, float y, float period, float phase, float n_x, float n_y);

#endif //MEDIPIX_HELPER_H
