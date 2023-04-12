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
 * @brief Homogeneous exposure of the Medipix
 * @param medipix
 * @param energy
 * @param number_of_photons
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted
 */
[[maybe_unused]] void homogeneous_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float exposure_time, float flux_density);

/**
 *
 * @param medipix
 * @param energy
 * @param m
 * @param c
 * @param number_of_photons
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted
 */
[[maybe_unused]] void edge_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float m, float c, float exposure_time, float flux_density);

/**
 *
 * @param medipix
 * @param energy
 * @param number_of_photons
 * @param period
 * @param phase
 * @param n_x
 * @param n_y
 * @param flux_density Flux density is used to calculate the time period in which number_of_photons are emitted
 */
[[maybe_unused]] void frequency_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float exposure_time, float period, float phase, float n_x, float n_y, float flux_density);

/**
 *
 * @param medipix
 * @param energy in keV
 * @param exposure_time in us
 * @param flux_density in photons / (s mm^2)
 * @param photon_interacting
 */
void exposure(const std::shared_ptr<Medipix>& medipix, float energy, float exposure_time, float flux_density, const std::function<bool (float, float)>& photon_interacting);

bool edge(float x, float y, float m, float c);

bool frequency(float x, float y, float period, float phase, float n_x, float n_y);

#endif //MEDIPIX_HELPER_H
