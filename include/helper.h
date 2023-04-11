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

[[maybe_unused]] void homogeneous_exposure(const std::shared_ptr<Medipix>& medipix, float energy, unsigned int number_of_photons, float duration);

[[maybe_unused]] void edge_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float m, float c, unsigned int number_of_photons, float duration);

[[maybe_unused]] [[maybe_unused]] void frequency_exposure(const std::shared_ptr<Medipix>& medipix,float energy, unsigned int number_of_photons, float period, float phase, float n_x, float n_y, float duration);

void exposure(const std::shared_ptr<Medipix>& medipix, float energy, unsigned int number_of_photons, float duration, const std::function<bool (float, float)>& photon_interacting);

bool edge(float x, float y, float m, float c);

bool frequency(float x, float y, float period, float phase, float n_x, float n_y);

#endif //MEDIPIX_HELPER_H
