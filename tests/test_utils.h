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

#ifndef MEDIPIX_TEST_UTILS_H
#define MEDIPIX_TEST_UTILS_H
#include "MedipixSPM.h"
#include "MedipixCSM.h"
#include "Medipix.h"


/**
 * A test class to expose protected methods of the different Medipix implementations.
 * @tparam T
 */
template<typename T>
class MedipixTest : public T {
public:
    MedipixTest(bool timed, unsigned int nx, unsigned int ny) : T(timed, nx, ny) {}
    float calculate_shared_energy(float x1, float y1, float energy, float x2, float y2) {
        return T::calculate_shared_energy(x1, y1, energy, x2, y2);
    }

};
#endif //MEDIPIX_TEST_UTILS_H
