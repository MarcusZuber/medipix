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

#include <gtest/gtest.h>
#include <memory>
#include "test_utils.h"

TEST(ChargeSharing, ChargeFractions) {
    /**
     * Checks charge sharing. A photon with a small psf is used.
     */

    MedipixTest<Medipix> m(false, 256, 256);
    m.set_psf_sigma(0.1);
    auto [x, y] = m.get_pixel_center(0, 0);
    float energy = 30.0f;
    float deposited_energy;

    // Check if center hit gives full energy.
    deposited_energy = m.calculate_shared_energy(x, y, energy, x, y);
    EXPECT_FLOAT_EQ(deposited_energy, energy);

    // Check if a hit in the edge center gives 1/2 energy.
    deposited_energy = m.calculate_shared_energy(x,
                                                 y,
                                                 energy,
                                                 x + m.get_pixel_pitch() / 2.f,
                                                 y);
    EXPECT_FLOAT_EQ(deposited_energy, energy / 2.f);

    // Check if a hit in a corner center gives 1/4 energy.
    deposited_energy = m.calculate_shared_energy(x,
                                                 y,
                                                 30.0f,
                                                 x + m.get_pixel_pitch() / 2.f,
                                                 y + m.get_pixel_pitch() / 2.f);
    EXPECT_FLOAT_EQ(deposited_energy, energy / 4.);
}

TEST(ChargeSharing, ChargeSum) {
    /**
     * Checks if overall charge is conserved.
     */

    MedipixTest<Medipix> m(false, 256, 256);
    m.set_psf_sigma(12);
    auto [x, y] = m.get_pixel_center(128, 128);
    float energy = 30.0f;
    float deposited_energy;

    // Check if center hit gives full energy.
    deposited_energy = 0.0f;
    for (int i = -10; i <= 10; ++i) {
        for (int j = -10; j <= 10; ++j) {
            auto [x_i, y_j] = m.get_pixel_center(128 + i, 128 + j);
            deposited_energy += m.calculate_shared_energy(x, y, energy, x_i, y_j);
        }
    }
    EXPECT_FLOAT_EQ(deposited_energy, energy);


    // Check edge center hit
    deposited_energy = 0.0f;
    for (int i = -10; i <= 10; ++i) {
        for (int j = -10; j <= 10; ++j) {
            auto [x_i, y_j] = m.get_pixel_center(128 + i, 128 + j);
            deposited_energy += m.calculate_shared_energy(x, y, energy, x_i + m.get_pixel_pitch() / 2.f, y_j);
        }
    }
    EXPECT_FLOAT_EQ(deposited_energy, energy);

    // Check corner hit
    deposited_energy = 0.0f;
    for (int i = -10; i <= 10; ++i) {
        for (int j = -10; j <= 10; ++j) {
            auto [x_i, y_j] = m.get_pixel_center(128 + i, 128 + j);
            deposited_energy += m.calculate_shared_energy(x, y, energy, x_i + m.get_pixel_pitch() / 2.f,
                                                          y_j + m.get_pixel_pitch() / 2.f);
        }
    }
    EXPECT_FLOAT_EQ(deposited_energy, energy);
}
