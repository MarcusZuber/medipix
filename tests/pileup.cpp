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

TEST(Pileup, SpmSinglePixel) {
    MedipixTest<MedipixSPM> m(true, 8, 8);
    m.set_psf_sigma(1.0f);
    float energy = 30.0f;
    m.set_th0(6.0f);
    auto [x, y] = m.get_pixel_center(4, 4);

    // Test single hit
    m.start_frame();
    m.add_photon(energy, x, y, 3, 10);
    m.finish_frame();
    EXPECT_EQ(m.get_pixel_value(4, 4), 1);

    // Test two hit with pileup (pixel response is in the order of 1us)
    m.start_frame();
    m.add_photon(energy, x, y, 3, 10.);
    m.add_photon(energy, x, y, 3, 10.5);
    m.finish_frame();
    EXPECT_EQ(m.get_pixel_value(4, 4), 1);

    // Test two hit without pileup (pixel response is in the order of 1us)
    m.start_frame();
    m.add_photon(energy, x, y, 3, 10.);
    m.add_photon(energy, x, y, 3, 12);
    m.finish_frame();
    EXPECT_EQ(m.get_pixel_value(4, 4), 2);

}

TEST(Pileup, SpmChargeSharing) {
    MedipixTest<MedipixSPM> m(true, 8, 8);
    m.set_psf_sigma(15.0f);
    float energy = 30.0f;
    m.set_th0(6.0f);

    // Two photons close to the border between pixel (4, 4) and (4, 5).
    auto [x, y] = m.get_pixel_center(4, 4);
    y += 0.48f * m.get_pixel_pitch();
    auto [x2, y2] = m.get_pixel_center(4, 5);
    y2 += -0.48f * m.get_pixel_pitch();

    // Test two hit from adjacent pixels with pileup (pixel response is in the order of 1us)
    m.start_frame();
    m.add_photon(energy, x2, y2, 3, 10.);
    m.add_photon(energy, x, y, 3, 10.);
    m.finish_frame();
    EXPECT_EQ(m.get_pixel_value(4, 4), 1);
    EXPECT_EQ(m.get_pixel_value(4, 5), 1);

    // Test two hit from adjacent pixels without pileup (pixel response is in the order of 1us)
    m.start_frame();
    m.add_photon(energy, x2, y2, 3, 10.);
    m.add_photon(energy, x, y, 3, 12.);
    m.finish_frame();
    EXPECT_EQ(m.get_pixel_value(4, 4), 2);
    EXPECT_EQ(m.get_pixel_value(4, 5), 2);

}