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

#include "MedipixSPM.h"
#include "helper.h"

/**
 * Simulation that generated images with different frequencies
 */
int main() {
    auto m = std::make_shared<MedipixSPM>(false);
    m->set_psf_sigma(15.f);
    m->random_threshold_dispersion(1.0f);
    m->start_frame();
    frequency_exposure(m, 59.6f, 0.01, 550.0f, 0.0f, 1.0f, 0.0f, 1E6);
    m->finish_frame();
    m->save_image("frequency_pattern_spm.raw");
}