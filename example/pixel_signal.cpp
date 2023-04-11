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

#include <Medipix.h>
#include <helper.h>
#include <iostream>
#include "MedipixSPM.h"
#include "MedipixCSM.h"
#include <fstream>

int main(){
    auto m = std::make_shared<MedipixSPM>(true, 2, 2);
    m->random_threshold_dispersion(1.0f);
    m->set_th0(8.0f);
    m->set_psf_sigma(0.1);
    m->start_frame();
    homogeneous_exposure(m, 59.6, 1E-4, 2E7);
    m->finish_frame();
    m->save_pixel_signals("pixel_signal_0_0.raw", 0, 0);
    m->save_pixel_signals("pixel_signal_0_1.raw", 0, 1);
    m->save_pixel_signals("pixel_signal_1_0.raw", 1, 0);
    m->save_pixel_signals("pixel_signal_1_1.raw", 1, 1);
    std::cout << m->get_real_photons() << std::endl;
    std::cout << m->get_total_counts() << std::endl;
}
