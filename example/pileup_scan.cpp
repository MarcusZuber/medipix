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

#include <fstream>
#include "helper.h"
#include "MedipixSPM.h"

/**
 * Simulation of flux dependent count-rate.
 */
int main(){
    auto m = std::make_shared<MedipixSPM>(true, 16, 16);
    m->random_threshold_dispersion(1.0f);
    m->set_th0(6.0f);
    m->set_psf_sigma(13.0);
    std::ofstream data_file;
    data_file.open("pileup_scan.txt");
    data_file << "# flux_density counts" << std::endl;
    double exposure_time = 1/10.f;
    double flux_scale = exposure_time * (m->get_pixel_pitch() * m->get_pixel_pitch() * 1E-3 * 1E-3 * m->get_num_pixels_y() * m->get_num_pixels_x());
    for(float flux_density=1E4; flux_density<1E8; flux_density*=1.5f){
        m->start_frame();
        homogeneous_exposure(m, 30., exposure_time, flux_density);
        m->finish_frame();

        data_file << flux_density << " " << m->get_real_photons() << " " << m->get_total_counts()/flux_scale << std::endl;

    }
    data_file.close();
}
