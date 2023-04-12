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

int main(){
    auto m = std::make_shared<MedipixSPM>(true, 8, 8);
    m->random_threshold_dispersion(1.0f);
    m->set_th0(8.0f);
    m->set_psf_sigma(13.0);
    std::ofstream data_file;
    data_file.open("pileup_scan.txt");
    data_file << "# flux_density counts" << std::endl;
    for(float flux_density=1E5; flux_density<1E9; flux_density+=5E5){
        m->start_frame();
        homogeneous_exposure(m, 59.6, 1/1000., flux_density);
        m->finish_frame();
        data_file << flux_density << " " << m->get_real_photons() << " " << m->get_total_counts() << std::endl;

    }
    data_file.close();
}
