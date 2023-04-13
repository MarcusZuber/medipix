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

/**
 * Performs a threshold scan in SPM and CSM
 */
int main(){
    auto m = std::make_shared<MedipixSPM>(false);
    m->set_psf_sigma(14.f);
    m->random_threshold_dispersion(1.5f);
    std::ofstream data_file;
    data_file.open("threshold_scan_spm.txt");
    data_file << "# th0 counts" << std::endl;
    for(float th = 6.0f; th < 50.f; th += 0.5f){
        m->set_th0(th);
        std::cout << "th0: " << th << std::endl;
        m->start_frame();
        homogeneous_exposure(m, 40.0f, 0.01, 1E6);
        m->finish_frame();
        data_file << th << ' ' << m->get_total_counts() << std::endl;
    }
    data_file.close();


    auto m2 = std::make_shared<MedipixCSM>(false);
    m2->set_psf_sigma(14.f);
    m2->random_threshold_dispersion(1.5f);
    std::ofstream data_file_csm;
    data_file_csm.open("threshold_scan_csm.txt");
    data_file_csm << "# th1 counts" << std::endl;
    m2->set_th0(5.0f);
    for(float th = 6.0f; th < 50.f; th += 0.5f){
        std::cout << "th1: " << th << std::endl;
        m2->set_th1(th);
        m2->start_frame();
        homogeneous_exposure(m2, 40.f, 0.01, 1E6);
        m2->finish_frame();
        data_file_csm << th << " " << (m2->get_total_counts()) << std::endl;
    }
    data_file_csm.close();
}

