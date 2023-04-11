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
#include <vector>
#include <iostream>

int main(){
    auto m = std::make_shared<MedipixSPM>(true, 256, 256);
    m->set_psf_sigma(12.0f);
    m->set_th0(6.0f);
    m->random_threshold_dispersion(1.0f);

    std::vector<float> flux{ 1E4, 5E4, 1E5, 5E5, 1E6, 5E6, 1E7, 5E7, 1E8, 5E8, 1E9};
    for (auto& f: flux){
        std::cout << f << std::endl;
        m->start_frame();
        homogeneous_exposure(m, 30.0f, 1/30.f, f);
        m->finish_frame();
        m->save_image("flux_" + std::to_string(f) + ".raw");
    }

}
