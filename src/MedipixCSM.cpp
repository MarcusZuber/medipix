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

#include <map>
#include <random>
#include "MedipixCSM.h"

MedipixCSM::MedipixCSM(bool timed, unsigned int nx, unsigned ny): Medipix(timed, nx, ny){
    th1_dispersion.resize(n_pixel_x * n_pixel_y);

}

MedipixCSM::MedipixCSM(): Medipix(false) {
    th0_dispersion.resize(n_pixel_x * n_pixel_y);
}

void MedipixCSM::add_photon(float energy, float position_x, float position_y, int radius, float time) {
    Medipix::add_photon(energy, position_x, position_y, radius, time);


    if (!timed) {
        // NOTE: We assume that the charge is only be shared between 4 pixel -> only one summing node is activated!
        // TODO: Implement a more general solution
        radius = 1;
        std::list<std::pair<unsigned int, unsigned int>> pixels;
        auto [center_position_x, center_position_y] = get_pixel_index(position_x, position_y);
        auto [pixel_center_index_x, pixel_center_index_y] = get_pixel_center(center_position_x, center_position_y);
        unsigned int x_shift = 0;
        if (position_x < pixel_center_index_x){
            x_shift = -1;
        }
        else{
            x_shift = 1;
        }
        unsigned int y_shift = 0;
        if (position_y < pixel_center_index_y){
            y_shift = -1;
        }
        else{
            y_shift = 1;
        }
        float summed_energy = 0;
        std::vector<unsigned int> index_i { center_position_x, center_position_x + x_shift};
        std::vector<unsigned int> index_j { center_position_y, center_position_y + y_shift};
        for(auto& i: index_i){
            for(auto& j: index_j){
                if (i >= 0 && i < n_pixel_x && j >=0 && j < n_pixel_y){
                    auto pixel_center = get_pixel_center(i, j);
                    float dep_energy = calculate_shared_energy(position_x, position_y, energy, pixel_center.first, pixel_center.second);
                    if(dep_energy > get_th0(i, j)){
                        summed_energy += dep_energy;
                    }
                }
            }
        }
        if (summed_energy > get_th1(center_position_x, center_position_y)){
            increase_counter(center_position_x, center_position_y);
        }


    }

    else{
        std::list<std::pair<unsigned int, unsigned int>> pixels;
        auto [center_position_x, center_position_y] = get_pixel_index(position_x, position_y);
        for (int i = -radius; i<radius; ++i){
            for (int j = -radius; j < radius; ++j){
                int x_index = int(center_position_x) + i;
                int y_index = int(center_position_y) + j;
                if (x_index >= 0 && x_index < n_pixel_x && y_index >=0 && y_index < n_pixel_y)
                    pixels.emplace_back(x_index, y_index);
            }
        }
        for (auto& pixel : pixels){
            unsigned int i = pixel.first;
            unsigned int j = pixel.second;
            auto pixel_center = get_pixel_center(i, j);
            float dep_energy = calculate_shared_energy(position_x, position_y, energy, pixel_center.first, pixel_center.second);
            Event event(time, dep_energy);
            events[i*n_pixel_y + j].emplace_back(event);
        }
    }

}

void MedipixCSM::finish_frame() {
    Medipix::finish_frame();
}

void MedipixCSM::random_threshold_dispersion(float sigma) {
    Medipix::random_threshold_dispersion(sigma);
    std::normal_distribution<float> distribution(0.f, sigma);
    unsigned int seed = time(nullptr);
    std::default_random_engine generator(seed);

#pragma omp parallel for default(none) shared(generator, distribution)
    for (unsigned int i = 0; i < n_pixel_y * n_pixel_x; ++i) {
        th1_dispersion[i] = distribution(generator);
    }
}

void MedipixCSM::set_th1(float t) {
    th1 = t;
}

float MedipixCSM::get_th1(unsigned int i, unsigned int j) {
    return th1 + th1_dispersion[i*n_pixel_y + j];
}

[[maybe_unused]] float MedipixCSM::get_th1() {
    return th1;
}
