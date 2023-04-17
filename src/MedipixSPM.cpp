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

#include <list>
#include "MedipixSPM.h"

MedipixSPM::MedipixSPM(bool timed, unsigned nx, unsigned int ny) : Medipix(timed, nx, ny) {

}

MedipixSPM::MedipixSPM() : Medipix(false) {

}

void MedipixSPM::add_photon(float energy, float position_x, float position_y, int radius, float time) {
    Medipix::add_photon(energy, position_x, position_y, radius, time);

    if (!timed) {
        std::list<std::pair<unsigned int, unsigned int>> pixels;
        auto [center_position_x, center_position_y] = get_pixel_index(position_x, position_y);
        for (int i = -radius; i < radius; ++i) {
            for (int j = -radius; j < radius; ++j) {
                int x_index = int(center_position_x) + i;
                int y_index = int(center_position_y) + j;
                if (x_index >= 0 && x_index < n_pixel_x && y_index >= 0 && y_index < n_pixel_y)
                    pixels.emplace_back(x_index, y_index);
            }
        }
        for (auto &pixel: pixels) {
            unsigned int i = pixel.first;
            unsigned int j = pixel.second;
            auto pixel_center = get_pixel_center(i, j);
            float dep_energy = calculate_shared_energy(position_x, position_y, energy, pixel_center.first,
                                                       pixel_center.second);
            if (dep_energy > get_th0(i, j)) {
                increase_counter(i, j);
            }
        }
    } else {
        std::list<std::pair<unsigned int, unsigned int>> pixels;
        auto [center_position_x, center_position_y] = get_pixel_index(position_x, position_y);
        for (int i = -radius; i < radius; ++i) {
            for (int j = -radius; j < radius; ++j) {
                int x_index = int(center_position_x) + i;
                int y_index = int(center_position_y) + j;
                if (x_index >= 0 && x_index < n_pixel_x && y_index >= 0 && y_index < n_pixel_y)
                    pixels.emplace_back(x_index, y_index);
            }
        }
        for (auto &pixel: pixels) {
            unsigned int i = pixel.first;
            unsigned int j = pixel.second;
            auto pixel_center = get_pixel_center(i, j);
            float dep_energy = calculate_shared_energy(position_x, position_y, energy, pixel_center.first,
                                                       pixel_center.second);
            Event event(time, dep_energy);
            events[i * n_pixel_y + j].emplace_back(event);
        }
    }
}

void MedipixSPM::finish_frame() {
    Medipix::finish_frame();
    if (timed) {
        std::lock_guard<std::mutex> lk(image_write_mutex);
        #pragma parallel for default(none) shared(events)
        for (unsigned int index = 0; index < n_pixel_x * n_pixel_y; ++index) {
            unsigned int i = index / n_pixel_y;
            unsigned int j = index % n_pixel_y;
            float threshold = get_th0(i, j);
            auto pixel_response = calculate_pixel_signal(i, j);
            for (unsigned int t = 1; t < pixel_response.size(); ++t) {
                if (pixel_response[t - 1] < threshold && pixel_response[t] > threshold) {
                    image[i * n_pixel_y + j] += 1;
                }
            }
        }
    }
}
