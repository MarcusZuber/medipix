//
// Created by marcus on 10.04.23.
//

#include <list>
#include "MedipixSPM.h"

MedipixSPM::MedipixSPM(bool timed, unsigned nx, unsigned int ny): Medipix(timed, nx, ny){

}

MedipixSPM::MedipixSPM(): Medipix(false) {

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

void MedipixSPM::finish_frame() {
    Medipix::finish_frame();
    if (timed){
        #pragma parallel for collapse(2) default(none) shared(events)
        for (unsigned int i = 0; i < n_pixel_x; ++i){
            for (unsigned int j = 0; j < n_pixel_y; ++j){
                auto pixel_response = calculate_pixel_signal(i, j);
                for(unsigned int t = 1; t < pixel_response.size() - 1; ++t){
                    if(pixel_response[t-1] < get_th0(i, j) && pixel_response[t+1] > get_th0(i, j)) {
                        increase_counter(i, j);
                    }
                }
            }
        }
    }
}
