//
// Created by marcus on 05.04.2023.
//

#include "Medipix.h"

#include <cmath>
#include <list>
#include <random>
#include <fstream>
#include <future>


Medipix::Medipix() {
    image.reset(new unsigned int[n_pixel_x * n_pixel_y]);
    for (unsigned int i = 0; i < n_pixel_y*n_pixel_x; ++i){
        image[i] = 0;
    }
}

[[maybe_unused]] void Medipix::reset_image() {
    image.reset(new unsigned int[n_pixel_x * n_pixel_y]);
    for (unsigned int i = 0; i < n_pixel_y*n_pixel_x; ++i){
        image[i] = 0;
    }
}

float Medipix::get_min_x() const {
    return  -pixel_pitch * float(n_pixel_x) / 2.f;
}

float Medipix::get_max_x() const {
    return pixel_pitch * float(n_pixel_x) / 2.f;
}

float Medipix::get_min_y() const {
    return  -pixel_pitch * float(n_pixel_y) / 2.f;
}

float Medipix::get_max_y() const {
    return pixel_pitch * float(n_pixel_y) / 2.f;
}

[[maybe_unused]] float Medipix::get_th0() const {
    return th0;
}

[[maybe_unused]] float Medipix::get_pixel_pitch() const {
    return pixel_pitch;
}

[[maybe_unused]] float Medipix::get_psf_sigma() const {
    return psf_sigma;
}

std::pair<unsigned int, unsigned int> Medipix::get_pixel_index(float position_x, float position_y) const {
    int i = int(position_x / pixel_pitch + float(n_pixel_x) / 2.f - 0.5f);
    int j =  int(position_y / pixel_pitch + float(n_pixel_y) / 2.f - 0.5f);
    return std::make_pair(i, j);
}

std::pair<float, float> Medipix::get_pixel_center(unsigned int i, unsigned int j) const {
    float x = pixel_pitch * (float(i) - float(n_pixel_x) / 2.f + 0.5f);
    float y = pixel_pitch * (float(j) - float(n_pixel_y) / 2.f + 0.5f);
    return std::make_pair(x, y);
}

float Medipix::calculate_shared_energy(float x, float y, float energy, float pixel_center_x, float pixel_center_y) const {
    float x_b = (pixel_center_x + pixel_pitch / 2.f);
    float x_a = (pixel_center_x - pixel_pitch / 2.f);
    float y_b = (pixel_center_y + pixel_pitch / 2.f);
    float y_a = (pixel_center_y - pixel_pitch / 2.f);

    float x_component = std::erf((x_b - x) / (psf_sigma * float(std::numbers::sqrt2))) - std::erf((x_a - x) / (psf_sigma * float(std::numbers::sqrt2)));
    float y_component = std::erf((y_b - y) / (psf_sigma * float(std::numbers::sqrt2))) - std::erf((y_a - y) / (psf_sigma * float(std::numbers::sqrt2)));

    return 0.25f * energy * x_component * y_component;

}

void Medipix::add_photon(float energy, float position_x, float position_y, int radius) {
    if (energy < th0)
        return;

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
        if (dep_energy > th0){
            increase_counter(i, j);
        }
    }

}

void Medipix::increase_counter(unsigned int x, unsigned int y) {
    std::lock_guard<std::mutex> lk(image_write_mutex);
    image[x * n_pixel_y + y] += 1;
}

void Medipix::homogeneous_exposure(float energy, unsigned int number_of_photons) {
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distribution_x(get_min_x(),get_max_x());
    std::uniform_real_distribution<float> distribution_y(get_min_y(),get_max_y());
    for (unsigned int i = 0; i < number_of_photons; ++i){
        float x = distribution_x(generator);
        float y = distribution_x(generator);
        add_photon(energy, x, y, 3);
    }

}

unsigned int Medipix::get_total_counts() {
    unsigned int counts = 0;
    for (unsigned int i=0; i < n_pixel_y*n_pixel_y; ++i){
        counts += image[i];
    }
    return counts;
}

void Medipix::save_image(std::string filename) const {
    std::ofstream image_file(filename, std::ios::out | std::ios::binary);
    for (unsigned int i=0; i < n_pixel_y*n_pixel_y; ++i){
        uint32_t pixel_value = image[i];
        image_file.write((char *) &pixel_value, sizeof(pixel_value));
    }
    image_file.close();

}

void Medipix::set_psf_sigma(float s) {
    psf_sigma = s;
}

void Medipix::set_th0(float s) {
    th0 = s;
}
