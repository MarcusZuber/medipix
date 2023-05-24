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

#include "Medipix.h"

#include <cmath>
#include <list>
#include <random>
#include <fstream>
#include <iostream>
#include <ctime>
#include <fftw3.h>

[[maybe_unused]] void Medipix::start_frame() {
    image.resize(static_cast<std::vector<unsigned int>::size_type>(n_pixel_x) * n_pixel_y);
    for (auto &pixel: image) {
        pixel = 0;
    }
    max_time = 0.0f;
    real_photons = 0;

    for (auto &event: events) {
        event.clear();
    }
    shutter_open = true;
}

float Medipix::get_min_x() const {
    return -pixel_pitch * float(n_pixel_x) / 2.f;
}

float Medipix::get_max_x() const {
    return pixel_pitch * float(n_pixel_x) / 2.f;
}

float Medipix::get_min_y() const {
    return -pixel_pitch * float(n_pixel_y) / 2.f;
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
    int j = int(position_y / pixel_pitch + float(n_pixel_y) / 2.f - 0.5f);
    return std::make_pair(i, j);
}

std::pair<float, float> Medipix::get_pixel_center(unsigned int i, unsigned int j) const {
    float x = pixel_pitch * (float(i) - float(n_pixel_x) / 2.f + 0.5f);
    float y = pixel_pitch * (float(j) - float(n_pixel_y) / 2.f + 0.5f);
    return std::make_pair(x, y);
}

float
Medipix::calculate_shared_energy(float x, float y, float energy, float pixel_center_x, float pixel_center_y) const {
    float x_b = (pixel_center_x + pixel_pitch / 2.f);
    float x_a = (pixel_center_x - pixel_pitch / 2.f);
    float y_b = (pixel_center_y + pixel_pitch / 2.f);
    float y_a = (pixel_center_y - pixel_pitch / 2.f);

    float x_component = std::erf((x_b - x) / (psf_sigma * float(std::numbers::sqrt2))) -
                        std::erf((x_a - x) / (psf_sigma * float(std::numbers::sqrt2)));
    float y_component = std::erf((y_b - y) / (psf_sigma * float(std::numbers::sqrt2))) -
                        std::erf((y_a - y) / (psf_sigma * float(std::numbers::sqrt2)));

    return 0.25f * energy * x_component * y_component;

}


void Medipix::increase_counter(unsigned int x, unsigned int y) {
    std::lock_guard<std::mutex> lk(image_write_mutex);
    image[x * n_pixel_y + y] += 1;
}

unsigned int Medipix::get_total_counts() {
    std::lock_guard<std::mutex> lk(image_write_mutex);

    unsigned int counts = 0;
    for (unsigned int i = 0; i < n_pixel_y * n_pixel_x; ++i) {
        counts += image[i];
    }
    return counts;
}

void Medipix::save_image(const std::string &filename) {
    std::ofstream image_file(filename, std::ios::out | std::ios::binary);
    std::lock_guard<std::mutex> lk(image_write_mutex);
    for (unsigned int i = 0; i < n_pixel_y * n_pixel_x; ++i) {
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


void Medipix::random_threshold_dispersion(float sigma) {
    std::normal_distribution<float> distribution(0.f, sigma);
    unsigned int seed = time(nullptr);
    std::default_random_engine generator(seed);

#pragma omp parallel for default(none) shared(generator, distribution)
    for (unsigned int i = 0; i < n_pixel_y * n_pixel_x; ++i) {
        th0_dispersion[i] = distribution(generator);
    }
}

Medipix::Medipix(bool timed, unsigned int nx, unsigned int ny) : timed(timed), n_pixel_x(nx), n_pixel_y(ny) {
    th0_dispersion.resize(static_cast<std::vector<float>::size_type>(n_pixel_x) * n_pixel_y);
    for (auto &pixel: th0_dispersion) {
        pixel = 0.0f;
    }

    if (timed) {
        events.resize(static_cast<std::vector<Event>::size_type>(n_pixel_x) * n_pixel_y);
        for (auto &pixel_events: events) {
            pixel_events.clear();
        }
    }
}

void Medipix::finish_frame() {
    if (timed) {
        build_i_krum_response(i_krum);
    }
    shutter_open = false;
}

void Medipix::add_photon(float energy, float position_x, float position_y, int radius, float time) {
    if (!shutter_open)
        throw std::logic_error("Shutter not open. Call start_frame() before.");
    if (timed)
        max_time = std::max(max_time, time);
    std::lock_guard<std::mutex> lk(image_write_mutex);
    real_photons++;
}

void Medipix::build_i_krum_response() {
    // We sample the response function at 100 points per us
    float max_resp_time = 2.f;
    unsigned int n_response_points = int(max_resp_time * float(samples_per_us));
    response_function.clear();
    response_function.resize(n_response_points, 0.f);

    // Extremely rough estimation for IKrum 20 setting
    // # TODO: make this configurable
    for (unsigned int i = 0; i < samples_per_us * 1; ++i) {
        response_function[i] = 1.f - float(i) / float(samples_per_us * 1);
    }
    for (unsigned int i = samples_per_us * 1; i < n_response_points; ++i) {
        response_function[i] = 0.f;
    }
}

std::vector<float> Medipix::calculate_pixel_signal(unsigned int i, unsigned int j) {
    auto pixel_events = events[i * n_pixel_y + j];
    std::vector<float> pixel_signal(int(max_time * float(samples_per_us)) + response_function.size(), 0.f);

    for (auto event: pixel_events) {
        unsigned int start_index = int(event.time * float(samples_per_us));
        for (unsigned int index = 0; index < response_function.size(); ++index) {
            if (start_index + index < pixel_signal.size()) {
                pixel_signal[start_index + index] += event.energy * response_function[index];
            }
        }
    }
    return pixel_signal;
}

unsigned int Medipix::get_num_pixels_x() const {
    return n_pixel_x;
}

unsigned int Medipix::get_num_pixels_y() const {
    return n_pixel_y;
}

unsigned int Medipix::get_real_photons() const {
    return real_photons;
}

void Medipix::save_pixel_signals(const std::string &filename, unsigned int i, unsigned int j) {
    if (shutter_open)
        throw std::logic_error("Shutter open. Call finish_frame() before.");
    auto pixel_signal = calculate_pixel_signal(i, j);
    std::ofstream signal_file(filename, std::ios::out | std::ios::binary);
    for (auto &value: pixel_signal) {
        signal_file.write((char *) &value, sizeof(value));
    }
    signal_file.close();
}

bool Medipix::get_timed() const {
    return timed;
}

unsigned int Medipix::get_pixel_value(unsigned int i, unsigned int j) const {
    if (shutter_open)
        throw std::logic_error("Shutter open. Call finish_frame() before.");
    return image[i * n_pixel_y + j];
}

bool Medipix::get_shutter_open() const {
    return shutter_open;
}

std::vector<float> Medipix::get_fourier_spectrum() {
    if (shutter_open)
        throw std::logic_error("Shutter open. Call finish_frame() before.");
    std::vector<fftw_complex> fourier_spectrum((size_t) n_pixel_x * n_pixel_y);
    std::vector<double> image_float((size_t) n_pixel_x * n_pixel_y);
    for (unsigned int k = 0; k < n_pixel_x * n_pixel_y; ++k) {
        image_float[k] = double(image[k]);
    }
    auto plan = fftw_plan_dft_r2c_2d(int(n_pixel_x), int(n_pixel_y), image_float.data(), fourier_spectrum.data(),
                                     FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    unsigned int n_k = std::min(n_pixel_x, n_pixel_y) / 2;
    std::vector<float> spectrum(n_k, 0.f);
    std::vector<unsigned int> spectrum_count(n_k, 0);
    float center[2] = {static_cast<float>(n_pixel_x / 2.0), static_cast<float>((n_pixel_y) / 2.0)};
    for (int i = 0; i < n_pixel_x; ++i) {
        for (int j = 0; j < n_pixel_y; ++j) {
            double x_distance, y_distance;
            if (i < int(center[0])) {
                x_distance = double(i);
            } else {
                x_distance = double(i) - center[0];
            }
            if (j < int(center[1])) {
                y_distance = double(j);
            } else {
                y_distance = double(j) - center[1];
            }
            x_distance = x_distance / n_pixel_x * n_k * 2;
            y_distance = y_distance / n_pixel_y * n_k * 2;
            unsigned int r = int(std::sqrt(x_distance * x_distance + y_distance * y_distance));
            if (r > n_k)
                continue;
            spectrum[r] += float(std::sqrt(fourier_spectrum[i * n_pixel_y + j][0] * fourier_spectrum[i * n_pixel_y + j][0] +
                                     fourier_spectrum[i * n_pixel_y + j][1] * fourier_spectrum[i * n_pixel_y + j][1]));
            spectrum_count[r] += 1;
        }
    }
    for (unsigned int i = 0; i < n_k; ++i) {
        if (spectrum_count[i] > 0) {
            spectrum[i] /= float(spectrum_count[i]);
        }
    }

    return spectrum;
}

void Medipix::save_fourier_spectrum(const std::string &filename) {
    if (shutter_open)
        throw std::logic_error("Shutter open. Call finish_frame() before.");
    auto spectrum = get_fourier_spectrum();
    std::ofstream image_file(filename, std::ios::out | std::ios::binary);

    for (auto &p: spectrum) {
        float pixel_value = p;
        image_file.write((char *) &pixel_value, sizeof(pixel_value));
    }
    image_file.close();
}


void Medipix::build_i_krum_response(int _i_krum) {
    if (_i_krum < 1 or _i_krum > 100)
        throw std::invalid_argument("i_krum must be between 1 and 100.");

    // Fitted values from the paper https://iopscience.iop.org/article/10.1088/1748-0221/10/01/C01047/
    // with the response function "response" below.
    std::vector<int> i_krum_model = {1, 10, 20, 50, 70, 100};
    std::vector<float> wn = {4.89624784, 5.86311405, 8.00195927, 16.50273734, 15.77001966, 17.26883039};
    std::vector<float> wd = {14.90772755, 11.88862317, 43.40896533, 7.22433423, 12.83117555, 19.05247813};
    std::vector<float> d = {2.9689961, 2.65441403, 1.54012187, 0.69959799, 0.66982567, 0.53208576};

    int index_l=0, index_r=0;
    for(int i = 0; i < (i_krum_model.size()-1); ++i){
        if(_i_krum >= i_krum_model[i] && _i_krum < i_krum_model[i + 1]){
            index_l = i;
            index_r = i+1;
            break;
        }
    }
    float t = float(_i_krum - i_krum_model[index_l]) / float(i_krum_model[index_r] - i_krum_model[index_l]);
    float wn_i = wn[index_l] + t * (wn[index_r] - wn[index_l]);
    float wd_i = wd[index_l] + t * (wd[index_r] - wd[index_l]);
    float d_i = d[index_l] + t * (d[index_r] - d[index_l]);

    auto response = [](float x_, float wn_, float wd_, float d_)
    {
        float y= 0;
        // https://www.tutorialspoint.com/control_systems/control_systems_response_second_order.htm
        if(0 < d_ && d_ < 1) {
            y = (wn_ * expf(-d_ * wn_ * x_)) / sqrtf(1 - powf(d_, 2)) * sinf(wd_ * x_);
        }
        if(d_ > 1) {
            y = (wn_ / (2 * sqrtf(powf(d_, 2) - 1))) * (
                    expf(-(d_ * wn_ - wn_ * sqrtf(powf(d_, 2) - 1)) * x_) - expf(-(d_ * wn_ + wn_ * sqrtf(powf(d_, 2) - 1)) * x_));
        }
        if(d_ == 1) {
            y = powf(wn_, 2) * x_ * expf(-wn_ * x_);
        }
        if(d_ <= 0) {
            y = wn_ * sinf(wn_ * x_);
        }
        return y;
    };

    float max_resp_time = 5.f;
    unsigned int n_response_points = int(max_resp_time * float(samples_per_us));
    response_function.clear();
    response_function.resize(n_response_points, 0.f);
    float max_response = 0.f;
    for(unsigned int i = 0; i < n_response_points; ++i){
        float x = float(i) / float(samples_per_us);
        response_function[i] = response(x, wn_i, wd_i, d_i);
        if(response_function[i] > max_response)
            max_response = response_function[i];
    }
    for(auto &r: response_function)
        r /= max_response;
}

int Medipix::get_i_krum() const {
    return i_krum;
}

void Medipix::set_i_krum(int value) {
    if(value < 1 or value > 100)
        throw std::invalid_argument("i_krum must be between 1 and 100.");
    i_krum = value;
}
