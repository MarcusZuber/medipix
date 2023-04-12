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

#ifndef MEDIPIX_MEDIPIX_H
#define MEDIPIX_MEDIPIX_H

#include <utility>
#include <memory>
#include <list>
#include <vector>


struct Event {
    /**
     * Time in us.
     */
    float time;

    /**
     * Energy in keV.
     */
    float energy;
};

class Medipix {
public:
    explicit Medipix(bool timed = false, unsigned int nx=256, unsigned int ny=256);

    /**
     * Resets the current image
     */
    [[maybe_unused]] void start_frame();

    /**
     * Finishes the current frame.
     * In case of timed mode, the actual events are processed.
     */
    virtual void finish_frame();

    /**
     * Simulates the interaction of a single photon.
     *
     * @param energy Energy of the interacting in keV
     * @param position_x x-component of the position of the interacting photon in um
     * @param position_y y-component of the position of the interacting photon in um
     * @param radius Radius in pixel, in which shared charge could be deposited
     * @param time interaction Time in us
     */
    virtual void add_photon(float energy, float position_x, float position_y, int radius, float time);

    /**
     * Getter for x dimension of the sensor
     */
    [[nodiscard]] float get_min_x() const;

    /**
     * Getter for x dimension of the sensor
     */
    [[nodiscard]] float get_max_x() const;

    /**
     * Getter for y dimension of the sensor
     */
    [[nodiscard]] float get_min_y() const;

    /**
     * Getter for y dimension of the sensor
     */
    [[nodiscard]] float get_max_y() const ;

    /**
     * Getter for the threshold zero
     * @return Energy threshold in keV
     */
    [[maybe_unused]] [[nodiscard]] float get_th0() const;

    /**
     * Getter for pixel pitch
     * @return Pixel pitch in um
     */
    [[maybe_unused]] [[nodiscard]] float get_pixel_pitch() const;

    /**
     * Getter for sigma of gaussian shaped psf
     * @return sigma in um
     */
    [[maybe_unused]] [[nodiscard]] float get_psf_sigma() const;

    /**
     * Setter for sigma of gaussian shaped psf
     * @param s sigma in um
     */
    [[maybe_unused]] [[maybe_unused]] void set_psf_sigma(float s);

    /**
     * Setter for threshold zero
     * @param t Threshold in keV
     */
    [[maybe_unused]] void set_th0(float t);

    /**
     * Gets pixel index of the pixel where (position_x, position_y) is located in
     * @param position_x
     * @param position_y
     * @return pair of unsigned int with the pixel indices
     */
    [[nodiscard]] std::pair<unsigned int, unsigned int> get_pixel_index(float position_x, float position_y) const;

    /**
     * Calculates the pixel center of a pixel with index (i, j)
     * @param i
     * @param j
     * @return pair of floats with the pixel center
     */
    [[nodiscard]] std::pair<float, float> get_pixel_center(unsigned int i, unsigned int j) const;

    /**
     * Saves the current image to a raw file
     * @param filename
     */
    void save_image(const std::string& filename);

    /**
     * Returns the total number of counts in the current image
     * @return
     */
    unsigned int get_total_counts();

    /**
     * Sets a normal distributed threshold dispersion
     *
     * @param sigma in keV
     */
    void random_threshold_dispersion(float sigma);


    [[nodiscard]] unsigned int get_num_pixels_x() const;
    [[nodiscard]] unsigned int get_num_pixels_y() const;
    [[nodiscard]] unsigned int get_real_photons() const;
    void save_pixel_signals(const std::string& filename, unsigned int i, unsigned int j);
protected:
    /**
     * Calculates the amount of energy in a single pixel
     *
     *
     *   \f[
     *   \frac{E}{2\pi\sigma} \int\limits_{\mbox{pixel_center_x}-\mbox{pixel_pitch}/2}^{\mbox{pixel_center_x}-\mbox{pixel_pitch}/2}
     *   \exp(\frac{x-x_0}{2\sigma^2}) dx
     *   \int\limits_{\mbox{pixel_center_y}-\mbox{pixel_pitch}/2}^{\mbox{pixel_center_y}-\mbox{pixel_pitch}/2}
     *   \exp(\frac{y-y_0}{2\sigma^2}) dy
     *   \f]
     *
     *   \f[
     *      \frac{E}{4} ((erf(\frac{\mbox{pixel_center_x}-\mbox{pixel_pitch}/2 - x_0}{\sqrt{2}\sigma}))
     *      -(erf(\frac{\mbox{pixel_center_x}+\mbox{pixel_pitch}/2 - x_0}{\sqrt{2}\sigma})))
     *      ((erf(\frac{\mbox{pixel_center_y}-\mbox{pixel_pitch}/2 - y_0}{\sqrt{2}\sigma}))
     *      -(erf(\frac{\mbox{pixel_center_y}+\mbox{pixel_pitch}/2 - y_0}{\sqrt{2}\sigma})))
     *   \f]
     *
     * @param x position of the photon interaction in um
     * @param y position of the photon interaction in um
     * @param energy energy of the photon in keV
     * @param pixel_center_x in um
     * @param pixel_center_y in um
     * @return
     */
    [[nodiscard]] float calculate_shared_energy( float x, float y, float energy, float pixel_center_x, float pixel_center_y) const;

    /**
     * Getter for pixel wise threshold (including threshold dispersion)
     * @param i pixel
     * @param j pixel
     * @return threshold in keV
     */
    float get_th0(unsigned int i, unsigned int j);

    /**
     * Threshold zero in keV
     */
    float th0 = 6.0;

    /**
     * pixel pitch in um
     */
    float pixel_pitch = 55.0;

    /**
     * sigma of gaussian shaped psf in um
     */
    float psf_sigma = 13.0;

    /**
     * x dimension of the sensor in pixel
     */
    unsigned int n_pixel_x;

    /**
     * y dimension of the sensor in pixel
     */
    unsigned int n_pixel_y;

    /**
     * Last timepoint of an interaction in the current exposure in us
     */
    float max_time = 0.f;

    /**
     * Current image
     */
    std::vector<unsigned int> image;

    /**
     * Threshold dispersion map
     */
    std::vector<float> th0_dispersion;

    /**
     * Increase the counter of the pixel (i, j) by one
     * @param x
     * @param y
     */
    void increase_counter(unsigned int x, unsigned int y);

    /**
     * Mutex for image write access
     */
    std::mutex image_write_mutex;

    /**
     * Specifies if the current exposure is timed. Times means that pulse-pileup is handled.
     */
    bool timed = false;

    /**
     * Map of lists storing the pixel-wise events
     */
    std::vector<std::list<Event>> events;

    /**
     * Response function of the preamplifier
     */
    std::vector<float> response_function;

    /**
     * Calculates the response function of the preamplifier
     */
    void build_i_krum_response();

    /**
     * Samples per microsecond for the preamplifier response and pileup simulation
     */
    unsigned int samples_per_us = 100;

    /**
     * Calculates the signal of a pixel (events convolved with the preamp response).
     * @param i pixel
     * @param j pixel
     * @return vector of floats with the signal
     */
    std::vector<float> calculate_pixel_signal(unsigned int i, unsigned int j);
    unsigned int real_photons = 0;
};


#endif //MEDIPIX_MEDIPIX_H
