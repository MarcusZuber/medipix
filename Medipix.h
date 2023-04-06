//
// Created by marcus on 05.04.2023.
//

#ifndef MEDIPIX_MEDIPIX_H
#define MEDIPIX_MEDIPIX_H

#include <utility>
#include <memory>

class Medipix {
public:
    Medipix();

    /**
     * Resets the current image
     */
    [[maybe_unused]] void reset_image();

    /**
     * Simulates the interaction of a single photon.
     *
     * @param energy Energy of the interacting in keV
     * @param position_x x-component of the position of the interacting photon in um
     * @param position_y y-component of the position of the interacting photon in um
     * @param radius Radius in pixel, in which shared charge could be deposited
     */
    void add_photon(float energy, float position_x, float position_y, int radius);

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
     * Simulates a homogeneous exposure of *number_of_photons* with the energy *energy*
     * @param energy Energy of the photons in keV
     * @param number_of_photons
     */
    void homogeneous_exposure(float energy, unsigned int number_of_photons);

    /**
     * Returns the total number of counts in the current image
     * @return
     */
    unsigned int get_total_counts();
private:
    /**
     * Calculates the amount of energy in a single pixel
     *
     *
     *   \f[
     *   \frac{E}{2\pi\sigma} \int\limits_{\text{pixel_center_x}-\text{pixel_pitch}/2}^{\text{pixel_center_x}-\text{pixel_pitch}/2}
     *   \exp(\frac{x-x_0}{2\sigma^2}) dx
     *   \int\limits_{\text{pixel_center_y}-\text{pixel_pitch}/2}^{\text{pixel_center_y}-\text{pixel_pitch}/2}
     *   \exp(\frac{y-y_0}{2\sigma^2}) dy
     *   \f]
     *
     *   \f[
     *      \frac{E}{4} ((erf(\frac{\text{pixel_center_x}-\text{pixel_pitch}/2 - x_0}{\sqrt{2}\sigma}))
     *      -(erf(\frac{\text{pixel_center_x}+\text{pixel_pitch}/2 - x_0}{\sqrt{2}\sigma})))
     *      ((erf(\frac{\text{pixel_center_y}-\text{pixel_pitch}/2 - y_0}{\sqrt{2}\sigma}))
     *      -(erf(\frac{\text{pixel_center_y}+\text{pixel_pitch}/2 - y_0}{\sqrt{2}\sigma})))
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


    float th0 = 6.0;
    float pixel_pitch = 55.0;
    float psf_sigma = 15.0;
    unsigned int n_pixel_x = 256;
    unsigned int n_pixel_y = 256;
    std::shared_ptr<unsigned int[]> image;
    void increase_counter(unsigned int x, unsigned int y);
    std::mutex image_write_mutex;
};


#endif //MEDIPIX_MEDIPIX_H
