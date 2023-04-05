//
// Created by marcus on 05.04.2023.
//

#ifndef MEDIPIX_MEDIPIX_H
#define MEDIPIX_MEDIPIX_H

#include <utility>
#include <memory>
// all lentghs in um
// all energies in keV

class Medipix {
public:
    Medipix();
    void reset_image();
    void add_photon(float energy, float position_x, float position_y, int radius);
    [[nodiscard]] float get_min_x() const;
    [[nodiscard]] float get_max_x() const;
    [[nodiscard]] float get_min_y() const;
    [[nodiscard]] float get_max_y() const ;
    [[nodiscard]] float get_th0() const;
    [[nodiscard]] float get_pixel_pitch() const;
    [[nodiscard]] float get_psf_sigma() const;
    [[nodiscard]] std::pair<unsigned int, unsigned int> get_pixel_index(float position_x, float position_y) const;
    [[nodiscard]] std::pair<float, float> get_pixel_center(unsigned int i, unsigned int j) const;
    void save_image(std::string filename) const;
    void homogeneous_exposure(float energy, unsigned int number_of_photons);
    unsigned int get_total_counts();
private:
    float calculate_shared_energy( float x, float y, float energy, float pixel_center_x, float pixel_center_y) const;
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
