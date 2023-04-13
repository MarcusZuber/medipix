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

#ifndef MEDIPIX_MEDIPIX_CSM_H
#define MEDIPIX_MEDIPIX_CSM_H



#include "Medipix.h"

class MedipixCSM: public Medipix {
public:
    /**
     * Constructor for the Medipix with charge summing.
     * @param timed If true, pile-up events are simulated. If false, the simulation is significantly faster.
     * @param nx Number of pixels in x direction
     * @param ny Number of pixels in y direction
     *
     */
    explicit MedipixCSM(bool timed, unsigned int nx=256, unsigned int ny=256);

    /**
     * Constructor for the Medipix with charge summing.
     */
    MedipixCSM();

    /**
     * Adds a photon to the current frame.
     * @param energy in keV
     * @param position_x in um
     * @param position_y in um
     * @param radius area in *pixel* in which the charge distribution is calculated
     * @param time interaction time in us. Only relevant for timed mode.
     */
    void add_photon(float energy, float position_x, float position_y, int radius, float time) override;

    /**
     * Finishes the current frame. In timed mode here the pile-up events are processed.
     * This can take a while.
     */
    void finish_frame() override;

    /**
     * Sets a normal distribution for the threshold dispersion.
     * @param sigma in keV
     */
    void random_threshold_dispersion(float sigma) override;

    /**
     * Getter for the th1
     * @return in keV
     */
    [[maybe_unused]] float get_th1();

    /**
     * Setter for the th1
     * @param t threshold in keV
     */
    void set_th1(float t);
protected:
    /**
     * Getter for the pixel wise th1 value (including the threshold dispersion)
     * @param i pixel index in x direction
     * @param j piexel index in y direction
     * @return in keV
     */
    float get_th1(unsigned int i, unsigned int j);

    /**
     * Value of the threshold 1 in keV
     */
    float th1 = 6.0;

    /**
     * Vector holding the threshold displacement for each pixel
     */
    std::vector<float> th1_dispersion;

};

#endif //MEDIPIX_MEDIPIX_CSM_H
