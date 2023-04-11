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
    explicit MedipixCSM(bool timed, unsigned int nx=256, unsigned int ny=256);
    MedipixCSM();
    void add_photon(float energy, float position_x, float position_y, int radius, float time) override;
    void finish_frame() override;
    void random_threshold_dispersion(float sigma) override;
    float get_th1();
    void set_th1(float t);
protected:
    float get_th1(unsigned int i, unsigned int j);
    float th1 = 6.0;
    std::vector<float> th1_dispersion;

};

#endif //MEDIPIX_MEDIPIX_CSM_H
