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

#include "MedipixCSM.h"

MedipixCSM::MedipixCSM(bool timed, unsigned int nx, unsigned ny): Medipix(timed, nx, ny){

}

MedipixCSM::MedipixCSM(): Medipix(false) {

}

void MedipixCSM::add_photon(float energy, float position_x, float position_y, int radius, float time) {
    Medipix::add_photon(energy, position_x, position_y, radius, time);

}

void MedipixCSM::finish_frame() {
    Medipix::finish_frame();
}
