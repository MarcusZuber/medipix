//
// Created by marcus on 10.04.23.
//

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
