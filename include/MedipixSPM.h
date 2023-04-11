//
// Created by marcus on 10.04.23.
//

#ifndef MEDIPIX_MEDIPIX_SPM_H
#define MEDIPIX_MEDIPIX_SPM_H


#include "Medipix.h"
#include <list>



class MedipixSPM: public Medipix {
public:
    explicit MedipixSPM(bool timed, unsigned int nx=256, unsigned int ny=256);
    MedipixSPM();
    void add_photon(float energy, float position_x, float position_y, int radius, float time) override;
    void finish_frame() override;
};


#endif //MEDIPIX_MEDIPIX_SPM_H
