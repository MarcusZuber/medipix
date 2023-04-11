//
// Created by marcus on 10.04.23.
//

#ifndef MEDIPIX_MEDIPIX_CSM_H
#define MEDIPIX_MEDIPIX_CSM_H



#include "Medipix.h"

class MedipixCSM: public Medipix {
public:
    explicit MedipixCSM(bool timed, unsigned int nx=256, unsigned int ny=256);
    MedipixCSM();
    void add_photon(float energy, float position_x, float position_y, int radius, float time) override;
    void finish_frame() override;
};

#endif //MEDIPIX_MEDIPIX_CSM_H
