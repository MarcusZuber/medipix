//
// Created by marcus on 11.04.23.
//

#include "MedipixSPM.h"
#include "helper.h"

int main() {
    auto m = std::make_shared<MedipixSPM>(true);
    m->set_psf_sigma(15.f);
    m->random_threshold_dispersion(1.0f);
    int num_photons = 10000000;
    m->start_frame();
    frequency_exposure(m, 59.6f, num_photons, 550.0f, 0.0f, 1.0f, 0.0f, 10.f);
    m->finish_frame();
    m->save_image("frequency_pattern_spm.raw");
}