//
// Created by marcus on 10.04.23.
//

#ifndef MEDIPIX_HELPER_H
#define MEDIPIX_HELPER_H

#include <memory>
#include <functional>

class Medipix;

[[maybe_unused]] void homogeneous_exposure(const std::shared_ptr<Medipix>& medipix, float energy, unsigned int number_of_photons, float duration);

[[maybe_unused]] void edge_exposure(const std::shared_ptr<Medipix>& medipix, float energy, float m, float c, unsigned int number_of_photons, float duration);

[[maybe_unused]] [[maybe_unused]] void frequency_exposure(const std::shared_ptr<Medipix>& medipix,float energy, unsigned int number_of_photons, float period, float phase, float n_x, float n_y, float duration);

void exposure(const std::shared_ptr<Medipix>& medipix, float energy, unsigned int number_of_photons, float duration, const std::function<bool (float, float)>& photon_interacting);

bool edge(float x, float y, float m, float c);

bool frequency(float x, float y, float period, float phase, float n_x, float n_y);

#endif //MEDIPIX_HELPER_H
