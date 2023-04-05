#include <iostream>
#include "Medipix.h"


int main() {
    Medipix m;
    int num_photons = 1000000;
    m.homogeneous_exposure(30.f, num_photons);
    std::cout << float(m.get_total_counts())/float(num_photons) << std::endl;
    m.save_image("test.raw");

    return 0;
}
