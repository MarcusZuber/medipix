#include <iostream>
#include "Medipix.h"


int main() {
    Medipix m;
    int num_photons = 10000000;
    m.homogeneous_exposure(30.f, num_photons);
    std::cout << float(m.get_total_counts())/float(num_photons) << std::endl;
    std::cout << float(m.get_total_counts())/(256*256) << std::endl;
    m.save_image("test.raw");

    return 0;
}
