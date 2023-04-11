//
// Created by marcus on 09.04.23.
//
#include <Medipix.h>
#include <helper.h>
#include <iostream>
#include "MedipixSPM.h"
#include "MedipixCSM.h"
#include <fstream>

int main(){
    auto m = std::make_shared<MedipixSPM>(false);
    //m->set_psf_sigma(15.f);
    m->random_threshold_dispersion(1.0f);
    int num_photons = 100000;
    std::ofstream data_file;
    data_file.open("threshold_scan_spm.txt");
    data_file << "# num_photons=" << num_photons << std::endl;
    data_file << "# th0 counts" << std::endl;
    for(float th = 6.0f; th < 70.f; th += 1.0f){
        m->set_th0(th);

        m->start_frame();
        homogeneous_exposure(m, 59.6f, num_photons, 1.0f);
        m->finish_frame();

        data_file << th << " " << (m->get_total_counts()) << std::endl;
    }
    data_file.close();


    auto m2 = std::make_shared<MedipixCSM>(false);
    m2->random_threshold_dispersion(3.0f);
    num_photons = 100000;
    std::ofstream data_file_csm;
    data_file_csm.open("threshold_scan_csm.txt");
    data_file_csm << "# num_photons=" << num_photons << std::endl;
    data_file_csm << "# th0 counts" << std::endl;
    for(float th = 6.0f; th < 40.f; th += 1.0f){
        m2->set_th0(th);

        m2->start_frame();
        homogeneous_exposure(m2, 30.f, num_photons, 1.0f);
        m2->finish_frame();

        data_file_csm << th << " " << (m2->get_total_counts()) << std::endl;
    }
    data_file_csm.close();
}

