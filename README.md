![Tests](https://github.com/MarcusZuber/medipix/actions/workflows/test.yml/badge.svg)

# Medipix
Simple Monte-Carlo simulation of photon counting detectors from the Medipix family

## Installation

### Requirements
* c++20
* cmake
* openMP
* gnuplot (for crating the plots of the examples)
* fftw3

## Simulation

### Detector properties for the Simulation

* Gaussian shaped point spread function (PSF) is defined by *psf_sigma*. This is mostly relevant for charge sharing.
* Threshold dispersion. This is relevant for pile-up simulation. 
* Preamplifier feedback (so far only a triangular response is implemented)

### Assumptions

* Perfect sensor:
  * Only a gaussian shaped charge distribution is considered.
  * No depth-dependence of the charge collection efficiency.
  * No fluorescence.

* Single pixel mode:
  * Simplified I_Krum response function for pile-up simulation.

* Charge summing mode:
  * No pile-up simulation.
  * Assuming, that charge is only shared in a 2x2 pixel area.


