#include "ElectronTemperatureAux.h"
#include "Constants.h"
#include <iostream>

double ElectronTemperatureAux::computeValue(const std::map<std::string, double>& u) const{
    if (u.count("electron_energy") == 0 || u.count("ion_number_density") == 0)
        throw std::out_of_range("ERROR: Required variable(s) not found.");
    return 2.0/3 * u.at("electron_energy") / (u.at("ion_number_density") * pconst::k_B);
}