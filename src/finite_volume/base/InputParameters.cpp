// #include "InputParameters.h"
// #include <exception>

// void InputParameters::addScalar(const std::string& key, double value){ _scalars[key] = value; }

// double InputParameters::getRequiredScalar(const std::string& key) const{
//     if (_scalars.count(key) == 0) throw std::out_of_range("ERROR: Required parameter not found");
//     return _scalars.at(key);
// }

// double InputParameters::getScalar(const std::string& key, double def) const noexcept{
//     if (_scalars.count(key) == 0) return def;
//     return _scalars.at(key);
// }

// void InputParameters::addVector(const std::string& key, const Eigen::Vector3d& value){ _vectors[key] = value; }

// Eigen::Vector3d InputParameters::getRequiredVector(const std::string& key) const{
//     if (_vectors.count(key) == 0) throw std::out_of_range("ERROR: Required parameter not found");
//     return _vectors.at(key);
// }

// Eigen::Vector3d InputParameters::getVector(const std::string& key, const Eigen::Vector3d& def) const noexcept{
//     if (_vectors.count(key) == 0) return def;
//     return _vectors.at(key);
// }

// void InputParameters::addNLVariable(const std::string& key, Eigen::Index id){
//     _nlVariables[key] = id;
// }

// Eigen::Index InputParameters::getNLVariableID(const std::string& key) const{
//     if (_nlVariables.count(key) == 0) throw std::out_of_range("ERROR: Required parameter not found");
//     return _nlVariables.at(key);    
// }

// void InputParameters::addAuxVariable(const std::string& key){ _auxVariables.emplace(key); }
