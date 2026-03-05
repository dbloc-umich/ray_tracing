// #ifndef INPUT_PARAMETERS_H
// #define INPUT_PARAMETERS_H

// #include <map>
// #include <set>
// #include <string>
// #include "Eigen/Dense"

// class InputParameters{
//     public:
//     template<typename T>
//     using Map = std::map<std::string, T>;
//     using Set = std::set<std::string>;
//     InputParameters(){}

//     void addScalar(const std::string& key, double value);
//     double getRequiredScalar(const std::string& key) const;
//     double getScalar(const std::string& key, double def) const noexcept;

//     void addVector(const std::string& key, const Eigen::Vector3d& value);
//     Eigen::Vector3d getRequiredVector(const std::string& key) const;
//     Eigen::Vector3d getVector(const std::string& key, const Eigen::Vector3d& def) const noexcept;

//     void addNLVariable(const std::string& key, Eigen::Index id);
//     Eigen::Index getNLVariableID(const std::string& key) const;

//     void addAuxVariable(const std::string& key);

//     protected:
//     Map<double> _scalars;
//     Map<Eigen::Vector3d> _vectors;
//     Map<Eigen::Index> _nlVariables;
//     Set _auxVariables;
// };

// #endif