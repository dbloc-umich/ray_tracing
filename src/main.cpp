#include <chrono>
#include <random>
#include <stack>
#include <vector>
#include <type_traits>

#include "Direction.h"
#include "Octree.h"
#include "RayTracing.h"
#include "Sphere.h"
using namespace std;

int main(){
    // double r1 = 2.5;
    // double r2 = 1.5;
    // double refrac = 1.3;
    // std::vector<Node> nodes;
    // Point p(-2*r1, r2, 0);
    // Direction dir(1, 0, 0);

    // double Sigma_t = 3.0e-6;
    // double pow = sqrt(sqrt(10));
    // for (unsigned i = 0; i < 27; i++){
    //     nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, Sigma_t, refrac));
    //     nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, Sigma_t, refrac));
    //     Octree tree(nodes);        
    //     double I1 = intensity(tree, p, dir, false, false);
    //     double I2 = intensity(tree, p, dir, true, true);
        
    //     cout << Sigma_t*r1 << ", " << I1 << ", " << I2 << endl;
    //     Sigma_t *= pow;
    // }

    {
        double r1 = 2.5;
        double r2 = 1.5;
        double Sigma_t = 0.3;
        double refrac = 1.3;
        
        std::vector<Node> nodes;
        nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, Sigma_t, refrac));
        nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, Sigma_t, refrac));
        Octree tree(nodes);
        Point p(-2*r1, 0.5*r1, 0);
        Direction dir(1, 0, 0);

        double y = 0.5*r1;
        double s = 2*sqrt(r1*r1-y*y) + 2*sqrt(r2*r2-y*y);
        double A = exp(-Sigma_t*s);
        cout << "Analytical expressions:" << endl;
        cout << "Intensity without refraction and reflection: " << A << endl;
        cout << "Calculated expressions:" << endl;
        double I1 = intensity(tree, p, dir, false, false);
        double I2 = intensity(tree, p, dir, true, true);
        cout << "Intensity without refraction and reflection: " << I1 << endl;
        cout << "Intensity with refraction and reflection: " << I2 << endl;
        cout << "Error: " << fabs(I2-I1)/I2 << endl;
        cout << endl;
    }

    return 0;
}