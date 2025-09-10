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
    unsigned int caseNum = 1;
    cout << "Case " << caseNum++ << ": Two adjacent spheres, refraction changes direction" << endl;
    double r1 = 2.5;
    double r2 = 1.5;
    double Sigma_t = 0.3;
    double refrac = 1.3;
    
    std::vector<Node> nodes;
    nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, Sigma_t, refrac));
    nodes.emplace_back(std::make_unique<Sphere>(2*(r1+r2), 0, 0, r2, Sigma_t, refrac));
    Octree tree(nodes);
    
    Direction dir(1, 0, 0);

    unsigned N = 100;
    // vector<double> y(N+1);
    for (unsigned i = 0; i <= N; i++){
        double y = r1 - 2*i*r1/N;
        Point p(-2*r1, y, 0);
        double I1 = intensity(tree, p, dir, false, false);
        double I2 = intensity(tree, p, dir, true, true);
        cout << y << "," << I1 << "," << I2 << endl;
    }

    return 0;
}