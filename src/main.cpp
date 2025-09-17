#include "Direction.h"
#include "Octree.h"
#include "RayTracing.h"
#include "Sphere.h"

#include <random>

using namespace std;

int main(){

    // Study the impact of refraction and reflection on a system of varying packing factor
    double L = 10;
    double Rmax = L/2;
    double Rmin = 0;

    std::default_random_engine rng(3);
    std::uniform_real_distribution<double> xDist(-L, L);
    std::uniform_real_distribution<double> yDist(-L, L);
    std::uniform_real_distribution<double> zDist(-L, L);
    std::uniform_real_distribution<double> rDist(Rmin, Rmax);
    vector<Node> nodes;


    unsigned N = 2;
    double V = 0;
    for (unsigned i = 0; i < N; i++){
        double x = xDist(rng);
        double y = yDist(rng);
        double z = zDist(rng);
        double r = rDist(rng);
        nodes.emplace_back(std::make_unique<Sphere>(x, y, z, r));
        V += nodes.back()->volume();
    }

    Octree tree(nodes);
    cout << "Packing factor = " << V/tree.root()->volume();

    return 0;
}