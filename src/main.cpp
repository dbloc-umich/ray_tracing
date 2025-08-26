#include <chrono>
#include <random>
#include <vector>
#include <type_traits>

#include "Direction.h"
#include "Octree.h"
#include "RayTracing.h"
#include "Sphere.h"
using namespace std;

int main(){
    // {
    //     cout << "Case 1: Single sphere, refraction doesn't change direction" << endl;
    //     double r = 2;
    //     double Sigma_t = 0.3;
    //     double refrac = 1.3;
    //     std::vector<Node> nodes;
    //     nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, refrac));
    //     Octree tree(nodes);

    //     Point p(-2*r, 0, 0);
    //     Direction dir(1, 0, 0);

    //     double A = exp(-2*Sigma_t*r);
    //     double R = (1-refrac)/(1+refrac);
    //     R *= R;
    //     double T = 1-R;

    //     cout << "Analytical expressions:" << endl;
    //     cout << "Intensity without refraction and reflection: " << A << endl;
    //     cout << "Intensity with refraction and reflection: " << R + T*T*A*(1+R*A)/(1-R*R*A*A) << endl;
    //     cout << "Calculated expressions:" << endl;
    //     double I1 = intensity(tree, p, dir, false, false);
    //     double I2 = intensity(tree, p, dir, true, true);
    //     cout << "Intensity without refraction and reflection: " << I1 << endl;
    //     cout << "Intensity with refraction and reflection: " << I2 << endl;
    //     cout << "Error: " << fabs(I2-I1)/I2 << endl;
    //     cout << endl;
    // }

    // {
    //     cout << "Case 2: Single sphere, refraction changes direction" << endl;
    //     double r = 2;
    //     double Sigma_t = 0.3;
    //     double refrac = 1.3;
    //     std::vector<Node> nodes;
    //     nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r, Sigma_t, refrac));
    //     Octree tree(nodes);

    //     Point p(-2*r, 0.5*r, 0);
    //     Direction dir(1, 0, 0);
        
    //     double A = exp(-sqrt(3)*Sigma_t*r);
    //     cout << "Analytical expressions:" << endl;
    //     cout << "Intensity without refraction and reflection: " << A << endl;
    //     cout << "Calculated expressions:" << endl;
    //     double I1 = intensity(tree, p, dir, false, false);
    //     double I2 = intensity(tree, p, dir, true, true);
    //     cout << "Intensity without refraction and reflection: " << I1 << endl;
    //     cout << "Intensity with refraction and reflection: " << I2 << endl;
    //     cout << "Error: " << fabs(I2-I1)/I2 << endl;
    //     cout << endl;       
    // }

    {
        cout << "Case 3: Group of randomly generated spheres" << endl;
        std::vector<Node> nodes;
        
        double r = 10;
        double Sigma_t = 0.3;
        double refrac = 1.3;
        std::mt19937 xGen(2), yGen(4), zGen(8);
        std::default_random_engine rGen(16);
        std::uniform_real_distribution<double> xDist(-r, r);
        std::uniform_real_distribution<double> rDist(0, 0.1*r);
        auto xBind = [&xDist, &xGen]() { return xDist(xGen); };
        auto yBind = [&xDist, &yGen]() { return xDist(yGen); };
        auto zBind = [&xDist, &zGen]() { return xDist(zGen); };
        auto rBind = [&rDist, &rGen]() { return rDist(rGen); };

        for (int i = 0; i < 10; i++){
            double x = xBind(), y = yBind(), z = zBind();
            r = rBind();
            nodes.emplace_back(std::make_unique<Sphere>(x, y, z, r, Sigma_t, refrac));
        }

        try{
            Octree tree(nodes);
            cout << tree << endl;
            Point p1(tree.xMin(), tree.yMin(), (tree.zMin()+tree.zMax())/2);
            Point p2(-3.58927, 1.95112, -1.55298);
            Direction dir(p1, p2);

            double s;
            auto nextNode = tree.nextNode(p1, dir, nullptr, s);
            if (nextNode) cout << *nextNode << endl;
            else cout << "nullptr" << endl;

            // cout << "Calculated expressions:" << endl;
            // double I1 = intensity(tree, p, dir, false, false);
            // double I2 = intensity(tree, p, dir, true, true);
            // cout << "Intensity without refraction and reflection: " << I1 << endl;
            // cout << "Intensity with refraction and reflection: " << I2 << endl;
            // cout << "Error: " << fabs(I2-I1)/I2 << endl;
            // cout << endl;

        } catch (const std::invalid_argument& ex){
            cerr << ex.what() << endl;
        }
    }

    // // Miscellaneous testing
    // {
    //     BoundingBox bbox(-8.07637, -7.48154, -9.88248, 9.45813, 8.1869, 8.5984);
    //     Sphere sph(-3.58927, 1.95112, -1.55298, sqrt(0.91219));
    //     cout << bbox.octant(sph) << endl;
    // }

    return 0;
}