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
    // {
    //     cout << "Case " << caseNum++ << ": Single sphere, refraction doesn't change direction" << endl;
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
    //     cout << "Case " << caseNum++ << ": Single sphere, refraction changes direction" << endl;
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

    // {
    //     cout << "Case " << caseNum++ << ": Two adjacent spheres, refraction doesn't change direction" << endl;
    //     double r1 = 2.5;
    //     double r2 = 1.5;
    //     double Sigma_t = 0.3;
    //     double refrac = 1.3;
        
    //     std::vector<Node> nodes;
    //     nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, Sigma_t, refrac));
    //     nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, Sigma_t, refrac));
    //     Octree tree(nodes);
    //     Point p(-2*r1, 0, 0);
    //     Direction dir(1, 0, 0);

    //     double A = exp(-2*Sigma_t*(r1+r2));
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

    {
        cout << "Case " << caseNum++ << ": Two adjacent spheres, refraction changes direction" << endl;
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
        double I2 = intensity(tree, p, dir, true, false);
        cout << "Intensity without refraction and reflection: " << I1 << endl;
        cout << "Intensity with refraction and reflection: " << I2 << endl;
        cout << "Error: " << fabs(I2-I1)/I2 << endl;
        cout << endl;
    }

    // {
    //     cout << "Case " << caseNum++ << ": Two adjacent spheres, point travels tangentially to both spheres" << endl;
    //     double r1 = 2.5;
    //     double r2 = 1.5;
    //     double Sigma_t = 0.3;
    //     double refrac = 1.3;
        
    //     std::vector<Node> nodes;
    //     nodes.emplace_back(std::make_unique<Sphere>(0, 0, 0, r1, Sigma_t, refrac));
    //     nodes.emplace_back(std::make_unique<Sphere>(r1+r2, 0, 0, r2, Sigma_t, refrac));
    //     Octree tree(nodes);
    //     Point p(r1, 0.5*r1, 0);
    //     Direction dir(0, -1, 0);

    //     cout << "Analytical expressions:" << endl;
    //     cout << "Intensity without refraction and reflection: " << 1.0 << endl;
    //     cout << "Intensity with refraction and reflection: " << 1.0 << endl;
    //     cout << "Calculated expressions:" << endl;
    //     double I1 = intensity(tree, p, dir, false, false);
    //     double I2 = intensity(tree, p, dir, true, true);
    //     cout << "Intensity without refraction and reflection: " << I1 << endl;
    //     cout << "Intensity with refraction and reflection: " << I2 << endl;
    //     cout << "Error: " << fabs(I2-I1)/I2 << endl;
    //     cout << endl;
    // }

    // {
    //     cout << "Case " << caseNum++ << ": Group of randomly generated spheres" << endl;
    //     std::vector<Node> nodes;
        
    //     double r = 10;
    //     double Sigma_t = 0.3;
    //     double refrac = 1.3;
    //     std::mt19937 xGen(2), yGen(4), zGen(8);
    //     std::default_random_engine rGen(16);
    //     std::uniform_real_distribution<double> xDist(-0.9*r, 0.9*r);
    //     std::uniform_real_distribution<double> rDist(0.05*r, 0.1*r);
    //     auto xBind = [&xDist, &xGen]() { return xDist(xGen); };
    //     auto yBind = [&xDist, &yGen]() { return xDist(yGen); };
    //     auto zBind = [&xDist, &zGen]() { return xDist(zGen); };
    //     auto rBind = [&rDist, &rGen]() { return rDist(rGen); };

    //     for (int i = 0; i < 10; i++){
    //         double x = xBind(), y = yBind(), z = zBind();
    //         r = rBind();
    //         nodes.emplace_back(std::make_unique<Sphere>(x, y, z, r, Sigma_t, refrac));
    //     }

    //     try{
    //         Octree tree(nodes);
    //         std::cout << tree;
    //         Point p(tree.xMin(), 0, 0);
    //         Direction dir(1, 0.5, -0.25);

    //         cout << "Calculated expressions:" << endl;
    //         double I1 = intensity(tree, p, dir, false, false);
    //         double I2 = intensity(tree, p, dir, true, true);
    //         cout << "Intensity without refraction and reflection: " << I1 << endl;
    //         cout << "Intensity with refraction and reflection: " << I2 << endl;
    //         cout << "Error: " << fabs(I2-I1)/I2 << endl;
    //         cout << endl;

    //     } catch (const std::invalid_argument& ex){
    //         cerr << ex.what() << endl;
    //     }
    // }

    return 0;
}