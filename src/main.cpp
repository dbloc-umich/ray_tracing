#include "Constants.h"
#include "Direction.h"
#include "kdTree.h"
#include "Octree.h"
#include "RayTracing.h"
#include "Sphere.h"

#include <random>

int main(){

    /*
    // Study the impact of refraction and reflection on a system of varying packing factor
    double L = 10;
    double Rmax = L/10;
    double Rmin = L/20;

    std::default_random_engine rng(3);
    std::uniform_real_distribution<double> xDist(-(L-Rmax), L-Rmax);
    std::uniform_real_distribution<double> yDist(-(L-Rmax), L-Rmax);
    std::uniform_real_distribution<double> zDist(-(L-Rmax), L-Rmax);
    std::uniform_real_distribution<double> rDist(Rmin, Rmax);
    std::vector<Node> nodes;

    const Direction dir(1, 0, 0);

    double V0 = 8*L*L*L;
    std::vector<double> target_pf{0.001, 0.01, 0.1};
    std::vector<double> actual_pf(3);
    std::vector<double> avg_error(3);
    std::vector<double> max_error(3);
    std::vector<Point> argmax_pt(3);
    double Ny = 50;

    for (unsigned i = 0; i < 3; i++){
        double V = 0;
        while (V < V0*target_pf[i]){
            Point r0(xDist(rng), yDist(rng), zDist(rng));
            double r = rDist(rng);

            // Adjust radius to match packing fraction
            Sphere sph = Sphere(r0, r);
            if (sph.volume() + V > V0*target_pf[i]){
                double V_diff = V0*target_pf[i] - V;
                r = std::cbrt(3*V_diff/(4*Constants::PI));
                sph.setRadius(r);
            }

            // Check for collision
            bool overlaps = false;
            for (auto& node: nodes){
                if (node->overlaps(sph)){
                    overlaps = true;
                    break;
                }
            }
            if (!overlaps){
                nodes.emplace_back(std::make_unique<Sphere>(r0, r, 0.3, 1.3));
                V += nodes.back()->volume();
            }
        }
        
        const kdTree tree(nodes);
        actual_pf[i] = V/tree.root()->volume();
        
        // Ray tracing
        double x = tree.root()->xMin() - L;
        double dy = (tree.root()->yMax() - tree.root()->yMin())/Ny;
        double dz = (tree.root()->zMax() - tree.root()->zMin())/Ny;
        double overallI1 = 0;
        double overallI2 = 0;

        for (int iy = 0; iy < Ny+1; iy++){
            double y = tree.root()->yMin() + dy*iy;
            for (int iz = 0; iz < Ny+1; iz++){
                double z = tree.root()->zMin() + dz*iz;
                Point p(x, y, z);
                double I1 = intensity(tree, p, dir, false, false);
                double I2 = intensity(tree, p, dir, true, true);
                overallI1 += I1;
                overallI2 += I2;

                double err = std::fabs(I2-I1)/I2;
                avg_error[i] += err;
                if (err > max_error[i]){
                    max_error[i] = err;
                    argmax_pt[i] = p;
                }
            }
        }
        avg_error[i] /= (Ny*Ny);
        std::cout << "Packing fraction: " << actual_pf[i] << ", ";
        std::cout << "integral error: " << std::fabs(overallI2-overallI1)/overallI2 << ", ";
        std::cout << "average error: " << avg_error[i] << ", ";
        std::cout << "maximum error: " << max_error[i] << " occuring at " << argmax_pt[i] << std::endl;
    }
    */
    return 0;
}