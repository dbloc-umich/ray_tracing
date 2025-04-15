#include <chrono>
#include <random>
#include <vector>
#include <type_traits>

#include "Direction.h"
#include "kDTree.h"
#include "RayTracing.h"
#include "Sphere.h"
using namespace std;

int main(){
    const double R = 2;
    const double Sigma_t = 0.5;
    const double refrac = 1.3;

    NodeList nodes;
    for (int i = 0; i < 8; i++){
        nodes.emplace_back(std::make_unique<Sphere>(0, 2*i*R, 0, R, Sigma_t, refrac));
        nodes.emplace_back(std::make_unique<Sphere>(4*R, 2*i*R, 0, R, Sigma_t, refrac));
    }
    kDTree tree(nodes);
    std::cout << tree << std::endl;

    auto test = [&tree](Point p, Direction dir) -> void
    {
        try{
            double s;
            auto node = tree.nextNode(p, dir, nullptr, s);
            if (s == 0.0){
                node = tree.nextNode(p, dir, node, s);
                if (s == 0.0) cerr << "Cannot find the right node." << endl;
            }
            if (node) cout << "p=" << p << ", dir=" << dir << ", s=" << s << ", next node=" << *node << endl;
        } catch(const exception& ex){
            cerr << ex.what() << endl;
        }
    };

    // Test for a point outside of the entire Box
    {
        Point p(-4, 0, 0);
        Direction dir(1, 0, 0);
        test(p, dir);
    }

    // Test for a Point inside of the Box
    {
        Point p(3, 3, 0);
        Direction dir(1, 0, 0);
        test(p, dir);
    }

    return 0;
}