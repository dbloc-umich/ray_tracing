#include <chrono>
#include <random>
#include <vector>
#include <type_traits>

#include "Box.h"
#include "Direction.h"
#include "kdTree.h"
#include "RayTracing.h"
#include "Sphere.h"

int main(){

    std::vector<std::unique_ptr<Shape>> ptrs;
    ptrs.push_back(std::make_unique<Sphere>(0, 0, 0, 2));
    ptrs.push_back(std::make_unique<Sphere>(4, 4, 4, 2));
    ptrs.push_back(std::make_unique<Sphere>(8, 8, 8, 2));
    
    kdTree tree(ptrs);
    tree.insert(std::make_unique<Sphere>(0, 4, 0, 2));
    std::cout << tree << std::endl;

}