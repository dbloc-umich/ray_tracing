# Ray Tracing in Stochastic Droplet Media

Currently under development, this package aims to model laser interaction within a media of randomly distributed liquid droplets through means of ray tracing.

Requires:
[Eigen](github.com/PX4/eigen) for numerical linear algebra.
[CoolProp]() for fluid properties. CoolProp is compiled to a static library using the following command:

cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_C_COMPILER=/ucrt64/bin/gcc.exe \
  -DCMAKE_CXX_COMPILER=/ucrt64/bin/g++.exe \
  -DCOOLPROP_STATIC_LIBRARY=ON \
  -DCOOLPROP_SHARED_LIBRARY=OFF \
  -DCOOLPROP_EES_MODULE=OFF

In CMakeLists.txt, I commented out lines 711 through 724 to avoid bitness maldetection.

[catch2](github.com/catchorg/Catch2) for testing purposes.