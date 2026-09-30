# Analysis tooling for the multi-dimensional NSCBC tests
# (Exec/RegTests/NSCBC-COVO, NSCBC-Acoustic, NSCBC-FlameOutflow).
#
#   cmake -S . -B build   -DAMReX_DIR=<amrex-2d-install>/lib/cmake/AMReX
#   cmake -S . -B build3d -DAMReX_DIR=<amrex-3d-install>/lib/cmake/AMReX
#   cmake --build build
#   ./build/fielddump <plotfile> pressure out.dat
#   python3 metrics.py circularity <p_amb> <c> out.dat     # 2-D
#   python3 metrics.py sphericity  <p_amb> <c> out.dat     # 3-D
#
# The dimensionality comes from the AMReX build it is pointed at; fielddump.cpp
# is written against AMREX_SPACEDIM and metrics.py detects 2-D or 3-D from the
# header fielddump writes.  Needs nothing from PeleC or PelePhysics.
cmake_minimum_required(VERSION 3.20)
project(nscbcfields CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(AMReX REQUIRED)
add_executable(fielddump fielddump.cpp)
if(TARGET AMReX::amrex_3d)
  target_link_libraries(fielddump PRIVATE AMReX::amrex_3d)
elseif(TARGET AMReX::amrex_2d)
  target_link_libraries(fielddump PRIVATE AMReX::amrex_2d)
else()
  target_link_libraries(fielddump PRIVATE AMReX::amrex_1d)
endif()
