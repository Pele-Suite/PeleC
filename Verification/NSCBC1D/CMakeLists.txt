# Standalone 1-D verification driver for Source/NSCBC.H.
#
# Deliberately independent of the PeleC application build: it compiles the
# production boundary-condition header against a minimal 1-D solver so that a
# failure can be attributed to the boundary condition rather than to the AMReX
# plumbing around it.  See the header comment in nscbc1d.cpp.
#
#   cmake -S . -B build -DAMReX_DIR=<amrex-install>/lib/cmake/AMReX \
#         -DPELEC_DIR=../.. -DPELE_MECHANISM=air
#   cmake --build build && ./build/nscbc1d sweep
cmake_minimum_required(VERSION 3.20)
project(nscbc1d CXX)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PELEC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../.." CACHE PATH "PeleC root")
set(PELE_MECHANISM "air" CACHE STRING "PelePhysics mechanism")
# Fuego (ideal gas, the default) or SRK (real gas): the kernel restricts
# itself to density-carrying EOS entry points precisely so that it compiles
# and runs under SRK, and this switch is how that claim is checked.
set(PELE_EOS "Fuego" CACHE STRING "EOS model: Fuego or SRK")
# Passive advected scalars ride the lambda_0 family through pack_ghost;
# a build with PELE_NUM_ADV > 0 is how that path is compiled and checked.
set(PELE_NUM_ADV 0 CACHE STRING "Number of passive advected scalars")
find_package(AMReX REQUIRED)

set(PP "${PELEC_DIR}/Submodules/PelePhysics")
add_executable(nscbc1d
  nscbc1d.cpp
  ${PP}/Mechanisms/${PELE_MECHANISM}/mechanism.cpp)
target_include_directories(nscbc1d PRIVATE
  ${PELEC_DIR}/Source
  ${PP}/Source
  ${PP}/Source/Eos
  ${PP}/Source/Transport
  ${PP}/Source/Utility
  ${PP}/Source/Utility/Utilities
  ${PP}/Source/Utility/Filter
  ${PP}/Source/Utility/BlackBoxFunction
  ${PP}/Source/Utility/PMF
  ${PP}/Source/Utility/TurbInflow
  ${PP}/Source/Utility/Diagnostics
  ${PP}/Mechanisms/${PELE_MECHANISM})
if (PELE_EOS STREQUAL "SRK")
  # SRK forbids constant transport (PelePhysicsConstraints.H); pair it with
  # Simple.  Nothing in the driver reads transport coefficients -- the C12
  # conduction lives in the mini solver -- so the choice only has to compile.
  set(PELE_EOS_DEF USE_SRK_EOS)
  set(PELE_TRANS_DEF USE_SIMPLE_TRANSPORT)
else ()
  set(PELE_EOS_DEF USE_FUEGO_EOS)
  set(PELE_TRANS_DEF USE_CONSTANT_TRANSPORT)
endif ()
target_compile_definitions(nscbc1d PRIVATE
  ${PELE_EOS_DEF} ${PELE_TRANS_DEF} NUM_ADV=${PELE_NUM_ADV} NUM_AUX=0)
target_link_libraries(nscbc1d PRIVATE AMReX::amrex_1d)
