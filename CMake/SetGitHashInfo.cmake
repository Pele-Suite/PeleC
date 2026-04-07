# Extract git hashes from each submodule at configure time
# and generate a header that PeleLMeX C++ code can use directly.

find_package(Git QUIET)

function(get_git_hash _dir _output_var)
  if(GIT_FOUND AND EXISTS "${_dir}")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --always --dirty --long
      WORKING_DIRECTORY "${_dir}"
      OUTPUT_VARIABLE _hash
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET
    )
  endif()
  if(NOT _hash)
    set(_hash "unknown")
  endif()
  set(${_output_var} "${_hash}" PARENT_SCOPE)
endfunction()

# --- Paths to each repository ---
set(PELELMEX_DIR "${CMAKE_SOURCE_DIR}")
set(AMREX_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Submodules/amrex")
set(PELEPHYSICS_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics")
set(SUNDIALS_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Submodules/sundials")

# --- Get hashes ---
get_git_hash("${PELELMEX_DIR}" PELELMEX_GIT_HASH)
get_git_hash("${AMREX_DIR}" AMREX_GIT_HASH)
get_git_hash("${PELEPHYSICS_DIR}" PELEPHYSICS_GIT_HASH)
get_git_hash("${SUNDIALS_DIR}" SUNDIALS_GIT_HASH)

message(STATUS "PeleLMeX    git hash: ${PELELMEX_GIT_HASH}")
message(STATUS "AMReX       git hash: ${AMREX_GIT_HASH}")
message(STATUS "PelePhysics git hash: ${PELEPHYSICS_GIT_HASH}")
message(STATUS "SUNDIALS    git hash: ${SUNDIALS_GIT_HASH}")

# --- Generate the header ---
configure_file(
  "${CMAKE_SOURCE_DIR}/CMake/PeleGitHashes.H.in"
  "${CMAKE_BINARY_DIR}/PeleGitHashes.H"
  @ONLY
)
