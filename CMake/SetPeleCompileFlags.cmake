list(APPEND PELE_CXX_FLAGS "-Wall" "-Wextra" "-pedantic" "-Wno-unused-function")
if(CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL 7.0)
    list(APPEND PELE_CXX_FLAGS "-faligned-new"
                                "-Wunreachable-code"
                                "-Wnull-dereference"
                                "-Wfloat-conversion"
                                "-Wshadow"
                                "-Woverloaded-virtual")
    if(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")
      list(APPEND PELE_CXX_FLAGS "-Wno-pass-failed" "-Wno-c++17-extensions")
      if((CMAKE_CXX_COMPILER_ID STREQUAL "Clang") AND PELE_ENABLE_FPE_TRAP_FOR_TESTS)
        list(APPEND PELE_CXX_FLAGS "-ffp-exception-behavior=maytrap")
      endif()
    endif()
  endif()
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  list(APPEND PELE_CXX_FLAGS "-diag-disable:11074,11076,10397,15335")
endif()
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
  list(APPEND PELE_CXX_FLAGS "-Wno-sign-compare"
                              "-Wno-unused-parameter"
                              "-Wno-unused-variable")
endif()

list(APPEND PELE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
separate_arguments(PELE_CXX_FLAGS)
target_compile_options(${pele_physics_lib_name} PUBLIC $<$<COMPILE_LANGUAGE:CXX>:${PELE_CXX_FLAGS}>)

# Use flags to avoid warnings from certain files
if((NOT PELE_ENABLE_CUDA) AND (CMAKE_CXX_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$"))
  set(NO_WARN_CXX_FLAGS "-w")
endif()

set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_Sundials.H PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_Sundials.cpp PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.cpp PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.H PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.cpp PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
set_source_files_properties(${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.H PROPERTIES COMPILE_OPTIONS "${NO_WARN_CXX_FLAGS}")
