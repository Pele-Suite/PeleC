function(build_pele_exe pele_exe_name pele_physics_lib_name)

  set(SRC_DIR ${CMAKE_SOURCE_DIR}/Source)

  add_executable(${pele_exe_name} "")

  if(CLANG_TIDY_EXE)
    set_target_properties(${pele_exe_name} PROPERTIES CXX_CLANG_TIDY
                          "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
  endif()

  target_sources(${pele_exe_name}
     PRIVATE
       prob_parm.H
       prob.H
       prob.cpp
  )

  target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_include_directories(${pele_exe_name} PRIVATE ${SRC_DIR})
  target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_BINARY_DIR})
  target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_SOURCE_DIR}/Source/Params/param_includes)
  #Adv and Aux Variables
  if (PELEC_NUM_ADV GREATER 0)
    target_compile_definitions(${pele_exe_name} PRIVATE NUM_ADV=${PELEC_NUM_ADV})
  endif()
  if (PELEC_NUM_AUX GREATER 0)
    target_compile_definitions(${pele_exe_name} PRIVATE NUM_AUX=${PELEC_NUM_AUX})
  endif()

  target_sources(${pele_exe_name}
     PRIVATE
       ${SRC_DIR}/Advance.cpp
       ${SRC_DIR}/BCfill.cpp
       ${SRC_DIR}/Bld.cpp
       ${SRC_DIR}/Constants.H
       ${SRC_DIR}/Derive.H
       ${SRC_DIR}/Derive.cpp
       ${SRC_DIR}/Diffterm.H
       ${SRC_DIR}/Diffusion.H
       ${SRC_DIR}/Diffusion.cpp
       ${SRC_DIR}/EB.H
       ${SRC_DIR}/EB.cpp
       ${SRC_DIR}/EBStencilTypes.H
       ${SRC_DIR}/External.cpp
       ${SRC_DIR}/Forcing.cpp
       ${SRC_DIR}/GradUtil.H
       ${SRC_DIR}/Hydro.H
       ${SRC_DIR}/Hydro.cpp
       ${SRC_DIR}/Geometry.H
       ${SRC_DIR}/Geometry.cpp
       ${SRC_DIR}/Godunov.H
       ${SRC_DIR}/Godunov.cpp
       ${SRC_DIR}/PLM.H
       ${SRC_DIR}/PPM.H
       ${SRC_DIR}/PPM.cpp
       ${SRC_DIR}/InitEB.cpp
       ${SRC_DIR}/IndexDefines.H
       ${SRC_DIR}/IO.H
       ${SRC_DIR}/IO.cpp
       ${SRC_DIR}/LES.H
       ${SRC_DIR}/LES.cpp
       ${SRC_DIR}/MOL.H
       ${SRC_DIR}/MOL.cpp
       ${SRC_DIR}/PeleC.H
       ${SRC_DIR}/PeleC.cpp
       ${SRC_DIR}/PeleCAmr.H
       ${SRC_DIR}/PeleCAmr.cpp
       ${SRC_DIR}/ProblemSpecificFunctions.H
       ${SRC_DIR}/React.cpp
       ${SRC_DIR}/Riemann.H
       ${SRC_DIR}/Setup.cpp
       ${SRC_DIR}/Sources.cpp
       ${SRC_DIR}/SparseData.H
       ${SRC_DIR}/SumIQ.cpp
       ${SRC_DIR}/SumUtils.cpp
       ${SRC_DIR}/Tagging.H
       ${SRC_DIR}/Tagging.cpp
       ${SRC_DIR}/Timestep.H
       ${SRC_DIR}/TransCoeff.H
       ${SRC_DIR}/PeleCUtilities.H
       ${SRC_DIR}/Utilities.cpp
       ${SRC_DIR}/WENO.H
       ${SRC_DIR}/RotSource.cpp
  )

  if(PELE_PHYSICS_ENABLE_SPRAY)
    target_sources(${pele_exe_name} PRIVATE
                   SprayParticlesInitInsert.cpp
                   ${SRC_DIR}/Particle.cpp)
  endif()

  if(PELE_PHYSICS_ENABLE_SOOT)
      target_sources(${pele_exe_name} PRIVATE
                   ${SRC_DIR}/Soot.cpp)
  endif()

  if(PELE_ENABLE_ASCENT)
    target_sources(${pele_exe_name}
      PRIVATE
        ${SRC_DIR}/PeleAscent.H
        ${SRC_DIR}/PeleAscent.cpp
    )
  endif()

  if(PELE_ENABLE_MASA)
    target_compile_definitions(${pele_exe_name} PRIVATE PELE_USE_MASA)
    target_sources(${pele_exe_name} PRIVATE ${SRC_DIR}/MMS.cpp)
    target_link_libraries(${pele_exe_name} PRIVATE MASA::MASA)
    if(PELE_ENABLE_FPE_TRAP_FOR_TESTS)
      set_source_files_properties(${SRC_DIR}/PeleC.cpp PROPERTIES COMPILE_DEFINITIONS PELE_ENABLE_FPE_TRAP)
    endif()
  endif()

  if(NOT "${pele_exe_name}" STREQUAL "${PROJECT_NAME}-UnitTests")
    target_sources(${pele_exe_name}
       PRIVATE
         ${CMAKE_SOURCE_DIR}/Source/main.cpp
    )
  endif()

  if(PELE_ENABLE_CUDA)
    set(pctargets "${pele_exe_name};${pele_physics_lib_name}")
    foreach(tgt IN LISTS pctargets)
      get_target_property(PELE_SOURCES ${tgt} SOURCES)
      list(FILTER PELE_SOURCES INCLUDE REGEX "\\.cpp")
      set_source_files_properties(${PELE_SOURCES} PROPERTIES LANGUAGE CUDA)
    endforeach()
    set_target_properties(${pele_exe_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${pele_exe_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:-Xptxas --disable-optimizer-constants>)
  endif()

  target_link_libraries(${pele_exe_name} PRIVATE ${pele_physics_lib_name} AMReX::amrex)

  install(TARGETS ${pele_exe_name}
          RUNTIME DESTINATION bin
          ARCHIVE DESTINATION lib
          LIBRARY DESTINATION lib)

endfunction()
