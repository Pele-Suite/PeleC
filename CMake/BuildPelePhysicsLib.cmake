function(build_pele_physics_lib pele_physics_lib_name)
  if (NOT (TARGET ${pele_physics_lib_name}))
    add_library(${pele_physics_lib_name} OBJECT)

    set(PELE_PHYSICS_SRC_DIR "${CMAKE_SOURCE_DIR}/Submodules/PelePhysics")
    set(PELE_PHYSICS_TRANSPORT_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Transport")
    set(PELE_PHYSICS_EOS_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Eos")
    set(PELE_PHYSICS_MECHANISM_DIR "${PELE_PHYSICS_SRC_DIR}/Mechanisms/${PELE_PHYSICS_CHEMISTRY_MODEL}")
    set(PELE_PHYSICS_UTILITY_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Utility")
    set(PELE_PHYSICS_REACTIONS_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Reactions")
    set(PELE_PHYSICS_SPRAY_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Spray")
    set(PELE_PHYSICS_SOOT_DIR "${PELE_PHYSICS_SRC_DIR}/Source/Soot")
    set(AMREX_SUNDIALS_DIR "${AMREX_SUBMOD_LOCATION}/Src/Extern/SUNDIALS")

    if(CLANG_TIDY_EXE)
      set_target_properties(${pele_physics_lib_name} PROPERTIES CXX_CLANG_TIDY
                            "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
    endif()

    include(SetPeleCompileFlags)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/TurbInflow/turbinflow.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/TurbInflow/turbinflow.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/TurbInflow)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagBase.H
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagBase.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagConditional.H
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagConditional.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagFilter.H
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagFilter.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagFramePlane.H
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagFramePlane.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagPDF.H
      ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics/DiagPDF.cpp)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/Diagnostics)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/PltFileManager/PltFileManager.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/PltFileManager/PltFileManager.H
      ${PELE_PHYSICS_UTILITY_DIR}/PltFileManager/PltFileManagerBCFill.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/PltFileManager)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/BlackBoxFunction.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/BlackBoxFunction.H
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/BlackBoxFunctionFactory.H
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/NeuralNetHomerolled.H
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/NeuralNetLayerDef.H
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/NeuralNetModelDef.H
      ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction/Table.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/BlackBoxFunction)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/Utilities/Utilities.H
      ${PELE_PHYSICS_UTILITY_DIR}/Utilities/UnitConversions.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/Utilities)

    target_sources(${pele_physics_lib_name}
      PRIVATE
      ${PELE_PHYSICS_UTILITY_DIR}/Filter/Filter.cpp
      ${PELE_PHYSICS_UTILITY_DIR}/Filter/Filter.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_UTILITY_DIR}/Filter)

    target_sources(${pele_physics_lib_name} PRIVATE ${AMREX_SUNDIALS_DIR}/AMReX_Sundials.H
                                             ${AMREX_SUNDIALS_DIR}/AMReX_Sundials_Core.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_Sundials_Core.H
                                             ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_NVector_MultiFab.H
                                             ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.cpp
                                             ${AMREX_SUNDIALS_DIR}/AMReX_SUNMemory.H)
    target_include_directories(${pele_physics_lib_name} SYSTEM PUBLIC ${AMREX_SUNDIALS_DIR})

    target_include_directories(${pele_physics_lib_name} PUBLIC "${PELE_PHYSICS_SRC_DIR}/Source")

    target_sources(${pele_physics_lib_name} PRIVATE
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Transport.cpp
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportParams.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/TransportTypes.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Constant.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Simple.H
                   ${PELE_PHYSICS_TRANSPORT_DIR}/Sutherland.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_TRANSPORT_DIR})
    if("${PELE_PHYSICS_TRANSPORT_MODEL}" STREQUAL "Constant")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_CONSTANT_TRANSPORT)
    endif()
    if("${PELE_PHYSICS_TRANSPORT_MODEL}" STREQUAL "Simple")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_SIMPLE_TRANSPORT)
    endif()
    if("${PELE_PHYSICS_TRANSPORT_MODEL}" STREQUAL "Sutherland")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_SUTHERLAND_TRANSPORT)
    endif()

    target_sources(${pele_physics_lib_name} PRIVATE
                   ${PELE_PHYSICS_EOS_DIR}/EOS.H
                   ${PELE_PHYSICS_EOS_DIR}/GammaLaw.H
                   ${PELE_PHYSICS_EOS_DIR}/Fuego.H
                   ${PELE_PHYSICS_EOS_DIR}/SRK.H)
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_EOS_DIR})
    if("${PELE_PHYSICS_EOS_MODEL}" STREQUAL "GammaLaw")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_GAMMALAW_EOS)
    endif()
    if("${PELE_PHYSICS_EOS_MODEL}" STREQUAL "Fuego")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_FUEGO_EOS)
    endif()
    if("${PELE_PHYSICS_EOS_MODEL}" STREQUAL "Soave-Redlich-Kwong")
      target_compile_definitions(${pele_physics_lib_name} PUBLIC USE_SRK_EOS)
    endif()

    target_sources(${pele_physics_lib_name} PRIVATE
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.cpp
                   ${PELE_PHYSICS_MECHANISM_DIR}/mechanism.H)
    target_include_directories(${pele_physics_lib_name} SYSTEM PUBLIC ${PELE_PHYSICS_MECHANISM_DIR})

    target_sources(${pele_physics_lib_name}
      PRIVATE
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorArkode.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorArkode.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorBase.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorBase.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvode.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvode.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeCustomLinSolver.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeCustomLinSolver.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeJacobian.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeJacobian.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodePreconditioner.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodePreconditioner.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeUtils.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorCvodeUtils.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorNull.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorNull.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorRK64.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorRK64.cpp
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorTypes.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorUtils.H
        ${PELE_PHYSICS_REACTIONS_DIR}/ReactorUtils.cpp
    )
    target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_REACTIONS_DIR})

    if(PELE_PHYSICS_ENABLE_SPRAY)
      target_compile_definitions(${pele_physics_lib_name} PUBLIC PELE_USE_SPRAY)
      target_compile_definitions(${pele_physics_lib_name} PUBLIC SPRAY_FUEL_NUM=${PELE_PHYSICS_SPRAY_FUEL_NUM})
      target_sources(${pele_physics_lib_name} PRIVATE
                     ${PELE_PHYSICS_SPRAY_DIR}/Drag.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayDerive.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayFuelData.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayIO.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayInjection.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayInterpolation.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayJet.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.H
                     ${PELE_PHYSICS_SPRAY_DIR}/SprayParticles.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/SpraySB.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/SpraySetup.cpp
                     ${PELE_PHYSICS_SPRAY_DIR}/WallFunctions.H
                     ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/AhamedSplash.H
                     ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/ReitzKHRT.H
                     ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/SBData.H
                     ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/TABBreakup.H
                     ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash/WallFilm.H
                     ${PELE_PHYSICS_SPRAY_DIR}/Distribution/DistBase.H
                     ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.H
                     ${PELE_PHYSICS_SPRAY_DIR}/Distribution/Distributions.cpp)
      target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR})
      target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/Distribution)
      target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_SPRAY_DIR}/BreakupSplash)
    endif()

    if(PELE_PHYSICS_ENABLE_SOOT)
      set(SOOT_MOMENTS_VALUES 3 6)
      if(NOT PELE_PHYSICS_NUM_SOOT_MOMENTS IN_LIST SOOT_MOMENTS_VALUES)
        message(FATAL_ERROR "NUM_SOOT_MOMENTS must be either 3 or 6")
      endif()
      target_compile_definitions(${pele_physics_lib_name} PUBLIC PELE_USE_SOOT)
      target_compile_definitions(${pele_physics_lib_name} PUBLIC NUM_SOOT_MOMENTS=${PELE_PHYSICS_NUM_SOOT_MOMENTS})
      target_sources(${pele_physics_lib_name} PRIVATE
                     ${PELE_PHYSICS_SOOT_DIR}/SootModel.H
                     ${PELE_PHYSICS_SOOT_DIR}/SootModel.cpp
                     ${PELE_PHYSICS_SOOT_DIR}/SootModel_react.cpp
                     ${PELE_PHYSICS_SOOT_DIR}/SootModel_derive.H
                     ${PELE_PHYSICS_SOOT_DIR}/SootModel_derive.cpp
                     ${PELE_PHYSICS_SOOT_DIR}/Constants_Soot.H
                     ${PELE_PHYSICS_SOOT_DIR}/SootData.H
                     ${PELE_PHYSICS_SOOT_DIR}/SootReactions.H)
      target_include_directories(${pele_physics_lib_name} PUBLIC ${PELE_PHYSICS_SOOT_DIR})
    endif()

    include(AMReXBuildInfo)
    generate_buildinfo(${pele_physics_lib_name} ${CMAKE_SOURCE_DIR})
    target_include_directories(${pele_physics_lib_name} SYSTEM PUBLIC ${AMREX_SUBMOD_LOCATION}/Tools/C_scripts)

    target_link_libraries(${pele_physics_lib_name} PUBLIC sundials_arkode sundials_cvode)

    if(PELE_ENABLE_CUDA)
      target_link_libraries(${pele_physics_lib_name} PUBLIC sundials_nveccuda sundials_sunlinsolcusolversp sundials_sunmatrixcusparse)
    elseif(PELE_ENABLE_HIP)
      target_link_libraries(${pele_physics_lib_name} PUBLIC sundials_nvechip)
    elseif(PELE_ENABLE_SYCL)
      target_link_libraries(${pele_physics_lib_name} PUBLIC sundials_nvecsycl)
    endif()

    if(PELE_ENABLE_MPI)
      target_link_libraries(${pele_physics_lib_name} PUBLIC $<$<BOOL:${MPI_CXX_FOUND}>:MPI::MPI_CXX>)
    endif()

    target_link_libraries(${pele_physics_lib_name} PUBLIC AMReX::amrex)

  endif()
endfunction()
