# target_link_libraries_system
#
# This function is similar to target_link_libraries but allows the includes
# determined from the library to be added as system includes to suppress
# warnings generated from those header files
#
# https://stackoverflow.com/questions/52135983/cmake-target-link-libraries-include-as-system-to-suppress-compiler-warnings/52136398#52136398
#
function(target_link_libraries_system target visibility)
  set(libs ${ARGN})
  foreach(lib ${libs})
    get_target_property(lib_include_dirs ${lib} INTERFACE_INCLUDE_DIRECTORIES)
    target_include_directories(${target} SYSTEM ${visibility} ${lib_include_dirs})
    target_link_libraries(${target} ${visibility} ${lib})
  endforeach(lib)
endfunction(target_link_libraries_system)

macro(init_code_checks)
  if(PELE_ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES "clang-tidy")
    if(CLANG_TIDY_EXE)
      message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
    else()
      message(WARNING "clang-tidy not found.")
    endif()
  endif()

  if(PELE_ENABLE_CPPCHECK)
    find_program(CPPCHECK_EXE NAMES "cppcheck")
    if(CPPCHECK_EXE)
      message(STATUS "cppcheck found: ${CPPCHECK_EXE}")
      include(ProcessorCount)
      ProcessorCount(NP)
      if(NP EQUAL 0)
        set(NP 1)
      endif()
      if("${CMAKE_SYSTEM_NAME}" MATCHES "Darwin")
        set(NP 1)
      endif()
      add_custom_target(cppcheck ALL
          COMMAND ${CMAKE_COMMAND} -E echo "Running cppcheck on project using ${NP} cores..."
          COMMAND ${CMAKE_COMMAND} -E make_directory cppcheck
          # cppcheck ignores -isystem directories, so we change them to regular -I include directories (with no spaces either)
          COMMAND sed "s/isystem /I/g" compile_commands.json > cppcheck/cppcheck_compile_commands.json
          COMMAND ${CPPCHECK_EXE} --version
          COMMAND ${CPPCHECK_EXE} --template=gcc --inline-suppr --suppress=unusedFunction --suppress=useStlAlgorithm --suppress=missingIncludeSystem --std=c++20 --language=c++ --enable=all --project=cppcheck/cppcheck_compile_commands.json --output-file=cppcheck/cppcheck-full-report.txt -j ${NP}
          COMMAND egrep "information:|error:|performance:|portability:|style:|warning:" cppcheck/cppcheck-full-report.txt | egrep -v "Submodules/AMReX|Submodules/sundials|Submodules/PelePhysics/Mechanisms" | sort | uniq > cppcheck/cppcheck-report.txt
          COMMAND wc -l cppcheck/cppcheck-report.txt
          COMMENT "Run cppcheck on project compile_commands.json"
          BYPRODUCTS cppcheck/cppcheck-report.txt
          WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
          VERBATIM USES_TERMINAL
      )
    else()
      message(WARNING "cppcheck not found.")
    endif()
  endif()
endmacro(init_code_checks)
