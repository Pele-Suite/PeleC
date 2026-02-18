PeleC: An adaptive mesh refinement solver for compressible reacting flows
-------------------------------------------------------------------------

`Documentation <https://amrex-combustion.github.io/PeleC/>`_ | `Nightly Test Results <https://my.cdash.org/index.php?project=Pele>`_ | `PeleC Citation <https://doi.org/10.1177/10943420221121151>`_ | `Pele Citation <https://doi.org/10.1137/1.9781611977967.2>`_

Getting Started
~~~~~~~~~~~~~~~

To compile and run `PeleC`, one needs a C++ compiler that supports the C++17 standard.  A hierarchical strategy for parallelism is supported, based on MPI, MPI + OpenMP, or MPI + GPU (CUDA/HIP/DPC++).  The code should work with all major MPI and OpenMP implementations.  PeleC should build and run with no modifications to the `make` system if using a Linux system with the GNU compilers, version 7 and above. The preferred build strategy for most users is to use GNU Make to build executables for individual cases as needed, from the directory containing the case files. CMake, although used mostly for testing, is also an option for building the code.

The best way to download `PeleC` and its necessary dependencies AMReX, PelePhysics, and SUNDIALS (as git submodules) is through a recursive
git clone. To reduce the download size, the ``--shallow-submodules`` and ``--single-branch`` flags can be added to the clone command to omit
extraneous parts of the git history. By default, you will download the latest version of PeleC. If you'd like to use a specific
`released version <https://github.com/AMReX-Combustion/PeleC/releases/>`_
of PeleC, you can add the option ``--branch=<version>``, e.g. ``--branch=v25.04``, to your clone command. ::

    git clone --recursive https://github.com/AMReX-Combustion/PeleC.git

To run `PeleC` for a sample 3D flame problem::

    cd PeleC/Exec/RegTests/PMF
    make TPLrealclean && make realclean && make TPL && make -j
    ./Pele3d.xxx.yyy.ex example.inp

1. In the exec line above, xxx.yyy is a tag identifying your compiler and various build options, and will vary across pltaform.  (Note that GNU compilers must be at least version 7, and MPI should be at least of standard version 3).

   a. To run with MPI support, ensure ``USE_MPI = TRUE`` in the ``GNUmakefile`` or specify it on the command line::

        make TPLrealclean && make realclean && make TPL USE_MPI=TRUE && make -j USE_MPI=TRUE
        mpiexec -n 4 ./PeleC3d.xxx.MPI.ex example.inp

   b. If the string ``.MPI.`` does not appear in your executable name, you have not successfully built PeleC with MPI support, and if you execute in parallel you will be running multiple instances of the solver in serial.

   c. See the ``GNUmakefile`` for other build options, including the compiler used (``COMP``), and certain model settings that must determined at compile time, such as the solver dimensionality (``DIM``) and the chemical mechanism used (``Chemistry_Model``).

2. The example is a 3D premixed flame, flowing vertically upward through the domain with no gravity. The lateral boundaries are periodic.  A detailed hydrogen model is used.  The solution is initialized with a wrinkled/perturbed 1D steady flame solution computed using Cantera (historically these used the PREMIX code).  Two levels of solution-adaptive refinement are automatically triggered by the presence of the flame intermediate, HO2.

3. In addition to informative output to the terminal, periodic plotfiles are written in the run folder.  These may be viewed with AMReX's `Amrvis <https://amrex-codes.github.io/amrex/docs_html/Visualization.html>`_, `VisIt <https://visit-dav.github.io/visit-website/>`_, or `ParaView <https://www.paraview.org>`_:

   a. In VisIt, direct the File->Open dialogue to select the file named "Header" that is inside each plotfile folder.

   b. In ParaView, navigate to the case directory, open the plotfile folder.

   c. With Amrvis, ``$ amrvis3d plt00030``, for example.

Dependencies
~~~~~~~~~~~~

`PeleC` is built on the `AMReX` and `PelePhysics` libraries. PeleC also requires the `SUNDIALS <https://github.com/LLNL/sundials>`_ ODE solver library.


Development model
~~~~~~~~~~~~~~~~~

To add a new feature to PeleC, the procedure is:

1. Create a branch for the new feature (locally): ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch: ::

    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                      # fix any identified conflicts between local and remote branches of "development"
    git checkout AmazingNewFeature
    git rebase development        # fix any identified conflicts between "development" and "AmazingNewFeature"

3. Build and run

   a. Build and run the full test suite using CMake and CTest (See the ``Build/`` directory for an example script). Please do not introduce warnings. PeleC is checked against ``clang-tidy`` and ``cppcheck`` in the CI. To use ``cppcheck`` and ``clang-tidy`` locally use these CMake options: ::

        -DPELE_ENABLE_CLANG_TIDY:BOOL=ON
        -DPELE_ENABLE_CPPCHECK:BOOL=ON

   b. Run ``clang-tidy`` by using an LLVM compiler and making sure ``clang-tidy`` is found during configure. Then ``make`` will run ``clang-tidy`` along with compilation. Once verifying ``cppcheck`` was found during configure, using the ``make cppcheck`` target should run its checks on the ``compile_commands.json`` database generated by CMake. More information on these checks can be seen in the CI files used for GitHub Actions in the ``.github/workflows/`` directory.

   c. To easily format all source files before commit, use the following command: ::

        find ./Source ./Exec \( -name "*.cpp" -o -name "*.H" -o -name "*.h" -o -name "*.C" \) -exec clang-format -i {} +

4. If you don't already have a fork of the PeleC repository, follow the `Github instructions <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo>`_ to create one. Then, push a feature branch to your forked PeleC repository: ::

    git remote add remotename git@github.com:remoteurl # add a remote pointing to the user's fork
    git push -u remotename AmazingNewFeature # Note: -u option required only for the first push of new branch

5. Submit a pull request through git@github.com:AMReX-Combustion/PeleC.git, and make sure you are requesting a merge against the development branch

6. Check the CI status on Github and make sure the tests passed for merge request

.. note::

   Github CI uses the CMake build system and CTest to test the core source files of PeleC. If you are adding source files, you will need to add them to the list of source files in the ``CMake`` directory for the tests to pass. Make sure to add them to the GNU make makefiles as well.


Test Status
~~~~~~~~~~~

Nightly test results for PeleC against multiple compilers and machines can be seen on its `CDash page <https://my.cdash.org/index.php?project=Pele>`_.

Documentation
~~~~~~~~~~~~~

The full documentation for Pele exists in the Docs directory; at present this is maintained inline using
Sphinx  `Sphinx <http://www.sphinx-doc.org>`_. With
Sphinx, documentation is written in *Restructured Text*. reST is a markup language
similar to Markdown, but with somewhat greater capabilities (and idiosyncrasies). There
are several `primers <http://thomas-cokelaer.info/tutorials/sphinx/rest_syntax.html>`_
available to get started. One gotcha is that indentation matters. To build ::

    $ cd Docs && mkdir build && cd build && sphinx-build -M html ../sphinx .


Citation
~~~~~~~~

To cite the PeleC software and refer to its computational performance, use the following journal articles for `PeleC <https://doi.org/10.1177/10943420221121151>`_ and the `Pele software suite <https://doi.org/10.1137/1.9781611977967.2>`_::

    @article{PeleC_IJHPCA,
      author = {Marc T {Henry de Frahan} and Jon S Rood and Marc S Day and Hariswaran Sitaraman and Shashank Yellapantula and Bruce A Perry and Ray W Grout and Ann Almgren and Weiqun Zhang and John B Bell and Jacqueline H Chen},
      title = {{PeleC: An adaptive mesh refinement solver for compressible reacting flows}},
      journal = {The International Journal of High Performance Computing Applications},
      volume = {37},
      number = {2},
      pages = {115-131},
      year = {2022},
      doi = {10.1177/10943420221121151},
      url = {https://doi.org/10.1177/10943420221121151}
    }

    @article{PeleSoftware,
      author = {Marc T. {Henry de Frahan} and Lucas Esclapez and Jon Rood and Nicholas T. Wimer and Paul Mullowney and Bruce A. Perry and Landon Owen and Hariswaran Sitaraman and Shashank Yellapantula and Malik Hassanaly and Mohammad J. Rahimi and Michael J. Martin and Olga A. Doronina and Sreejith N. A. and Martin Rieth and Wenjun Ge and Ramanan Sankaran and Ann S. Almgren and Weiqun Zhang and John B. Bell and Ray Grout and Marc S. Day and Jacqueline H. Chen},
      title = {The Pele Simulation Suite for Reacting Flows at Exascale},
      booktitle = {Proceedings of the 2024 SIAM Conference on Parallel Processing for Scientific Computing},
      journal = {Proceedings of the 2024 SIAM Conference on Parallel Processing for Scientific Computing},
      chapter = {},
      pages = {13-25},
      doi = {10.1137/1.9781611977967.2},
      URL = {https://epubs.siam.org/doi/abs/10.1137/1.9781611977967.2},
      eprint = {https://epubs.siam.org/doi/pdf/10.1137/1.9781611977967.2},
      year = {2024},
      publisher = {Proceedings of the 2024 SIAM Conference on Parallel Processing for Scientific Computing}
    }


Additionally, to cite the application of PeleC to compressible reacting flows, use the following `Combustion and Flame journal article <https://doi.org/10.1016/j.combustflame.2021.111531>`_::

  @article{Sitaraman2021,
    author = {Hariswaran Sitaraman and Shashank Yellapantula and Marc T. {Henry de Frahan} and Bruce Perry and Jon Rood and Ray Grout and Marc Day},
    title = {Adaptive mesh based combustion simulations of direct fuel injection effects in a supersonic cavity flame-holder},
    journal = {Combustion and Flame},
    volume = {232},
    pages = {111531},
    year = {2021},
    issn = {0010-2180},
    doi = {https://doi.org/10.1016/j.combustflame.2021.111531},
    url = {https://www.sciencedirect.com/science/article/pii/S0010218021002741},
  }

A full list of publications documenting the development of the Pele suite and its
application to various reacting flow and other simulations is available on the main
`Pele suite page <https://amrex-combustion.github.io/pubs.html>`_. After publication,
if you'd like your work to be included on that list, you can request to have it added
`here <https://github.com/AMReX-Combustion/AMReX-Combustion.github.io/discussions/3>`_.

Versioning
~~~~~~~~~~

PeleC now uses a type of semantic versioning to help users navigate different versions of the code,
which are labeled with `GitHub tags <https://github.com/AMReX-Combustion/PeleC/releases/>`_. These tagged versions are not exhaustive, but they adhere to
the following convention. Given a version number MAJOR.MINOR.PATCH:
1. MAJOR version for changes to key aspects of the solver affecting input/source files for all cases, when a key model is changed to significantly affect results of simulations, when a major new capability is added
2. MINOR version for when a significant feature is added (in a backward compatible manner), accumulation of smaller features, or changes to input file compatibility for less central aspects of the solver (e.g., post-processing) or aspects not affecting all cases
3. PATCH version for backward compatible bug fixes and minor features

PeleC previously used YY.MM formatting for versions. These should be interpreted as version 0 subversions,
e.g. v25.04 is equivalent to v0.25.04.


Acknowledgment
~~~~~~~~~~~~~~

This research was supported by the Exascale Computing Project (ECP), Project
Number: 17-SC-20-SC, a collaborative effort of two DOE organizations -- the
Office of Science and the National Nuclear Security Administration --
responsible for the planning and preparation of a capable exascale ecosystem --
including software, applications, hardware, advanced system engineering, and
early testbed platforms -- to support the nation's exascale computing
imperative.
