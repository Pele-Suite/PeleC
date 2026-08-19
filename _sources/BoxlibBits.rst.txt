Development Reference
=====================

#Function Listing for PeleC
#--------------------------
#
#When built, the full doxygen documentation for PeleC can be found 
#`here <../../../doxygen_output/html/index.html>`_.



AMReX functions useful for Pele development
-------------------------------------------

Pele is built on AMReX (available at `https://github.com/AMReX-Codes/amrex <https://github.com/AMReX-Codes/amrex>`_), an adaptive mesh refinement software framework, which provides the underlying software infrastructure for block structured AMR operations. Below is a quick reference list with links to many of the AMReX tools used to build up PeleC. The full AMReX documentation can be found `here <https://amrex-codes.github.io/AMReXUsersGuide.pdf>`_. 


Solution environment
~~~~~~~~~~~~~~~~~~~~

* `amrex::BoxArray <https://amrex-codes.github.io/amrex/docs_html/Basics.html#boxarray>`_
* `amrex::StateData <https://amrex-codes.github.io/amrex/docs_html/AmrLevel.html?highlight=statedata#statedata>`_

  * amrex::StateData::copyOld
  * amrex::StateData::setTimeLevel
  * amrex::StateData::hasOldData
  * amrex::StateData::swapTimeLevels

Data structures
~~~~~~~~~~~~~~~

* `amrex::Array <https://amrex-codes.github.io/amrex/docs_html/Basics.html#vector-and-array>`_
* `amrex::FArrayBox <https://amrex-codes.github.io/amrex/docs_html/Basics.html#basefab-farraybox-and-iarraybox>`_

  * BaseFab< T >::dataPtr
  * BaseFab< T >::loVect
  * BaseFab< T >::hiVect
  * BaseFab< T >::nComp
  * BaseFab< T >::setVal(Real)
  * amrex::FArrayBox::resize
  * amrex::FArrayBox::plus

* `TagBox <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=tagbox#tagbox-and-cluster>`_
* `Box <https://amrex-codes.github.io/amrex/docs_html/Basics.html#box-intvect-and-indextype>`_

  * amrex::Box::surroundingNodes(const Box&, int)
  * amrex::grow(const Box&, int)

* `iMultiFab, amrex::MultiFab<FArraybox> <https://amrex-codes.github.io/amrex/docs_html/Basics.html#fabarray-multifab-and-imultifab>`_
   * amrex::MultiFab::SumBoundary(int, int, const Periodicity&)
   * amrex::MultiFab::setVal(value_type)
   * amrex::MultiFab::Copy
   * amrex::MultiFab::Add
   * amrex::MultiFab::define(const BoxArray&, int, int, FabAlloc, const IntVect&, ParallelDescriptor::Color)
   * amrex::MultiFab::sum
   * amrex::MultiFab::clear
   * amrex::MultiFab::FillBoundary
   * amrex::MultiFab::Saxpy
   * amrex::RealBox

PeleC Implementation 
~~~~~~~~~~~~~~~~~~~~

* `amrex::AmrLevel <https://amrex-codes.github.io/amrex/docs_html/AmrLevel.html#amrlevel-class>`_
* StateDescriptor
* `MFIter <https://amrex-codes.github.io/amrex/docs_html/Basics.html#mfiter-and-tiling>`_
* `FillPatch(AmrLevel&, MultiFab&, int, Real, int, int, int, int) <https://amrex-codes.github.io/amrex/docs_html/AsyncIter.html?highlight=fillpatch>`_



Multilevel tools
~~~~~~~~~~~~~~~~
* `amrex::average_down(MultiFab&, MultiFab&, const Geometry&, const Geometry&, int, int, const int) <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=averagedown>`_
* `FluxRegister <https://amrex-codes.github.io/amrex/docs_html/AmrCore.html?highlight=fluxregister#using-fluxregisters>`_
   * FluxRegister::FineAdd(const MultiFab&, int, int, int, int, Real)
   * FluxRegister::CrseInit(const MultiFab&, int, int, int, int, Real, FrOp)
 


Contributing to PeleC
---------------------

Development Model
~~~~~~~~~~~~~~~~~

To add a new feature to PeleC, the procedure is:

1. Create a branch for the new feature (locally) ::

    git checkout -b AmazingNewFeature

2. Develop the feature, merging changes often from the development branch into your AmazingNewFeature branch ::
   
    git commit -m "Developed AmazingNewFeature"
    git checkout development
    git pull                     [fix any identified conflicts between local and remote branches of "development"]
    git checkout AmazingNewFeature
    git merge development        [fix any identified conflicts between "development" and "AmazingNewFeature"]

3. Push feature branch to PeleC repository ::

    git push -u origin AmazingNewFeature [Note: -u option required only for the first push of new branch]

    Check that Travis CI build passed. Currently to reduce load on the CI framework this is just a test that one of the regression tests will compile.

4.  Submit a merge request through git@github.com:Pele-Suite/PeleC.git - be sure you are requesting to merge your branch to the development branch.




Building Regression Test Suite Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. LocalTesting:

It is possibly---and desirable---to create the regression testing framework locally and ensure that all tests pass successfully. This is also a good way to ensure all is properly installed on a new machine to be used for PeleC calculations. The initial setup is somewhat tedious but is worth the effort. What needs to be done is: (1) make a scratch area on the local machine where you manually run the regression tests.  As part of the process, a set of "gold" solutions will be generated using code from the current versions of PeleC, PelePhysics and amrex.  Regression tests afterwards will compare to those solutions and indicate binary compatibility. (2) Clone the required repositories into this area and set required environment variables that point to where everything is. (3) Run the tests to generate the "gold" benchmark data.


1. Make scratch area ::

     mkdir ~/REG_TEST_AREA; cd ~/REG_TEST_AREA
   
2. Clone repositories (amrex, PeleC, PelePhysics, regression_testing (AMReX's driver scripts) and PeleRegressionTesting (Pele-specific stuff)) ::

     git clone git@github.com:AMReX-Combustion/PeleRegressionTesting.git
     cd PeleRegressionTesting; git checkout development
     mkdir -p TestData/PeleC  # this is where the test results will be written
     mkdir Repositories   # this is where the src code to be tested is put
     cd Repositories
     export PELE_HOME=`pwd`/PeleC; git clone git@github.com:Pele-Suite/PeleC.git $PELE_HOME
     cd $PELE_HOME; git checkout development; cd ..
     export PELE_PHYSICS_HOME=`pwd`/PelePhysics; git clone git@github.com:Pele-Suite/PelePhysics.git $PELE_PHYSICS_HOME
     cd $PELE_PHYSICS_HOME; git checkout development; cd ..
     export AMREX_HOME=`pwd`/amrex; git clone git@github.com:AMReX-Codes/amrex.git $AMREX_HOME
     cd $AMREX_HOME; git checkout development; cd ..
     export AMREX_REGTEST_HOME=`pwd`/regression_testing; git clone git@github.com:AMReX-Codes/regression_testing.git $AMREX_REGTEST_HOME
     cd ..

3. Run the script to execute the tests to generate benchmarks. After it finishes building and running (12 as of March 2019) tests, it will archive the pltfiles that result from each into a folder in the TestData/PeleC folder in the PeleRegressionTesting folder.  Once the benchmarks exist, any changes to the repositories in the PeleRegressionTesting/Repositories can be tested to diff clean against these benchmarks by running the script (again from within the PeleRegressionTesting folder) ::

     ./Scripts/genbenchPC.sh

4. Run the tests to execute the tests again, using local modifications (the scripts use a switch for the regtest.py function that says do NOT pull the latest copies off the web for all the repositories in the local testing area) ::

     ./Scripts/runtestsPC.sh

