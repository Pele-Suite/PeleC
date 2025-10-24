#include <AMReX_CArena.H>
#include <AMReX_REAL.H>
#include <AMReX_Amr.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_EB2.H>

// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
#include <AMReX_Sundials.H>
#endif

#include "PeleC.H"
#include "PeleCAmr.H"

using namespace amrex::mpidatatypes;

std::string inputs_name;

void initialize_EB2(
  const amrex::Geometry& geom,
  const int eb_max_level,
  const int max_level,
  const int coarsening,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  const amrex::IntVect& max_grid_size);

amrex::LevelBld* getLevelBld();

void
override_default_parameters()
{
  amrex::ParmParse pp("eb2");
  if (not pp.contains("geom_type")) {
    std::string geom_type("all_regular");
    pp.add("geom_type", geom_type);
  }
}

int
main(int argc, char* argv[])
{
  if (argc <= 1) {
    amrex::Abort("Error: no inputs file provided on command line.");
  }

  // check to see if it contains --describe
  if (argc >= 2) {
    for (auto i = 1; i < argc; i++) {
      if (std::string(argv[i]) == "--describe") {
        PeleC::writeBuildInfo(std::cout);
        return 0;
      }
    }
  }

#ifdef AMREX_USE_HIP
  // Explicitly initialize the HIP runtime before MPI/AMReX to avoid lazy init
  // race conditions at scale.
  hipError_t herr = hipInit(0);
  if (herr != hipSuccess) {
    fprintf(stderr, "hipInit failed: %s\n", hipGetErrorString(herr));
    return 1;
  }
#endif

  amrex::Initialize(
    argc, argv, true, MPI_COMM_WORLD, override_default_parameters);
// Defined and initialized when in gnumake, but not defined in cmake and
// initialization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Initialize(amrex::OpenMP::get_max_threads());
#endif

  // Save the inputs file name for later.
  if (strchr(argv[1], '=') == nullptr) {
    inputs_name = argv[1];
  }

  BL_PROFILE_VAR("main()", pmain);

  amrex::Real dRunTime1 = amrex::ParallelDescriptor::second();

  amrex::Print() << std::setprecision(10);

  int max_step{-1};
  amrex::Real strt_time{0.0};
  amrex::Real stop_time{-1.0};
  amrex::Real max_wall_time{-1.0};
  amrex::ParmParse pp;

  bool pause_for_debug = false;
  pp.query("pause_for_debug", pause_for_debug);
  if (pause_for_debug) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Enter any string to continue" << std::endl;
      std::string text;
      std::cin >> text;
    }
    amrex::ParallelDescriptor::Barrier();
  }

  pp.query("max_step", max_step);
  pp.query("strt_time", strt_time);
  pp.query("stop_time", stop_time);
  pp.query("max_wall_time", max_wall_time);

  if (strt_time < 0.0) {
    amrex::Abort("MUST SPECIFY a non-negative strt_time");
  }

  if (max_step < 0 && stop_time < 0.0 && max_wall_time < 0.0) {
    amrex::Abort(
      "Exiting because neither max_step nor stop_time nor "
      "max_wall_time is non-negative.");
  }

  // Print the current date and time
  time_t time_type;
  struct tm time_now;
  time(&time_type);
  gmtime_r(&time_type, &time_now);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << std::setfill('0') << "\nStarting run at " << std::setw(2)
                   << time_now.tm_hour << ":" << std::setw(2) << time_now.tm_min
                   << ":" << std::setw(2) << time_now.tm_sec << " UTC on "
                   << time_now.tm_year + 1900 << "-" << std::setw(2)
                   << time_now.tm_mon + 1 << "-" << std::setw(2)
                   << time_now.tm_mday << "." << std::endl;
  }

  bool raise_failure = false;
  amrex::Real dRunTime2 = 0.0;
  {
    std::unique_ptr<PeleCAmr> amrptr =
      std::make_unique<PeleCAmr>(getLevelBld());
    amrex::AmrLevel::SetEBSupportLevel(
      amrex::EBSupport::full); // need both area and volume fractions
    amrex::AmrLevel::SetEBMaxGrowCells(
      PeleC::numGrow() + 1, PeleC::numGrow() + 1, PeleC::numGrow() + 1);

    initialize_EB2(
      amrptr->Geom(PeleC::getEBMaxLevel()), PeleC::getEBMaxLevel(),
      amrptr->maxLevel(), PeleC::getEBCoarsening(), amrptr->refRatio(),
      amrptr->maxGridSize(amrptr->maxLevel()));

    amrptr->init(strt_time, stop_time);

#ifdef AMREX_USE_ASCENT
    amrptr->doInSituViz(amrptr->levelSteps(0));
#endif

    // If we set the regrid_on_restart flag and if we are *not* going to take
    // a time step then we want to go ahead and regrid here.
    if (
      amrptr->RegridOnRestart() && ((amrptr->levelSteps(0) >= max_step) ||
                                    (amrptr->cumTime() >= stop_time))) {
      // Regrid only!
      amrptr->RegridOnly(amrptr->cumTime());
    }

    dRunTime2 = amrex::ParallelDescriptor::second();
    amrex::Real wall_time_elapsed{0.0};

    while (
      (amrptr->okToContinue() != 0) &&
      (amrptr->levelSteps(0) < max_step || max_step < 0) &&
      (amrptr->cumTime() < stop_time || stop_time < 0.0) &&
      (wall_time_elapsed < (max_wall_time * 3600.0) || max_wall_time < 0.0)) {
      // Do a timestep
      amrptr->coarseTimeStep(stop_time);
#ifdef AMREX_USE_ASCENT
      amrptr->doInSituViz(amrptr->levelSteps(0));
#endif
      // Get the elapsed time
      wall_time_elapsed = amrex::ParallelDescriptor::second() - dRunTime1;
      amrex::ParallelDescriptor::ReduceRealMax(wall_time_elapsed);
    }

    // Write final checkpoint
    if (amrptr->stepOfLastCheckPoint() < amrptr->levelSteps(0)) {
      amrptr->checkPoint();
    }

    // Write final plotfile
    if (amrptr->stepOfLastPlotFile() < amrptr->levelSteps(0)) {
      amrptr->writePlotFile();
    }

    time(&time_type);
    gmtime_r(&time_type, &time_now);

    if (amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << std::setfill('0') << "\nEnding run at " << std::setw(2)
                     << time_now.tm_hour << ":" << std::setw(2)
                     << time_now.tm_min << ":" << std::setw(2)
                     << time_now.tm_sec << " UTC on " << time_now.tm_year + 1900
                     << "-" << std::setw(2) << time_now.tm_mon + 1 << "-"
                     << std::setw(2) << time_now.tm_mday << "." << std::endl;
    }

    raise_failure = (amrptr->okToContinue() == 0);
  }

  const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();

  amrex::Real dRunTime3 = amrex::ParallelDescriptor::second();

  amrex::Real runtime_total = dRunTime3 - dRunTime1;
  amrex::Real runtime_timestep = dRunTime3 - dRunTime2;

  amrex::ParallelDescriptor::ReduceRealMax(runtime_total, IOProc);
  amrex::ParallelDescriptor::ReduceRealMax(runtime_timestep, IOProc);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Run time = " << runtime_total << std::endl;
    amrex::Print() << "Run time w/o init = " << runtime_timestep << std::endl;
  }

  if (auto* arena = dynamic_cast<amrex::CArena*>(amrex::The_Arena())) {
    // A barrier to make sure our output follows that of RunStats.
    amrex::ParallelDescriptor::Barrier();
    // We're using a CArena -- output some FAB memory stats.
    // This'll output total # of bytes of heap space in the Arena.
    // It's actually the high water mark of heap space required by FABs.
    char buf[256];

    snprintf(
      buf, 256, "CPU(%d): Heap Space (bytes) used by Coalescing FAB Arena: %zu",
      amrex::ParallelDescriptor::MyProc(), arena->heap_space_used());

    amrex::Print() << buf << std::endl;
  }

  BL_PROFILE_VAR_STOP(pmain);
  BL_PROFILE_SET_RUN_TIME(dRunTime2);

  if (raise_failure) {
    amrex::Abort("Simulation stopped early");
  }

// Defined and finalized when in gnumake, but not defined in cmake and
// finalization done manually
#ifndef AMREX_USE_SUNDIALS
  amrex::sundials::Finalize();
#endif
  amrex::Finalize();

  return 0;
}
