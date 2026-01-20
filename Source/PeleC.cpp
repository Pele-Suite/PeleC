#include <memory>
#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_Vector.H>
#include <AMReX_TagBox.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EBAmrUtil.H>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#ifdef PELE_USE_SPRAY
#include "SprayParticles.H"
#endif
#endif

#ifdef PELE_USE_SOOT
#include "SootModel.H"
#endif

#ifdef PELE_USE_MASA
#include <masa.h>
using namespace MASA;
#endif

#include "PeleC.H"
#include "Derive.H"
#include "prob.H"
#include "Timestep.H"
#include "PeleCUtilities.H"
#include "Tagging.H"
#include "IndexDefines.H"

#ifdef PELE_ENABLE_FPE_TRAP
#if defined(__linux__)
#include <cfenv>
#elif defined(__APPLE__)
#include <fenv.h>
#endif
#endif

bool PeleC::signalStopJob = false;
int PeleC::verbose = 0;
int PeleC::radius_grow = 1;
amrex::BCRec PeleC::phys_bc;
amrex::Real PeleC::frac_change = std::numeric_limits<amrex::Real>::max();
int PeleC::Density = -1;
int PeleC::Eden = -1;
int PeleC::Eint = -1;
int PeleC::Temp = -1;
int PeleC::Xmom = -1;
int PeleC::Ymom = -1;
int PeleC::Zmom = -1;
int PeleC::FirstSpec = -1;
int PeleC::FirstAux = -1;
int PeleC::FirstAdv = -1;
int PeleC::FirstLin = -1;
int PeleC::NumSootVars = 0;
int PeleC::FirstSootVar = -1;

#include "pelec_defaults.H"

bool PeleC::do_diffuse = false;

#ifdef PELE_USE_SPRAY
bool PeleC::do_spray_particles = true;
#else
bool PeleC::do_spray_particles = false;
#endif

#ifdef PELE_USE_SOOT
bool PeleC::add_soot_src = true;
bool PeleC::plot_soot = true;
#else
bool PeleC::add_soot_src = false;
#endif

#ifdef PELE_USE_MASA
bool PeleC::mms_initialized = false;
#endif

int PeleC::les_model = 0;
int PeleC::les_filter_type = no_filter;
int PeleC::les_filter_fgr = 1;
int PeleC::les_test_filter_type = box_3pt_optimized_approx;
int PeleC::les_test_filter_fgr = 2;

bool PeleC::eb_in_domain = false;
bool PeleC::eb_initialized = false;
int PeleC::eb_max_lvl_gen = -1;
bool PeleC::body_state_set = false;
amrex::GpuArray<amrex::Real, NVAR> PeleC::body_state;

bool PeleC::do_react_load_balance = false;
bool PeleC::do_mol_load_balance = false;

amrex::Vector<std::string> PeleC::spec_names;
amrex::Vector<std::string> PeleC::adv_names;
amrex::Vector<std::string> PeleC::aux_names;

pele::physics::PeleParams<pele::physics::transport::TransParm<
  pele::physics::PhysicsType::eos_type,
  pele::physics::PhysicsType::transport_type>>
  PeleC::trans_parms;

pele::physics::turbinflow::TurbInflow PeleC::turb_inflow;
amrex::Vector<std::string> PeleC::m_diagVars;

amrex::Vector<int> PeleC::src_list;

// this will be reset upon restart
amrex::Real PeleC::previousCPUTimeUsed = 0.0;
amrex::Real PeleC::startCPUTime = 0.0;
int PeleC::num_state_type = 0;

amrex::Real PeleC::typical_rhoY_val_min = 1.e-10;
bool PeleC::use_typical_vals_chem = false;
bool PeleC::use_typical_vals_chem_usr = false;
int PeleC::reset_typical_vals_int = 10;
amrex::Vector<amrex::Real> PeleC::typical_values_chem_usr;

static bool eb_initialized = false;

bool
ebInitialized()
{
  return eb_initialized;
}

void
ebInitialized(bool eb_init_val)
{
  eb_initialized = eb_init_val;
}

int
PeleC::getEBMaxLevel()
{
  // Look into amr PP
  amrex::ParmParse ppa("amr");
  int max_eb_level = -1;

  // Default to amr.max_level
  ppa.query("max_level", max_eb_level);

  // Get the level in the restart file if present
  std::string restart_file;
  ppa.query("restart", restart_file);
  if (!restart_file.empty()) {
    std::string FullPathEBLevelFile = restart_file;
    FullPathEBLevelFile += "/EBMaxLevel";
    std::ifstream EBLevelFile(FullPathEBLevelFile);
    if (!EBLevelFile.fail()) {
      EBLevelFile >> max_eb_level;
      EBLevelFile.close();
    }
  }

  // Allow manual overwrite
  amrex::ParmParse ppeb("eb2");
  ppeb.query("max_level_generation", max_eb_level);

  return max_eb_level;
}

int
PeleC::getEBCoarsening()
{
  amrex::ParmParse ppa("amr");
  std::string restart_file;
  ppa.query("restart", restart_file);

  if (restart_file.empty() and !(init_pltfile.empty())) {
    return init_pltfile_coarse_levels;
  }
  return 0;
}

void
PeleC::read_params()
{
  static bool read_params_done = false;

  if (read_params_done) {
    return;
  }

  read_params_done = true;

  amrex::ParmParse pp("pelec");

#include "pelec_queries.H"

  pp.query("v", verbose);

  // Get boundary conditions
  amrex::Vector<std::string> lo_bc_char(AMREX_SPACEDIM);
  amrex::Vector<std::string> hi_bc_char(AMREX_SPACEDIM);
  pp.getarr("lo_bc", lo_bc_char, 0, AMREX_SPACEDIM);
  pp.getarr("hi_bc", hi_bc_char, 0, AMREX_SPACEDIM);

  amrex::Vector<int> lo_bc(AMREX_SPACEDIM);
  amrex::Vector<int> hi_bc(AMREX_SPACEDIM);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (lo_bc_char[dir] == "Interior") {
      lo_bc[dir] = 0;
    } else if (lo_bc_char[dir] == "Hard") {
      lo_bc[dir] = 1;
    } else if (lo_bc_char[dir] == "FOExtrap") {
      lo_bc[dir] = 2;
    } else if (lo_bc_char[dir] == "Symmetry") {
      lo_bc[dir] = 3;
    } else if (lo_bc_char[dir] == "SlipWall") {
      lo_bc[dir] = 4;
    } else if (lo_bc_char[dir] == "NoSlipWall") {
      lo_bc[dir] = 5;
    } else if (lo_bc_char[dir] == "UserBC") {
      lo_bc[dir] = 6;
    } else {
      amrex::Abort(
        "Wrong boundary condition word in lo_bc, please use: "
        "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }

    if (hi_bc_char[dir] == "Interior") {
      hi_bc[dir] = 0;
    } else if (hi_bc_char[dir] == "Hard") {
      hi_bc[dir] = 1;
    } else if (hi_bc_char[dir] == "FOExtrap") {
      hi_bc[dir] = 2;
    } else if (hi_bc_char[dir] == "Symmetry") {
      hi_bc[dir] = 3;
    } else if (hi_bc_char[dir] == "SlipWall") {
      hi_bc[dir] = 4;
    } else if (hi_bc_char[dir] == "NoSlipWall") {
      hi_bc[dir] = 5;
    } else if (hi_bc_char[dir] == "UserBC") {
      hi_bc[dir] = 6;
    } else {
      amrex::Abort(
        "Wrong boundary condition word in hi_bc, please use: "
        "Interior, UserBC, Symmetry, SlipWall, NoSlipWall");
    }
  }

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    phys_bc.setLo(dir, lo_bc[dir]);
    phys_bc.setHi(dir, hi_bc[dir]);
  }

  typical_values_chem_usr.resize(NUM_SPECIES + 1, 1.0e-10);
  pp.query("use_typ_vals_chem", use_typical_vals_chem);
  pp.query("use_typ_vals_chem_usr", use_typical_vals_chem_usr);

  if (use_typical_vals_chem_usr) {
    use_typical_vals_chem = true;
  }

  pp.query("typical_rhoY_val_min", typical_rhoY_val_min);
  pp.query("reset_typical_vals_int", reset_typical_vals_int);
  pp.queryarr(
    "typical_values_chem", typical_values_chem_usr, 0, NUM_SPECIES + 1);

  // Check phys_bc against possible periodic geometry
  // if periodic, must have internal BC marked.
  // Check, periodic means interior in those directions.
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (amrex::DefaultGeometry().isPeriodic(dir)) {
      if (
        lo_bc[dir] != PCPhysBCType::interior &&
        amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:periodic in direction " << dir
                  << " but low BC is not Interior\n";
        amrex::Error();
      }
      if (
        hi_bc[dir] != PCPhysBCType::interior &&
        amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:periodic in direction " << dir
                  << " but high BC is not Interior\n";
        amrex::Error();
      }
    } else {
      // If not periodic, should not be interior.
      if (
        lo_bc[dir] == PCPhysBCType::interior &&
        amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
      if (
        hi_bc[dir] == PCPhysBCType::interior &&
        amrex::ParallelDescriptor::IOProcessor()) {
        std::cerr << "PeleC::read_params:interior bc in direction " << dir
                  << " but not periodic\n";
        amrex::Error();
      }
    }
  }
  if (do_isothermal_walls) {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (
        domlo_isothermal_temp[dir] > 0.0 &&
        lo_bc[dir] != PCPhysBCType::slip_wall &&
        lo_bc[dir] != PCPhysBCType::no_slip_wall &&
        lo_bc[dir] != PCPhysBCType::user_bc &&
        lo_bc[dir] != PCPhysBCType::inflow) {
        amrex::Abort(
          "Cannot have isothermal wall on a BC that isn't a wall or "
          "user defined BC");
      }
      if (
        domhi_isothermal_temp[dir] > 0.0 &&
        hi_bc[dir] != PCPhysBCType::slip_wall &&
        hi_bc[dir] != PCPhysBCType::no_slip_wall &&
        hi_bc[dir] != PCPhysBCType::user_bc &&
        hi_bc[dir] != PCPhysBCType::inflow) {
        amrex::Abort(
          "Cannot have isothermal wall on a BC that isn't a wall or "
          "user defined BC");
      }
    }
  }

  if (amrex::DefaultGeometry().IsRZ() && (lo_bc[0] != PCPhysBCType::symmetry)) {
    amrex::Error(
      "PeleC::read_params: must set r=0 boundary condition to "
      "Symmetry for r-z");
  }

  // TODO: Any reason to support spherical in PeleC?
  if (amrex::DefaultGeometry().IsRZ()) {
    amrex::Abort("We don't support cylindrical coordinate systems in 3D");
  } else if (amrex::DefaultGeometry().IsSPHERICAL()) {
    amrex::Abort("We don't support spherical coordinate systems in 3D");
  }

  do_diffuse = (diffuse_temp || diffuse_enth || diffuse_spec || diffuse_vel);

  // sanity checks
  if (cfl <= 0.0 || cfl > 1.0) {
    amrex::Error("Invalid CFL factor; must be between zero and one.");
  }
  if (do_hydro && do_mol && (cfl > 0.3)) {
    amrex::Print() << "WARNING -- CFL should be <= 0.3 when using MOL hydro."
                   << std::endl;
  }

  if ((do_les || use_explicit_filter) && (AMREX_SPACEDIM != 3)) {
    amrex::Abort("Using LES/filtering currently requires 3d.");
  }

  if (do_les) {
    pp.query("les_model", les_model);
    pp.query("les_test_filter_type", les_test_filter_type);
    pp.query("les_test_filter_fgr", les_test_filter_fgr);
  }

  if (use_explicit_filter) {
    pp.query("les_filter_type", les_filter_type);
    pp.query("les_filter_fgr", les_filter_fgr);
  }

  // Check on PPM type
  if (do_hydro && (!do_mol)) {
    if (ppm_type != 0 && ppm_type != 1) {
      amrex::Error("PeleC::ppm_type must be 0 (PLM) or 1 (PPM)");
    }
    if (eb_problem_state) {
      amrex::Error("We do not support filling the EB state with do_mol = 0");
    }
  }

  if (ppm_type != 1 && use_hybrid_weno) {
    amrex::Error("PeleC::ppm_type must be 1 (PPM) to use WENO method");
  }

  if (do_hydro) {
    if (do_mol) {
      if ((mol_iorder != 1) && (mol_iorder != 2)) {
        amrex::Error("PeleC::mol_iorder must be 1, or 2.");
      }
    } else if (ppm_type == 0) {
      if ((plm_iorder != 1) && (plm_iorder != 2) && (plm_iorder != 4)) {
        amrex::Error("PeleC::plm_iorder must be 1, 2, or 4");
      }
    }
  }

  // for the moment, ppm_type = 0 does not support ppm_trace_sources --
  // we need to add the momentum sources to the states (and not
  // add it in trans_3d
  if (ppm_type == 0 && ppm_trace_sources) {
    amrex::Print()
      << "WARNING: ppm_trace_sources = 1 not implemented for ppm_type = 0"
      << std::endl;
    ppm_trace_sources = false;
    pp.add("ppm_trace_sources", ppm_trace_sources);
  }

  if (max_dt < fixed_dt) {
    amrex::Error("Cannot have max_dt < fixed_dt");
  }

#ifdef PELE_USE_SPRAY
  readSprayParams();
#endif

#ifdef PELE_USE_SOOT
  pp.query("add_soot_src", add_soot_src);
  pp.query("plot_soot", plot_soot);
  soot_model.readSootParams();
#endif

  // TODO: What is this?
  amrex::StateDescriptor::setBndryFuncThreadSafety(
    static_cast<int>(bndry_func_thread_safe));

  // Get some useful amr inputs
  amrex::ParmParse ppa("amr");
  // This turns on the lb stuff inside Amr, but we use our own flag to signal
  // whether to gather data
  ppa.query("loadbalance_with_workestimates", do_mol_load_balance);
  ppa.query("loadbalance_with_workestimates", do_react_load_balance);

  if (eb_in_domain) {
    int local_bf;
    if (static_cast<bool>(ppa.query("blocking_factor", local_bf))) {
      if (local_bf < 8) {
        amrex::Error("Blocking factor must be at least 8 for EB");
      }
    }
  }
}

PeleC::PeleC()
  : old_sources(num_src),
    new_sources(num_src)
#ifdef PELE_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
  nGrowF = 0;
}

PeleC::PeleC(
  amrex::Amr& papa,
  int lev,
  const amrex::Geometry& level_geom,
  const amrex::BoxArray& bl,
  const amrex::DistributionMapping& dm,
  amrex::Real time)
  : AmrLevel(papa, lev, level_geom, bl, dm, time),
    old_sources(num_src),
    new_sources(num_src)
#ifdef PELE_USE_MASA
    ,
    mms_src_evaluated(false)
#endif
{
  buildMetrics();

  init_eb();

  const amrex::MultiFab& S_new = get_new_data(State_Type);

  for (int src : src_list) {
    int oldGrow = numGrow();
    int newGrow = S_new.nGrow();
    old_sources[src] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, oldGrow, amrex::MFInfo(), Factory());
    new_sources[src] = std::make_unique<amrex::MultiFab>(
      grids, dmap, NVAR, newGrow, amrex::MFInfo(), Factory());
  }

  int nGrowS = numGrow();
#ifdef PELE_USE_SPRAY
  if (do_spray_particles) {
    if (level > 0) {
      nGrowS = amrex::max<int>(
        nGrowS, sprayStateGhosts(parent->MaxRefRatio(level - 1)));
      defineSpraySource(parent->MaxRefRatio(level - 1));
    } else {
      defineSpraySource(1);
    }
  }
#endif
  if (do_hydro || do_diffuse || do_spray_particles) {
    Sborder.define(grids, dmap, NVAR, nGrowS, amrex::MFInfo(), Factory());
  }

  if (!do_mol) {
    if (do_hydro) {
      hydro_source.define(
        grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());

      // This array holds the sum of all source terms that affect the
      // hydrodynamics. If we are doing the source term predictor, we'll also
      // use this after the hydro update to store the sum of the new-time
      // sources, so that we can compute the time derivative of the source
      // terms.
      sources_for_hydro.define(
        grids, dmap, NVAR, numGrow(), amrex::MFInfo(), Factory());
    }
  } else {
    Sborder.define(grids, dmap, NVAR, nGrowS, amrex::MFInfo(), Factory());
  }

  if (do_reflux && level > 0) {
    flux_reg = std::make_unique<amrex::EBFluxRegister>(
      bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
      level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, NVAR);

    if (!amrex::DefaultGeometry().IsCartesian()) {
      // pres_reg.define(
      // bl, papa.boxArray(level - 1), dm, papa.DistributionMap(level - 1),
      // level_geom, papa.Geom(level - 1), papa.refRatio(level - 1), level, 1);
      amrex::Abort("We don't do rz.");
    }
  }

  get_new_data(Reactions_Type).setVal(0.0);

  // Don't need this in pure C++?
  // initialize the Godunov state array used in hydro -- we wait
  // until here so that ngroups is defined (if needed) in
  // rad_params_module
  // if (do_hydro)
  //{
  //  init_godunov_indices();
  //}

  // Initialize the reactor
  if (do_react) {
    init_reactor();
  }

  // initialize LES variables
  if (do_les) {
    init_les();
  }

  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter) {
    init_filters();
  }

  // initialize diagnostics (only level 0 calls them)
  init_diagnostics();
}

PeleC::~PeleC()
{
  if (do_react) {
    close_reactor();
  }
}

void
PeleC::buildMetrics()
{
  const int ngrd = static_cast<int>(grids.size());

  radius.resize(ngrd);

  const amrex::Real* dx = geom.CellSize();

  for (int i = 0; i < ngrd; i++) {
    const amrex::Box& b = grids[i];
    int ilo = b.smallEnd(0) - radius_grow;
    int ihi = b.bigEnd(0) + radius_grow;
    int len = ihi - ilo + 1;

    radius[i].resize(len);

    amrex::Real* rad = radius[i].dataPtr();

    if (amrex::DefaultGeometry().IsCartesian()) {
      for (int j = 0; j < len; j++) {
        rad[j] = 1.0;
      }
    } else {
      amrex::RealBox gridloc =
        amrex::RealBox(grids[i], geom.CellSize(), geom.ProbLo());

      const amrex::Real xlo = gridloc.lo(0) + (0.5 - radius_grow) * dx[0];

      for (int j = 0; j < len; j++) {
        rad[j] = xlo + j * dx[0];
      }
    }
  }

  volume.clear();
  volume.define(
    grids, dmap, 1, numGrow(), amrex::MFInfo(), amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, numGrow(), amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }

  vfrac.clear();
  vfrac.define(grids, dmap, 1, numGrow(), amrex::MFInfo(), Factory());
  const auto& ebfactory =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
  amrex::MultiFab::Copy(vfrac, ebfactory.getVolFrac(), 0, 0, 1, numGrow());

  level_mask.clear();
  level_mask.define(grids, dmap, 1, 3);
  level_mask.BuildMask(
    geom.Domain(), geom.periodicity(), constants::level_mask_covered(),
    constants::level_mask_notcovered(), constants::level_mask_physbnd(),
    constants::level_mask_interior());
}

void
PeleC::setTimeLevel(amrex::Real time, amrex::Real dt_old, amrex::Real dt_new)
{
  AmrLevel::setTimeLevel(time, dt_old, dt_new);
}

/*
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_prob_close()
{
}
*/

void
PeleC::initData()
{
  BL_PROFILE("PeleC::initData()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::h_prob_parm_device,
    PeleC::h_prob_parm_device + 1, PeleC::d_prob_parm_device);

  amrex::MultiFab& S_new = get_new_data(State_Type);

  S_new.setVal(0.0);

#if AMREX_SPACEDIM > 1
  // make sure dx = dy = dz -- that's all we guarantee to support
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::Real small = 1.e-13;
  const bool unequal_dx =
    amrex::max<amrex::Real>(AMREX_D_DECL(
      static_cast<amrex::Real>(0.0),
      static_cast<amrex::Real>(std::abs(dx[0] - dx[1])),
      static_cast<amrex::Real>(std::abs(dx[0] - dx[2])))) > small * dx[0];
  if (eb_in_domain && unequal_dx) {
    amrex::Abort("dx != dy != dz not supported with EB");
  }
#endif

  if (verbose != 0) {
    amrex::Print() << "Initializing the data at level " << level << std::endl;
  }

  get_new_data(Reactions_Type).setVal(0.0);

  if (do_mol_load_balance || do_react_load_balance) {
    get_new_data(Work_Estimate_Type).setVal(1.0);
  }

  if (init_pltfile.empty()) {
    const auto geomdata = geom.data();
    const ProbParmDevice* lprobparm = d_prob_parm_device;
    auto sarrs = S_new.arrays();
    amrex::ParallelFor(
      S_new, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        pc_initdata(i, j, k, sarrs[nbx], geomdata, *lprobparm);
        // Verify that the sum of (rho Y)_i = rho at every cell
        pc_check_initial_species(i, j, k, sarrs[nbx]);
      });
    amrex::Gpu::synchronize();
  } else {
    initLevelDataFromPlt(level, init_pltfile, S_new);
  }

  enforce_consistent_e(S_new);

  set_body_state(S_new);
  amrex::Real cur_time = state[State_Type].curTime();
  const amrex::StateDescriptor* desc = state[State_Type].descriptor();
  const auto& bcs = desc->getBCs();
  InitialRedistribution(cur_time, bcs, S_new);

#ifdef PELE_USE_SPRAY
  if (level == 0) {
    initParticles();
  } else {
    particle_redistribute(level - 1);
  }
#endif

  if (verbose != 0) {
    amrex::Print() << "Done initializing level " << level << " data "
                   << std::endl;
  }
}

void
PeleC::init(AmrLevel& old)
{
  BL_PROFILE("PeleC::init(old)");

  auto* oldlev = (PeleC*)&old;

  // Create new grid data by fillpatching from old.
  amrex::Real dt_new = parent->dtLevel(level);
  amrex::Real cur_time = oldlev->state[State_Type].curTime();
  amrex::Real prev_time = oldlev->state[State_Type].prevTime();
  amrex::Real dt_old = cur_time - prev_time;
  setTimeLevel(cur_time, dt_old, dt_new);

  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillPatch(old, S_new, 0, cur_time, State_Type, 0, NVAR);

  amrex::MultiFab& React_new = get_new_data(Reactions_Type);

  if (do_react) {
    FillPatch(
      old, React_new, 0, cur_time, Reactions_Type, 0, React_new.nComp());
  } else {
    React_new.setVal(0);
  }

  if (do_mol_load_balance || do_react_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    FillPatch(
      old, work_estimate_new, 0, cur_time, Work_Estimate_Type, 0,
      work_estimate_new.nComp());
  }
}

void
PeleC::init()
{
  // This version inits the data on a new level that did not
  // exist before regridding.
  BL_PROFILE("PeleC::init()");

  amrex::Real dt = parent->dtLevel(level);
  amrex::Real cur_time = getLevel(level - 1).state[State_Type].curTime();
  amrex::Real prev_time = getLevel(level - 1).state[State_Type].prevTime();

  amrex::Real dt_old =
    (cur_time - prev_time) / (amrex::Real)parent->MaxRefRatio(level - 1);

  setTimeLevel(cur_time, dt_old, dt);
  amrex::MultiFab& S_new = get_new_data(State_Type);
  FillCoarsePatch(S_new, 0, cur_time, State_Type, 0, NVAR);

  if (do_mol_load_balance || do_react_load_balance) {
    amrex::MultiFab& work_estimate_new = get_new_data(Work_Estimate_Type);
    int ncomp = work_estimate_new.nComp();
    FillCoarsePatch(
      work_estimate_new, 0, cur_time, Work_Estimate_Type, 0, ncomp);
  }
}

amrex::Real
PeleC::initialTimeStep()
{
  BL_PROFILE("PeleC::initialTimeStep()");

  amrex::Real init_dt = 0.0;

  if (initial_dt > 0.0) {
    init_dt = initial_dt;
  } else {
    const amrex::Real dummy_dt = 0.0;
    init_dt = init_shrink * estTimeStep(dummy_dt);
  }

  return init_dt;
}

amrex::Real
PeleC::estTimeStep(amrex::Real /*dt_old*/)
{
  BL_PROFILE("PeleC::estTimeStep()");

  if (fixed_dt > 0.0) {
    return fixed_dt;
  }

  amrex::Real estdt = max_dt;

  const amrex::MultiFab& stateMF = get_new_data(State_Type);

  const amrex::Real* dx = geom.CellSize();

  std::string limiter = "pelec.max_dt";

  // Start the hydro with the max_dt value, but divide by CFL
  // to account for the fact that we multiply by it at the end.
  // This ensures that if max_dt is more restrictive than the hydro
  // criterion, we will get exactly max_dt for a timestep.

  const amrex::Real max_dt_over_cfl = max_dt / cfl;
  amrex::Real estdt_hydro = max_dt_over_cfl;
  amrex::Real estdt_vdif = max_dt_over_cfl;
  amrex::Real estdt_tdif = max_dt_over_cfl;
  amrex::Real estdt_edif = max_dt_over_cfl;
  if (do_hydro || do_mol || diffuse_vel || diffuse_temp || diffuse_enth) {

    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(stateMF.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

    amrex::Real AMREX_D_DECL(dx1 = dx[0], dx2 = dx[1], dx3 = dx[2]);

    if (do_hydro) {
      amrex::Real dt = amrex::ReduceMin(
        stateMF, flags, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr)
          -> amrex::Real {
          return pc_estdt_hydro(
            bx, fab_arr, flag_arr, AMREX_D_DECL(dx1, dx2, dx3));
        });
      estdt_hydro = amrex::min<amrex::Real>(estdt_hydro, dt);
    }

    if (diffuse_vel) {
      auto const& geomdata = geom.data();
      auto const* ltransparm = trans_parms.device_parm();
      const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
      amrex::Real dt = amrex::ReduceMin(
        stateMF, flags, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr)
          -> amrex::Real {
          return pc_estdt_veldif(
            bx, fab_arr, flag_arr, geomdata, ltransparm, *lprobparm);
        });
      estdt_vdif = amrex::min<amrex::Real>(estdt_vdif, dt);
    }

    if (diffuse_temp) {
      auto const& geomdata = geom.data();
      auto const* ltransparm = trans_parms.device_parm();
      const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
      amrex::Real dt = amrex::ReduceMin(
        stateMF, flags, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr)
          -> amrex::Real {
          return pc_estdt_tempdif(
            bx, fab_arr, flag_arr, geomdata, ltransparm, *lprobparm);
        });
      estdt_tdif = amrex::min<amrex::Real>(estdt_tdif, dt);
    }

    if (diffuse_enth) {
      auto const& geomdata = geom.data();
      auto const* ltransparm = trans_parms.device_parm();
      const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
      amrex::Real dt = amrex::ReduceMin(
        stateMF, flags, 0,
        [=] AMREX_GPU_HOST_DEVICE(
          amrex::Box const& bx, const amrex::Array4<const amrex::Real>& fab_arr,
          const amrex::Array4<const amrex::EBCellFlag>& flag_arr)
          -> amrex::Real {
          return pc_estdt_enthdif(
            bx, fab_arr, flag_arr, geomdata, ltransparm, *lprobparm);
        });
      estdt_edif = amrex::min<amrex::Real>(estdt_edif, dt);
    }

    estdt_hydro = amrex::min<amrex::Real>(
      estdt_hydro,
      amrex::min<amrex::Real>(
        estdt_vdif, amrex::min<amrex::Real>(estdt_tdif, estdt_edif)));

    amrex::ParallelDescriptor::ReduceRealMin(estdt_hydro);
    estdt_hydro *= cfl;

    if (verbose != 0) {
      amrex::Print() << "...estimated hydro-limited timestep at level " << level
                     << ": " << estdt_hydro << std::endl;
    }

    // Determine if this is more restrictive than the maximum timestep limiting
    if (estdt_hydro < estdt) {
      limiter = "hydro";
      estdt = estdt_hydro;
    }
  }

#ifdef PELE_USE_SPRAY
  amrex::Real estdt_particle = max_dt;
  if (do_spray_particles) {
    estTimeStepParticles(estdt_particle);
    if (estdt_particle < estdt) {
      limiter = "particles";
      estdt = estdt_particle;
    }
  }
#endif

  if (verbose != 0) {
    amrex::Print() << "PeleC::estTimeStep (" << limiter << "-limited) at level "
                   << level << ":  estdt = " << estdt << '\n';
  }

  return estdt;
}

void
PeleC::computeNewDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_min,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time,
  int post_regrid_flag)
{
  BL_PROFILE("PeleC::computeNewDt()");

  // We are at the start of a coarse grid timecycle.
  // Compute the timesteps for the next iteration.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = std::numeric_limits<amrex::Real>::max();
  int n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    PeleC& adv_level = getLevel(i);
    dt_min[i] = adv_level.estTimeStep(dt_level[i]);
  }

  if (fixed_dt <= 0.0) {
    if (post_regrid_flag == 1) {
      // Limit dt's by pre-regrid dt
      for (int i = 0; i <= finest_level; i++) {
        dt_min[i] = amrex::min<amrex::Real>(dt_min[i], dt_level[i]);
      }
    } else {
      // Limit dt's by change_max * old dt
      for (int i = 0; i <= finest_level; i++) {
        if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
          if (dt_min[i] > change_max * dt_level[i]) {
            amrex::Print() << "PeleC::compute_new_dt : limiting dt at level "
                           << i << '\n';
            amrex::Print() << " ... new dt computed: " << dt_min[i] << '\n';
            amrex::Print() << " ... but limiting to: "
                           << change_max * dt_level[i] << " = " << change_max
                           << " * " << dt_level[i] << '\n';
          }
        }
        dt_min[i] =
          amrex::min<amrex::Real>(dt_min[i], change_max * dt_level[i]);
      }
    }
  }

  // Find the minimum over all levels
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_0 = amrex::min<amrex::Real>(dt_0, n_factor * dt_min[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
PeleC::computeInitialDt(
  int finest_level,
  int /*sub_cycle*/,
  amrex::Vector<int>& n_cycle,
  const amrex::Vector<amrex::IntVect>& /*ref_ratio*/,
  amrex::Vector<amrex::Real>& dt_level,
  amrex::Real stop_time)
{
  BL_PROFILE("PeleC::computeInitialDt()");

  // Grids have been constructed, compute dt for all levels.
  if (level > 0) {
    return;
  }

  amrex::Real dt_0 = std::numeric_limits<amrex::Real>::max();
  int n_factor = 1;
  // TODO: This will need to change for optimal subcycling.
  for (int i = 0; i <= finest_level; i++) {
    dt_level[i] = getLevel(i).initialTimeStep();
    n_factor *= n_cycle[i];
    dt_0 = amrex::min<amrex::Real>(dt_0, n_factor * dt_level[i]);
  }

  // Limit dt's by the value of stop_time.
  const amrex::Real dt_eps = 0.001 * dt_0;
  amrex::Real cur_time = state[State_Type].curTime();
  if (stop_time >= 0.0) {
    if ((cur_time + dt_0) > (stop_time - dt_eps)) {
      dt_0 = stop_time - cur_time;
    }
  }

  n_factor = 1;
  for (int i = 0; i <= finest_level; i++) {
    n_factor *= n_cycle[i];
    dt_level[i] = dt_0 / n_factor;
  }
}

void
PeleC::post_timestep(int iteration)
{
  BL_PROFILE("PeleC::post_timestep()");

  const int finest_level = parent->finestLevel();

#ifdef PELE_USE_SPRAY
  postTimeStepParticles(iteration);
#else
  amrex::ignore_unused(iteration);
#endif

  if (do_reflux && level < finest_level) {
    reflux();
  }

  if (level < finest_level) {
    // We need to do this before anything else because refluxing changes the
    // values of coarse cells underneath fine grids with the assumption they'll
    // be over-written by averaging down
    avgDown();

    // fillpatcher on level+1 needs to be reset because data on this
    // level have changed.
    getLevel(level + 1).resetFillPatcher();
  }

  // Re-compute temperature after all the other updates.
  amrex::MultiFab& S_new = get_new_data(State_Type);
  int ng_pts = 0;
  computeTemp(S_new, ng_pts);
  set_body_state(S_new);

#ifdef PELE_USE_SOOT
  clipSootMoments(S_new, ng_pts);
#endif

  problem_post_timestep();

  if (level == 0) {
    int nstep = parent->levelSteps(0);
    amrex::Real dtlev = parent->dtLevel(0);
    amrex::Real cumtime = parent->cumTime() + dtlev;

    bool sum_int_test = (sum_interval > 0 && nstep % sum_interval == 0);

    bool sum_per_test = false;

    if (sum_per > 0.0) {
      const int num_per_old =
        static_cast<int>(std::floor((cumtime - dtlev) / sum_per));
      const int num_per_new = static_cast<int>(std::floor((cumtime) / sum_per));

      if (num_per_old != num_per_new) {
        sum_per_test = true;
      }
    }

    if (sum_int_test || sum_per_test) {
      sum_integrated_quantities();
      if (track_extrema) {
        monitor_extrema();
      }
    }
  }

  if (
    do_react && use_typical_vals_chem &&
    parent->levelSteps(0) % reset_typical_vals_int == 0) {
    set_typical_values_chem();
  }

  do_diagnostics();
}

void
PeleC::post_restart()
{
  BL_PROFILE("PeleC::post_restart()");

  // Copy problem parameter structs to device
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::h_prob_parm_device,
    PeleC::h_prob_parm_device + 1, PeleC::d_prob_parm_device);

#ifdef PELE_USE_SPRAY
  postRestartParticles(parent->theRestartFile());
#endif
  // Initialize the reactor
  if (do_react) {
    init_reactor();
    if (use_typical_vals_chem) {
      set_typical_values_chem();
    }
  }

  // initialize LES variables
  if (do_les) {
    init_les();
  }

  // initialize filters and variables
  nGrowF = 0;
  if (use_explicit_filter) {
    init_filters();
  }

  // initialize diagnostics
  init_diagnostics();

  problem_post_restart();
}

void
PeleC::postCoarseTimeStep(amrex::Real cumtime)
{
  BL_PROFILE("PeleC::postCoarseTimeStep()");
  AmrLevel::postCoarseTimeStep(cumtime);
}

void
PeleC::post_regrid(int lbase, int /*new_finest*/)
{
  BL_PROFILE("PeleC::post_regrid()");
  fine_mask.clear();

#ifdef PELE_USE_SPRAY
  if (lbase == level) {
    particle_redistribute(lbase);
  }
#else
  amrex::ignore_unused(lbase);
#endif

  if ((do_react) && (use_typical_vals_chem)) {
    set_typical_values_chem();
  }
}

void
PeleC::post_init(amrex::Real /*stop_time*/)
{
  BL_PROFILE("PeleC::post_init()");

  amrex::Real dtlev = parent->dtLevel(level);
  amrex::Real cumtime = parent->cumTime();

  // Fill Reactions_Type data based on initial dt
  if (do_react) {

    bool react_init = true;
    if (use_typical_vals_chem) {
      set_typical_values_chem();
    }

    react_state(cumtime, dtlev, react_init);
  }

  if (level > 0) {
    return;
  }

  // Average data down from finer levels
  // so that conserved data is consistent between levels.
  if (do_avg_down) {
    int finest_level = parent->finestLevel();
    for (int k = finest_level - 1; k >= 0; k--) {
      getLevel(k).avgDown();
    }
  }

  // Allow the user to define their own post_init functions.
  problem_post_init();

#ifdef PELE_USE_SPRAY
  postInitParticles();
#endif

  int nstep = parent->levelSteps(0);
  if (cumtime != 0.0) {
    cumtime += dtlev;
  }

  bool sum_int_test = false;

  if (sum_interval > 0) {
    if (nstep % sum_interval == 0) {
      sum_int_test = true;
    }
  }

  bool sum_per_test = false;

  if (sum_per > 0.0) {
    const int num_per_old =
      static_cast<int>(std::floor((cumtime - dtlev) / sum_per));
    const int num_per_new = static_cast<int>(std::floor((cumtime) / sum_per));

    if (num_per_old != num_per_new) {
      sum_per_test = true;
    }
  }

  if (sum_int_test || sum_per_test) {
    sum_integrated_quantities();
    if (track_extrema) {
      monitor_extrema();
    }
  }

  do_diagnostics();
}

int
PeleC::okToContinue()
{
  if (level > 0) {
    return 1;
  }

  int test = 1;

  if (signalStopJob) {
    test = 0;

    amrex::Print()
      << " Signalling a stop of the run due to signalStopJob = true."
      << std::endl;
  } else if (parent->dtLevel(0) < dt_cutoff) {
    test = 0;

    amrex::Print() << " Signalling a stop of the run because dt < dt_cutoff."
                   << std::endl;
  }

  if (
    get_new_data(State_Type).contains_nan() ||
    get_new_data(State_Type).contains_inf()) {
    test = 0;
    amrex::Print() << " Signalling a stop because NaNs detected in the Solution"
                   << std::endl;
  }

  return test;
}

void
PeleC::update_flux_registers(
  const amrex::Real dt,
  const amrex::MFIter& mfi,
  const amrex::FabType& typ,
  const std::array<amrex::FArrayBox const*, AMREX_SPACEDIM>& flux,
  const amrex::FArrayBox& dm_as_fine)
{
  BL_PROFILE("PeleC::update_flux_registers()");

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::Abort("Flux registers not r-z compatible yet");
  }

  if (!do_reflux) {
    return;
  }

  amrex::EBFluxRegister* fr_as_crse =
    (level < parent->finestLevel()) ? &getFluxReg(level + 1) : nullptr;
  amrex::EBFluxRegister* fr_as_fine =
    (level > 0) ? &getFluxReg(level) : nullptr;

  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  const amrex::Real dx1 = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxD = {
    {AMREX_D_DECL(dx1, dx1, dx1)}};

  if (typ == amrex::FabType::singlevalued) {
    if (fr_as_crse != nullptr) {
      fr_as_crse->CrseAdd(
        mfi, flux, dxD.data(), dt, vfrac[mfi],
        {AMREX_D_DECL(
          &((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
          &((*areafrac[2])[mfi]))},
        amrex::RunOn::Device);
    }

    if (fr_as_fine != nullptr) {
      fr_as_fine->FineAdd(
        mfi, flux, dxD.data(), dt, vfrac[mfi],
        {AMREX_D_DECL(
          &((*areafrac[0])[mfi]), &((*areafrac[1])[mfi]),
          &((*areafrac[2])[mfi]))},
        dm_as_fine, amrex::RunOn::Device);
    }

  } else if (typ == amrex::FabType::regular) {
    if ((level < parent->finestLevel()) && (fr_as_crse != nullptr)) {
      fr_as_crse->CrseAdd(mfi, flux, dxD.data(), dt, amrex::RunOn::Device);
    }

    if ((level > 0) && (fr_as_fine != nullptr)) {
      fr_as_fine->FineAdd(mfi, flux, dxD.data(), dt, amrex::RunOn::Device);
    }
  } else if (typ == amrex::FabType::multivalued) {
    amrex::Abort("multi-valued EB flux register update to be implemented");
  }
}

void
PeleC::reflux()
{
  BL_PROFILE("PeleC::reflux()");

  AMREX_ASSERT(level < parent->finestLevel());

  const amrex::Real strt = amrex::ParallelDescriptor::second();

  PeleC& fine_level = getLevel(level + 1);
  amrex::MultiFab& S_crse = get_new_data(State_Type);
  amrex::MultiFab& S_fine = fine_level.get_new_data(State_Type);
  getFluxReg(level + 1).Reflux(S_crse, vfrac, S_fine, fine_level.vfrac);

  if (!amrex::DefaultGeometry().IsCartesian() && eb_in_domain) {
    amrex::Abort("rz not yet compatible with EB");
  }

  if (!amrex::DefaultGeometry().IsCartesian()) {
    amrex::MultiFab dr(
      volume.boxArray(), volume.DistributionMap(), 1, volume.nGrow(),
      amrex::MFInfo(), amrex::FArrayBoxFactory());
    dr.setVal(geom.CellSize(0));
    amrex::Abort("PeleC reflux not yet ready for r-z");
  }

  if (verbose != 0) {
    const int IOProc = amrex::ParallelDescriptor::IOProcessorNumber();
    amrex::Real end = amrex::ParallelDescriptor::second() - strt;

    amrex::ParallelDescriptor::ReduceRealMax(end, IOProc);

    amrex::Print() << "PeleC::reflux() at level " << level
                   << " : time = " << end << std::endl;
  }
}

void
PeleC::avgDown()
{
  BL_PROFILE("PeleC::avgDown()");

  if (level == parent->finestLevel()) {
    return;
  }

  avgDown(State_Type);
  avgDown(Reactions_Type);
}

void
PeleC::enforce_consistent_e(amrex::MultiFab& S)
{

  auto sarrs = S.arrays();
  amrex::ParallelFor(
    S, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      auto sarr = sarrs[nbx];
      const amrex::Real rhoInv = 1.0 / sarr(i, j, k, URHO);
      const amrex::Real u = sarr(i, j, k, UMX) * rhoInv;
      const amrex::Real v = sarr(i, j, k, UMY) * rhoInv;
      const amrex::Real w = sarr(i, j, k, UMZ) * rhoInv;
      sarr(i, j, k, UEDEN) = sarr(i, j, k, UEINT) + 0.5 * sarr(i, j, k, URHO) *
                                                      (u * u + v * v + w * w);
    });

  amrex::Gpu::synchronize();

  if (do_rf) {
    const auto geomdata = geom.data();
    int axis = rf_axis;
    amrex::Real omega = rf_omega;
    amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
      AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};

    amrex::ParallelFor(
      S, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        auto sarr = sarrs[nbx];
        amrex::IntVect iv(AMREX_D_DECL(i, j, k));
        amrex::Real rotenrg =
          get_rot_energy(iv, omega, axis, axis_loc, geomdata);
        sarr(i, j, k, UEDEN) =
          sarr(i, j, k, UEINT) - sarr(i, j, k, URHO) * rotenrg;
      });
  }
}

void
PeleC::avgDown(int state_indx)
{
  BL_PROFILE("PeleC::avgDown(state_indx)");

  if (level == parent->finestLevel()) {
    return;
  }

  amrex::MultiFab& S_crse = get_new_data(state_indx);
  const amrex::MultiFab& S_fine = getLevel(level + 1).get_new_data(state_indx);

  if (eb_in_domain) {
    PeleC& fine_lev = getLevel(level + 1);

    amrex::EB_average_down(
      S_fine, S_crse, fine_lev.Volume(), fine_lev.volFrac(), 0, S_fine.nComp(),
      fine_ratio);

    if (state_indx == State_Type) {
      set_body_state(S_crse); // TODO: Is this necessary?
    }
  } else {
    const amrex::Geometry& fgeom = getLevel(level + 1).geom;
    const amrex::Geometry& cgeom = geom;

    amrex::average_down(
      S_fine, S_crse, fgeom, cgeom, 0, S_fine.nComp(), fine_ratio);
  }
}

void
PeleC::allocOldData()
{
  for (int k = 0; k < num_state_type; k++) {
    state[k].allocOldData();
  }
}

void
PeleC::removeOldData()
{
  AmrLevel::removeOldData();
}

void
PeleC::errorEst(
  amrex::TagBoxArray& tags,
  int /*clearval*/,
  int /*tagval*/,
  amrex::Real time,
  int /*n_error_buf*/,
  int /*ngrow*/)
{
  BL_PROFILE("PeleC::errorEst()");

  amrex::MultiFab S_data(
    get_new_data(State_Type).boxArray(),
    get_new_data(State_Type).DistributionMap(), NVAR, 1, amrex::MFInfo(),
    Factory());
  const amrex::Real cur_time = state[State_Type].curTime();
  FillPatch(
    *this, S_data, S_data.nGrow(), cur_time, State_Type, Density, NVAR, 0);

  amrex::Vector<amrex::BCRec> bcs(NVAR);
  const char tagval = amrex::TagBox::SET;

  // Tag EB
  if (eb_in_domain) {
    if (
      ((tagging_parm->eb_refine_type == "static") &&
       (level < tagging_parm->max_eb_refine_lev)) ||
      ((tagging_parm->eb_refine_type == "adaptive") &&
       (level < tagging_parm->adapt_eb_refined_lev))) {
      amrex::TagCutCells(tags, S_data);
    }
  }

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_data.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& tilebox = mfi.tilebox();
      const auto Sfab = S_data.array(mfi);
      auto tag_arr = tags.array(mfi);
      const auto datbox = amrex::grow(tilebox, 1);
      const auto vfrac_arr = vfrac.array(mfi);

      amrex::FArrayBox S_derData(datbox, 1, amrex::The_Async_Arena());
      auto S_derarr = S_derData.array();
      const int ncp = S_derData.nComp();
      const int* bc = bcs[0].data();

      const auto& flag_arr = flags.const_array(mfi);

      // Tagging density
      if (level < tagging_parm->max_denerr_lev) {
        const amrex::Real captured_denerr = tagging_parm->denerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(
              i, j, k, flag_arr, tag_arr, Sfab, captured_denerr, tagval);
          });
      }
      if (level < tagging_parm->max_dengrad_lev) {
        const amrex::Real captured_dengrad = tagging_parm->dengrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, Sfab, captured_dengrad, tagval);
          });
      }
      if (level < tagging_parm->max_denratio_lev) {
        const amrex::Real captured_denratio = tagging_parm->denratio;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_ratioerror(
              i, j, k, flag_arr, tag_arr, Sfab, captured_denratio, tagval);
          });
      }

      // Tagging pressure
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_derpres(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_presserr_lev) {
        const amrex::Real captured_presserr = tagging_parm->presserr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_presserr, tagval);
          });
      }
      if (level < tagging_parm->max_pressgrad_lev) {
        const amrex::Real captured_pressgrad = tagging_parm->pressgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_pressgrad, tagval);
          });
      }

      // Tagging vel_x
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_derdividebyrho<UMX, 1>(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging vel_y
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_derdividebyrho<UMY, 1>(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging vel_z
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_derdividebyrho<UMZ, 1>(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_velerr_lev) {
        const amrex::Real captured_velerr = tagging_parm->velerr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velerr, tagval);
          });
      }
      if (level < tagging_parm->max_velgrad_lev) {
        const amrex::Real captured_velgrad = tagging_parm->velgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_velgrad, tagval);
          });
      }

      // Tagging magnitude of vorticity
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_dermagvort(
        tilebox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_vorterr_lev) {
        const amrex::Real vorterr =
          tagging_parm->vorterr * std::pow(2.0, level);
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_abserror(i, j, k, flag_arr, tag_arr, S_derarr, vorterr, tagval);
          });
      }

      // Tagging temperature
      S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
      pc_dertemp(
        datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
        level);
      if (level < tagging_parm->max_temperr_lev) {
        const amrex::Real captured_temperr = tagging_parm->temperr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_error(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_temperr, tagval);
          });
      }
      if (level < tagging_parm->max_lotemperr_lev) {
        const amrex::Real captured_lotemperr = tagging_parm->lotemperr;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_loerror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_lotemperr, tagval);
          });
      }
      if (level < tagging_parm->max_tempgrad_lev) {
        const amrex::Real captured_tempgrad = tagging_parm->tempgrad;
        amrex::ParallelFor(
          tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            tag_graderror(
              i, j, k, flag_arr, tag_arr, S_derarr, captured_tempgrad, tagval);
          });
      }

      // Tagging flame tracer
      if (!flame_trac_name.empty()) {
        int idx = find_position(spec_names, flame_trac_name);

        if (idx >= 0) {
          S_derData.setVal<amrex::RunOn::Device>(0.0, datbox);
          pc_derspectrac(
            datbox, S_derData, ncp, Sfab.nComp(), S_data[mfi], geom, time, bc,
            level, idx);

          if (level < tagging_parm->max_ftracerr_lev) {
            const amrex::Real captured_ftracerr = tagging_parm->ftracerr;
            amrex::ParallelFor(
              tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tag_error(
                  i, j, k, flag_arr, tag_arr, S_derarr, captured_ftracerr,
                  tagval);
              });
          }
          if (level < tagging_parm->max_ftracgrad_lev) {
            const amrex::Real captured_ftracgrad = tagging_parm->ftracgrad;
            amrex::ParallelFor(
              tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                tag_graderror(
                  i, j, k, flag_arr, tag_arr, S_derarr, captured_ftracgrad,
                  tagval);
              });
          }

        } else {
          amrex::Abort("Unknown species identified as flame_trac_name");
        }
      }

      if (eb_in_domain) {
        // Tagging volume fraction
        if (level < tagging_parm->max_vfracerr_lev) {
          amrex::ParallelFor(
            tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              tag_error_bounds(
                i, j, k, flag_arr, tag_arr, vfrac_arr, 0.0, 1.0, tagval);
            });

          const int local_i = mfi.LocalIndex();
          const auto Nebg = sv_eb_bndry_geom[local_i].size();
          EBBndryGeom* ebg = sv_eb_bndry_geom[local_i].data();
          amrex::ParallelFor(Nebg, [=] AMREX_GPU_DEVICE(int L) {
            const auto& iv = ebg[L].iv;
            if (tilebox.contains(iv)) {
              tag_arr(iv) = tagval;
            }
          });
        }
      }
    }
  }

  // amrex tagging utils
  for (const auto& err_tag : tagging_parm->err_tags) {
    std::unique_ptr<amrex::MultiFab> mf;
    if (!err_tag.Field().empty()) {
      mf = derive(err_tag.Field(), time, err_tag.NGrow());
    }
    err_tag(
      tags, mf.get(), amrex::TagBox::CLEAR, amrex::TagBox::SET, time, level,
      geom);
  }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(S_data, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& tilebox = mfi.tilebox();
      const auto Sfab = S_data.array(mfi);
      auto tag_arr = tags.array(mfi);
      const auto& flag_arr = flags.const_array(mfi);

      // Problem specific tagging
      const ProbParmDevice* lprobparm = d_prob_parm_device;
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
        geom.CellSizeArray();
      const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
        geom.ProbLoArray();
      const auto captured_level = level;
      amrex::ParallelFor(
        tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          ProblemSpecificFunctions::set_problem_tags(
            i, j, k, flag_arr, tag_arr, Sfab, tagval, dx, prob_lo, time,
            captured_level, *lprobparm);
        });
    }
  }

  // Untag cell close to EB, do this last
  if (
    eb_in_domain && (tagging_parm->eb_refine_type == "static") &&
    (level >= tagging_parm->max_eb_refine_lev)) {
    // Get distance function at current level
    const auto& ebfactory =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    amrex::MultiFab signDist(grids, dmap, 1, 0, amrex::MFInfo(), ebfactory);
    eb_distance(level, signDist);

    // Estimate how far I need to derefine
    const amrex::Real safetyFac = tagging_parm->detag_eb_factor;
    const auto& dx =
      parent->Geom(tagging_parm->max_eb_refine_lev).CellSizeArray();
    const amrex::Real dx_max = *std::max_element(dx.begin(), dx.end());
    amrex::Real clearTagDist = dx_max *
                               static_cast<amrex::Real>(parent->nErrorBuf(
                                 tagging_parm->max_eb_refine_lev)) *
                               safetyFac;
    const int finest_level = parent->finestLevel();
    for (int ilev = tagging_parm->max_eb_refine_lev + 1; ilev <= finest_level;
         ++ilev) {
      clearTagDist +=
        static_cast<amrex::Real>(parent->nErrorBuf(ilev)) * dx_max * safetyFac;
    }

    // Untag cells too close to EB
    const auto& dists = signDist.const_arrays();
    const auto& tagarrs = tags.arrays();
    amrex::ParallelFor(
      tags, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        if (dists[nbx](i, j, k) < clearTagDist) {
          tagarrs[nbx](i, j, k) = amrex::TagBox::CLEAR;
        }
      });
    amrex::Gpu::synchronize();
  }
}

std::unique_ptr<amrex::MultiFab>
PeleC::derive(const std::string& name, amrex::Real time, int ngrow)
{
  if ((do_les) && ((name == "C_s2") || (name == "C_I") || (name == "Pr_T"))) {
    AMREX_ASSERT(ngrow <= LES_Coeffs.nGrow());
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, ngrow));
    if (name == "C_s2") {
      amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_Cs2, 0, 1, ngrow);
    } else if (name == "C_I") {
      amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_CI, 0, 1, ngrow);
    } else if ((les_model != 1) && (name == "Pr_T")) {
      amrex::MultiFab::Copy(*derive_dat, LES_Coeffs, comp_PrT, 0, 1, ngrow);
    } else { // Pr_T, les_model==1
      amrex::MultiFab::Copy(
        *derive_dat, LES_Coeffs, comp_Cs2ovPrT, 0, 1, ngrow);
      amrex::MultiFab::Divide(*derive_dat, LES_Coeffs, comp_Cs2, 0, 1, ngrow);
    }
    return derive_dat;
  }

  if (name == "vfrac") {
    std::unique_ptr<amrex::MultiFab> mf(
      new amrex::MultiFab(grids, dmap, 1, ngrow, amrex::MFInfo()));
    if (ngrow > 0) {
      mf->setBndry(0);
    }
    amrex::MultiFab::Copy(*mf, vfrac, 0, 0, 1, 0);
    return mf;
  }

  // Can't use the AmrLevel derive for state variables with ghost cells:
  // That will fillpatch them individually, but we require the full state
  // for physBCs. Therefore, instead fillpatch the full state and then grab
  // the component we need.
  int index, scomp;
  if (isStateVariable(name, index, scomp) && (ngrow > 0)) {
    amrex::MultiFab S_data(
      get_new_data(State_Type).boxArray(),
      get_new_data(State_Type).DistributionMap(), NVAR, ngrow, amrex::MFInfo(),
      Factory());
    FillPatch(
      *this, S_data, S_data.nGrow(), time, State_Type, Density, NVAR, 0);
    std::unique_ptr<amrex::MultiFab> derive_dat(
      new amrex::MultiFab(grids, dmap, 1, ngrow, amrex::MFInfo(), Factory()));
    amrex::MultiFab::Copy(*derive_dat, S_data, scomp, 0, 1, ngrow);
    return derive_dat;
  }

  // For those using GrowBoxByOne we need this
  if ((name == "enstrophy") || (name == "magvort") || (name == "divu")) {
    ngrow += 1;
  }

  return AmrLevel::derive(name, time, ngrow);
}

void
PeleC::derive(
  const std::string& name, amrex::Real time, amrex::MultiFab& mf, int dcomp)
{
  if (name == "vfrac") {
    amrex::MultiFab::Copy(mf, vfrac, 0, dcomp, 1, 0);
  } else {
    AmrLevel::derive(name, time, mf, dcomp);
  }
}

void
PeleC::clear_prob()
{
  pc_prob_close();
}

void
PeleC::init_reactor()
{
  reactor = pele::physics::reactions::ReactorBase::create(chem_integrator);
  if (do_react && (chem_integrator == "ReactorNull")) {
    amrex::Print() << "WARNING: turning on reactions while using ReactorNull. "
                      "Make sure this is intended."
                   << std::endl;
  }
  reactor->init(1, 1);
}

void
PeleC::close_reactor()
{
  reactor->close();
}

void
PeleC::init_les()
{
  // Fill with default coefficient values
  LES_Coeffs.define(grids, dmap, nCompC, 1, amrex::MFInfo(), Factory());
  LES_Coeffs.setVal(0.0);
  if (les_model == 2) {
    LES_Coeffs.setVal(Cw * Cw, comp_Cs2, 1, LES_Coeffs.nGrow());
  } else if (les_model == 3) {
    LES_Coeffs.setVal(2.5 * Cs * Cs, comp_Cs2, 1, LES_Coeffs.nGrow());
  } else {
    LES_Coeffs.setVal(Cs * Cs, comp_Cs2, 1, LES_Coeffs.nGrow());
  }
  LES_Coeffs.setVal(CI, comp_CI, 1, LES_Coeffs.nGrow());
  if (les_model == 1) {
    LES_Coeffs.setVal(Cs * Cs * PrT, comp_Cs2ovPrT, 1, LES_Coeffs.nGrow());
  } else {
    LES_Coeffs.setVal(PrT, comp_PrT, 1, LES_Coeffs.nGrow());
  }

  amrex::Print() << "WARNING: LES with Fuego assumes Cp is a weak function of T"
                 << std::endl;
#if NUM_SPECIES > 2
  amrex::Abort("LES is not supported for multi-component systems");
#elif NUM_SPECIES == 2
  amrex::Print() << "WARNING: LES is not supported for multi-component systems"
                 << std::endl;
#endif
  if (std::is_same_v<
        pele::physics::PhysicsType::eos_type, pele::physics::eos::SRK>) {
    amrex::Abort("LES is not supported for non-ideal equations of state");
  }
}

void
PeleC::init_filters()
{
  if (level > 0) {
    amrex::IntVect ref_ratio = parent->refRatio(level - 1);
    les_filter = Filter(
      les_filter_type,
      static_cast<int>(les_filter_fgr * std::pow(ref_ratio[0], level)));
  } else {
    les_filter = Filter(les_filter_type, les_filter_fgr);
  }

  nGrowF = les_filter.get_filter_ngrow();

  // Add grow cells necessary for explicit filtering of source terms
  if (do_hydro) {
    Sborder.define(
      grids, dmap, NVAR, Sborder.nGrow() + nGrowF, amrex::MFInfo(), Factory());
    hydro_source.define(
      grids, dmap, NVAR, hydro_source.nGrow() + nGrowF, amrex::MFInfo(),
      Factory());
    sources_for_hydro.define(
      grids, dmap, NVAR, sources_for_hydro.nGrow() + nGrowF, amrex::MFInfo(),
      Factory());
  }

  volume.clear();
  volume.define(
    grids, dmap, 1, numGrow() + nGrowF, amrex::MFInfo(),
    amrex::FArrayBoxFactory());
  geom.GetVolume(volume);

  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    area[dir].clear();
    area[dir].define(
      getEdgeBoxArray(dir), dmap, 1, numGrow() + nGrowF, amrex::MFInfo(),
      amrex::FArrayBoxFactory());
    geom.GetFaceArea(area[dir], dir);
  }
}

void
PeleC::init_diagnostics()
{
  if (level == 0) {
    // initialize diagnostics (only level 0 calls them)
    const std::string pele_prefix = "pelec";
    amrex::ParmParse pp(pele_prefix);
    const int n_diags = pp.countval("diagnostics");
    amrex::Vector<std::string> diags(n_diags);
    for (int n = 0; n < n_diags; ++n) {
      pp.get("diagnostics", diags[n], n);
      const std::string diag_prefix = pele_prefix + "." + diags[n];
      amrex::ParmParse ppd(diag_prefix);
      std::string diag_type;
      ppd.get("type", diag_type);
      m_diagnostics.emplace_back(DiagBase::create(diag_type));
      m_diagnostics[n]->init(diag_prefix, diags[n]);
      m_diagnostics[n]->addVars(m_diagVars);
    }
    // Remove duplicates from m_diagVars and check that all the variables exist
    std::sort(m_diagVars.begin(), m_diagVars.end());
    auto last = std::unique(m_diagVars.begin(), m_diagVars.end());
    m_diagVars.erase(last, m_diagVars.end());
    int index = 0;
    int scomp = 0;
    for (auto& v : m_diagVars) {
      const bool itexists = derive_lst.canDerive(v) ||
                            isStateVariable(v, index, scomp) || v == "vfrac";
      if (!itexists) {
        amrex::Abort(
          "PeleC::init_diagnostics(): Field " + v + " is not available");
      }
    }
  }
}

void
PeleC::do_diagnostics()
{

  // Deal with Diagnostics
  if (level == 0) {

    // Timing info
    const int nstep = parent->levelSteps(0);
    const int finest_level = parent->finestLevel();
    const amrex::Real dtlev = parent->dtLevel(0);
    const amrex::Real cumtime = parent->cumTime() + dtlev;

    bool do_diags = false;
    for (const auto& m_diagnostic : m_diagnostics) {
      do_diags = do_diags || m_diagnostic->doDiag(cumtime, nstep);
    }

    if (do_diags) {

      // Need to update some internal data as the grid changes
      amrex::Vector<amrex::Geometry> geomAll(finest_level + 1);
      amrex::Vector<amrex::BoxArray> gridAll(finest_level + 1);
      amrex::Vector<amrex::DistributionMapping> dmapAll(finest_level + 1);
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto& amrlevel = parent->getLevel(lev);
        geomAll[lev] = amrlevel.Geom();
        gridAll[lev] = amrlevel.boxArray();
        dmapAll[lev] = amrlevel.DistributionMap();
      }
      for (const auto& m_diagnostic : m_diagnostics) {
        m_diagnostic->prepare(
          finest_level + 1, geomAll, gridAll, dmapAll, m_diagVars);
      }

      // Assemble a vector of MF containing the requested data
      amrex::Vector<std::unique_ptr<amrex::MultiFab>> diagMFVec(
        finest_level + 1);
      for (int lev = 0; lev <= finest_level; ++lev) {
        auto& amrlevel = parent->getLevel(lev);
        amrex::MultiFab S_data(
          amrlevel.get_new_data(State_Type).boxArray(),
          amrlevel.get_new_data(State_Type).DistributionMap(), NVAR, 1,
          amrex::MFInfo(), amrlevel.Factory());
        amrex::MultiFab R_data(
          amrlevel.get_new_data(Reactions_Type).boxArray(),
          amrlevel.get_new_data(Reactions_Type).DistributionMap(),
          NUM_SPECIES + 2, 1, amrex::MFInfo(), amrlevel.Factory());
        FillPatch(
          amrlevel, S_data, S_data.nGrow(), cumtime, State_Type, Density, NVAR,
          0);
        FillPatch(
          amrlevel, R_data, R_data.nGrow(), cumtime, Reactions_Type, 0,
          NUM_SPECIES + 2, 0);

        diagMFVec[lev] = std::make_unique<amrex::MultiFab>(
          amrlevel.boxArray(), amrlevel.DistributionMap(), m_diagVars.size(),
          1);
        for (int v{0}; v < m_diagVars.size(); ++v) {
          // Already tested: either a derive or a state variable
          if (derive_lst.canDerive(m_diagVars[v]) || m_diagVars[v] == "vfrac") {
            auto mf = amrlevel.derive(m_diagVars[v], cumtime, 1);
            int varIdx{0};
            if (derive_lst.canDerive(m_diagVars[v])) {
              const amrex::DeriveRec* rec = derive_lst.get(m_diagVars[v]);
              for (int vd{0}; vd < rec->numDerive(); ++vd) {
                if (m_diagVars[v] == rec->variableName(vd)) {
                  varIdx = vd;
                  break;
                }
              }
            }
            amrex::MultiFab::Copy(*diagMFVec[lev], *mf, varIdx, v, 1, 1);
          } else {
            int StIndex = 0;
            int scomp = 0;
            isStateVariable(m_diagVars[v], StIndex, scomp);
            if (StIndex == State_Type) {
              amrex::MultiFab::Copy(*diagMFVec[lev], S_data, scomp, v, 1, 1);
            } else if (StIndex == Reactions_Type) {
              amrex::MultiFab::Copy(*diagMFVec[lev], R_data, scomp, v, 1, 1);
            }
          }
        }
      }

      for (const auto& m_diagnostic : m_diagnostics) {
        if (m_diagnostic->doDiag(cumtime, nstep)) {
          m_diagnostic->processDiag(
            nstep, cumtime, amrex::GetVecOfConstPtrs(diagMFVec), m_diagVars);
        }
      }
    }
  }
}

#ifdef PELE_USE_MASA
void
PeleC::init_mms()
{
  if (!mms_initialized) {
    if (verbose && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "Initializing MMS" << std::endl;
    }
// Shut of FPE for MASA initialization because it has FPEs
#ifdef PELE_ENABLE_FPE_TRAP
#if defined(__linux__)
    unsigned int prev_fpe_excepts = fegetexcept();
    fedisableexcept(prev_fpe_excepts);
#elif defined(__APPLE__)
    static fenv_t prev_fpe_excepts;
    fegetenv(&prev_fpe_excepts);
    static fenv_t new_fpe_excepts;
    new_fpe_excepts.__control |= FE_ALL_EXCEPT;
    new_fpe_excepts.__mxcsr |= FE_ALL_EXCEPT << 7;
    fesetenv(&new_fpe_excepts);
#endif
#endif
    masa_init("mms", masa_solution_name.c_str());
    masa_set_param("Cs", PeleC::Cs);
    masa_set_param("CI", PeleC::CI);
    masa_set_param("PrT", PeleC::PrT);
    mms_initialized = true;
#ifdef PELE_ENABLE_FPE_TRAP
#if defined(__linux__)
    if (prev_fpe_excepts != 0) {
      feenableexcept(prev_fpe_excepts);
    }
#elif defined(__APPLE__)
    fesetenv(&prev_fpe_excepts);
#endif
#endif
  }
}
#endif

void
PeleC::reset_internal_energy(amrex::MultiFab& S_new, int ng)
{
#ifndef AMREX_USE_GPU
  amrex::Real sum0 = 0.0;
  if (parent->finestLevel() == 0 && print_energy_diagnostics) {
    // Pass in the multifab and the component
    sum0 = volWgtSumMF(S_new, Eden, true);
  }
#endif

  // Ensure (rho e) isn't too small or negative
  {
    const auto captured_allow_small_energy = allow_small_energy;
    const auto captured_allow_negative_energy = allow_negative_energy;
    const auto captured_dual_energy_update_E_from_e =
      dual_energy_update_E_from_e;
    const auto captured_dual_energy_eta2 = dual_energy_eta2;
    auto sarrs = S_new.arrays();
    const amrex::IntVect ngs(ng);

    // for Rotational frames
    if (do_rf) {
      const auto geomdata = geom.data();
      int axis = rf_axis;
      amrex::Real omega = rf_omega;
      amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
        AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
      amrex::ParallelFor(
        S_new, ngs,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
          amrex::IntVect iv(AMREX_D_DECL(i, j, k));
          amrex::Real rad = get_rotaxis_dist(iv, axis, axis_loc, geomdata);
          pc_rst_int_e(
            i, j, k, sarrs[nbx], captured_allow_small_energy,
            captured_allow_negative_energy,
            captured_dual_energy_update_E_from_e, captured_dual_energy_eta2,
            omega, rad);
        });
      amrex::Gpu::synchronize();
    } else {
      amrex::ParallelFor(
        S_new, ngs,
        [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
          pc_rst_int_e(
            i, j, k, sarrs[nbx], captured_allow_small_energy,
            captured_allow_negative_energy,
            captured_dual_energy_update_E_from_e, captured_dual_energy_eta2);
        });
      amrex::Gpu::synchronize();
    }
  }

#ifndef AMREX_USE_GPU
  if (parent->finestLevel() == 0 && print_energy_diagnostics) {
    // Pass in the multifab and the component
    amrex::Real sum = volWgtSumMF(S_new, Eden, true);
    amrex::ParallelDescriptor::ReduceRealSum(sum0);
    amrex::ParallelDescriptor::ReduceRealSum(sum);
    if (amrex::ParallelDescriptor::IOProcessor() && std::abs(sum - sum0) > 0) {
      amrex::Print() << "(rho E) added from reset terms                 : "
                     << sum - sum0 << " out of " << sum0 << std::endl;
    }
  }
#endif
}

void
PeleC::computeTemp(amrex::MultiFab& S, int ng)
{
  reset_internal_energy(S, ng);

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  auto const& sarrs = S.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::ParallelFor(
    S, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      if (!flagarrs[nbx](i, j, k).isCovered()) {
        pc_cmpTemp(i, j, k, sarrs[nbx]);
      }
    });
  amrex::Gpu::synchronize();
}

amrex::Real
PeleC::getCPUTime()
{
  int numCores = amrex::ParallelDescriptor::NProcs();
#ifdef AMREX_USE_OMP
  numCores = numCores * omp_get_max_threads();
#endif

  amrex::Real T =
    numCores * (amrex::ParallelDescriptor::second() - startCPUTime) +
    previousCPUTimeUsed;

  return T;
}

amrex::MultiFab&
PeleC::build_fine_mask()
{
  // Mask for zeroing covered cells
  AMREX_ASSERT(level > 0);

  if (!fine_mask.empty()) {
    return fine_mask;
  }

  const amrex::BoxArray& cba = parent->boxArray(level - 1);
  const amrex::DistributionMapping& cdm = parent->DistributionMap(level - 1);

  fine_mask.define(cba, cdm, 1, 0, amrex::MFInfo(), amrex::FArrayBoxFactory());
  fine_mask.setVal(1.0);

  amrex::BoxArray fba = parent->boxArray(level);
  amrex::iMultiFab ifine_mask = makeFineMask(cba, cdm, fba, crse_ratio, 1, 0);

  const auto& arrs = fine_mask.arrays();
  const auto& iarrs = ifine_mask.const_arrays();
  amrex::ParallelFor(
    fine_mask, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
#ifdef AMREX_USE_OMP
#pragma omp atomic write
#endif
      arrs[nbx](i, j, k) = iarrs[nbx](i, j, k);
    });
  amrex::Gpu::synchronize();
  return fine_mask;
}

const amrex::iMultiFab*
PeleC::build_interior_boundary_mask(int ng)
{
  for (const auto& ibm : ib_mask) {
    if (ibm->nGrow() == ng) {
      return ibm.get();
    }
  }

  // If we got here, we need to build a new one
  if (ib_mask.empty()) {
    ib_mask.resize(0);
  }

  ib_mask.push_back(
    std::make_unique<amrex::iMultiFab>(
      grids, dmap, 1, ng, amrex::MFInfo(),
      amrex::DefaultFabFactory<amrex::IArrayBox>()));

  amrex::iMultiFab* imf = ib_mask.back().get();
  int ghost_covered_by_valid = 0;
  int other_cells =
    1; // uncovered ghost, valid, and outside domain cells are set to 1

  imf->BuildMask(
    geom.Domain(), geom.periodicity(), ghost_covered_by_valid, other_cells,
    other_cells, other_cells);

  return imf;
}

amrex::Real
PeleC::clean_state(const amrex::MultiFab& /*S*/)
{
  // In the past, we enforced a minimum density and normalization of species.
  return 0.0;
}

amrex::Real
PeleC::clean_state(const amrex::MultiFab& /*S*/, amrex::MultiFab& /*S_old*/)
{
  // In the past, we enforced a minimum density and normalization of species.
  return 0.0;
}
