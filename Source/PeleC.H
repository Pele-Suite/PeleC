#ifndef PELE_H
#define PELE_H

#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrLevel.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBFluxRegister.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_MultiCutFab.H>

#ifdef AMREX_PARTICLES
#include <AMReX_Particles.H>
#endif

#include "Filter.H"
#include "PeleCUtilities.H"
#include "Tagging.H"
#include "IndexDefines.H"
#include "prob_parm.H"
#include "PelePhysics.H"
#include "ReactorBase.H"
#include "turbinflow.H"
#include "SparseData.H"
#include "EBStencilTypes.H"
#include "DiagBase.H"

enum StateType { State_Type = 0, Reactions_Type, Work_Estimate_Type };

// Create storage for all source terms.

enum sources {
  ext_src = 0,
  forcing_src,
  spray_src,
  diff_src,
  les_src,
  soot_src,
  rot_src,
#ifdef PELE_USE_MASA
  mms_src,
#endif
  num_src
};

// These match AMReX_BC_TYPES.H, except user_bc is added
namespace PCPhysBCType {
enum phys_bc_type {
  interior = 0,
  inflow,
  outflow,
  symmetry,
  slip_wall,
  no_slip_wall,
  user_bc
};
}

// Forward declarations
#ifdef PELE_USE_SOOT
class SootModel;
#endif
#ifdef PELE_USE_SPRAY
class SprayParticleContainer;
#endif

AMREX_FORCE_INLINE static bool
ebInDomain()
{
  amrex::ParmParse pp("eb2");
  std::string geom_type("all_regular");
  pp.query("geom_type", geom_type);
  return (geom_type != "all_regular");
}

// AmrLevel-derived class for hyperbolic conservation equations for stellar
// media

class PeleC : public amrex::AmrLevel
{
public:
  // Default constructor.  Builds invalid object.
  PeleC();

  PeleC(
    amrex::Amr& papa,
    int lev,
    const amrex::Geometry& level_geom,
    const amrex::BoxArray& bl,
    const amrex::DistributionMapping& dm,
    amrex::Real time);

  ~PeleC() override;

  // Restart from a checkpoint file.
  void restart(
    amrex::Amr& papa, std::istream& is, bool bReadSpecial = false) override;

  // This is called only when we restart from an old checkpoint.
  void
  set_state_in_checkpoint(amrex::Vector<int>& state_in_checkpoint) override;
  bool check_state_in_checkpoint(const StateType state_type);

  // Call AmrLevel::checkPoint and then add radiation info
  void checkPoint(
    const std::string& dir,
    std::ostream& os,
    amrex::VisMF::How how,
    bool dump_old) override;

  void setPlotVariables() override;

  // Write a plotfile to specified directory.
  void writePlotFilePost(const std::string& dir, std::ostream& os) override;
  void writeJobInfo(const std::string& dir);
  static void writeBuildInfo(std::ostream& os);

  // Initialize level from plt file
  void initLevelDataFromPlt(
    const int lev, const std::string& pltFile, amrex::MultiFab& S_new);

  // Define data descriptors.
  static void variableSetUp();

  // Cleanup data descriptors at end of run.
  static void variableCleanUp();

  // Initialize grid data at problem start-up.
  void initData() override;

#ifdef PELE_USE_SPRAY

  // Read spray related inputs
  static void readSprayParams();

  // Define particle related variables
  static void defineParticles();

  // Initialize particle locations and velocities (and strengths if relevant)
  void initParticles();

  // Create relevant particle data
  void createDataParticles();

  // How to initialize at restart
  void postRestartParticles(const std::string& restart_file);

  // Post timestep particle functions
  void postTimeStepParticles(int iteration);

  // Post init for spray particles
  void postInitParticles();

  // Redistribute
  void particle_redistribute(int lbase = 0, bool init_part = false) override;

  // Setup ghost and virtual particles and do moveKickDrift
  void particleMKD(
    amrex::Real time,
    amrex::Real dt,
    int sub_iteration,
    int sub_ncycle,
    int amr_ncycle);

  // Do the moveKick for active and ghost particles
  void particleMK(
    amrex::Real time,
    amrex::Real dt,
    int sub_iteration,
    int sub_ncycle,
    int amr_ncycle);

  // Setup virtual particles if necessary
  static void setupVirtualParticles(const int level, const int finest_level);

  // Remove virtual particles if necessary
  static void removeVirtualParticles(const int level);

  // Setup ghost particles (for finer levels) if necessary
  static void
  setupGhostParticles(const int level, const int finest_level, const int ngrow);

  // Remove ghost particles (for this level) if necessary
  static void removeGhostParticles(const int level);

  // Time step control based on particles
  void estTimeStepParticles(amrex::Real& est_dt);

  // Initialize the temporary spray source
  void defineSpraySource(int amr_ncycle);
  amrex::MultiFab tmp_spray_source;

  // Determine the number of ghost cells for the state MF
  int sprayStateGhosts(int amr_ncycle);

  // Number of source terms in temporary spray data
  static int num_spray_src;

  static std::unique_ptr<SprayParticleContainer> SprayPC;
  static std::unique_ptr<SprayParticleContainer> VirtPC;
  static std::unique_ptr<SprayParticleContainer> GhostPC;
#endif

  // Get problem metrics.
  amrex::MultiFab* Area();
  amrex::MultiFab& Area(int dir);

  // The volume of the multifab.
  amrex::MultiFab& Volume();

  // Set time levels of state data.
  void setTimeLevel(
    amrex::Real time, amrex::Real dt_old, amrex::Real dt_new) override;

  // Initialize data on this level from another PeleC (during regrid).
  void init(amrex::AmrLevel& old) override;

  // Initialize data on this level after regridding if old level did not
  // previously exist
  void init() override;

  // Initialize EB geometry for finest_level and level grids for
  // other levels for the Amr class to do timed load balances.
  int WorkEstType() override { return Work_Estimate_Type; }

  const amrex::MultiFab& volFrac() const { return vfrac; }

  void init_eb();

  void initialize_eb2_structs();

  void define_body_state();

  void set_body_state(amrex::MultiFab& S);

  void initialize_signed_distance();

  void eb_distance(const int lev, amrex::MultiFab& signDistLev);

  void
  extend_signed_distance(amrex::MultiFab* signDist, amrex::Real extendFactor);

  void set_typical_values_chem();

  // Proceed with next timestep?
  int okToContinue() override;

  // Advance grids at this level in time.
  amrex::Real
  advance(amrex::Real time, amrex::Real dt, int iteration, int ncycle) override;

  amrex::Real do_mol_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  amrex::Real do_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  void initialize_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  static void finalize_sdc_advance(
    amrex::Real time, amrex::Real dt, int amr_iteration, int amr_ncycle);

  amrex::Real do_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  static void initialize_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sdc_iteration,
    int sdc_ncycle);

  static void finalize_sdc_iteration(
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sdc_iteration,
    int sdc_ncycle);

  void construct_Snew(
    amrex::MultiFab& S_new, const amrex::MultiFab& S_old, amrex::Real dt);

  void construct_hydro_source(
    const amrex::MultiFab& S,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void construct_old_source(
    int src,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void construct_new_source(
    int src,
    amrex::Real time,
    amrex::Real dt,
    int amr_iteration,
    int amr_ncycle,
    int sub_iteration,
    int sub_ncycle);

  void sum_of_sources(amrex::MultiFab& source);

  void construct_old_ext_source(amrex::Real time, amrex::Real dt);

  void construct_new_ext_source(amrex::Real time, amrex::Real dt);

  void fill_ext_source(
    amrex::Real time,
    amrex::Real dt,
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& state_new,
    amrex::MultiFab& ext_src,
    int ng);

  void construct_old_forcing_source(amrex::Real time, amrex::Real dt);

  void construct_new_forcing_source(amrex::Real time, amrex::Real dt);

  static void fill_forcing_source(
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& state_new,
    amrex::MultiFab& forcing_src,
    int ng);

  static bool add_soot_src;
#ifdef PELE_USE_SOOT
  static void setSootIndx();

  static void addSootDerivePlotVars(
    amrex::DeriveList& derive_lst, const amrex::DescriptorList& desc_lst);

  void construct_old_soot_source(amrex::Real time, amrex::Real dt);

  void construct_new_soot_source(amrex::Real time, amrex::Real dt);

  void fill_soot_source(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& state,
    amrex::MultiFab& soot_src,
    int ng);

  static void clipSootMoments(amrex::MultiFab& S_new, const int ng);

  static bool plot_soot;
#endif

#ifdef PELE_USE_MASA
  void construct_old_mms_source(amrex::Real time);

  void construct_new_mms_source(amrex::Real time);

  void fill_mms_source(
    amrex::Real time,
    const amrex::MultiFab& S,
    amrex::MultiFab& mms_src,
    int ng);
#endif

  void construct_old_rot_source(amrex::Real time, amrex::Real dt);

  void construct_new_rot_source(amrex::Real time, amrex::Real dt);

  void fill_rot_source(
    const amrex::MultiFab& state_old,
    const amrex::MultiFab& state_new,
    amrex::MultiFab& rot_src,
    int ng);

  void getLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void updateFluxRegistersLES(
    amrex::Real reflux_factor,
    const amrex::MultiFab& LESTerm,
    amrex::Real dt,
    const amrex::MFIter& mfi,
    amrex::FabType typ,
    const std::array<amrex::FArrayBox const*, AMREX_SPACEDIM>& flux_ec);

  void getSmagorinskyLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void getDynamicSmagorinskyLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void getWALELESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void getVremanLESTerm(
    amrex::Real time,
    amrex::Real dt,
    amrex::MultiFab& LESTerm,
    amrex::Real flux_factor);

  void construct_old_les_source(
    amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

  void construct_new_les_source(
    amrex::Real time, amrex::Real dt, int sub_iteration, int sub_ncycle);

  static void set_active_sources();

  // Estimate time step.
  amrex::Real estTimeStep(amrex::Real dt_old);

  // Compute initial time step.
  amrex::Real initialTimeStep();

  // Compute initial `dt'.
  void computeInitialDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time) override;

  // Compute new `dt'.
  void computeNewDt(
    int finest_level,
    int sub_cycle,
    amrex::Vector<int>& n_cycle,
    const amrex::Vector<amrex::IntVect>& ref_ratio,
    amrex::Vector<amrex::Real>& dt_min,
    amrex::Vector<amrex::Real>& dt_level,
    amrex::Real stop_time,
    int post_regrid_flag) override;

  // Allocate data at old time.
  void allocOldData() override;

  // Remove data at old time.
  void removeOldData() override;

  // Print information about energy budget.
  void do_energy_diagnostics();

  // Do work after timestep().
  void post_timestep(int iteration) override;

  // Contains operations to be done only after a full coarse timestep.
  void postCoarseTimeStep(amrex::Real cumtime) override;

  // Do work after regrid().
  void post_regrid(int lbase, int new_finest) override;

  // Do work after a restart().
  void post_restart() override;

  // Do work after init().
  void post_init(amrex::Real stop_time) override;

  // Error estimation for regridding.
  void errorEst(
    amrex::TagBoxArray& tags,
    int clearval,
    int tagval,
    amrex::Real time,
    int n_error_buf = 0,
    int ngrow = 0) override;

  // Returns a MultiFab containing the derived data for this level.
  // The user is responsible for deleting this pointer when done
  // with it.  If ngrow>0 the MultiFab is built on the appropriately
  // grown BoxArray.
  std::unique_ptr<amrex::MultiFab>
  derive(const std::string& name, amrex::Real time, int ngrow) override;

  // This version of derive() fills the dcomp'th component of mf with the
  // derived quantity.
  void derive(
    const std::string& name,
    amrex::Real time,
    amrex::MultiFab& mf,
    int dcomp) override;

  static int numGrow();

  void react_state(
    amrex::Real time,
    amrex::Real dt,
    bool init = false,
    amrex::MultiFab* aux_src = nullptr);

  void reset_internal_energy(amrex::MultiFab& S_new, int ng);

  void computeTemp(amrex::MultiFab& State, int ng);

  void getMOLSrcTerm(
    const amrex::MultiFab& S,
    amrex::MultiFab& MOLSrcTerm,
    amrex::Real time,
    amrex::Real dt,
    amrex::Real flux_factor);

  void enforce_consistent_e(amrex::MultiFab& S);

  amrex::Real volWgtSum(
    const std::string& name,
    amrex::Real time,
    bool local = false,
    bool finemask = true);
  amrex::Real volWgtSquaredSum(
    const std::string& name, amrex::Real time, bool local = false);
  amrex::Real volWgtSumMF(
    const amrex::MultiFab& mf,
    int comp,
    bool local = false,
    bool finemask = false);
  amrex::Real
  volWgtSquaredSumDiff(int comp, amrex::Real time, bool local = false);
  amrex::Real
  maxDerive(const std::string& name, amrex::Real time, bool local = false);
  amrex::Real
  minDerive(const std::string& name, amrex::Real time, bool local = false);

  // derives that need variables part of this class (e.g. trans_parm)
  static void pc_derviscosity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derbulkviscosity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derconductivity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_derdiffusivity(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  template <int dir>
  static void pc_dervel_if(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_dermagvel_if(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);
  static void pc_dereint1(
    const amrex::Box& bx,
    amrex::FArrayBox& derfab,
    int dcomp,
    int ncomp,
    const amrex::FArrayBox& datfab,
    const amrex::Geometry& geomdata,
    amrex::Real time,
    const int* bcrec,
    const int level);

  static int Density, Xmom, Ymom, Zmom, Eden, Eint, Temp;

  static int FirstAdv;
  static int FirstSpec;
  static int FirstAux;
  static int FirstLin;

  /// Index locations for soot modeling components
  static int FirstSootVar, NumSootVars;

  // This MultiFab is on the coarser level.  This is useful for the coarser
  // level to mask out the finer level.  We only build this when it is needed.
  // This coarse MultiFab has to live on the fine level because it must be
  // updated even when only the fine level changes.
  amrex::MultiFab fine_mask;
  amrex::MultiFab& build_fine_mask();

  static bool eb_in_domain;

  std::unique_ptr<amrex::EBFluxRegister> flux_reg;
  amrex::EBFluxRegister& getFluxReg() const;
  amrex::EBFluxRegister& getFluxReg(int lev);

  static bool ebInitialized();
  static int getEBMaxLevel();
  static int getEBCoarsening();

  void InitialRedistribution(
    const amrex::Real time,
    const amrex::Vector<amrex::BCRec> bcs,
    amrex::MultiFab& S_new);

  void avgDown();
  void avgDown(int state_indx);

  static ProbParmDevice* h_prob_parm_device;
  static ProbParmDevice* d_prob_parm_device;
  static ProbParmHost* prob_parm_host;
  static TaggingParm* tagging_parm;
  static pele::physics::PeleParams<pele::physics::transport::TransParm<
    pele::physics::PhysicsType::eos_type,
    pele::physics::PhysicsType::transport_type>>
    trans_parms;
  static pele::physics::turbinflow::TurbInflow turb_inflow;

  // A set of runtime diagnostics from PelePhysics lib
  amrex::Vector<std::unique_ptr<DiagBase>> m_diagnostics;
  static amrex::Vector<std::string> m_diagVars;

#ifdef PELE_USE_SOOT
  static SootModel soot_model;
#endif

protected:
  amrex::iMultiFab level_mask;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> areafrac;
  std::array<const amrex::MultiCutFab*, AMREX_SPACEDIM> facecent;

  // Build a mask that ghost cells overlapping with interior cells in the same
  // multifab are set to 0, whereas others are set to 1.
  amrex::Vector<std::unique_ptr<amrex::iMultiFab>> ib_mask;
  const amrex::iMultiFab* build_interior_boundary_mask(int ng);

  // A state array with ghost zones.
  amrex::MultiFab Sborder;

  // Source terms to the hydrodynamics solve.
  amrex::MultiFab sources_for_hydro;

  // Source term representing hydrodynamics update.
  amrex::MultiFab hydro_source;

  // Non-hydro source terms.
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> old_sources;
  amrex::Vector<std::unique_ptr<amrex::MultiFab>> new_sources;

  std::unique_ptr<pele::physics::reactions::ReactorBase> reactor;
  void init_reactor();
  void close_reactor();

  void init_les();
  void init_filters();

  void init_diagnostics();
  void do_diagnostics();

#ifdef PELE_USE_MASA
  static void init_mms();
#endif

  static void clear_prob();

  static void read_params();

  static void check_params();

  static void read_tagging_params();

  PeleC& getLevel(int lev);

  void update_flux_registers(
    const amrex::Real dt,
    const amrex::MFIter& mfi,
    const amrex::FabType& typ,
    const std::array<amrex::FArrayBox const*, AMREX_SPACEDIM>& flux,
    const amrex::FArrayBox& dm_as_fine);

  void reflux();

  static amrex::Real clean_state(const amrex::MultiFab& S);

  static amrex::Real
  clean_state(const amrex::MultiFab& S, amrex::MultiFab& S_old);

  void buildMetrics();

  // integrate derived quantities over domain

  amrex::Real
  sumDerive(const std::string& name, amrex::Real time, bool local = false);

  int find_datalog_index(const std::string& logname);

  void sum_integrated_quantities();

  void monitor_extrema();

  void write_info();

  static void stopJob();

  // The data.
  amrex::MultiFab volume;
  amrex::MultiFab area[AMREX_SPACEDIM];
  amrex::MultiFab dLogArea[1];
  amrex::Vector<amrex::Vector<amrex::Real>> radius;

  // Static data members.
#include "pelec_params.H"

  static bool signalStopJob;
  static int radius_grow;
  static int verbose;
  static amrex::BCRec phys_bc;

  // for keeping track of mass changes from negative density resets
  static amrex::Real frac_change;

  // for keeping track of the amount of CPU time used -- this will persist
  // after restarts
  static amrex::Real previousCPUTimeUsed;
  static amrex::Real startCPUTime;

  static amrex::Real getCPUTime();

  static bool do_diffuse;

  static int les_model;
  static int les_filter_type;
  static int les_filter_fgr;
  Filter les_filter;
  int nGrowF;
  static int les_test_filter_type;
  static int les_test_filter_fgr;
  amrex::MultiFab LES_Coeffs;
  amrex::MultiFab filtered_les_source;

#ifdef PELE_USE_MASA
  static bool mms_initialized;
  bool mms_src_evaluated;
  amrex::MultiFab mms_source;
#endif

  // There can be only one Diffusion object, it covers all levels:
  static class Diffusion* diffusion;

  // Do spray particles
  static bool do_spray_particles;

  static int num_state_type;

  static amrex::Vector<std::string> spec_names;
  static amrex::Vector<std::string> adv_names;
  static amrex::Vector<std::string> aux_names;

  static amrex::Vector<int> src_list;

  static bool use_typical_vals_chem;
  static bool use_typical_vals_chem_usr;
  static amrex::Real typical_rhoY_val_min;
  static int reset_typical_vals_int;
  static amrex::Vector<amrex::Real> typical_values_chem_usr;

  void problem_post_timestep();
  void problem_post_restart();
  void problem_post_init();

  static amrex::GpuArray<amrex::Real, NVAR> body_state;
  static bool body_state_set;
  static bool eb_initialized;
  static int eb_max_lvl_gen;

  amrex::MultiFab vfrac;

  amrex::Vector<amrex::Gpu::DeviceVector<EBBndryGeom>> sv_eb_bndry_geom;
  amrex::Vector<amrex::Gpu::DeviceVector<EBBndrySten>> sv_eb_bndry_grad_stencil;
  amrex::
    GpuArray<amrex::Vector<amrex::Gpu::DeviceVector<FaceSten>>, AMREX_SPACEDIM>
      flux_interp_stencil;

  amrex::Vector<SparseData<amrex::Real, EBBndrySten>> sv_eb_flux;
  amrex::Vector<SparseData<amrex::Real, EBBndrySten>> sv_eb_bcval;

  amrex::MultiFab signed_dist_0;
  static bool do_react_load_balance;
  static bool do_mol_load_balance;
};

void pc_bcfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

void pc_reactfill_hyp(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

void pc_nullfill(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  const int numcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec>& bcr,
  const int bcomp,
  const int scomp);

// Inlines

AMREX_FORCE_INLINE
int
PeleC::numGrow()
{
  int ng = 4;
  if (eb_in_domain) {
    // IF MOL     + FRD then 5
    // IF Godunov + FRD then 6
    // IF MOL     + SRD then 6
    // IF Godunov + SRD then 7
    if (do_mol) {
      if (redistribution_type == "FluxRedist") {
        ng += 1;
      } else if (redistribution_type == "StateRedist") {
        ng += 2;
      }
    } else {
      if (redistribution_type == "FluxRedist") {
        ng += 2;
      } else if (redistribution_type == "StateRedist") {
        ng += 3;
      }
    }
  }
  return ng;
}

AMREX_FORCE_INLINE
amrex::MultiFab*
PeleC::Area()
{
  return area;
}

AMREX_FORCE_INLINE
amrex::MultiFab&
PeleC::Area(int dir)
{
  return area[dir];
}

AMREX_FORCE_INLINE
amrex::MultiFab&
PeleC::Volume()
{
  return volume;
}

AMREX_FORCE_INLINE
PeleC&
PeleC::getLevel(int lev)
{
  return *(PeleC*)&parent->getLevel(lev);
}

AMREX_FORCE_INLINE
void
PeleC::stopJob()
{
  signalStopJob = true;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_check_initial_species(
  const int i, const int j, const int k, amrex::Array4<amrex::Real> const& Sfab)
{
  amrex::Real spec_sum = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    spec_sum = spec_sum + Sfab(i, j, k, UFS + n);
  }
  if (std::abs(Sfab(i, j, k, URHO) - spec_sum) > 1.e-8 * Sfab(i, j, k, URHO)) {
    // print *,'Sum of (rho X)_i vs rho at (i,j,k):
    // ',i,j,k,spec_sum,state(i,j,k,URHO)
    amrex::Abort("Error:: Failed check of initial species summing to 1");
  }
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getFluxReg() const
{
  return *flux_reg;
}

AMREX_FORCE_INLINE
amrex::EBFluxRegister&
PeleC::getFluxReg(int lev)
{
  return getLevel(lev).getFluxReg();
}

#endif
