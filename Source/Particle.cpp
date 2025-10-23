#include "PeleC.H"
#include "SprayParticles.H"

namespace {
bool virtual_particles_set = false;

// Indices for spray source MultiFab
int sprayRhoSrcIndx = 0;
int sprayMomSrcIndx = 1;
int sprayEngSrcIndx = 1 + AMREX_SPACEDIM;
int spraySpecSrcIndx = 2 + AMREX_SPACEDIM;

int particle_verbose = 0;
} // namespace

std::unique_ptr<SprayParticleContainer> PeleC::SprayPC = nullptr;
std::unique_ptr<SprayParticleContainer> PeleC::VirtPC = nullptr;
std::unique_ptr<SprayParticleContainer> PeleC::GhostPC = nullptr;

// momentum + density + fuel species + energy
int PeleC::num_spray_src = AMREX_SPACEDIM + 2 + SPRAY_FUEL_NUM;

void
PeleC::estTimeStepParticles(amrex::Real& est_dt)
{
  if (!do_spray_particles) {
    return;
  }
  BL_PROFILE("PeleC::particleEstTimeStep()");
  amrex::Real est_dt_particle = SprayPC->estTimestep(level);

  if (est_dt_particle > 0) {
    est_dt = amrex::min<amrex::Real>(est_dt, est_dt_particle);
  }

  if (particle_verbose >= 2 && amrex::ParallelDescriptor::IOProcessor()) {
    if (est_dt_particle > 0.) {
      amrex::Print() << "...estdt from particles at level " << level << ": "
                     << est_dt_particle << '\n';
    } else {
      amrex::Print() << "...there are no particles at level " << level << '\n';
    }
  }
}

void
PeleC::readSprayParams()
{
  amrex::ParmParse pp("pelec");

  pp.query("do_spray_particles", do_spray_particles);
  if (do_spray_particles) {
    SprayParticleContainer::readSprayParams(particle_verbose);
  }
}

void
PeleC::defineParticles()
{
  if (!do_spray_particles) {
    return;
  }
  // There must be at least as many fuel species in the spray as
  // there are species in the fluid
  if (SPRAY_FUEL_NUM > NUM_SPECIES) {
    amrex::Abort("Cannot have more spray fuel species than fluid species");
  }
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ext_force = {0.0};
  for (int i = 0; i < static_cast<int>(external_forcing.size()); i++) {
    ext_force[i] = external_forcing[i];
  }
  SprayParticleContainer::spraySetup(ext_force.data());
  SprayComps scomps;
  scomps.rhoIndx = PeleC::Density;
  scomps.momIndx = PeleC::Xmom;
  scomps.engIndx = PeleC::Eden;
  scomps.utempIndx = PeleC::Temp;
  scomps.specIndx = PeleC::FirstSpec;
  scomps.rhoSrcIndx = sprayRhoSrcIndx;
  scomps.momSrcIndx = sprayMomSrcIndx;
  scomps.specSrcIndx = spraySpecSrcIndx;
  scomps.engSrcIndx = sprayEngSrcIndx;
  SprayParticleContainer::AssignSprayComps(scomps);
}

int
PeleC::sprayStateGhosts(const int amr_ncycle)
{
  if (!do_spray_particles) {
    return 0;
  }
  int finest_level = parent->finestLevel();
  return SprayParticleContainer::getStateGhostCells(
    level, finest_level, amr_ncycle);
}

void
PeleC::defineSpraySource(const int amr_ncycle)
{
  int finest_level = parent->finestLevel();
  int spray_source_ghosts = SprayParticleContainer::getSourceGhostCells(
    level, finest_level, amr_ncycle);
  // We must make a temporary spray source term to ensure number of ghost
  // cells are correct
  tmp_spray_source.define(
    grids, dmap, num_spray_src, spray_source_ghosts, amrex::MFInfo(),
    Factory());
}

void
PeleC::setupVirtualParticles(const int level, const int finest_level)
{
  BL_PROFILE("PeleC::setupVirtualParticles()");
  if (SprayPC != nullptr && !virtual_particles_set) {
    if (level < finest_level) {
      SprayParticleContainer::AoS virts;
      setupVirtualParticles(level + 1, finest_level);
      VirtPC->CreateVirtualParticles(level + 1, virts);
      VirtPC->AddParticlesAtLevel(virts, level);

      SprayPC->CreateVirtualParticles(level + 1, virts);
      VirtPC->AddParticlesAtLevel(virts, level);
    }
    virtual_particles_set = true;
  }
}

void
PeleC::removeVirtualParticles(const int level)
{
  if (VirtPC != nullptr) {
    VirtPC->RemoveParticlesAtLevel(level);
  }
  virtual_particles_set = false;
}

void
PeleC::setupGhostParticles(
  const int level, const int finest_level, const int ngrow)
{
  BL_PROFILE("PeleC::setupGhostParticles()");
  if (SprayPC != nullptr && level < finest_level) {
    SprayParticleContainer::AoS ghosts;
    SprayPC->CreateGhostParticles(level, ngrow, ghosts);
    GhostPC->AddParticlesAtLevel(ghosts, level + 1, ngrow);
  }
}

void
PeleC::removeGhostParticles(const int level)
{
  if (GhostPC != nullptr) {
    GhostPC->RemoveParticlesAtLevel(level);
  }
}

// Create new particle data
void
PeleC::createDataParticles()
{
  SprayPC = std::make_unique<SprayParticleContainer>(parent, phys_bc);
  SprayPC->SetVerbose(particle_verbose);
  VirtPC = std::make_unique<SprayParticleContainer>(parent, phys_bc);
  GhostPC = std::make_unique<SprayParticleContainer>(parent, phys_bc);
}

// Initialize the particles on the grid at level 0
void
PeleC::initParticles()
{
  BL_PROFILE("PeleC::initParticles()");

  if (level > 0) {
    return;
  }

  if (do_spray_particles) {
    AMREX_ASSERT(SprayPC == nullptr);
    createDataParticles();

    SprayPC->SprayInitialize();
  }
}

void
PeleC::postRestartParticles(const std::string& restart_file)
{
  if (level > 0) {
    defineSpraySource(parent->MaxRefRatio(level - 1));
    return;
  }

  if (do_spray_particles) {
    defineSpraySource(1);
    AMREX_ASSERT(SprayPC == nullptr);
    createDataParticles();

    SprayPC->SprayInitialize(restart_file);
    amrex::Gpu::Device::streamSynchronize();
  }
}

void
PeleC::particleMKD(
  amrex::Real time,
  amrex::Real dt,
  int sub_iteration,
  int /*sub_ncycle*/,
  int amr_ncycle)
{
  if (sub_iteration != 0) {
    return;
  }
  old_sources[spray_src]->setVal(0.);
  tmp_spray_source.setVal(0.);
  // Setup ghost particles for use in finer levels. Note that ghost
  // particles that will be used by this level have already been created,
  // the particles being set here are only used by finer levels.
  const int finest_level = parent->finestLevel();
  if (level < finest_level) {
    // Setup the virtual particles that represent particles on finer levels
    setupVirtualParticles(level, finest_level);
    // Make a copy of the particles on this level into ghost particles
    // for the finer level
    const int finer_ref = parent->MaxRefRatio(level);
    // Determine the number of ghost cells on the next level we need for making
    // ghost particles
    const int ghost_width = SprayParticleContainer::getGhostPartCells(
      level + 1, finest_level, finer_ref);
    setupGhostParticles(level, finest_level, ghost_width);
  }

  // Advance the particle velocities to the half-time and the positions to
  // the new time
  if (particle_verbose >= 1) {
    amrex::Print()
      << "moveKickDrift ... updating particle positions and velocity\n";
  }
  const int spray_state_ghosts = sprayStateGhosts(amr_ncycle);
  const int spray_source_ghosts = tmp_spray_source.nGrow();
  auto const* ltransparm = PeleC::trans_parms.device_parm();
  // We will make a temporary copy of the source term array inside
  // moveKickDrift and we are only going to use the spray force out to one ghost
  // cell so we need only define spray_force_old with one ghost cell

  AMREX_ASSERT(old_sources[spray_src]->nGrow() >= 1);

  // Do the valid particles themselves
  SprayPC->moveKickDrift(
    Sborder, tmp_spray_source, level, dt, time, false, false,
    spray_state_ghosts, spray_source_ghosts, true, ltransparm);

  // Only need the coarsest virtual particles here.
  if (level < finest_level) {
    VirtPC->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, true, false,
      spray_state_ghosts, spray_source_ghosts, true, ltransparm);
  }

  // Miiiight need all Ghosts
  if (GhostPC != nullptr && level != 0) {
    GhostPC->moveKickDrift(
      Sborder, tmp_spray_source, level, dt, time, false, true,
      spray_state_ghosts, spray_source_ghosts, true, ltransparm);
  }
  // Must call transfer source after moveKick and moveKickDrift
  // on all particle types
  SprayPC->transferSource(
    spray_source_ghosts, level, tmp_spray_source, *old_sources[spray_src]);
}

void
PeleC::particleMK(
  amrex::Real time,
  amrex::Real dt,
  int sub_iteration,
  int sub_ncycle,
  int amr_ncycle)
{
  if (sub_iteration != sub_ncycle - 1 && sub_ncycle != 0) {
    return;
  }
  const int spray_state_ghosts = sprayStateGhosts(amr_ncycle);
  const int spray_source_ghosts = tmp_spray_source.nGrow();
  new_sources[spray_src]->setVal(0.);
  auto const* ltransparm = PeleC::trans_parms.device_parm();
  if (particle_verbose >= 1) {
    amrex::Print() << "moveKick ... updating velocity only\n";
  }
  SprayPC->moveKick(
    Sborder, tmp_spray_source, level, dt, time, false, false,
    spray_state_ghosts, spray_source_ghosts, ltransparm);

  if (level < parent->finestLevel()) {
    VirtPC->moveKick(
      Sborder, tmp_spray_source, level, dt, time, true, false,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }

  if (GhostPC != nullptr && level != 0) {
    GhostPC->moveKick(
      Sborder, tmp_spray_source, level, dt, time, false, true,
      spray_state_ghosts, spray_source_ghosts, ltransparm);
  }
  SprayPC->transferSource(
    spray_source_ghosts, level, tmp_spray_source, *new_sources[spray_src]);
}

void
PeleC::postTimeStepParticles(int iteration)
{
  const int finest_level = parent->finestLevel();
  const int ncycle = parent->nCycle(level);
  if (do_spray_particles) {
    // Remove virtual particles at this level if we have any.
    if (VirtPC != nullptr) {
      removeVirtualParticles(level);
    }

    // Remove Ghost particles on the final iteration
    if (iteration == ncycle) {
      removeGhostParticles(level);
    }

    // Do particle injection
    const int nstep = parent->levelSteps(0);
    const amrex::Real dtlev = parent->dtLevel(0);
    const amrex::Real cumtime = parent->cumTime() + dtlev;
    BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
    const bool injectParts =
      SprayPC->injectParticles(cumtime, dtlev, nstep, level, finest_level);
    BL_PROFILE_VAR_STOP(INJECT_SPRAY);
    // Sync up if we're level 0 or if we have particles that may have moved
    // off the next finest level and need to be added to our own level, or
    // if we injected particles
    if (
      (iteration < ncycle && level < finest_level) || level == 0 ||
      injectParts) {
      // TODO: Determine how many ghost cells to use here
      const int nGrow = iteration;
      SprayPC->Redistribute(level, SprayPC->finestLevel(), nGrow);
    }
  }
}

void
PeleC::postInitParticles()
{
  const amrex::Real dtlev = parent->dtLevel(level);
  const amrex::Real cumtime = parent->cumTime();
  if (do_spray_particles) {
    BL_PROFILE_VAR("SprayParticles::injectParticles()", INJECT_SPRAY);
    const bool injectParts =
      SprayPC->injectParticles(cumtime, dtlev, 0, level, parent->finestLevel());
    BL_PROFILE_VAR_STOP(INJECT_SPRAY);
    if (injectParts) {
      SprayPC->Redistribute(level, SprayPC->finestLevel(), 0);
    }
  }
}

void
PeleC::particle_redistribute(int lbase, bool init_part)
{
  if (!do_spray_particles) {
    return;
  }
  BL_PROFILE("PeleC::particle_redistribute()");
  int flev = parent->finestLevel();
  if (SprayPC != nullptr) {
    // If we are calling with init_part = true, then we want to force the
    // redistribute without checking whether the grids have changed.
    if (init_part) {
      SprayPC->Redistribute(lbase);
      return;
    }

    // These are usually the BoxArray and DMap from the last regridding.
    static amrex::Vector<amrex::BoxArray> ba;
    static amrex::Vector<amrex::DistributionMapping> dm;

    bool changed = false;

    while (parent->getAmrLevels()[flev] == nullptr) {
      flev--;
    }

    if (ba.size() != flev + 1) {
      ba.resize(flev + 1);
      dm.resize(flev + 1);
      changed = true;
    } else {
      for (int i = 0; i <= flev && !changed; i++) {
        // Check if BoxArrays have changed during regridding
        if (ba[i] != parent->boxArray(i)) {
          changed = true;
        }

        if (!changed) {
          // Check DistributionMaps have changed during regridding
          if (dm[i] != parent->getLevel(i).get_new_data(0).DistributionMap()) {
            changed = true;
          }
        }
      }
    }

    if (changed) {
      // We only need to call Redistribute if the BoxArrays or DistMaps have
      // changed.
      // We also only call it for particles >= lbase. This is
      // because if we called redistribute during a subcycle, there may be
      // particles not in the proper position on coarser levels.
      if (particle_verbose >= 2 && amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print() << "Calling redistribute because grid has changed "
                       << '\n';
      }
      if (flev == 0) {
        // Do a local redistribute
        SprayPC->Redistribute(lbase, SprayPC->finestLevel(), 0, 1);
      } else {
        SprayPC->Redistribute(lbase, SprayPC->finestLevel(), 0);
      }

      // Use the new BoxArray and DistMap to define ba and dm for next time.
      for (int i = 0; i <= flev; i++) {
        ba[i] = parent->boxArray(i);
        dm[i] = parent->getLevel(i).get_new_data(0).DistributionMap();
      }
    } else {
      if (particle_verbose >= 2 && amrex::ParallelDescriptor::IOProcessor()) {
        amrex::Print()
          << "NOT calling redistribute because grid has NOT changed " << '\n';
      }
    }
  }
}
