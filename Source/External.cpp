#include "PeleC.H"
#include "IndexDefines.H"
#include "prob.H"

void
PeleC::construct_old_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);

  int ng = 0; // None filled

  old_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_old, *old_sources[ext_src], ng);

  old_sources[ext_src]->FillBoundary(geom.periodicity());
}

void
PeleC::construct_new_ext_source(amrex::Real time, amrex::Real dt)
{
  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);

  int ng = 0;

  new_sources[ext_src]->setVal(0.0);

  if (!add_ext_src) {
    return;
  }

  fill_ext_source(time, dt, S_old, S_new, *new_sources[ext_src], ng);
}

void
PeleC::fill_ext_source(
  amrex::Real time,
  amrex::Real dt,
  const amrex::MultiFab& state_old,
  const amrex::MultiFab& state_new,
  amrex::MultiFab& ext_src,
  int ng)
{
  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(state_old.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  auto const& Sns = state_new.const_arrays();
  auto const& Farrs = ext_src.arrays();
  auto const& flagarrs = flags.const_arrays();
  const amrex::IntVect ngs(ng);
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> ext_force = {0.0};
  for (int i = 0; i < static_cast<int>(external_forcing.size()); i++) {
    ext_force[i] = external_forcing[i];
  }
  amrex::ParallelFor(
    ext_src, ngs, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      if (!flagarrs[nbx](i, j, k).isCovered()) {
        amrex::Real e_force = 0.0;
        for (int idir = 0; idir < AMREX_SPACEDIM; idir++) {
          Farrs[nbx](i, j, k, UMX + idir) = ext_force[idir];
          e_force += Sns[nbx](i, j, k, UMX + idir) / Sns[nbx](i, j, k, URHO) *
                     ext_force[idir];
        }
        Farrs[nbx](i, j, k, UEDEN) = e_force;
      }
    });
  amrex::Gpu::synchronize();

  // problem_modify_ext_sources runs on the host, so give it the host copy of
  // the problem parameters. The device copy is not readable from the host.
  const ProbParmDevice* lprobparm = PeleC::h_prob_parm_device;
  const auto geomdata = geom.data();
  ProblemSpecificFunctions::problem_modify_ext_sources(
    time, dt, state_old, state_new, ext_src, ng, geomdata, *lprobparm);
}
