#include "PeleC.H"

namespace {
// A derived quantity can be NaN or Inf in a covered cell, for example when it
// divides by a density that was zeroed in the body. Multiplying by a zero
// volume fraction does not remove those values, so zero them explicitly.
void
zero_covered_cells(amrex::MultiFab& mf, const amrex::MultiFab& vfrac)
{
  auto const& mf_arrs = mf.arrays();
  auto const& vf_arrs = vfrac.const_arrays();
  amrex::ParallelFor(
    mf, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      if (vf_arrs[nbx](i, j, k) <= 0.0) {
        mf_arrs[nbx](i, j, k) = 0.0;
      }
    });
  amrex::Gpu::synchronize();
}
} // namespace

amrex::Real
PeleC::sumDerive(const std::string& name, amrex::Real time, bool local)
{
  if (eb_in_domain) {
    amrex::Abort("sumDerive undefined for EB");
  }

  auto mf = derive(name, time, 0);

  AMREX_ASSERT(!(mf == nullptr));

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  return mf->sum(0, local);
}

amrex::Real
PeleC::volWgtSum(
  const std::string& name, amrex::Real time, bool local, bool finemask)
{
  BL_PROFILE("PeleC::volWgtSum()");

  auto mf = derive(name, time, 0);

  AMREX_ASSERT(mf != nullptr);

  if (level < parent->finestLevel() && finemask) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  if (eb_in_domain) {
    zero_covered_cells(*mf, vfrac);
    amrex::MultiFab::Multiply(*mf, vfrac, 0, 0, 1, 0);
  }

  // Dot does the parallel reduction unless local is true
  return amrex::MultiFab::Dot(*mf, 0, volume, 0, 1, 0, local);
}

amrex::Real
PeleC::volWgtSquaredSum(const std::string& name, amrex::Real time, bool local)
{
  BL_PROFILE("PeleC::volWgtSquaredSum()");

  auto mf = derive(name, time, 0);

  AMREX_ASSERT(mf != nullptr);

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(*mf, mask, 0, 0, 1, 0);
  }

  if (eb_in_domain) {
    zero_covered_cells(*mf, vfrac);
  }

  amrex::MultiFab::Multiply(*mf, *mf, 0, 0, 1, 0);

  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, volume, 0, 0, 1, 0);
  if (eb_in_domain) {
    amrex::MultiFab::Multiply(vol, vfrac, 0, 0, 1, 0);
  }
  return amrex::MultiFab::Dot(*mf, 0, vol, 0, 1, 0, local);
}

amrex::Real
PeleC ::volWgtSquaredSumDiff(int comp, amrex::Real /*time*/, bool local)
{
  // Calculate volume weighted sum of the square of the difference
  // between the old and new quantity

  const amrex::MultiFab& S_old = get_old_data(State_Type);
  const amrex::MultiFab& S_new = get_new_data(State_Type);
  amrex::MultiFab diff(grids, dmap, 1, 0);

  // Calculate the difference between the states
  amrex::MultiFab::Copy(diff, S_old, comp, 0, 1, 0);
  amrex::MultiFab::Subtract(diff, S_new, comp, 0, 1, 0);

  if (level < parent->finestLevel()) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(diff, mask, 0, 0, 1, 0);
  }

  amrex::MultiFab::Multiply(diff, diff, 0, 0, 1, 0);

  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, volume, 0, 0, 1, 0);
  if (eb_in_domain) {
    amrex::MultiFab::Multiply(vol, vfrac, 0, 0, 1, 0);
  }
  return amrex::MultiFab::Dot(diff, 0, vol, 0, 1, 0, local);
}

amrex::Real
PeleC::volWgtSumMF(
  const amrex::MultiFab& mf, int comp, bool local, bool finemask)
{
  BL_PROFILE("PeleC::volWgtSumMF()");

  amrex::MultiFab vol(grids, dmap, 1, 0);
  amrex::MultiFab::Copy(vol, mf, comp, 0, 1, 0);

  if (level < parent->finestLevel() && finemask) {
    const amrex::MultiFab& mask = getLevel(level + 1).build_fine_mask();
    amrex::MultiFab::Multiply(vol, mask, 0, 0, 1, 0);
  }

  if (eb_in_domain) {
    zero_covered_cells(vol, vfrac);
    amrex::MultiFab::Multiply(vol, vfrac, 0, 0, 1, 0);
  }

  return amrex::MultiFab::Dot(vol, 0, volume, 0, 1, 0, local);
}

amrex::Real
PeleC::maxDerive(const std::string& name, amrex::Real time, bool local)
{
  // Note: Includes all cells in MF, even those covered by finer grids or EB
  auto mf = derive(name, time, 0);

  BL_ASSERT(!(mf == nullptr));

  return mf->max(0, 0, local);
}

amrex::Real
PeleC::minDerive(const std::string& name, amrex::Real time, bool local)
{
  // Note: Includes all cells in MF, even those covered by finer grids or EB
  auto mf = derive(name, time, 0);

  BL_ASSERT(!(mf == nullptr));

  return mf->min(0, 0, local);
}

int
PeleC::find_datalog_index(const std::string& logname)
{
  for (int ii = 0; ii < parent->NumDataLogs(); ii++) {
    if (logname == parent->DataLogName(ii)) {
      return ii;
    }
  }
  return -1; // Requested log not found
}
