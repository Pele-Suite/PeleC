#include <memory>

#include "AMReX_EB_Redistribution.H"
#include "EB.H"
#include "prob.H"
#include "PeleCUtilities.H"
#include "Geometry.H"

inline bool
PeleC::ebInitialized()
{
  return eb_initialized;
}

void
PeleC::init_eb()
{
  if (!eb_in_domain) {
    return;
  }

  // Build the geometry information; this is done for each new set of grids
  initialize_eb2_structs();

  // Set up CC signed distance container to control EB refinement
  initialize_signed_distance();
}

// Set up PeleC EB Datastructures from AMReX EB2 constructs
// At the end of this routine, the following structures are populated:
//  - MultiFAB vfrac
//  - sv_eb_bndry_geom

void
PeleC::initialize_eb2_structs()
{
  BL_PROFILE("PeleC::initialize_eb2_structs()");
  amrex::Print() << "Initializing EB2 structs" << std::endl;

  static_assert(
    std::is_standard_layout_v<EBBndryGeom>,
    "EBBndryGeom is not standard layout");

  const auto& ebfactory =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());

  // These are the data sources
  amrex::MultiFab::Copy(vfrac, ebfactory.getVolFrac(), 0, 0, 1, numGrow());
  const auto* bndrycent = &(ebfactory.getBndryCent());
  areafrac = ebfactory.getAreaFrac();
  facecent = ebfactory.getFaceCent();

  // First pass over fabs to fill sparse per cut-cell ebg structures
  sv_eb_bndry_geom.resize(vfrac.local_size());
  sv_eb_bndry_grad_stencil.resize(vfrac.local_size());
  sv_eb_flux.resize(vfrac.local_size());
  sv_eb_bcval.resize(vfrac.local_size());

  auto const& flags = ebfactory.getMultiEBCellFlagFab();

  // Boundary stencil option: 0 = original, 1 = amrex way, 2 = least squares
  amrex::ParmParse pp("ebd");

  int bgs = 0;
  pp.query("boundary_grad_stencil_type", bgs);

  if (bgs == 0) {
    amrex::Print() << "Using quadratic stencil for the EB gradient\n";
  } else if (bgs == 1) {
    amrex::Print() << "Using least-squares stencil for the EB gradient\n";
  } else {
    amrex::Print() << "Unknown or unspecified EB gradient stencil type:" << bgs
                   << std::endl;
    amrex::Abort();
  }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
    const amrex::Box tbox = mfi.growntilebox();
    const amrex::EBCellFlagFab& flagfab = flags[mfi];
    const amrex::FabType typ = flagfab.getType(tbox);
    const int iLocal = mfi.LocalIndex();

    if ((typ == amrex::FabType::regular) || (typ == amrex::FabType::covered)) {
      // do nothing
    } else if (typ == amrex::FabType::singlevalued) {
      auto const& flag_arr = flags.const_array(mfi);

      const auto nallcells = static_cast<int>(tbox.numPts());
      amrex::Gpu::DeviceVector<int> cutcell_offset(nallcells, 0);
      auto* d_cutcell_offset = cutcell_offset.data();
      const auto ncutcells = amrex::Scan::PrefixSum<int>(
        nallcells,
        [=] AMREX_GPU_DEVICE(int icell) -> int {
          const auto iv = tbox.atOffset(icell);
          return static_cast<int>(flag_arr(iv).isSingleValued());
        },
        [=] AMREX_GPU_DEVICE(int icell, int const& x) {
          d_cutcell_offset[icell] = x;
        },
        amrex::Scan::Type::exclusive, amrex::Scan::retSum);

      AMREX_ASSERT(ncutcells == flagfab.getNumCutCells(tbox));

      sv_eb_bndry_geom[iLocal].resize(ncutcells);
      if (ncutcells > 0) {
        auto* d_sv_eb_bndry_geom = sv_eb_bndry_geom[iLocal].data();
        amrex::ParallelFor(
          tbox,
          [=] AMREX_GPU_DEVICE(int i, int j, int AMREX_D_PICK(, , k)) noexcept {
            const amrex::IntVect iv(amrex::IntVect(AMREX_D_DECL(i, j, k)));
            if (flag_arr(iv).isSingleValued()) {
              const auto icell = tbox.index(iv);
              const auto idx = d_cutcell_offset[icell];
              d_sv_eb_bndry_geom[idx].iv = iv;
            }
          });
      }

      // Now fill the sv_eb_bndry_geom
      auto const& vfrac_arr = vfrac.const_array(mfi);
      auto const& bndrycent_arr = bndrycent->const_array(mfi);
      AMREX_D_TERM(
        auto const& apx = areafrac[0]->const_array(mfi);
        , auto const& apy = areafrac[1]->const_array(mfi);
        , auto const& apz = areafrac[2]->const_array(mfi);)
      pc_fill_sv_ebg(
        tbox, ncutcells, vfrac_arr, bndrycent_arr, AMREX_D_DECL(apx, apy, apz),
        sv_eb_bndry_geom[iLocal].data());

      // Fill in boundary gradient for cut cells in this grown tile
      sv_eb_bndry_grad_stencil[iLocal].resize(ncutcells);
      const auto& dx = geom.CellSizeArray();
      if (bgs == 0) {
        pc_fill_bndry_grad_stencil_quadratic(
          tbox, dx, ncutcells, sv_eb_bndry_geom[iLocal].data(), ncutcells,
          sv_eb_bndry_grad_stencil[iLocal].data());
      } else if (bgs == 1) {
        pc_fill_bndry_grad_stencil_ls(
          tbox, dx, ncutcells, sv_eb_bndry_geom[iLocal].data(), ncutcells,
          flags.array(mfi), sv_eb_bndry_grad_stencil[iLocal].data());
      } else {
        amrex::Print()
          << "Unknown or unspecified boundary gradient stencil type:" << bgs
          << std::endl;
        amrex::Abort();
      }

      if (eb_noslip or eb_isothermal) {
        amrex::Box sbox = amrex::grow(tbox, -3);
        pc_check_bndry_grad_stencil(
          sbox, ncutcells, flags.array(mfi),
          sv_eb_bndry_grad_stencil[iLocal].data());
      }

      sv_eb_flux[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], NVAR);
      sv_eb_bcval[iLocal].define(sv_eb_bndry_grad_stencil[iLocal], QVAR);

      if (eb_isothermal && (diffuse_temp || diffuse_enth)) {
        sv_eb_bcval[iLocal].setVal(eb_boundary_T, QTEMP);
      }
      if (eb_noslip && diffuse_vel) {

        if (do_rf) {
#if AMREX_SPACEDIM > 2
          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxlev =
            geom.CellSizeArray();
          const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo =
            geom.ProbLoArray();
          int axis = rf_axis;
          amrex::Real omega = rf_omega;
          amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
            AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
          amrex::Real rfdist2 = rf_rad * rf_rad;

          auto* u_bcval = sv_eb_bcval[iLocal].dataPtr(QU);
          auto* v_bcval = sv_eb_bcval[iLocal].dataPtr(QV);
          auto* w_bcval = sv_eb_bcval[iLocal].dataPtr(QW);
          auto* sten = sv_eb_bndry_geom[iLocal].data();
          if (ncutcells > 0) {
            amrex::ParallelFor(ncutcells, [=] AMREX_GPU_DEVICE(int L) {
              const auto& iv = sten[L].iv;

              amrex::Real xloc =
                prob_lo[0] + (iv[0] + 0.5 + sten[L].eb_centroid[0]) * dxlev[0];
              amrex::Real yloc =
                prob_lo[1] + (iv[1] + 0.5 + sten[L].eb_centroid[1]) * dxlev[1];
              amrex::Real zloc =
                prob_lo[2] + (iv[2] + 0.5 + sten[L].eb_centroid[2]) * dxlev[2];

              amrex::RealVect r(AMREX_D_DECL(0.0, 0.0, 0.0));
              amrex::RealVect w(AMREX_D_DECL(0.0, 0.0, 0.0));

              r[0] = xloc - axis_loc[0];
              r[1] = yloc - axis_loc[1];
              r[2] = zloc - axis_loc[2];

              amrex::Real rmag2 = r.radSquared();
              rmag2 -= r[axis] * r[axis]; // only in plane

              if (rmag2 > rfdist2) {
                w[axis] = omega;
                amrex::RealVect w_cross_r = w.crossProduct(r);

                u_bcval[L] = -w_cross_r[0];
                v_bcval[L] = -w_cross_r[1];
                w_bcval[L] = -w_cross_r[2];
              } else {
                u_bcval[L] = 0.0;
                v_bcval[L] = 0.0;
                w_bcval[L] = 0.0;
              }
            });
          }
#endif
        } else {
          sv_eb_bcval[iLocal].setVal(0, QU, AMREX_SPACEDIM);
        }
      }

    } else {
      amrex::Print() << "unknown (or multivalued) fab type" << std::endl;
      amrex::Abort();
    }
  }

  // Second pass over dirs and fabs to fill flux interpolation stencils
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    flux_interp_stencil[dir].resize(vfrac.local_size());

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(vfrac, false); mfi.isValid(); ++mfi) {
      const amrex::Box tbox = mfi.growntilebox(numGrow());
      const auto& flagfab = flags[mfi];
      const amrex::FabType typ = flagfab.getType(tbox);
      const int iLocal = mfi.LocalIndex();

      if (typ == amrex::FabType::singlevalued) {
        auto const& flag_arr = flagfab.const_array();
        const auto afrac_arr = (*areafrac[dir])[mfi].const_array();
        const auto facecent_arr = (*facecent[dir])[mfi].const_array();
        const auto fbox = amrex::surroundingNodes(tbox, dir);
        const auto nallfaces = static_cast<int>(fbox.numPts());
        amrex::Gpu::DeviceVector<int> cutfaces_offset(nallfaces, 0);
        auto* d_cutface_offset = cutfaces_offset.data();
        const auto ncutfaces = amrex::Scan::PrefixSum<int>(
          nallfaces,
          [=] AMREX_GPU_DEVICE(int iface) -> int {
            const auto iv = fbox.atOffset(iface);
            const bool covered_cells =
              flag_arr(iv).isCovered() &&
              flag_arr(iv - amrex::IntVect::TheDimensionVector(dir))
                .isCovered();
            return static_cast<int>((afrac_arr(iv) < 1.0) && (!covered_cells));
          },
          [=] AMREX_GPU_DEVICE(int iface, int const& x) {
            d_cutface_offset[iface] = x;
          },
          amrex::Scan::Type::exclusive, amrex::Scan::retSum);

        if (ncutfaces > 0) {
          flux_interp_stencil[dir][iLocal].resize(ncutfaces);
          auto* d_flux_interp_stencil = flux_interp_stencil[dir][iLocal].data();
          amrex::ParallelFor(
            fbox, [=] AMREX_GPU_DEVICE(
                    int i, int j, int AMREX_D_PICK(, , k)) noexcept {
              const amrex::IntVect iv(amrex::IntVect(AMREX_D_DECL(i, j, k)));
              const bool covered_cells =
                flag_arr(iv).isCovered() &&
                flag_arr(iv - amrex::IntVect::TheDimensionVector(dir))
                  .isCovered();
              if ((afrac_arr(iv) < 1.0) && (!covered_cells)) {
                const auto iface = fbox.index(iv);
                const auto idx = d_cutface_offset[iface];
                d_flux_interp_stencil[idx].iv = iv;
              }
            });

          pc_fill_flux_interp_stencil(
            tbox, ncutfaces, facecent_arr, afrac_arr, d_flux_interp_stencil);
        }
      } else if (
        (typ != amrex::FabType::regular) && (typ != amrex::FabType::covered)) {
        amrex::Abort("multi-valued flux interp stencil to be implemented");
      }
    }
  }
}

void
PeleC::define_body_state()
{
  BL_PROFILE("PeleC::define_body_state()");

  if (!eb_in_domain) {
    return;
  }

  if (eb_zero_body_state) {
    for (int n = 0; n < NVAR; ++n) {
      body_state[n] = 0.0;
    }
    body_state_set = true;
    return;
  }

  // Scan over data and find a point in the fluid to use to
  // set computable values in cells outside the domain
  // We look for the lexicographically first valid point
  if (!body_state_set) {
    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();
    auto flag_arrays = flags.const_arrays();

    amrex::MultiFab minIdxTmp(grids, dmap, 1, 0, amrex::MFInfo(), Factory());
    auto tmp_arrays = minIdxTmp.arrays();

    const auto& dbox = geom.Domain();
    const long max_idx = dbox.numPts();

    // function minimized at the first uncovered cell
    amrex::ParallelFor(
      minIdxTmp, [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        amrex::ignore_unused(k);
        const amrex::IntVect iv(amrex::IntVect(AMREX_D_DECL(i, j, k)));
        const long idx_masked =
          flag_arrays[nbx](iv).isRegular() ? dbox.index(iv) : max_idx;
        tmp_arrays[nbx](iv) = static_cast<amrex::Real>(idx_masked);
      });
    amrex::Gpu::synchronize();

    // select the data at the first uncovered cell
    const amrex::IntVect idx_min = minIdxTmp.minIndex(0);
    const amrex::Box tgt_box{idx_min, idx_min};
    const amrex::MultiFab& S = get_new_data(State_Type);
    for (int n = 0; n < NVAR; ++n) {
      body_state[n] = S.min(tgt_box, n);
    }
    body_state_set = true;
  }
}

void
PeleC::set_body_state(amrex::MultiFab& S)
{
  BL_PROFILE("PeleC::set_body_state()");

  if (!eb_in_domain) {
    return;
  }

  if (!body_state_set) {
    define_body_state();
  }

  auto const& fact = dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();

  auto const& sarrs = S.arrays();
  auto const& flagarrs = flags.const_arrays();
  auto const captured_body_state = body_state;
  const amrex::IntVect ngs(0);
  amrex::ParallelFor(
    S, ngs, NVAR,
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k, int n) noexcept {
      pc_set_body_state(
        i, j, k, n, flagarrs[nbx], captured_body_state, sarrs[nbx]);
    });
  amrex::Gpu::synchronize();
}

// Sets up implicit function using EB2 infrastructure
void
initialize_EB2(
  const amrex::Geometry& geom,
  const int eb_max_level,
  const int max_level,
  const int coarsening,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  const amrex::IntVect& max_grid_size)
{
  BL_PROFILE("PeleC::initialize_EB2()");

  amrex::Print() << "Initializing EB2" << std::endl;
  amrex::ParmParse ppeb2("eb2");

  std::string geom_type("all_regular");
  ppeb2.query("geom_type", geom_type);

  int max_coarsening_level = 0;
  for (int lev = 0; lev < max_level; ++lev) {
    // Since EB always coarsens by a factor of 2
    if (ref_ratio[lev] == 2) {
      max_coarsening_level += 1;
    } else if (ref_ratio[lev] == 4) {
      max_coarsening_level += 2;
    } else {
      amrex::Abort("initalize_EB2: Refinement ratio must be 2 or 4");
    }
  }

  // Custom types defined here - all_regular, plane, sphere, etc, will get
  // picked up by default (see AMReX_EB2.cpp around L100 )
  amrex::Vector<std::string> amrex_defaults(
    {"all_regular", "box", "cylinder", "plane", "sphere", "torus", "parser",
     "stl"});
  if (!(std::find(amrex_defaults.begin(), amrex_defaults.end(), geom_type) !=
        amrex_defaults.end())) {
    std::unique_ptr<pele::pelec::Geometry> geometry(
      pele::pelec::Geometry::create(geom_type));
    geometry->build(geom, max_coarsening_level + coarsening);
  } else {
    amrex::EB2::Build(
      geom, max_coarsening_level + coarsening,
      max_coarsening_level + coarsening);
  }

  // Add finer level, might be inconsistent with the coarser level created
  // above.
  // EY: This condition is not acceptable in AMReX with stl format
  if ((geom_type != "chkfile") && (geom_type != "stl")) {
    amrex::EB2::addFineLevels(max_level - eb_max_level);
  } else {
    // The AMReX implementation for these does not support addFineLevels
    AMREX_ALWAYS_ASSERT(max_level == eb_max_level);
  }

  bool write_chk_geom = false;
  ppeb2.query("write_chk_geom", write_chk_geom);
  if (write_chk_geom) {
    const auto& is = amrex::EB2::IndexSpace::top();
    const auto& eb_level = is.getLevel(geom);
    std::string chkfile = "chk_geom";
    ppeb2.query("chkfile", chkfile);

    eb_level.write_to_chkpt_file(
      chkfile, amrex::EB2::ExtendDomainFace(), max_grid_size[0]);
  }
}

void
PeleC::initialize_signed_distance()
{
  BL_PROFILE("PeleC::initialize_signed_distance()");
  if (level == 0) {
    const auto& ebfactory =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(Factory());
    signed_dist_0.define(grids, dmap, 1, 1, amrex::MFInfo(), ebfactory);

    // Estimate the maximum distance we need in terms of level 0 dx:
    auto extentFactor = static_cast<amrex::Real>(parent->nErrorBuf(0));
    for (int ilev = 1; ilev <= parent->maxLevel(); ++ilev) {
      extentFactor += static_cast<amrex::Real>(parent->nErrorBuf(ilev)) /
                      std::pow(
                        static_cast<amrex::Real>(parent->refRatio(ilev - 1)[0]),
                        static_cast<amrex::Real>(ilev));
    }
    extentFactor *= tagging_parm->detag_eb_factor;

    amrex::MultiFab signDist(
      convert(grids, amrex::IntVect::TheUnitVector()), dmap, 1, 1,
      amrex::MFInfo(), ebfactory);
    amrex::FillSignedDistance(signDist, true);

    const auto& sd_ccs = signed_dist_0.arrays();
    const auto& sd_nds = signDist.const_arrays();
    const amrex::IntVect ngs(signed_dist_0.nGrow());
    amrex::ParallelFor(
      signed_dist_0, ngs,
      [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
        const auto& sd_cc = sd_ccs[nbx];
        const auto& sd_nd = sd_nds[nbx];
        const amrex::Real fac = AMREX_D_PICK(0.5, 0.25, 0.125);
        sd_cc(i, j, k) = AMREX_D_TERM(
          sd_nd(i, j, k) + sd_nd(i + 1, j, k),
          +sd_nd(i, j + 1, k) + sd_nd(i + 1, j + 1, k),
          +sd_nd(i, j, k + 1) + sd_nd(i + 1, j, k + 1) +
            sd_nd(i, j + 1, k + 1) + sd_nd(i + 1, j + 1, k + 1));
        sd_cc(i, j, k) *= fac;
      });
    amrex::Gpu::synchronize();

    signed_dist_0.FillBoundary(parent->Geom(0).periodicity());
    extend_signed_distance(&signed_dist_0, extentFactor);
  }
}

void
PeleC::eb_distance(const int lev, amrex::MultiFab& signDistLev)
{
  BL_PROFILE("PeleC::eb_distance()");
  if (lev == 0) {
    amrex::MultiFab::Copy(signDistLev, signed_dist_0, 0, 0, 1, 0);
    return;
  }

  // A pair of MF to hold crse & fine dist data
  amrex::Array<std::unique_ptr<amrex::MultiFab>, 2> MFpair;

  // dummy bcs
  amrex::Vector<amrex::BCRec> bcrec_dummy(1);
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    bcrec_dummy[0].setLo(dir, amrex::BCType::int_dir);
    bcrec_dummy[0].setHi(dir, amrex::BCType::int_dir);
  }

  // Interpolate on successive levels up to lev
  const auto& pc_0lev = getLevel(0);
  for (int ilev = 1; ilev <= lev; ++ilev) {

    // Use MF EB interp
    auto& interpolater = amrex::eb_mf_lincc_interp;

    // Get signDist on coarsen fineBA
    const auto& pc_ilev = getLevel(ilev);
    const auto& grids_ilev = pc_ilev.grids;
    const auto& dmap_ilev = pc_ilev.DistributionMap();
    amrex::BoxArray coarsenBA(grids_ilev.size());
    for (int j = 0, N = static_cast<int>(coarsenBA.size()); j < N; ++j) {
      coarsenBA.set(
        j, interpolater.CoarseBox(grids_ilev[j], parent->refRatio(ilev - 1)));
    }
    amrex::MultiFab coarsenSignDist(coarsenBA, dmap_ilev, 1, 0);
    coarsenSignDist.setVal(0.0);
    const amrex::MultiFab* crseSignDist =
      (ilev == 1) ? &pc_0lev.signed_dist_0 : MFpair[0].get();
    coarsenSignDist.ParallelCopy(*crseSignDist, 0, 0, 1);

    // Interpolate on current ilev
    amrex::MultiFab* currentSignDist;
    if (ilev < lev) {
      const auto& ebfactory =
        dynamic_cast<amrex::EBFArrayBoxFactory const&>(pc_ilev.Factory());
      MFpair[1] = std::make_unique<amrex::MultiFab>(
        grids_ilev, dmap_ilev, 1, 0, amrex::MFInfo(), ebfactory);
    }
    currentSignDist = (ilev == lev) ? &signDistLev : MFpair[1].get();

    interpolater.interp(
      coarsenSignDist, 0, *currentSignDist, 0, 1, amrex::IntVect(0),
      parent->Geom(ilev - 1), parent->Geom(ilev), parent->Geom(ilev).Domain(),
      parent->refRatio(ilev - 1), {bcrec_dummy}, 0);

    // Swap MFpair
    if (ilev < lev) {
      std::swap(MFpair[0], MFpair[1]);
    }
  }
}

// Extend the cell-centered based signed distance function
void
PeleC::extend_signed_distance(
  amrex::MultiFab* signDist, amrex::Real extendFactor)
{
  // This is a not-so-pretty piece of code that'll take AMReX cell-averaged
  // signed distance and propagates it manually up to the point where we need to
  // have it for derefining.
  BL_PROFILE("PeleC::extend_signed_distance()");
  amrex::Real maxSignedDist = signDist->max(0);
  const auto& ebfactory =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(signDist->Factory());
  const auto& flags = ebfactory.getMultiEBCellFlagFab();
  int nGrowFac = flags.nGrow() + 1;
  const auto& dx = parent->Geom(0).CellSizeArray();
  const amrex::Real dx_max = *std::max_element(dx.begin(), dx.end());

  // First set the region far away at the max value we need
  auto const& sd_ccs = signDist->arrays();
  const amrex::IntVect ngs(signDist->nGrow());
  amrex::ParallelFor(
    *signDist, ngs,
    [=] AMREX_GPU_DEVICE(int nbx, int i, int j, int k) noexcept {
      const auto& sd_cc = sd_ccs[nbx];
      if (sd_cc(i, j, k) >= maxSignedDist - 1e-12) {
        sd_cc(i, j, k) = nGrowFac * dx_max * extendFactor;
      }
    });
  amrex::Gpu::synchronize();

  // Iteratively compute the distance function in boxes, propagating across
  // boxes using ghost cells If needed, increase the number of loop to extend
  // the reach of the distance function
  const int nMaxLoop = 4;
  for (int dloop = 1; dloop <= nMaxLoop; dloop++) {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(*signDist, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& bx = mfi.tilebox();
      const amrex::Box& gbx = grow(bx, 1);
      if (flags[mfi].getType(gbx) == amrex::FabType::covered) {
        continue;
      }
      auto const& sd_cc = signDist->array(mfi);
      ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        const auto glo = amrex::lbound(gbx);
        const auto ghi = amrex::ubound(gbx);
        const amrex::Real extendedDist = dx_max * extendFactor;
        if (sd_cc(i, j, k) >= maxSignedDist - 1e-12) {
          amrex::Real closestEBDist = 1e12;
          for (int kk = glo.z; kk <= ghi.z; ++kk) {
            for (int jj = glo.y; jj <= ghi.y; ++jj) {
              for (int ii = glo.x; ii <= ghi.x; ++ii) {
                if ((i != ii) || (j != jj) || (k != kk)) {
                  if (sd_cc(ii, jj, kk) > 0.0) {
                    const amrex::Real distToCell = std::sqrt(AMREX_D_TERM(
                      ((i - ii) * dx[0] * (i - ii) * dx[0]),
                      +((j - jj) * dx[1] * (j - jj) * dx[1]),
                      +((k - kk) * dx[2] * (k - kk) * dx[2])));
                    const amrex::Real distToEB = distToCell + sd_cc(ii, jj, kk);
                    closestEBDist = amrex::min(distToEB, closestEBDist);
                  }
                }
              }
            }
          }
          if (closestEBDist < 1e10) {
            sd_cc(i, j, k) = closestEBDist;
          } else {
            sd_cc(i, j, k) = extendedDist;
          }
        }
      });
    }
    signDist->FillBoundary(parent->Geom(0).periodicity());
  }
}

void
PeleC::InitialRedistribution(
  const amrex::Real time,
  const amrex::Vector<amrex::BCRec> bcs,
  amrex::MultiFab& S_new)
{
  BL_PROFILE("PeleC::InitialRedistribution()");

  // Don't redistribute if there is no EB or if the redistribution type is
  // anything other than StateRedist
  if ((!eb_in_domain) || (redistribution_type != "StateRedist")) {
    return;
  }

  if (verbose != 0) {
    amrex::Print() << "Doing initial redistribution... " << std::endl;
  }

  // Initial data are set at new time step
  amrex::MultiFab tmp(
    grids, dmap, S_new.nComp(), numGrow(), amrex::MFInfo(), Factory());

  amrex::MultiFab::Copy(tmp, S_new, 0, 0, S_new.nComp(), S_new.nGrow());
  FillPatch(*this, tmp, numGrow(), time, State_Type, 0, S_new.nComp());
  EB_set_covered(tmp, 0.0);

  amrex::Gpu::DeviceVector<amrex::BCRec> d_bcs(bcs.size());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, bcs.begin(), bcs.end(), d_bcs.begin());

  for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
       ++mfi) {
    const amrex::Box& bx = mfi.validbox();

    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(S_new.Factory());

    auto const& flags = fact.getMultiEBCellFlagFab()[mfi];
    amrex::Array4<const amrex::EBCellFlag> const& flag_arr =
      flags.const_array();

    if (
      (flags.getType(amrex::grow(bx, 1)) != amrex::FabType::covered) &&
      (flags.getType(amrex::grow(bx, 1)) != amrex::FabType::regular)) {
      amrex::Array4<const amrex::Real> AMREX_D_DECL(fcx, fcy, fcz), ccc,
        AMREX_D_DECL(apx, apy, apz);

      AMREX_D_TERM(
        fcx = facecent[0]->const_array(mfi);
        , fcy = facecent[1]->const_array(mfi);
        , fcz = facecent[2]->const_array(mfi););

      ccc = fact.getCentroid().const_array(mfi);

      AMREX_D_TERM(
        apx = areafrac[0]->const_array(mfi);
        , apy = areafrac[1]->const_array(mfi);
        , apz = areafrac[2]->const_array(mfi););

      const auto& sarr = S_new.array(mfi);
      const auto& tarr = tmp.array(mfi);
      ApplyInitialRedistribution(
        bx, NVAR, sarr, tarr, flag_arr, AMREX_D_DECL(apx, apy, apz),
        vfrac.const_array(mfi), AMREX_D_DECL(fcx, fcy, fcz), ccc,
        d_bcs.dataPtr(), geom, redistribution_type, eb_srd_max_order);

      // Make sure rho is same as sum rhoY after redistribution
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          amrex::Real drhoYsum = 0.0;
          for (int n = 0; n < NUM_SPECIES; n++) {
            drhoYsum += sarr(i, j, k, UFS + n) - tarr(i, j, k, UFS + n);
          }
          sarr(i, j, k, URHO) = tarr(i, j, k, URHO) + drhoYsum;
        });
      amrex::ParallelFor(
        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          pc_check_initial_species(i, j, k, sarr);
        });
    }
  }
  set_body_state(S_new);
}
