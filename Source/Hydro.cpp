#include "Hydro.H"

// Set up the source terms to go into the hydro.
void
PeleC::construct_hydro_source(
  const amrex::MultiFab& S,
  amrex::Real time,
  amrex::Real dt,
  int /*amr_iteration*/,
  int /*amr_ncycle*/,
  int sub_iteration,
  int sub_ncycle)
{
  if (do_mol) {
    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Zeroing Godunov-based hydro advance" << std::endl;
    }
    hydro_source.setVal(0);
  } else {
    BL_PROFILE("PeleC::construct_hydro_source()");

    if ((verbose != 0) && amrex::ParallelDescriptor::IOProcessor()) {
      amrex::Print() << "... Computing Godunov hydro advance" << std::endl;
    }

    AMREX_ASSERT(S.nGrow() >= numGrow() + nGrowF);

    const int ng = 0;
    sources_for_hydro.setVal(0.0);
    for (int src : src_list) {
      amrex::MultiFab::Saxpy(
        sources_for_hydro, 0.5, *new_sources[src], 0, 0, NVAR, ng);
      amrex::MultiFab::Saxpy(
        sources_for_hydro, 0.5, *old_sources[src], 0, 0, NVAR, ng);
    }

    // Add I_R terms to advective forcing
    if (do_react) {
      amrex::MultiFab::Add(
        sources_for_hydro, get_new_data(Reactions_Type), 0, FirstSpec,
        NUM_SPECIES, ng);
      amrex::MultiFab::Add(
        sources_for_hydro, get_new_data(Reactions_Type), NUM_SPECIES, Eden, 1,
        ng);
    }
    sources_for_hydro.FillBoundary(geom.periodicity());
    hydro_source.setVal(0);

    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx =
      geom.CellSizeArray();

    amrex::Real courno = std::numeric_limits<amrex::Real>::lowest();

    const amrex::MultiFab& S_new = get_new_data(State_Type);

    auto const& fact =
      dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
    auto const& flags = fact.getMultiEBCellFlagFab();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion()) \
  reduction(max : courno)
#endif
    {
      amrex::Real cflLoc = std::numeric_limits<amrex::Real>::lowest();

      const int* domain_lo = geom.Domain().loVect();
      const int* domain_hi = geom.Domain().hiVect();

      for (amrex::MFIter mfi(S_new, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        const amrex::Box& qbx = amrex::grow(bx, numGrow() + nGrowF);
        const amrex::Box& fbx = amrex::grow(bx, nGrowF);

        const auto& flag_fab = flags[mfi];
        const amrex::Array4<amrex::EBCellFlag const>& flag_arr =
          flag_fab.const_array();
        if (flag_fab.getType(bx) == amrex::FabType::covered) {
          continue;
        }

        amrex::GpuArray<amrex::FArrayBox, AMREX_SPACEDIM> flux;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          const amrex::Box& efbx = surroundingNodes(fbx, dir);
          flux[dir].resize(efbx, NVAR, amrex::The_Async_Arena());
          flux[dir].setVal<amrex::RunOn::Device>(0.0);
        }

        auto const& sarr = S.array(mfi);
        auto const& hyd_src = hydro_source.array(mfi);

        // Temporary Fabs
        amrex::FArrayBox q(qbx, QVAR, amrex::The_Async_Arena());
        amrex::FArrayBox qaux(qbx, NQAUX, amrex::The_Async_Arena());
        amrex::FArrayBox src_q(qbx, QVAR, amrex::The_Async_Arena());

        // Get Arrays to pass to the gpu.
        auto const& qarr = q.array();
        auto const& qauxar = qaux.array();
        auto const& srcqarr = src_q.array();

        {
          BL_PROFILE("PeleC::ctoprim()");
          amrex::ParallelFor(
            qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              if (!flag_arr(i, j, k).isCovered()) {
                pc_ctoprim(i, j, k, sarr, qarr, qauxar);
              } else {
                for (int n = 0; n < QVAR; n++) {
                  qarr(i, j, k, n) = 0.0;
                }
              }
            });
        }

        // TODO GPUize NCSCBC
        // Imposing Ghost-Cells Navier-Stokes Characteristic BCs if "UserBC" are
        // used For the theory, see Motheau et al. AIAA J. Vol. 55, No. 10 : pp.
        // 3399-3408, 2017.
        //
        // The user should provide a bcnormal routine in bc_fill_module with
        // additional optional arguments to temporary fill ghost-cells for
        // amrex::BCType::ext_dir and to provide target BC values. See the
        // examples.

        // Allocate fabs for bcMask. Note that we grow in the opposite direction
        // because the Riemann solver wants a face value in a ghost-cell
        /*
              for (int dir = 0; dir < AMREX_SPACEDIM ; dir++)  {
                const Box& bxtmp = amrex::surroundingNodes(fbx,dir);
                Box TestBox(bxtmp);
                for(int d=0; d<AMREX_SPACEDIM; ++d) {
                  if (dir!=d) TestBox.grow(d,1);
                }
                bcMask[dir].resize(TestBox,1, amrex::The_Async_Arena());
                bcMask[dir].setVal(0);
              }

              // Because bcMask is read in the Riemann solver in any case,
              // here we put physbc values in the appropriate faces for the
           non-nscbc case set_bc_mask(lo, hi, domain_lo, domain_hi,
                          AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                                 AMREX_TO_FORTRAN(bcMask[1]),
                                 AMREX_TO_FORTRAN(bcMask[2])));

              if (nscbc_adv == 1)
              {
                impose_NSCBC(lo, hi, domain_lo, domain_hi,
                             AMREX_TO_FORTRAN(*statein),
                             AMREX_TO_FORTRAN(q.fab()),
                             AMREX_TO_FORTRAN(qaux.fab()),
                             AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
                                    AMREX_TO_FORTRAN(bcMask[1]),
                                    AMREX_TO_FORTRAN(bcMask[2])),
                             &flag_nscbc_isAnyPerio, flag_nscbc_perio,
                             &time, dx, &dt);
              }
        */
        {
          BL_PROFILE("PeleC::srctoprim()");
          const auto& src_in = sources_for_hydro.array(mfi);
          amrex::ParallelFor(
            qbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
              if (!flag_arr(i, j, k).isCovered()) {
                pc_srctoprim(i, j, k, qarr, qauxar, src_in, srcqarr);
              } else {
                for (int n = 0; n < QVAR; n++) {
                  srcqarr(i, j, k, n) = 0.0;
                }
              }
            });
        }
        const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          flx_arr{
            {AMREX_D_DECL(flux[0].array(), flux[1].array(), flux[2].array())}};
        const amrex::GpuArray<
          const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
          a{{AMREX_D_DECL(
            area[0].const_array(mfi), area[1].const_array(mfi),
            area[2].const_array(mfi))}};

        amrex::FArrayBox dm_as_fine(
          amrex::Box::TheUnitBox(), hydro_source.nComp(),
          amrex::The_Async_Arena());
        amrex::FArrayBox fab_drho_as_crse(
          amrex::Box::TheUnitBox(), hydro_source.nComp(),
          amrex::The_Async_Arena());
        amrex::IArrayBox fab_rrflag_as_crse(
          amrex::Box::TheUnitBox(), 1, amrex::The_Async_Arena());

        const int ngrow_bx = (redistribution_type == "StateRedist") ? 3 : 2;
        const amrex::Box& fbxg_i = grow(fbx, ngrow_bx);
        if (flag_fab.getType(fbxg_i) == amrex::FabType::singlevalued) {

          BL_PROFILE("PeleC::umdrv_eb()")

          auto const& vfrac_arr = vfrac.const_array(mfi);

          amrex::EBFluxRegister* fr_as_crse =
            (do_reflux && (level < parent->finestLevel()))
              ? &getFluxReg(level + 1)
              : nullptr;
          amrex::EBFluxRegister* fr_as_fine =
            (do_reflux && (level > 0)) ? &getFluxReg(level) : nullptr;

          const int as_crse = static_cast<int>(fr_as_crse != nullptr);
          const int as_fine = static_cast<int>(fr_as_fine != nullptr);

          auto* p_drho_as_crse = (fr_as_crse != nullptr)
                                   ? fr_as_crse->getCrseData(mfi)
                                   : &fab_drho_as_crse;
          const auto* p_rrflag_as_crse = (fr_as_crse != nullptr)
                                           ? fr_as_crse->getCrseFlag(mfi)
                                           : &fab_rrflag_as_crse;

          if (fr_as_fine != nullptr) {
            const amrex::Box dbox1 = geom.growPeriodicDomain(1);
            const amrex::Box bx_for_dm(amrex::grow(fbx, 1) & dbox1);
            dm_as_fine.resize(bx_for_dm, hydro_source.nComp());
            dm_as_fine.setVal<amrex::RunOn::Device>(0.0);
          }

          const amrex::StateDescriptor* desc = state[State_Type].descriptor();
          const auto& bcs = desc->getBCs();
          amrex::Gpu::DeviceVector<amrex::BCRec> bcs_d(desc->nComp());
          amrex::Gpu::copy(
            amrex::Gpu::hostToDevice, bcs.begin(), bcs.end(), bcs_d.begin());

          const auto& dxInv = geom.InvCellSizeArray();

          pc_umdrv_eb(
            fbx, fbxg_i, mfi, geom, &fact, phys_bc.lo(), phys_bc.hi(), sarr,
            hyd_src, qarr, qauxar, srcqarr, vfrac_arr, flag_arr, dx, dxInv,
            flx_arr, as_crse, p_drho_as_crse->array(),
            p_rrflag_as_crse->array(), as_fine, dm_as_fine.array(),
            level_mask.const_array(mfi), dt, ppm_type, plm_iorder,
            use_flattening, difmag, bcs_d.data(), redistribution_type,
            eb_weights_type, eb_srd_max_order, eb_clean_massfrac,
            eb_clean_massfrac_threshold, cflLoc);

        } else if (flag_fab.getType(fbxg_i) == amrex::FabType::regular) {
          BL_PROFILE("PeleC::umdrv()");
          pc_umdrv(
            time, fbx, domain_lo, domain_hi, phys_bc.lo(), phys_bc.hi(), sarr,
            hyd_src, qarr, qauxar, srcqarr, dx, dt, ppm_type, plm_iorder,
            use_flattening, use_hybrid_weno, weno_scheme, difmag, flx_arr, a,
            volume.array(mfi), cflLoc);
        } else if (flag_fab.getType(fbxg_i) == amrex::FabType::multivalued) {
          amrex::Abort("multi-valued cells are not supported");
        }

        courno = amrex::max<amrex::Real>(courno, cflLoc);

        // Filter hydro source and fluxes here
        if (use_explicit_filter) {
          BL_PROFILE("PeleC::apply_filter()");
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            const amrex::Box& bxtmp = amrex::surroundingNodes(bx, dir);
            amrex::FArrayBox filtered_flux(
              bxtmp, NVAR, amrex::The_Async_Arena());
            les_filter.apply_filter(
              bxtmp, flux[dir], filtered_flux, Density, NVAR);

            setV(bxtmp, flux[dir].nComp(), flx_arr[dir], 0.0);
            copy_array4(
              bxtmp, flux[dir].nComp(), filtered_flux.array(), flx_arr[dir]);
          }

          amrex::FArrayBox filtered_source_out(
            bx, NVAR, amrex::The_Async_Arena());
          les_filter.apply_filter(
            bx, hydro_source[mfi], filtered_source_out, Density, NVAR);

          setV(bx, hyd_src.nComp(), hyd_src, 0.0);
          copy_array4(
            bx, hyd_src.nComp(), filtered_source_out.array(), hyd_src);
        }

        // Refluxing
        if (do_reflux && sub_iteration == sub_ncycle - 1) {
          const amrex::FabType gtyp = flag_fab.getType(amrex::grow(bx, 1));
          update_flux_registers(
            dt, mfi, gtyp,
            {{AMREX_D_DECL(flux.data(), &(flux[1]), &(flux[2]))}}, dm_as_fine);
        }
      }
    }

    if (courno > 1.0) {
      amrex::Print() << "WARNING -- EFFECTIVE CFL AT THIS LEVEL " << level
                     << " IS " << courno << '\n';
      if (hard_cfl_limit) {
        amrex::Abort(
          "CFL is too high at this level -- go back to a checkpoint "
          "and restart with lower cfl number");
      }
    }
  }
}

void
pc_umdrv(
  const amrex::Real /*time*/,
  amrex::Box const& bx,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  amrex::Array4<const amrex::Real> const& uin,
  amrex::Array4<amrex::Real> const& uout,
  amrex::Array4<const amrex::Real> const& q,
  amrex::Array4<const amrex::Real> const& qaux,
  amrex::Array4<const amrex::Real> const&
    src_q, // amrex::IArrayBox const& bcMask,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const bool use_hybrid_weno,
  const int weno_scheme,
  const amrex::Real difmag,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    a,
  amrex::Array4<amrex::Real> const& vol,
  amrex::Real /*cflLoc*/)
{
  // Set Up for Hydro Flux Calculations
  auto const& bxg2 = grow(bx, 2);
  amrex::FArrayBox qec[AMREX_SPACEDIM];
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const amrex::Box eboxes = amrex::surroundingNodes(bxg2, dir);
    qec[dir].resize(eboxes, NGDNV, amrex::The_Async_Arena());
  }
  amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM> qec_arr{
    {AMREX_D_DECL(qec[0].array(), qec[1].array(), qec[2].array())}};

  // Temporary FArrayBoxes
  amrex::FArrayBox divu(bxg2, 1, amrex::The_Async_Arena());
  amrex::FArrayBox pdivu(bx, 1, amrex::The_Async_Arena());
  auto const& divuarr = divu.array();
  auto const& pdivuarr = pdivu.array();

  {
    BL_PROFILE("PeleC::umeth()");
#if AMREX_SPACEDIM == 1
    amrex::Abort("PLM isn't implemented in 1D.");
    // pc_umeth_1D(
    //  bx, bclo, bchi, domlo, domhi, q, qaux, src_q, // bcMask,
    //  flx[0], qec_arr[0], a[0], pdivuarr, vol, dx, dt);
#elif AMREX_SPACEDIM == 2
    pc_umeth_2D(
      bx, bclo, bchi, domlo, domhi, q, qaux, src_q, flx, qec_arr, a, pdivuarr,
      vol, dx, dt, ppm_type, plm_iorder, use_flattening, use_hybrid_weno,
      weno_scheme);
#elif AMREX_SPACEDIM == 3
    pc_umeth_3D(
      bx, bclo, bchi, domlo, domhi, q, qaux, src_q, flx, qec_arr, a, pdivuarr,
      vol, dx, dt, ppm_type, plm_iorder, use_flattening, use_hybrid_weno,
      weno_scheme);
#endif
  }

  // divu
  amrex::ParallelFor(bxg2, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_divu(i, j, k, q, AMREX_D_DECL(dx[0], dx[1], dx[2]), divuarr);
  });

  pc_adjust_fluxes(
    bx, uin, flx, a, divuarr, dx, domlo, domhi, bclo, bchi, difmag);

  // consup
  pc_consup(bx, uout, flx, vol, pdivuarr);
}

void
pc_adjust_fluxes(
  const amrex::Box& bx,
  const amrex::Array4<const amrex::Real>& u,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>&
    a,
  const amrex::Array4<const amrex::Real>& divu,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& del,
  const int* domlo,
  const int* domhi,
  const int* bclo,
  const int* bchi,
  const amrex::Real difmag)
{
  BL_PROFILE("PeleC::pc_adjust_fluxes()");
  // Flux alterations
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    amrex::Box const& fbx = surroundingNodes(bx, dir);
    const amrex::Real dx = del[dir];
    const int domlo_dir = domlo[dir];
    const int domhi_dir = domhi[dir];
    const int bclo_dir = bclo[dir];
    const int bchi_dir = bchi[dir];
    amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      pc_artif_visc(
        AMREX_D_DECL(i, j, k), flx[dir], divu, u, dx, difmag, dir, domlo_dir,
        domhi_dir, bclo_dir, bchi_dir);
      // Normalize Species Flux
      pc_norm_spec_flx(i, j, k, flx[dir]);
      // Make flux extensive
      pc_ext_flx(i, j, k, flx[dir], a[dir]);
    });
  }
}

void
pc_adjust_fluxes_eb(
  const amrex::Box& /*bx*/,
  amrex::Array4<const amrex::Real> const& q_arr,
  amrex::Array4<const amrex::Real> const& u_arr,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  amrex::Array4<const amrex::Real> const& vfrac,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flux,
  const amrex::Geometry& geom,
  const int* bclo,
  const int* bchi,
  amrex::Real difmag)
{
  BL_PROFILE("PeleC::pc_adjust_fluxes_eb()");
  amrex::Real areafac;

  // These are the fluxes on face centroids -- they are defined in
  // eb_compute_div
  //    and are the fluxes that go into the flux registers
  AMREX_D_TERM(
    auto const& fx_arr = flux[0];, auto const& fy_arr = flux[1];
    , auto const& fz_arr = flux[2];);

  // Compute divu to be used in artificial viscosity
  auto bx_divu = amrex::Box(q_arr);
  bx_divu.grow(-1);
  bx_divu.surroundingNodes();
  amrex::FArrayBox divu(bx_divu, NVAR, amrex::The_Async_Arena());
  divu.setVal<amrex::RunOn::Device>(0.0);

  amrex::Box nddom =
    amrex::convert(geom.growPeriodicDomain(16), amrex::IntVect(1));
  bx_divu &= nddom;

  auto const& divu_arr = divu.array();
  amrex::ParallelFor(
    bx_divu, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      eb_divu(i, j, k, q_arr, divu_arr, vfrac, dxinv);
    });

  // Flux alterations
  amrex::Box const& fbx = amrex::Box(fx_arr);
#if (AMREX_SPACEDIM == 2)
  areafac = dx[1];
#elif (AMREX_SPACEDIM == 3)
  areafac = dx[1] * dx[2];
#endif

  auto const* domlo = geom.Domain().loVect();
  auto const* domhi = geom.Domain().hiVect();

  int domlo_dir = domlo[0];
  int domhi_dir = domhi[0];
  int bclo_dir = bclo[0];
  int bchi_dir = bchi[0];

  amrex::ParallelFor(fbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_artif_visc(
      AMREX_D_DECL(i, j, k), fx_arr, divu_arr, u_arr, dx[0], difmag, 0,
      domlo_dir, domhi_dir, bclo_dir, bchi_dir);

    // Normalize Species Flux
    if (apx(i, j, k) > 0.) {
      pc_norm_spec_flx(i, j, k, fx_arr);
    }

    // Make flux extensive
    pc_ext_flx_eb(i, j, k, fx_arr, areafac, apx);
  });

  // Flux alterations
  amrex::Box const& fby = amrex::Box(fy_arr);
#if (AMREX_SPACEDIM == 2)
  areafac = dx[0];
#elif (AMREX_SPACEDIM == 3)
  areafac = dx[0] * dx[2];
#endif

  domlo_dir = domlo[1];
  domhi_dir = domhi[1];
  bclo_dir = bclo[1];
  bchi_dir = bchi[1];

  amrex::ParallelFor(fby, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_artif_visc(
      AMREX_D_DECL(i, j, k), fy_arr, divu_arr, u_arr, dx[1], difmag, 1,
      domlo_dir, domhi_dir, bclo_dir, bchi_dir);

    // Normalize Species Flux
    if (apy(i, j, k) > 0.) {
      pc_norm_spec_flx(i, j, k, fy_arr);
    }

    // Make flux extensive
    pc_ext_flx_eb(i, j, k, fy_arr, areafac, apy);
  });

#if (AMREX_SPACEDIM == 3)
  // Flux alterations
  amrex::Box const& fbz = amrex::Box(fz_arr);
  areafac = dx[1] * dx[0];

  domlo_dir = domlo[2];
  domhi_dir = domhi[2];
  bclo_dir = bclo[2];
  bchi_dir = bchi[2];

  amrex::ParallelFor(fbz, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_artif_visc(
      i, j, k, fz_arr, divu_arr, u_arr, dx[2], difmag, 2, domlo_dir, domhi_dir,
      bclo_dir, bchi_dir);

    // Normalize Species Flux
    if (apz(i, j, k) > 0.) {
      pc_norm_spec_flx(i, j, k, fz_arr);
    }

    // Make flux extensive
    pc_ext_flx_eb(i, j, k, fz_arr, areafac, apz);
  });

#endif
}

void
pc_consup(
  const amrex::Box& bx,
  const amrex::Array4<amrex::Real>& update,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const amrex::Array4<const amrex::Real>& vol,
  const amrex::Array4<const amrex::Real>& pdivu)
{
  BL_PROFILE("PeleC::pc_consup()");
  // Combine for Hydro Sources
  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    pc_update(i, j, k, update, flx, vol, pdivu);
  });
}

void
pc_consup_eb(
  const amrex::Box& bx,
  amrex::Array4<const amrex::Real> const& q_arr,
  amrex::Array4<const amrex::Real> const& qaux_arr,
  amrex::Array4<amrex::Real> const& divc_arr,
  amrex::Array4<amrex::Real> const& redistwgt_arr,
  AMREX_D_DECL(
    amrex::Array4<amrex::Real> const& q1,
    amrex::Array4<amrex::Real> const& q2,
    amrex::Array4<amrex::Real> const& q3),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& apx,
    amrex::Array4<amrex::Real const> const& apy,
    amrex::Array4<amrex::Real const> const& apz),
  AMREX_D_DECL(
    amrex::Array4<amrex::Real const> const& fcx,
    amrex::Array4<amrex::Real const> const& fcy,
    amrex::Array4<amrex::Real const> const& fcz),
  amrex::Array4<const amrex::Real> const& vfrac,
  amrex::Array4<amrex::EBCellFlag const> const& flag,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>&
    flux_tmp,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flux,

  const int eb_weights_type)
{
  const amrex::Box& bxg_i = amrex::Box(divc_arr);

  // These are the fluxes we computed in MOL_umeth_eb and modified in
  // adjust_fluxes_eb
  //  -- they live at face centers
  AMREX_D_TERM(
    auto const& fx_in_arr = flux_tmp[0];, auto const& fy_in_arr = flux_tmp[1];
    , auto const& fz_in_arr = flux_tmp[2];);

  // These are the fluxes on face centroids -- they are defined in
  // eb_compute_div
  //    and are the fluxes that go into the flux registers
  AMREX_D_TERM(
    auto const& fx_out_arr = flux[0];, auto const& fy_out_arr = flux[1];
    , auto const& fz_out_arr = flux[2];);

  auto const& blo = bx.smallEnd();
  auto const& bhi = bx.bigEnd();

  amrex::ParallelFor(
    bxg_i, NVAR, [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
      // This does the divergence but not the redistribution -- we will do that
      // later We do compute the weights here though
      if (!flag(i, j, k).isCovered()) {
        eb_compute_div(
          i, j, k, n, blo, bhi, q_arr, qaux_arr, divc_arr,
          AMREX_D_DECL(fx_in_arr, fy_in_arr, fz_in_arr),
          AMREX_D_DECL(fx_out_arr, fy_out_arr, fz_out_arr), flag, vfrac,
          redistwgt_arr, AMREX_D_DECL(apx, apy, apz),
          AMREX_D_DECL(fcx, fcy, fcz), dxinv, eb_weights_type);
      }
    });

  // Update UEINT update with pdivu term
  amrex::ParallelFor(bxg_i, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    if (!flag(i, j, k).isCovered()) {
      pc_eb_pdivu(
        i, j, k, q_arr, AMREX_D_DECL(q1, q2, q3), divc_arr, flag, vfrac,
        AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), dxinv);
    }
  });
}

void
pc_umdrv_eb(
  const amrex::Box& bx,
  const amrex::Box& bxg_i,
  const amrex::MFIter& mfi,
  const amrex::Geometry& geom,
  const amrex::EBFArrayBoxFactory* fact,
  const int* bclo,
  const int* bchi,
  const amrex::Array4<const amrex::Real>& uin,
  const amrex::Array4<amrex::Real>& uout,
  const amrex::Array4<const amrex::Real>& q,
  const amrex::Array4<const amrex::Real>& qaux,
  const amrex::Array4<const amrex::Real>& src_q,
  const amrex::Array4<const amrex::Real>& vf,
  const amrex::Array4<amrex::EBCellFlag const>& flag,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dx,
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM>& dxInv,
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>& flx,
  const int as_crse,
  const amrex::Array4<amrex::Real>& drho_as_crse,
  const amrex::Array4<int const>& rrflag_as_crse,
  const int as_fine,
  const amrex::Array4<amrex::Real>& dm_as_fine,
  const amrex::Array4<int const>& lev_mask,
  const amrex::Real dt,
  const int ppm_type,
  const int plm_iorder,
  const bool use_flattening,
  const amrex::Real difmag,
  amrex::BCRec const* bcs_d_ptr,
  const std::string& redistribution_type,
  const int eb_weights_type,
  const int eb_srd_max_order,
  const bool eb_clean_massfrac,
  const amrex::Real eb_clean_massfrac_threshold,
  amrex::Real /*cflLoc*/)
{
  BL_PROFILE("PeleC::pc_umdrv_eb()");

  const amrex::Box& bxg_ii = grow(bxg_i, 1);

  amrex::Array4<amrex::Real const> AMREX_D_DECL(fcx, fcy, fcz),
    AMREX_D_DECL(apx, apy, apz), ccc;
  AMREX_D_TERM(
    fcx = fact->getFaceCent()[0]->const_array(mfi);
    , fcy = fact->getFaceCent()[1]->const_array(mfi);
    , fcz = fact->getFaceCent()[2]->const_array(mfi););
  AMREX_D_TERM(
    apx = fact->getAreaFrac()[0]->const_array(mfi);
    , apy = fact->getAreaFrac()[1]->const_array(mfi);
    , apz = fact->getAreaFrac()[2]->const_array(mfi););
  ccc = fact->getCentroid().const_array(mfi);

  const auto* domlo = geom.Domain().loVect();
  const auto* domhi = geom.Domain().hiVect();

  // Temporary FArrayBoxes
  amrex::FArrayBox divu(bxg_ii, 1, amrex::The_Async_Arena());
  auto const& divuarr = divu.array();

  amrex::FArrayBox qec[AMREX_SPACEDIM];
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    const amrex::Box eboxes = amrex::surroundingNodes(grow(bxg_ii, 1), dir);
    qec[dir].resize(eboxes, NGDNV, amrex::The_Async_Arena());
  }
  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
    qec_arr{{AMREX_D_DECL(qec[0].array(), qec[1].array(), qec[2].array())}};

  // Quantities for redistribution
  amrex::FArrayBox divc, redistwgt;
  if (redistribution_type == "StateRedist") {
    divc.resize(
      bxg_i, NVAR,
      amrex::The_Async_Arena()); // This will hold "dUdt" before redistribution
    redistwgt.resize(
      bxg_i, NVAR, amrex::The_Async_Arena()); // This will be "scratch" which
                                              // holds "Uold + dt*dUdt"
  } else {
    divc.resize(
      bxg_i, NVAR,
      amrex::The_Async_Arena()); // This will hold "dUdt" before redistribution
    redistwgt.resize(
      bxg_i, 1, amrex::The_Async_Arena()); // This will hold the weights used in
                                           // flux redistribution
  }
  divc.setVal<amrex::RunOn::Device>(0.0);
  redistwgt.setVal<amrex::RunOn::Device>(0.0);

  // Because we are going to redistribute, we put the divergence into divc
  //    rather than directly into dsdt_arr
  auto const& divc_arr = divc.array();
  auto const& redistwgt_arr = redistwgt.array();

  // We need one extra in flux_tmp so we can tangentially interpolate fluxes
  amrex::FArrayBox flux_tmp[AMREX_SPACEDIM];
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    flux_tmp[idim].resize(
      amrex::surroundingNodes(bxg_ii, idim), NVAR, amrex::The_Async_Arena());
    flux_tmp[idim].setVal<amrex::RunOn::Device>(0.);
  }

  const amrex::GpuArray<const amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
    flux_tmp_arr{{AMREX_D_DECL(
      flux_tmp[0].array(), flux_tmp[1].array(), flux_tmp[2].array())}};
  const amrex::GpuArray<const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
    ap{{AMREX_D_DECL(apx, apy, apz)}};

  // Define divc, the update before redistribution
  // Also construct the redistribution weights for flux redistribution if
  // necessary
#if AMREX_SPACEDIM == 1
  amrex::Abort("PLM isn't implemented in 1D.");
#elif AMREX_SPACEDIM == 2
  pc_umeth_eb_2D(
    amrex::Box(divc_arr), bclo, bchi, domlo, domhi, q, qaux, src_q,
    flux_tmp_arr, qec_arr, ap, flag, dx, dt, ppm_type, use_flattening,
    plm_iorder);
#elif AMREX_SPACEDIM == 3
  pc_umeth_eb_3D(
    amrex::Box(divc_arr), bclo, bchi, domlo, domhi, q, qaux, src_q,
    flux_tmp_arr, qec_arr, ap, flag, dx, dt, ppm_type, use_flattening,
    plm_iorder);
#endif

  // Construct divu
  AMREX_D_TERM(
    const amrex::Real dx0 = dx[0];, const amrex::Real dx1 = dx[1];
    , const amrex::Real dx2 = dx[2];);

  amrex::ParallelFor(
    bxg_ii, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      if (flag(i, j, k).isRegular()) {
        pc_divu(i, j, k, q, AMREX_D_DECL(dx0, dx1, dx2), divuarr);
      } else {
        divuarr(i, j, k) = amrex::Real(0.0);
      }
    });

  pc_adjust_fluxes_eb(
    bx, q, uin, AMREX_D_DECL(apx, apy, apz), vf, dx, dxInv, flux_tmp_arr, geom,
    bclo, bchi, difmag);

  pc_consup_eb(
    bx, q, qaux, divc_arr, redistwgt_arr,
    AMREX_D_DECL(qec_arr[0], qec_arr[1], qec_arr[2]),
    AMREX_D_DECL(apx, apy, apz), AMREX_D_DECL(fcx, fcy, fcz), vf, flag, dxInv,
    flux_tmp_arr, flx, eb_weights_type);

  const int l_ncomp = uout.nComp();
  const int level_mask_not_covered = constants::level_mask_notcovered();
  const bool use_wts_in_divnc = false;

  const amrex::Real fac_for_redist = 1.0;
  {
    BL_PROFILE("ApplyMLRedistribution()");
    ApplyMLRedistribution(
      bx, l_ncomp, uout, divc_arr, uin, redistwgt_arr, flag,
      AMREX_D_DECL(apx, apy, apz), vf, AMREX_D_DECL(fcx, fcy, fcz), ccc,
      bcs_d_ptr, geom, dt, redistribution_type, as_crse, drho_as_crse,
      rrflag_as_crse, as_fine, dm_as_fine, lev_mask, level_mask_not_covered,
      fac_for_redist, use_wts_in_divnc, 0, eb_srd_max_order);
  }
  auto const& flags = fact->getMultiEBCellFlagFab();
  const auto& flag_fab = flags[mfi];
  amrex::FabType typ = flag_fab.getType(bx);
  pc_post_eb_redistribution(
    bx, dt, eb_clean_massfrac, eb_clean_massfrac_threshold, uin, typ, flag,
    redistwgt_arr, uout);
}
