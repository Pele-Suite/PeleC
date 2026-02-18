#include "Diffusion.H"
#include "prob.H"

void
PeleC::getMOLSrcTerm(
  const amrex::MultiFab& S,
  amrex::MultiFab& MOLSrcTerm,
  const amrex::Real /*time*/,
  const amrex::Real dt,
  const amrex::Real reflux_factor)
{
  BL_PROFILE("PeleC::getMOLSrcTerm()");
  if (
    (!diffuse_temp) && (!diffuse_enth) && (!diffuse_spec) && (!diffuse_vel) &&
    (!do_hydro)) {
    MOLSrcTerm.setVal(0, 0, NVAR, MOLSrcTerm.nGrow());
    return;
  }

  bool using_rf = do_rf;
  amrex::Real omega = rf_omega;
  int axis = rf_axis;
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> axis_loc = {
    AMREX_D_DECL(rf_axis_x, rf_axis_y, rf_axis_z)};
  const auto dx = geom.CellSizeArray();

  /*
     Across all conserved state components, compute the method of lines rhs
     = -Div(Flux). The input state, S, contained the conserved variables, and
     is "fill patched" in the usual AMReX way, where values at Dirichlet
     boundaries actually are assumed to live on the inflow face.

     1. Convert S to Q, primitive variables (since the transport coefficients
     typically depend on mass fractions and temperature).

     2. Evaluate transport coefficients

     3. Evaluate the diffusion operator over all components

     a. Compute face-centered transport coefficients
     b. Evaluate face-centered diffusion fluxes
     c. Zero fluxes and divergence after the fact if subset of components have
     diffuse shut off (this allows that T be diffused alone)

     4. (optional) evaluate face-centered MOL hydro fluxes

     5. Compute divergence of face-centered fluxes with hybrid divergence
     operator, include wall fluxes for diffusion and hydro,

     6. Perform weighted redistribution at the EB

     Extra notes:

     A. The face-based transport coefficients that are computed with face-based
     Fill-Patched data at Dirichlet boundaries.

     B. Within the routine that computes diffusion fluxes, there is a need for
     computing species enthalpies at cell faces.  In the EGLib model, species
     enthalpies are a function of temperature.  At the moment, in diffterm, we
     evaluate enthalpies at cell centers and then take the face values to be the
     arithmetic average of cell values on either side.  Similarly, mass
     fractions at cell faces are needed to compute the barodiffusion and
     correction velocity expressions. Arithmetic averages are used there as
     well.  Thus, these face values are thermodynamically inconsistent.  Not
     sure what are the consequences of that.
  */

  const int nCompTr = dComp_lambda + 1;
  const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dxinv =
    geom.InvCellSizeArray();

  auto const& fact =
    dynamic_cast<amrex::EBFArrayBoxFactory const&>(S.Factory());
  auto const& flags = fact.getMultiEBCellFlagFab();
  amrex::MultiFab* cost = nullptr;

  if (do_mol_load_balance) {
    cost = &(get_new_data(Work_Estimate_Type));
  }

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
  {
    for (amrex::MFIter mfi(MOLSrcTerm, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box vbox = mfi.tilebox();
      int ng = numGrow();
      const amrex::Box gbox = amrex::grow(vbox, ng);
      const amrex::Box cbox = amrex::grow(vbox, ng - 1);
      auto const& MOLSrc = MOLSrcTerm.array(mfi);

      amrex::Real wt = amrex::ParallelDescriptor::second();
      const auto& flag_fab = flags[mfi];
      amrex::FabType typ = flag_fab.getType(vbox);
      if (typ == amrex::FabType::covered) {
        setV(vbox, NVAR, MOLSrc, 0);
        if (do_mol_load_balance && (cost != nullptr)) {
          wt = (amrex::ParallelDescriptor::second() - wt) / vbox.d_numPts();
          (*cost)[mfi].plus<amrex::RunOn::Device>(wt, vbox);
        }
        continue;
      }
      // Note on typ: if interior cells (vbox) are all covered, no need to
      // do anything. But otherwise, we need to do EB stuff if there are any
      // cut cells within 1 grow cell (cbox) due to EB redistribute
      typ = flag_fab.getType(cbox);

      const amrex::Box ebfluxbox = amrex::grow(vbox, 3);

      const int local_i = mfi.LocalIndex();
      const auto Ncut =
        (!eb_in_domain)
          ? 0
          : static_cast<int>(sv_eb_bndry_grad_stencil[local_i].size());
      SparseData<amrex::Real, EBBndrySten> eb_flux_thdlocal;
      if (Ncut > 0) {
        eb_flux_thdlocal.define(sv_eb_bndry_grad_stencil[local_i], NVAR);
      }
      auto* d_sv_eb_bndry_geom =
        (Ncut > 0 ? sv_eb_bndry_geom[local_i].data() : nullptr);

      const int nqaux = NQAUX > 0 ? NQAUX : 1;
      amrex::FArrayBox q(gbox, QVAR, amrex::The_Async_Arena());
      amrex::FArrayBox qaux(gbox, nqaux, amrex::The_Async_Arena());
      amrex::FArrayBox coeff_cc(gbox, nCompTr, amrex::The_Async_Arena());
      auto const& sar = S.array(mfi);
      auto const& qar = q.array();
      auto const& qauxar = qaux.array();

      // Get primitives, Q, including (Y, T, p, rho) from conserved state
      {
        const auto geomdata = geom.data();
        BL_PROFILE("PeleC::ctoprim()");
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            if (using_rf) {
              amrex::IntVect iv(AMREX_D_DECL(i, j, k));
              amrex::Real rad = get_rotaxis_dist(iv, axis, axis_loc, geomdata);
              pc_ctoprim(i, j, k, sar, qar, qauxar, omega, rad);
            } else {
              pc_ctoprim(i, j, k, sar, qar, qauxar);
            }
          });
      }
      // TODO deal with NSCBC
      /*
         for (int dir = 0; dir < AMREX_SPACEDIM ; dir++)  {
         const amrex::Box& bxtmp = amrex::surroundingNodes(vbox,dir);
         amrex::Box TestBox(bxtmp);
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

      if (nscbc_diff == 1)
      {
      impose_NSCBC(lo, hi, domain_lo, domain_hi,
      AMREX_TO_FORTRAN(Sfab),
      AMREX_TO_FORTRAN(q.fab()),
      AMREX_TO_FORTRAN(qaux.fab()),
      AMREX_D_DECL(AMREX_TO_FORTRAN(bcMask[0]),
      AMREX_TO_FORTRAN(bcMask[1]),
      AMREX_TO_FORTRAN(bcMask[2])),
      &flag_nscbc_isAnyPerio, flag_nscbc_perio,
      &time, dx, &dt);
      }
      */

      // Compute transport coefficients, coincident with Q
      auto const& coe_cc = coeff_cc.array();
      {
        auto const& qar_yin = q.array(QFS);
        auto const& qar_Tin = q.array(QTEMP);
        auto const& qar_rhoin = q.array(QRHO);
        auto const& coe_rhoD = coeff_cc.array(dComp_rhoD);
        auto const& coe_mu = coeff_cc.array(dComp_mu);
        auto const& coe_xi = coeff_cc.array(dComp_xi);
        auto const& coe_lambda = coeff_cc.array(dComp_lambda);
        BL_PROFILE("PeleC::get_transport_coeffs()");
        auto const* ltransparm = trans_parms.device_parm();
        auto const& geomdata = geom.data();
        const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
        const bool get_xi = true, get_mu = true, get_lam = true,
                   get_Ddiag = true, get_chi = false;
        amrex::ParallelFor(
          gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            amrex::Real muloc, xiloc, lamloc;
            amrex::Real Ddiag[NUM_SPECIES], Y[NUM_SPECIES] = {0.0};
            amrex::Real* chi_mix = nullptr;
            amrex::Real T = qar_Tin(i, j, k);
            amrex::Real rho = qar_rhoin(i, j, k);
            for (int n = 0; n < NUM_SPECIES; ++n) {
              Y[n] = qar_yin(i, j, k, n);
            }

            const amrex::RealVect x =
              pc_cmp_loc({AMREX_D_DECL(i, j, k)}, geomdata);
            pc_transcoeff(
              get_xi, get_mu, get_lam, get_Ddiag, get_chi, T, rho, Y, Ddiag,
              chi_mix, muloc, xiloc, lamloc, ltransparm, *lprobparm, x);

            for (int n = 0; n < NUM_SPECIES; ++n) {
              coe_rhoD(i, j, k, n) = Ddiag[n];
            }
            coe_mu(i, j, k) = muloc;
            coe_xi(i, j, k) = xiloc;
            coe_lambda(i, j, k) = lamloc;
          });
      }

      amrex::FArrayBox flux_ec[AMREX_SPACEDIM];
      const amrex::Box eboxes[AMREX_SPACEDIM] = {AMREX_D_DECL(
        amrex::surroundingNodes(cbox, 0), amrex::surroundingNodes(cbox, 1),
        amrex::surroundingNodes(cbox, 2))};
      amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM> flx;
      const amrex::GpuArray<
        const amrex::Array4<const amrex::Real>, AMREX_SPACEDIM>
        area_arr{{AMREX_D_DECL(
          area[0].array(mfi), area[1].array(mfi), area[2].array(mfi))}};
      for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
        flux_ec[dir].resize(eboxes[dir], NVAR, amrex::The_Async_Arena());
        flx[dir] = flux_ec[dir].array();
        setV(eboxes[dir], NVAR, flx[dir], 0);
      }

      amrex::FArrayBox Dfab(cbox, NVAR, amrex::The_Async_Arena());
      auto const& Dterm = Dfab.array();
      setV(cbox, NVAR, Dterm, 0.0);
      auto flag_arr = flags.const_array(mfi);

      {
        // Compute Extensive diffusion fluxes for X, Y, Z
        BL_PROFILE("PeleC::diffusion_flux()");
        const bool l_transport_harmonic_mean = transport_harmonic_mean;
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          if (
            (typ == amrex::FabType::singlevalued) ||
            (typ == amrex::FabType::regular)) {
            amrex::ParallelFor(
              eboxes[dir], [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                amrex::GpuArray<amrex::Real, dComp_lambda + 1> cf = {0.0};
                if (
                  flag_arr(i, j, k).isRegular() ||
                  flag_arr(i, j, k).isSingleValued()) {
                  for (int n = 0; n < static_cast<int>(cf.size()); n++) {
                    pc_move_transcoefs_to_ec(
                      AMREX_D_DECL(i, j, k), n, coe_cc, cf.data(), dir,
                      l_transport_harmonic_mean);
                  }
                }
                if (typ == amrex::FabType::singlevalued) {
                  pc_diffusion_flux_eb(
                    i, j, k, qar, cf, flag_arr, area_arr[dir], flx[dir], dxinv,
                    dir);
                } else if (typ == amrex::FabType::regular) {
                  pc_diffusion_flux(
                    i, j, k, qar, cf, area_arr[dir], flx[dir], dxinv, dir);
                }
              });
          } else if (typ == amrex::FabType::multivalued) {
            amrex::Abort("multi-valued cells are not supported");
          }
        }
      }

      if (do_isothermal_walls) {
        // Compute extensive diffusion flux at domain boundaries
        BL_PROFILE("PeleC::isothermal_wall_fluxes()");
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          if (
            (typ == amrex::FabType::singlevalued) ||
            (typ == amrex::FabType::regular)) {
            int normalarr[2] = {-1, 1};
            amrex::Real bc_temp_arr[2] = {
              domlo_isothermal_temp[dir], domhi_isothermal_temp[dir]};
            for (int inorm = 0; inorm < 2; inorm++) {
              int normal = normalarr[inorm];
              amrex::Real bc_temp = bc_temp_arr[inorm];
              if (bc_temp > 0.0) {
                amrex::Box bbox = surroundingNodes(vbox, dir);
                if (normal == -1) {
                  bbox.setBig(dir, geom.Domain().smallEnd(dir));
                } else {
                  bbox.setSmall(dir, geom.Domain().bigEnd(dir) + 1);
                }
                if (bbox.ok()) {
                  amrex::FArrayBox tmpfabtemp(
                    bbox, 1, amrex::The_Async_Arena());
                  amrex::Array4<amrex::Real> temp_arr = tmpfabtemp.array();
                  const ProbParmDevice* lprobparm = PeleC::d_prob_parm_device;
                  auto const* ltransparm = trans_parms.device_parm();
                  const auto geomdata = geom.data();
                  amrex::ParallelFor(
                    bbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                      ProblemSpecificFunctions::set_isothermal_wall_temperature(
                        i, j, k, dir, normal, bc_temp, geomdata, *lprobparm,
                        qar, temp_arr);
                      pc_isothermal_wall_fluxes(
                        i, j, k, dir, normal, qar, temp_arr, flag_arr,
                        area_arr[dir], flx[dir], geomdata, ltransparm,
                        *lprobparm);
                    });
                }
              }
            }
          }
        }
      }

      // Shut off unwanted diffusion after the fact.
      //      Under normal conditions, you either have diffusion on all or
      //      none, so this shouldn't be done this way.  However, the regression
      //      test for diffusion works by diffusing only temperature through
      //      this process.  Ideally, we'd redo that test to diffuse a passive
      //      scalar instead....
      if ((!diffuse_temp) && (!diffuse_enth)) {
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], Eden, Eint, flx[dir], 0.0);
        }
      }
      if (!diffuse_spec) {
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], FirstSpec, FirstSpec + NUM_SPECIES, flx[dir], 0.0);
        }
      }
      if (!diffuse_vel) {
        for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
          setC(eboxes[dir], Xmom, Xmom + 3, flx[dir], 0.0);
        }
      }

      // Compute and add in the hydro fluxes.
      if (do_hydro && do_mol) {
        // amrex::FArrayBox flatn(cbox, 1, amrex::The_Async_Arena());
        // flatn.setVal(1.0); // Set flattening to 1.0

        // If filtering, save off the diffusion fluxes (don't want to filter
        // these)
        amrex::FArrayBox diffusion_flux[AMREX_SPACEDIM];
        amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
          diffusion_flux_arr;
        if (use_explicit_filter) {
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            diffusion_flux[dir].resize(
              flux_ec[dir].box(), NVAR, amrex::The_Async_Arena());
            diffusion_flux_arr[dir] = diffusion_flux[dir].array();
            copy_array4(
              flux_ec[dir].box(), flux_ec[dir].nComp(), flx[dir],
              diffusion_flux_arr[dir]);
          }
        }

        { // Get face-centered hyperbolic fluxes
          BL_PROFILE("PeleC::pc_hyp_mol_flux()");
          pc_compute_hyp_mol_flux(
            cbox, qar, qauxar, flx, area_arr, plm_iorder, use_laxf_flux,
            flags.array(mfi), geom, axis_loc, omega, axis, using_rf);
        }

        // Filter hydro fluxes
        if (use_explicit_filter) {
          // Get the hydro term
          amrex::FArrayBox hydro_flux[AMREX_SPACEDIM];
          amrex::GpuArray<amrex::Array4<amrex::Real>, AMREX_SPACEDIM>
            hydro_flux_arr;
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            hydro_flux[dir].resize(
              flux_ec[dir].box(), NVAR, amrex::The_Async_Arena());
            hydro_flux_arr[dir] = hydro_flux[dir].array();
            lincomb_array4(
              flux_ec[dir].box(), Density, NVAR, flx[dir],
              diffusion_flux_arr[dir], 1.0, -1.0, hydro_flux_arr[dir]);
          }

          // Filter
          const amrex::Box fbox = amrex::grow(cbox, -nGrowF);
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            const amrex::Box& bxtmp = amrex::surroundingNodes(fbox, dir);
            amrex::FArrayBox filtered_hydro_flux(
              bxtmp, NVAR, amrex::The_Async_Arena());
            les_filter.apply_filter(
              bxtmp, hydro_flux[dir], filtered_hydro_flux, Density, NVAR);

            setV(bxtmp, hydro_flux[dir].nComp(), hydro_flux_arr[dir], 0.0);
            copy_array4(
              bxtmp, hydro_flux[dir].nComp(), filtered_hydro_flux.array(),
              hydro_flux_arr[dir]);
          }

          // Combine with diffusion
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            lincomb_array4(
              diffusion_flux[dir].box(), Density, NVAR, diffusion_flux_arr[dir],
              hydro_flux_arr[dir], 1.0, 1.0, flx[dir]);
          }
        }
      }

      // Compute divergence and refluxing
      if (typ == amrex::FabType::singlevalued) {

        // Set extensive diffusive flux at embedded boundary, potentially
        // non-zero only for heat flux on isothermal boundaries,
        // and momentum fluxes at no-slip walls
        const auto nFlux =
          sv_eb_flux.empty() ? 0 : sv_eb_flux[local_i].numPts();
        if (Ncut > 0) {
          eb_flux_thdlocal.setVal(0); // Default to Neumann for all fields

          const auto Nvals = sv_eb_bcval[local_i].numPts();

          AMREX_ASSERT(Nvals == Ncut);
          AMREX_ASSERT(nFlux == Ncut);

          if (eb_isothermal && (diffuse_temp || diffuse_enth)) {
            {
              BL_PROFILE("PeleC::pc_apply_eb_boundry_flux_stencil()");
              pc_apply_eb_boundry_flux_stencil(
                ebfluxbox, sv_eb_bndry_grad_stencil[local_i].data(), Ncut, qar,
                QTEMP, coe_cc, dComp_lambda,
                sv_eb_bcval[local_i].dataPtr(QTEMP), Nvals,
                eb_flux_thdlocal.dataPtr(Eden), nFlux, 1);
            }
          }
          // Compute momentum transfer at no-slip EB wall
          if (eb_noslip && diffuse_vel) {
            {
              BL_PROFILE("PeleC::pc_apply_eb_boundry_visc_flux_stencil()");
              pc_apply_eb_boundry_visc_flux_stencil(
                ebfluxbox, sv_eb_bndry_grad_stencil[local_i].data(), Ncut,
                d_sv_eb_bndry_geom, Ncut, qar, coe_cc,
                sv_eb_bcval[local_i].dataPtr(QU), Nvals,
                eb_flux_thdlocal.dataPtr(Xmom), nFlux);
            }
          }
          if (do_hydro && do_mol) {
            { // Get hyp flux at EB wall
              BL_PROFILE("PeleC::pc_hyp_mol_flux_eb()");
              amrex::Real* d_eb_flux_thdlocal =
                (nFlux > 0 ? eb_flux_thdlocal.dataPtr() : nullptr);
              pc_compute_hyp_mol_flux_eb(
                geom, cbox, qar, qauxar, dx, use_laxf_flux, eb_problem_state,
                vfrac.array(mfi), d_sv_eb_bndry_geom, Ncut, d_eb_flux_thdlocal,
                nFlux);
            }
          }
        }

        amrex::Gpu::DeviceVector<int> v_eb_tile_mask(Ncut, 0);
        int* eb_tile_mask = v_eb_tile_mask.dataPtr();
        amrex::ParallelFor(Ncut, [=] AMREX_GPU_DEVICE(int icut) {
          if (ebfluxbox.contains(d_sv_eb_bndry_geom[icut].iv)) {
            eb_tile_mask[icut] = 1;
          }
        });
        if (typ == amrex::FabType::singlevalued && Ncut > 0) {
          sv_eb_flux[local_i].merge(eb_flux_thdlocal, 0, NVAR, v_eb_tile_mask);
        }

        // Interpolate fluxes from face centers to face centroids
        // Note that hybrid divergence and redistribution algorithms require
        // that we be able to compute the conservative divergence on 2 grow
        // cells, so we need interpolated fluxes on 2 grow cells, and
        // therefore we need face centered fluxes on 3.
        {
          BL_PROFILE("PeleC::pc_apply_face_stencil()");
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            const auto Nsten =
              static_cast<int>(flux_interp_stencil[dir][local_i].size());
            const amrex::Box valid_interped_flux_box =
              amrex::Box(ebfluxbox).surroundingNodes(dir);
            if (Nsten > 0) {
              pc_apply_face_stencil(
                valid_interped_flux_box, stencil_volume_box,
                flux_interp_stencil[dir][local_i].data(), Nsten, dir, NVAR,
                flx[dir]);
            }
          }
          amrex::Gpu::Device::streamSynchronize();
        }

        // Compute flux divergence (1/Vol).Div(F.A)
        {
          BL_PROFILE("PeleC::pc_flux_div()");
          auto const& vol = volume.array(mfi);
          amrex::ParallelFor(
            cbox, NVAR,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              pc_flux_div(
                i, j, k, n, AMREX_D_DECL(flx[0], flx[1], flx[2]), vol, Dterm);
            });
        }

        // Get "hybrid flux divergence"
        //
        // This operation takes as input centroid-centered fluxes and a
        // corresponding
        //  divergence on three grid cells.  Actually, we assume that
        //  div=(1/VOL)Div(flux) (VOL = volume of the full cells), and that
        //  flux is EXTENSIVE, weighted with the full face areas.
        //
        // Upon return:
        // div = kappa.(1/Vol) Div(FluxC.Area)  Vol = kappa.VOL,
        // Area=aperture.AREA, defined over the valid box

        // TODO: Rework this for r-z, if applicable
        amrex::Real vol = 1;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          vol *= geom.CellSize()[dir];
        }

        if (Ncut > 0) {
          BL_PROFILE("PeleC::pc_eb_div()");
          pc_eb_div(
            vbox, vol, NVAR, d_sv_eb_bndry_geom, Ncut,
            AMREX_D_DECL(flx[0], flx[1], flx[2]), sv_eb_flux[local_i].dataPtr(),
            vfrac.array(mfi), Dterm);
        }
      } else if (typ == amrex::FabType::regular) {
        // Compute flux divergence (1/Vol).Div(F.A)
        {
          BL_PROFILE("PeleC::pc_flux_div()");
          auto const& vol = volume.array(mfi);
          amrex::ParallelFor(
            cbox, NVAR,
            [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
              pc_flux_div(
                i, j, k, n, AMREX_D_DECL(flx[0], flx[1], flx[2]), vol, Dterm);
            });
        }
      } else if (typ == amrex::FabType::multivalued) {
        amrex::Abort("multi-valued eb boundary fluxes to be implemented");
      }

      // Extrapolate to GhostCells
      if (MOLSrcTerm.nGrow() > 0) {
        BL_PROFILE("PeleC::diffextrap()");
        const int mg = MOLSrcTerm.nGrow();
        const auto* low = vbox.loVect();
        const auto* high = vbox.hiVect();
        auto dlo = amrex::lbound(Dterm);
        auto dhi = amrex::ubound(Dterm);
        const int AMREX_D_DECL(lx = low[0], ly = low[1], lz = low[2]);
        const int AMREX_D_DECL(hx = high[0], hy = high[1], hz = high[2]);
        amrex::ParallelFor(
          vbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            pc_diffextrap(
              i, j, k, Dterm, mg, UMX, UMZ + 1, AMREX_D_DECL(lx, ly, lz),
              AMREX_D_DECL(hx, hy, hz), dlo, dhi);
            pc_diffextrap(
              i, j, k, Dterm, mg, UFS, UFS + NUM_SPECIES,
              AMREX_D_DECL(lx, ly, lz), AMREX_D_DECL(hx, hy, hz), dlo, dhi);
            pc_diffextrap(
              i, j, k, Dterm, mg, UEDEN, UEDEN + 1, AMREX_D_DECL(lx, ly, lz),
              AMREX_D_DECL(hx, hy, hz), dlo, dhi);
          });
      }

      // EB redistribution
      amrex::FArrayBox dm_as_fine(
        amrex::Box::TheUnitBox(), MOLSrcTerm.nComp(), amrex::The_Async_Arena());
      if (eb_in_domain && (typ != amrex::FabType::regular)) {
        AMREX_D_TERM(
          auto apx = areafrac[0]->const_array(mfi);
          , auto apy = areafrac[1]->const_array(mfi);
          , auto apz = areafrac[2]->const_array(mfi););
        AMREX_D_TERM(
          auto fcx = facecent[0]->const_array(mfi);
          , auto fcy = facecent[1]->const_array(mfi);
          , auto fcz = facecent[2]->const_array(mfi););
        auto ccc = fact.getCentroid().const_array(mfi);

        amrex::FArrayBox tmpfab(
          Dfab.box(), S.nComp(), amrex::The_Async_Arena());
        if (redistribution_type == "FluxRedist") {
          tmpfab.setVal<amrex::RunOn::Device>(1.0, tmpfab.box());
        }
        amrex::Array4<amrex::Real> scratch = tmpfab.array();

        amrex::FArrayBox Dterm_tmpfab(
          Dfab.box(), S.nComp(), amrex::The_Async_Arena());
        amrex::Array4<amrex::Real> Dterm_tmp = Dterm_tmpfab.array();
        copy_array4(Dfab.box(), NVAR, Dterm, Dterm_tmp);

        const amrex::StateDescriptor* desc = state[State_Type].descriptor();
        const auto& bcs = desc->getBCs();
        amrex::Gpu::DeviceVector<amrex::BCRec> d_bcs(desc->nComp());
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice, bcs.begin(), bcs.end(), d_bcs.begin());

        amrex::EBFluxRegister* fr_as_crse =
          (do_reflux && (level < parent->finestLevel()))
            ? &getFluxReg(level + 1)
            : nullptr;
        amrex::EBFluxRegister* fr_as_fine =
          (do_reflux && (level > 0)) ? &getFluxReg(level) : nullptr;

        const int as_crse = static_cast<int>(fr_as_crse != nullptr);
        const int as_fine = static_cast<int>(fr_as_fine != nullptr);

        amrex::FArrayBox fab_drho_as_crse(
          amrex::Box::TheUnitBox(), MOLSrcTerm.nComp(),
          amrex::The_Async_Arena());
        amrex::IArrayBox fab_rrflag_as_crse(
          amrex::Box::TheUnitBox(), 1, amrex::The_Async_Arena());

        auto* drho_as_crse = (fr_as_crse != nullptr)
                               ? fr_as_crse->getCrseData(mfi)
                               : &fab_drho_as_crse;
        const auto* rrflag_as_crse = (fr_as_crse != nullptr)
                                       ? fr_as_crse->getCrseFlag(mfi)
                                       : &fab_rrflag_as_crse;

        if (fr_as_fine != nullptr) {
          const amrex::Box dbox1 = geom.growPeriodicDomain(1);
          const amrex::Box bx_for_dm(amrex::grow(vbox, 1) & dbox1);
          dm_as_fine.resize(bx_for_dm, MOLSrcTerm.nComp());
          dm_as_fine.setVal<amrex::RunOn::Device>(0.0);
        }

        const bool use_wts_in_divnc = false;

        const int level_mask_not_covered = constants::level_mask_notcovered();

        {
          BL_PROFILE("ApplyMLRedistribution()");
          const amrex::Real fac_for_redist = (do_mol) ? 0.5 : 1.0;
          ApplyMLRedistribution(
            vbox, S.nComp(), Dterm, Dterm_tmp, S.const_array(mfi), scratch,
            flag_arr, AMREX_D_DECL(apx, apy, apz), vfrac.const_array(mfi),
            AMREX_D_DECL(fcx, fcy, fcz), ccc, d_bcs.dataPtr(), geom, dt,
            redistribution_type, as_crse, drho_as_crse->array(),
            rrflag_as_crse->array(), as_fine, dm_as_fine.array(),
            level_mask.const_array(mfi), level_mask_not_covered, fac_for_redist,
            use_wts_in_divnc, 0, eb_srd_max_order);
        }

        pc_post_eb_redistribution(
          vbox, dt, eb_clean_massfrac, eb_clean_massfrac_threshold,
          S.const_array(mfi), typ, flag_arr, scratch, Dterm);
      }

      // Refluxing
      if (do_reflux && reflux_factor != 0) {
        if (typ == amrex::FabType::singlevalued) {
          for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
            const auto& ap = areafrac[dir]->const_array(mfi);
            amrex::ParallelFor(
              eboxes[dir], NVAR,
              [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
                if (ap(i, j, k) > 0.0) {
                  flx[dir](i, j, k, n) /= ap(i, j, k);
                }
              });
          }
        }

        update_flux_registers(
          reflux_factor * dt, mfi, typ,
          {AMREX_D_DECL(&flux_ec[0], &flux_ec[1], &flux_ec[2])}, dm_as_fine);
      }

      copy_array4(vbox, NVAR, Dterm, MOLSrc);

      if (do_mol_load_balance && (cost != nullptr)) {
        amrex::Gpu::streamSynchronize();
        wt = (amrex::ParallelDescriptor::second() - wt) / vbox.d_numPts();
        (*cost)[mfi].plus<amrex::RunOn::Device>(wt, vbox);
      }
    }
  }
}
