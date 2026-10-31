// Problem-specific derives
#ifndef PROBLEMDERIVE_H
#define PROBLEMDERIVE_H

#include <AMReX_Derive.H>
#include <AMReX_FArrayBox.H>

#include "PeleC.H"
#include "IndexDefines.H"
#include "NSCBC.H"

struct DefaultProblemSpecificFunctions
{
  static void add_problem_derive(
    amrex::DeriveList& /*derive_lst*/, amrex::DescriptorList& /*desc_lst*/)
  {
  }

  static void pc_varname(
    const amrex::Box& /*bx*/,
    amrex::FArrayBox& /*derfab*/,
    int /*dcomp*/,
    int /*ncomp*/,
    const amrex::FArrayBox& /*datfab*/,
    const amrex::Geometry& /*geomdata*/,
    amrex::Real /*time*/,
    const int* /*bcrec*/,
    int /*level*/)
  {
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_problem_tags(
    const int /*i*/,
    const int /*j*/,
    const int /*k*/,
    amrex::Array4<amrex::EBCellFlag const> const& /*flags*/,
    amrex::Array4<char> const& /*tag*/,
    amrex::Array4<amrex::Real const> const& /*field*/,
    char /*tagval*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*dx*/,
    const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> /*prob_lo*/,
    const amrex::Real /*time*/,
    const int /*level*/,
    ProbParmDevice const& /*d_prob_parm_device*/) noexcept
  {
  }

  static void set_aux_names(amrex::Vector<std::string>& a_aux_names)
  {
    a_aux_names.resize(NUM_AUX);
#if NUM_AUX > 0
    for (int n = 0; n < NUM_AUX; ++n) {
      a_aux_names[n] = "aux_" + std::to_string(n);
    }
#endif
  }

  static void set_adv_names(amrex::Vector<std::string>& a_adv_names)
  {
    a_adv_names.resize(NUM_ADV);
#if NUM_ADV > 0
    for (int n = 0; n < NUM_ADV; ++n) {
      a_adv_names[n] = "adv_" + std::to_string(n);
    }
#endif
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static void set_isothermal_wall_temperature(
    const int i,
    const int j,
    const int k,
    const int /*dir*/,
    const int /*normal*/,
    const amrex::Real Twall_in,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& /*prob_parm*/,
    const amrex::Array4<const amrex::Real>& /*q*/,
    const amrex::Array4<amrex::Real>& Twall)
  {

    /*
      if desired: compute x,y,z and make wall temp an arbitrary function of
      space

      const amrex::Real* prob_lo = geom.ProbLo();
      const amrex::Real* dx = geom.CellSize();
      const amrex::Real x[AMREX_SPACEDIM] = {AMREX_D_DECL(
        prob_lo[0] + static_cast<amrex::Real>(iv[0] + 0.5) * dx[0],
        prob_lo[1] + static_cast<amrex::Real>(iv[1] + 0.5) * dx[1],
        prob_lo[2] + static_cast<amrex::Real>(iv[2] + 0.5) * dx[2])};

      Twall(i,j,k) = function of x, dir, normal, prob_parm
                     or set it using the fluid values (q) for adabatic
     */
    Twall(i, j, k) = Twall_in;
  }

  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static bool problem_eb_state(
    amrex::GeometryData const& /*geomdata*/,
    const amrex::Real /*vf*/,
    const amrex::IntVect /*iv*/,
    AMREX_D_DECL(
      const amrex::Real /*nx*/,
      const amrex::Real /*ny*/,
      const amrex::Real /*nz*/),
    const amrex::Real* /*ql[5 + NUM_SPECIES]*/,
    const amrex::Real* /*spl[NUM_SPECIES]*/,
    const amrex::Real /*rhoe_l*/,
    const amrex::Real /*gamc_l*/,
    amrex::Real* /*qr[5 + NUM_SPECIES]*/,
    amrex::Real* /*spr[NUM_SPECIES]*/,
    amrex::Real& /*rhoe_r*/,
    amrex::Real& /*gamc_r*/,
    ProbParmDevice const* /*prob_parm*/)
  {
    return false;
  }

  static void problem_modify_ext_sources(
    amrex::Real /*time*/,
    amrex::Real /*dt*/,
    const amrex::MultiFab& /*state_old*/,
    const amrex::MultiFab& /*state_new*/,
    amrex::MultiFab& /*ext_src*/,
    int /*ng*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& /*prob_parm*/)
  {
    /* Notes: ext_src contains sources from velocity forcing coming in
       This function should add to rather than overwrite ext_src.
    */
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  static void problem_modify_transport_coeffs(
    const bool /*get_xi*/,
    const bool /*get_mu*/,
    const bool /*get_lam*/,
    const bool /*get_Ddiag*/,
    const bool /*get_chi*/,
    const amrex::Real /*Tloc*/,
    const amrex::Real /*rholoc*/,
    amrex::Real* /*Yloc*/,
    amrex::Real* /*Ddiag*/,
    amrex::Real* /*chi_mix*/,
    amrex::Real& /*mu*/,
    amrex::Real& /*xi*/,
    amrex::Real& /*lam*/,
    pele::physics::transport::TransParm<
      pele::physics::PhysicsType::eos_type,
      pele::physics::PhysicsType::transport_type> const* /*tparm*/,
    ProbParmDevice const& /*prob_parm*/,
    const amrex::RealVect& /*x*/)
  {
  }

  // Characteristic (NSCBC) boundary treatment.
  //
  // Return the target state for one boundary POINT.  The default returns
  // Type::off, which leaves that point to the ordinary bcnormal() path, so
  // pelec.bc_nscbc = 1 is a no-op for a problem that does not override this.
  // Because the decision is per point, a single face may mix an inflow, an
  // outflow and a wall.
  //
  // Only the fields relevant to the returned type are read:
  //     Type::outflow -> p
  //     Type::inflow  -> u[], T, Y[]
  // There is deliberately no relaxation coefficient for composition; see
  // Source/NSCBC.H.
  //
  // s_int is the state in the first interior cell along the boundary normal,
  // for problems whose target depends on the interior state (a back-pressure
  // that tracks a measured plenum, say).  Most problems will not need it.
  AMREX_GPU_DEVICE
  AMREX_FORCE_INLINE
  static pc_nscbc::Target bcnormal_nscbc(
    const amrex::Real* /*x[AMREX_SPACEDIM]*/,
    const amrex::Real* /*s_int[NVAR]*/,
    const int /*idir*/,
    const int /*sgn*/,
    const amrex::Real /*time*/,
    amrex::GeometryData const& /*geomdata*/,
    ProbParmDevice const& /*prob_parm*/)
  {
    return pc_nscbc::Target{};
  }

  /* TODO
  static void something_with_ADV_AUX_transport()
  {
  }
  */
};

#endif
