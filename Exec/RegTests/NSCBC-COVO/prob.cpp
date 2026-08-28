#include "prob.H"

void
pc_prob_close()
{
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* problo,
  const amrex::Real* probhi)
{
  auto* pp_d = PeleC::h_prob_parm_device;
  {
    amrex::ParmParse pp("prob");
    pp.query("probtype", pp_d->probtype);
    pp.query("nscbc_inflow", pp_d->nscbc_inflow);
    pp.query("p_amb", pp_d->p_amb);
    pp.query("T_amb", pp_d->T_amb);
    pp.query("mach", pp_d->mach);
    pp.query("vortex_strength", pp_d->vortex_strength);
    pp.query("vortex_radius", pp_d->vortex_radius);
    pp.query("pulse_amp", pp_d->pulse_amp);
    pp.query("pulse_width", pp_d->pulse_width);
    amrex::Vector<amrex::Real> c(AMREX_SPACEDIM);
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      c[d] = pp_d->centre[d];
    }
    pp.queryarr("centre", c, 0, AMREX_SPACEDIM);
    for (int d = 0; d < AMREX_SPACEDIM; d++) {
      pp_d->centre[d] = c[d];
    }
  }

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;
  amrex::Real rho = 0.0, eint = 0.0, cs = 0.0, gam = 0.0;
  eos.PYT2RE(pp_d->p_amb, massfrac, pp_d->T_amb, rho, eint);
  eos.RTY2Cs(rho, pp_d->T_amb, massfrac, cs);
  eos.RTY2G(rho, pp_d->T_amb, massfrac, gam);
  pp_d->rho_amb = rho;
  pp_d->cs_amb = cs;
  pp_d->gamma = gam;
  pp_d->u_mean = (pp_d->probtype == 0) ? pp_d->mach * cs : 0.0;

  for (int d = 0; d < AMREX_SPACEDIM; d++) {
    pp_d->Lbox[d] = probhi[d] - problo[d];
    pp_d->xc[d] = problo[d] + pp_d->centre[d] * pp_d->Lbox[d];
  }

  amrex::Print() << "  NSCBC-COVO: probtype = " << pp_d->probtype
                 << (pp_d->probtype == 0 ? " (convected vortex)"
                                         : " (circular pulse)")
                 << "\n"
                 << "              rho_amb = " << rho << " g/cc,  c = " << cs
                 << " cm/s,  gamma = " << gam << "\n"
                 << "              u_mean  = " << pp_d->u_mean
                 << " cm/s (M = " << (pp_d->probtype == 0 ? pp_d->mach : 0.0)
                 << ")\n";
}
}

void
PeleC::problem_post_timestep()
{
}

void
PeleC::problem_post_init()
{
}

void
PeleC::problem_post_restart()
{
}
