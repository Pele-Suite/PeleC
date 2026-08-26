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
    pp.query("p_amb", pp_d->p_amb);
    pp.query("T_amb", pp_d->T_amb);
    pp.query("amp", pp_d->amp);
    pp.query("width", pp_d->width);
    pp.query("pulse_type", pp_d->pulse_type);
    pp.query("u0", pp_d->u0);
    pp.query("force_amp", pp_d->force_amp);
    pp.query("force_freq", pp_d->force_freq);
    pp.query("force_feedforward", pp_d->force_feedforward);

    // x0 accepts either a single value -- the historical spelling, which sets
    // the streamwise position of the planar pulse -- or one per direction.
    amrex::Vector<amrex::Real> x0;
    if (pp.queryarr("x0", x0) != 0) {
      if (x0.size() == 1) {
        pp_d->x0[0] = x0[0];
      } else if (x0.size() == AMREX_SPACEDIM) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
          pp_d->x0[d] = x0[d];
        }
      } else {
        amrex::Abort("prob.x0 must have 1 or AMREX_SPACEDIM entries");
      }
    }
  }

  const amrex::Real L[AMREX_SPACEDIM] = {AMREX_D_DECL(
    probhi[0] - problo[0], probhi[1] - problo[1], probhi[2] - problo[2])};

  // The pulse width is a single physical length, so that a radial pulse is a
  // sphere and not an ellipsoid on an anisotropic mesh.  It is measured
  // against the streamwise extent for the planar pulse and against the
  // SMALLEST extent for the radial one, so that the pulse fits in the box
  // whatever the aspect ratio.
  amrex::Real Lref = L[0];
  if (pp_d->pulse_type != 0) {
    AMREX_D_TERM(
      Lref = L[0];, Lref = amrex::min(Lref, L[1]);
      , Lref = amrex::min(Lref, L[2]);)
  }
  pp_d->w = pp_d->width * Lref;
  AMREX_D_EXPR(
    pp_d->xc[0] = problo[0] + pp_d->x0[0] * L[0],
    pp_d->xc[1] = problo[1] + pp_d->x0[1] * L[1],
    pp_d->xc[2] = problo[2] + pp_d->x0[2] * L[2]);

  auto eos = pele::physics::PhysicsType::eos();
  amrex::Real massfrac[NUM_SPECIES] = {0.0};
  massfrac[0] = 1.0;
  amrex::Real rho = 0.0, eint = 0.0, cs = 0.0;
  eos.PYT2RE(pp_d->p_amb, massfrac, pp_d->T_amb, rho, eint);
  eos.RTY2Cs(rho, pp_d->T_amb, massfrac, cs);
  pp_d->rho_amb = rho;
  pp_d->cs_amb = cs;

  amrex::Print() << "\n  NSCBC-Acoustic (" << AMREX_SPACEDIM
                 << "D): rho_amb = " << rho << " g/cc, c_amb = " << cs
                 << " cm/s\n";
  if (pp_d->pulse_type == 2) {
    // Quarter-wave frequency of the duct, Doppler-corrected -- the resonance
    // t3/t5 straddle.  Printed so the launch script never re-derives it.
    const amrex::Real f0 =
      (cs * cs - pp_d->u0 * pp_d->u0) / (4.0 * L[0] * cs);
    amrex::Print() << "     DUCT: u0 = " << pp_d->u0 << " cm/s, forcing "
                   << pp_d->force_amp << " cm/s at " << pp_d->force_freq
                   << " Hz\n"
                   << "     quarter-wave f0 = " << f0
                   << " Hz  (t_a = " << 2.0 * L[0] / cs << " s)\n\n";
  } else if (pp_d->pulse_type == 0) {
    amrex::Print() << "     PLANAR pulse at x = " << pp_d->xc[0] << " cm, "
                   << "width " << pp_d->w << " cm, running toward x-hi\n"
                   << "     transit to the outflow = "
                   << (probhi[0] - pp_d->xc[0]) / cs << " s\n\n";
  } else {
    // Distances the front has to travel to reach a face, an edge and a corner.
    // Quoted because they are what the isotropy metric is measured against:
    // the front is a sphere only if all three arrive on the same clock.
    amrex::Real dmin = L[0];
    AMREX_D_TERM(
      dmin = 0.5 * L[0];, dmin = amrex::min(dmin, 0.5 * L[1]);
      , dmin = amrex::min(dmin, 0.5 * L[2]);)
    amrex::Print() << "     RADIAL pulse at ("
                   << AMREX_D_TERM(
                        pp_d->xc[0], << ", " << pp_d->xc[1],
                                     << ", " << pp_d->xc[2])
                   << ") cm, width " << pp_d->w << " cm, at rest -- it splits\n"
                   << "     every face is a characteristic outflow; the front "
                      "meets the faces at normal\n"
                   << "     incidence, the edges at 45 deg and the corners "
                      "along the body diagonal\n"
                   << "     nearest face at " << dmin
                   << " cm, i.e. t = " << dmin / cs << " s\n\n";
  }
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
