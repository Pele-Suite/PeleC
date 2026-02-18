#include "prob.H"

void
pc_prob_close()
{
}

void
parse_params(ProbParmDevice* prob_parm_device)
{
  // Parse params
  amrex::ParmParse pp("prob");
  pp.query("d_inlet", prob_parm_device->d_inlet);
  pp.query("l_inlet", prob_parm_device->l_inlet);
  pp.query("M_inlet", prob_parm_device->M_inlet);
  pp.query("T_inlet", prob_parm_device->T_inlet);
  pp.query("p_inlet", prob_parm_device->p_inlet);
  pp.query("l_nozzle", prob_parm_device->l_nozzle);
  pp.query("d_exit", prob_parm_device->d_exit);
  pp.query("p_exit", prob_parm_device->p_exit);
  pp.query("p0", prob_parm_device->p0);
  pp.query("T0", prob_parm_device->T0);
  pp.query("inlet_type", prob_parm_device->inlet_type);
  pp.query("outlet_type", prob_parm_device->outlet_type);
  pp.query("do_sponge_zones", prob_parm_device->do_sponge_zones);
  pp.query("turb_inflow_type", prob_parm_device->turb_inflow_type);
}

extern "C" {
void
amrex_probinit(
  const int* /*init*/,
  const int* /*name*/,
  const int* /*namelen*/,
  const amrex::Real* /*problo*/,
  const amrex::Real* /*probhi*/)
{
  parse_params(PeleC::h_prob_parm_device);
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

void
EBConvergingNozzle::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  parse_params(PeleC::h_prob_parm_device);
  ProbParmDevice const* pp = PeleC::h_prob_parm_device;

  amrex::EB2::CylinderIF main(
    0.5 * pp->d_inlet, 0,
    {AMREX_D_DECL(static_cast<amrex::Real>(0.5 * pp->l_inlet), 0, 0)}, true);

  amrex::Real slope_nozzle =
    (0.5 * pp->d_inlet - 0.5 * pp->d_exit) / pp->l_nozzle;
  amrex::Real norm = -1.0 / slope_nozzle;
  amrex::Real nmag = std::sqrt(1 + 1 / (norm * norm));
  amrex::EB2::PlaneIF nozzle_plane(
    {AMREX_D_DECL(0, 0, 0)},
    {AMREX_D_DECL(
      static_cast<amrex::Real>(1.0 / nmag), slope_nozzle / nmag, 0.0)},
    true);
  auto nozzle = amrex::EB2::translate(
    amrex::EB2::rotate(amrex::EB2::lathe(nozzle_plane), 90 * M_PI / 180, 1),
    {AMREX_D_DECL(
      pp->l_inlet + static_cast<amrex::Real>(0.5 * pp->d_inlet / slope_nozzle),
      0, 0)});

  amrex::EB2::CylinderIF exit(
    0.5 * pp->d_exit, 0, {AMREX_D_DECL(pp->l_inlet + pp->l_nozzle, 0, 0)},
    true);
  auto nozzle_exit = amrex::EB2::makeIntersection(nozzle, exit);

  auto polys = amrex::EB2::makeUnion(main, nozzle_exit);
  auto gshop = amrex::EB2::makeShop(polys);
  amrex::EB2::Build(
    gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
}
