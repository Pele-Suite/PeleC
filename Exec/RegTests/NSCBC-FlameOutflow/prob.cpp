#include "prob.H"

std::string
read_pmf_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

bool
checkQuotes(const std::string& str)
{
  int count = 0;
  for (char c : str) {
    if (c == '"') {
      count++;
    }
  }
  return (count % 2) == 0;
}

void
read_pmf(const std::string& myfile)
{
  std::string firstline;
  std::string secondline;
  std::string remaininglines;
  unsigned int pos1;
  unsigned int pos2;
  int variable_count;
  int line_count;

  std::ifstream infile(myfile);
  if (!infile.good()) {
    // Without this, a missing file hands the quote parser an empty first
    // line, and `pos1 < firstline.length() - 1` underflows to SIZE_MAX: the
    // run hangs forever inside amrex_probinit with no message.  (The
    // original in Exec/RegTests/PMF has the same trap.)
    amrex::Abort("read_pmf: cannot open prob.pmf_datafile = " + myfile);
  }
  const std::string memfile = read_pmf_file(infile);
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  if (!checkQuotes(firstline)) {
    amrex::Abort("PMF file variable quotes unbalanced");
  }
  std::getline(iss, secondline);
  pos1 = 0;
  pos2 = 0;
  variable_count = 0;
  while ((pos1 < firstline.length() - 1) && (pos2 < firstline.length() - 1)) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    variable_count++;
    pos1 = pos2 + 1;
  }

  pos1 = 0;
  for (int i = 0; i < variable_count; i++) {
    pos1 = firstline.find('"', pos1);
    pos2 = firstline.find('"', pos1 + 1);
    pos1 = pos2 + 1;
  }

  amrex::Print() << variable_count << " variables found in PMF file"
                 << std::endl;
  // for (int i = 0; i < variable_count; i++)
  //  amrex::Print() << "Variable found: " << pmf_names[i] <<
  //  std::endl;

  line_count = 0;
  while (std::getline(iss, remaininglines)) {
    line_count++;
  }
  amrex::Print() << line_count << " data lines found in PMF file" << std::endl;

  PeleC::h_prob_parm_device->pmf_N = line_count;
  PeleC::h_prob_parm_device->pmf_M = variable_count - 1;
  PeleC::prob_parm_host->h_pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->pmf_X.resize(PeleC::h_prob_parm_device->pmf_N);
  PeleC::prob_parm_host->h_pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);
  PeleC::prob_parm_host->pmf_Y.resize(
    static_cast<long>(PeleC::h_prob_parm_device->pmf_N) *
    PeleC::h_prob_parm_device->pmf_M);

  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline);
  std::getline(iss, secondline);
  for (int i = 0; i < PeleC::h_prob_parm_device->pmf_N; i++) {
    std::getline(iss, remaininglines);
    std::istringstream sinput(remaininglines);
    sinput >> PeleC::prob_parm_host->h_pmf_X[i];
    for (int j = 0; j < PeleC::h_prob_parm_device->pmf_M; j++) {
      sinput >> PeleC::prob_parm_host
                  ->h_pmf_Y[j * PeleC::h_prob_parm_device->pmf_N + i];
    }
  }

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_X.begin(),
    PeleC::prob_parm_host->h_pmf_X.end(), PeleC::prob_parm_host->pmf_X.begin());
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, PeleC::prob_parm_host->h_pmf_Y.begin(),
    PeleC::prob_parm_host->h_pmf_Y.end(), PeleC::prob_parm_host->pmf_Y.begin());
  PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
  PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
}

void
init_bc()
{
  amrex::Real vt;
  amrex::Real ek;
  amrex::Real T;
  amrex::Real rho;
  amrex::Real e;
  amrex::Real molefrac[NUM_SPECIES];
  amrex::Real massfrac[NUM_SPECIES];
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {{0.0}};

  if (PeleC::h_prob_parm_device->phi_in < 0) {
    const amrex::Real yl = 0.0;
    const amrex::Real yr = 0.0;
    // Use host pointers for host call to pmf()
    PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->h_pmf_X.data();
    PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->h_pmf_Y.data();
    pmf(yl, yr, pmf_vals, *PeleC::h_prob_parm_device);
    // Switch back to device pointers
    PeleC::h_prob_parm_device->d_pmf_X = PeleC::prob_parm_host->pmf_X.data();
    PeleC::h_prob_parm_device->d_pmf_Y = PeleC::prob_parm_host->pmf_Y.data();
    amrex::Real mysum = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      molefrac[n] = amrex::max<amrex::Real>(0.0, pmf_vals[3 + n]);
      mysum += molefrac[n];
    }
    molefrac[N2_ID] = 1.0 - (mysum - molefrac[N2_ID]);
    T = pmf_vals[0];
    PeleC::h_prob_parm_device->vn_in = pmf_vals[1];
  } else {
    const amrex::Real a = 0.5;
    for (amrex::Real& n : molefrac) {
      n = 0.0;
    }
    molefrac[O2_ID] =
      1.0 / (1.0 + PeleC::h_prob_parm_device->phi_in / a + 0.79 / 0.21);
    molefrac[H2_ID] = PeleC::h_prob_parm_device->phi_in * molefrac[O2_ID] / a;
    molefrac[N2_ID] = 1.0 - molefrac[H2_ID] - molefrac[O2_ID];
    T = PeleC::h_prob_parm_device->T_in;
  }
  const amrex::Real p = PeleC::h_prob_parm_device->pamb;

  auto eos = pele::physics::PhysicsType::eos();
  eos.X2Y(molefrac, massfrac);
  eos.PYT2RE(p, massfrac, T, rho, e);

  vt = PeleC::h_prob_parm_device->vn_in;
  ek = 0.5 * (vt * vt);

  PeleC::h_prob_parm_device->fuel_state[URHO] = rho;
  PeleC::h_prob_parm_device->fuel_state[UMX] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UMY] = rho * vt;
  PeleC::h_prob_parm_device->fuel_state[UMZ] = 0.0;
  PeleC::h_prob_parm_device->fuel_state[UEINT] = rho * e;
  PeleC::h_prob_parm_device->fuel_state[UEDEN] = rho * (e + ek);
  PeleC::h_prob_parm_device->fuel_state[UTEMP] = T;
  for (int n = 0; n < NUM_SPECIES; n++) {
    PeleC::h_prob_parm_device->fuel_state[UFS + n - 1] = rho * massfrac[n];
  }
}

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
  std::string pmf_datafile;

  amrex::ParmParse pp("prob");
  pp.query("pamb", PeleC::h_prob_parm_device->pamb);
  pp.query("phi_in", PeleC::h_prob_parm_device->phi_in);
  pp.query("T_in", PeleC::h_prob_parm_device->T_in);
  pp.query("vn_in", PeleC::h_prob_parm_device->vn_in);
  pp.query("pertmag", PeleC::h_prob_parm_device->pertmag);
  pp.query("nscbc_inflow", PeleC::h_prob_parm_device->nscbc_inflow);
  pp.query("u_ratio", PeleC::h_prob_parm_device->u_ratio);
  pp.query("sheet_x", PeleC::h_prob_parm_device->sheet_x);
  pp.query("wrinkle_amp", PeleC::h_prob_parm_device->wrinkle_amp);
  pp.query("wrinkle_k", PeleC::h_prob_parm_device->wrinkle_k);
  pp.query("pmf_flame_loc", PeleC::h_prob_parm_device->pmf_flame_loc);
  pp.query("pmf_datafile", pmf_datafile);
  amrex::Vector<amrex::Real> local_L(AMREX_SPACEDIM, -1);
  pp.queryarr("L", local_L, 0, AMREX_SPACEDIM);
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    PeleC::h_prob_parm_device->L[i] =
      (local_L[i] == -1.0) ? probhi[i] - problo[i] : local_L[i];
  }

  read_pmf(pmf_datafile);

  // The laminar flame speed is the inlet velocity of the 1-D profile, which is
  // by construction the speed at which that flame consumes fresh gas.
  auto* pp_d = PeleC::h_prob_parm_device;
  // Column 1 of the profile is velocity; point 0 is the fresh inlet.
  pp_d->s_L = PeleC::prob_parm_host->h_pmf_Y[1 * pp_d->pmf_N + 0];
  pp_d->u_in = pp_d->u_ratio * pp_d->s_L;
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    pp_d->L[i] = probhi[i] - problo[i];
  }
  pp_d->sheet_x0 = problo[0] + pp_d->sheet_x * pp_d->L[0];
  pp_d->ylo = problo[1];

  amrex::Print() << "\n  NSCBC-FlameOutflow: S_L = " << pp_d->s_L
                 << " cm/s,  U = " << pp_d->u_in << " cm/s (" << pp_d->u_ratio
                 << " S_L)\n";

  const amrex::Real kk = 2.0 * constants::PI() * pp_d->wrinkle_k / pp_d->L[1];
  const amrex::Real slope = pp_d->wrinkle_amp * kk;
  const amrex::Real xmin = pp_d->sheet_x0 - pp_d->wrinkle_amp;
  const amrex::Real xmax = pp_d->sheet_x0 + pp_d->wrinkle_amp;
  const amrex::Real xout = probhi[0];
  amrex::Print() << "     wrinkled sheet at x_f = " << pp_d->sheet_x0 << " + "
                 << pp_d->wrinkle_amp << " cos(2 pi " << pp_d->wrinkle_k
                 << " y / " << pp_d->L[1] << "),  x_f in [" << xmin << ", "
                 << xmax << "]\n"
                 << "     max tilt of the flame normal off the boundary "
                    "normal = "
                 << std::atan(slope) * 180.0 / constants::PI() << " deg\n";
  if ((xmin < xout) && (xmax > xout)) {
    // Fraction of the transverse extent for which the boundary is burnt,
    // i.e. x_f(y) < x_out:  cos(k y) < (x_out - x0)/A.
    const amrex::Real c0 = (xout - pp_d->sheet_x0) / pp_d->wrinkle_amp;
    const amrex::Real frac = std::acos(c0) / constants::PI();
    amrex::Print()
      << "     STRADDLES the outflow: " << 100.0 * frac
      << "% of the boundary is burnt, " << 100.0 * (1.0 - frac)
      << "% fresh, with the reaction zone crossing at 2 points\n";
  } else if (amrex::ParallelDescriptor::IOProcessor()) {
    // amrex::Warning prints on ALL ranks (Print::AllProcs to the error
    // stream), so gate it here -- probinit runs everywhere.
    amrex::Warning(
      "NSCBC-FlameOutflow: the sheet does not straddle the outflow at x = " +
      std::to_string(xout) +
      "; the reaction zone is not in the boundary cells and beta_s has "
      "nothing to correct");
  }
  amrex::Print() << "     drift speed U - S_L = " << pp_d->u_in - pp_d->s_L
                 << " cm/s (unanchored: the sheet leaves on its own)\n\n";

  init_bc();
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
