#include "prob.H"

#include <AMReX_EB2_IF_Box.H>
#include <AMReX_EB2_IF_Union.H>

// The chamber as interior EB: three slabs (bottom, top, back) forming a
// rectangular cavity open at its x-hi end, plus optionally the Sydney baffle
// -- two plates across the cavity leaving a central gap.  All dimensions are
// prob.* inputs so the two runnable variants (box, box+baffle) share one
// geometry.  Everything is strictly inside the domain: an EB body cutting a
// domain face is the documented NaN limitation.
void
NSCBCChamberBox::build(
  const amrex::Geometry& geom, const int max_coarsening_level)
{
  amrex::ParmParse pp("prob");
  amrex::Real x0 = 0.2;      // chamber interior back wall
  amrex::Real len = 1.2;     // chamber interior length (vent at x0 + len)
  amrex::Real h = 0.3;       // chamber interior height, centred in y
  amrex::Real wall_t = 0.05; // slab thickness
  int baffle = 0;
  amrex::Real baffle_x = 0.8;   // plate upstream face, absolute
  amrex::Real baffle_w = 0.05;  // plate thickness
  amrex::Real baffle_gap = 0.1; // central opening height
  pp.query("chamber_x0", x0);
  pp.query("chamber_len", len);
  pp.query("chamber_h", h);
  pp.query("wall_t", wall_t);
  pp.query("baffle", baffle);
  pp.query("baffle_x", baffle_x);
  pp.query("baffle_w", baffle_w);
  pp.query("baffle_gap", baffle_gap);

  // Chamber centreline: defaults to the domain centre, overridable so the
  // EB surfaces can be placed OFF the grid lines.  Keep them off: a surface
  // sitting exactly on a cell face makes degenerate cut cells, and the first
  // grid-aligned build of this geometry NaN'd at the re-entrant interior
  // corner cell under both hydro schemes.
  amrex::Real yc = 0.5 * (geom.ProbLo(1) + geom.ProbHi(1));
  pp.query("chamber_yc", yc);
  const amrex::Real ylo = yc - 0.5 * h, yhi = yc + 0.5 * h;

  // bottom and top walls (extended back over the closed end's thickness),
  // and the closed end itself
  amrex::EB2::BoxIF bot(
    {AMREX_D_DECL(x0 - wall_t, ylo - wall_t, 0.0)},
    {AMREX_D_DECL(x0 + len, ylo, 1.0)}, false);
  amrex::EB2::BoxIF top(
    {AMREX_D_DECL(x0 - wall_t, yhi, 0.0)},
    {AMREX_D_DECL(x0 + len, yhi + wall_t, 1.0)}, false);
  amrex::EB2::BoxIF back(
    {AMREX_D_DECL(x0 - wall_t, ylo, 0.0)},
    {AMREX_D_DECL(x0, yhi, 1.0)}, false);
  auto walls = amrex::EB2::makeUnion(bot, top, back);

  if (baffle != 0) {
    const amrex::Real gl = yc - 0.5 * baffle_gap, gh = yc + 0.5 * baffle_gap;
    amrex::EB2::BoxIF plate_lo(
      {AMREX_D_DECL(baffle_x, ylo, 0.0)},
      {AMREX_D_DECL(baffle_x + baffle_w, gl, 1.0)}, false);
    amrex::EB2::BoxIF plate_hi(
      {AMREX_D_DECL(baffle_x, gh, 0.0)},
      {AMREX_D_DECL(baffle_x + baffle_w, yhi, 1.0)}, false);
    auto gshop = amrex::EB2::makeShop(
      amrex::EB2::makeUnion(walls, plate_lo, plate_hi));
    amrex::EB2::Build(
      gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
  } else {
    auto gshop = amrex::EB2::makeShop(walls);
    amrex::EB2::Build(
      gshop, geom, max_coarsening_level, max_coarsening_level, 4, false);
  }

  amrex::Print() << "\n  NSCBC-Chamber BOX: interior [" << x0 << ", "
                 << x0 + len << "] x [" << ylo << ", " << yhi
                 << "], vent at x = " << x0 + len
                 << (baffle != 0 ? ", baffle at x = " : "")
                 << (baffle != 0 ? std::to_string(baffle_x) : std::string())
                 << "\n";
}

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
  pp.query("kernel_r", PeleC::h_prob_parm_device->kernel_r);
  pp.query("kernel_x", PeleC::h_prob_parm_device->kernel_x);
  pp.query("kernel_y", PeleC::h_prob_parm_device->kernel_y);

  // Box variants: clip the ignition structure to the chamber interior (the
  // closed-end slab is thinner than the kernel radius, and the preheat tail
  // is thicker than the side walls -- see ProbParmDevice::kernel_clip).
  {
    std::string geom_type;
    amrex::ParmParse ppeb("eb2");
    ppeb.query("geom_type", geom_type);
    if (geom_type == "nscbc-chamber-box") {
      amrex::Real x0 = 0.2, len = 1.2, h = 0.3;
      amrex::Real yc = 0.5 * (problo[1] + probhi[1]);
      pp.query("chamber_x0", x0);
      pp.query("chamber_len", len);
      pp.query("chamber_h", h);
      pp.query("chamber_yc", yc);
      auto* pd = PeleC::h_prob_parm_device;
      pd->kernel_clip = 1;
      pd->clip_lo[0] = x0;
      pd->clip_hi[0] = x0 + len;
      pd->clip_lo[1] = yc - 0.5 * h;
      pd->clip_hi[1] = yc + 0.5 * h;
      amrex::Print() << "     ignition clipped to the chamber interior ["
                     << pd->clip_lo[0] << ", " << pd->clip_hi[0] << "] x ["
                     << pd->clip_lo[1] << ", " << pd->clip_hi[1] << "]\n";
    }
  }
  pp.query("pmf_flame_loc", PeleC::h_prob_parm_device->pmf_flame_loc);
  pp.query("pmf_datafile", pmf_datafile);

  read_pmf(pmf_datafile);

  auto* pp_d = PeleC::h_prob_parm_device;
  for (int i = 0; i < AMREX_SPACEDIM; i++) {
    pp_d->L[i] = probhi[i] - problo[i];
  }
  // Column 1 of the profile is velocity; point 0 is the fresh inlet, and its
  // value is the laminar flame speed by construction.
  pp_d->s_L = PeleC::prob_parm_host->h_pmf_Y[1 * pp_d->pmf_N + 0];

  amrex::Print() << "\n  NSCBC-Chamber (mini-SydGex, Lesson 9): S_L = "
                 << pp_d->s_L << " cm/s\n"
                 << "     chamber " << pp_d->L[0] << " x " << pp_d->L[1]
                 << " cm, kernel r = " << pp_d->kernel_r
                 << " cm on the closed end at y = "
                 << problo[1] + pp_d->kernel_y * pp_d->L[1] << " cm\n"
                 << "     quiescent charge at p = " << pp_d->pamb
                 << "; the open end is at x = " << probhi[0] << "\n\n";
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
