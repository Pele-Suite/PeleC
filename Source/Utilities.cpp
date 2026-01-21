#include "PeleCUtilities.H"

// -----------------------------------------------------------
// Read a binary file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_binary(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  const size_t ncol,
  amrex::Vector<double>& data) // needs to be double
{
  std::ifstream infile(iname, std::ios::in | std::ios::binary);
  if (!infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }

  for (size_t i = 0; i < nx * ny * nz * ncol; i++) {
    infile.read(reinterpret_cast<char*>(&data[i]), sizeof(data[i]));
  }
  infile.close();
}

// -----------------------------------------------------------
// Read a csv file
// INPUTS/OUTPUTS:
// iname => filename
// nx    => input resolution
// ny    => input resolution
// nz    => input resolution
// data  <= output data
// -----------------------------------------------------------
void
read_csv(
  const std::string& iname,
  const size_t nx,
  const size_t ny,
  const size_t nz,
  amrex::Vector<double>& data) // needs to be double
{
  std::ifstream infile(iname, std::ios::in);
  const std::string memfile = read_file(infile);
  if (!infile.is_open()) {
    amrex::Abort("Unable to open input file " + iname);
  }
  infile.close();
  std::istringstream iss(memfile);

  // Read the file
  size_t nlines = 0;
  std::string firstline;
  std::string line;
  std::getline(iss, firstline); // skip header
  while (getline(iss, line)) {
    ++nlines;
  }

  // Quick sanity check
  if (nlines != nx * ny * nz) {
    amrex::Abort(
      "Number of lines in the input file (= " + std::to_string(nlines) +
      ") does not match the input resolution (=" + std::to_string(nx) + ")");
  }

  // Read the data from the file
  iss.clear();
  iss.seekg(0, std::ios::beg);
  std::getline(iss, firstline); // skip header
  int cnt = 0;
  while (std::getline(iss, line)) {
    std::istringstream linestream(line);
    std::string value;
    while (getline(linestream, value, ',')) {
      std::istringstream sinput(value);
      sinput >> data[cnt];
      cnt++;
    }
  }
}

std::string
convertIntGG(int number)
{
  std::stringstream ss; // create a stringstream
  ss << number;         // add number to the stream
  return ss.str();      // return a string with the contents of the stream
}

void
clean_massfrac(
  const amrex::Box& bx,
  const amrex::Real threshold,
  amrex::Array4<const int> const& mask,
  amrex::Array4<amrex::Real> const& state)
{
  auto const& rho = amrex::Array4<amrex::Real>(state, URHO, 1);
  auto const& rhoU = amrex::Array4<amrex::Real>(state, UMX, AMREX_SPACEDIM);
  auto const& rhoY = amrex::Array4<amrex::Real>(state, UFS, NUM_SPECIES);
  auto const& rhoe = amrex::Array4<amrex::Real>(state, UEINT, 1);
  auto const& rhoE = amrex::Array4<amrex::Real>(state, UEDEN, 1);

  amrex::ParallelFor(
    bx, [=] AMREX_GPU_DEVICE(
          int i, int j, AMREX_D_PICK(int /*k*/, int /*k*/, int k)) noexcept {
      const amrex::IntVect iv{AMREX_D_DECL(i, j, k)};
      if (mask(iv) != 0) {
        const amrex::Real rhoOld = rho(iv);
        const amrex::Real rhoOld_inv = 1.0 / rhoOld;

        // Check for OOB mass fraction
        bool clean = false;
        for (int n = 0; n < NUM_SPECIES; n++) {
          const auto mf = rhoY(iv, n) * rhoOld_inv;
          if ((mf < -threshold) || ((1.0 + threshold) < mf)) {
            clean = true;
          }
        }

        if (clean) {
          // Clip species rhoYs and get new rho
          amrex::Real rhoNew = 0.0;
          for (int n = 0; n < NUM_SPECIES; n++) {
            rhoY(iv, n) = amrex::min<amrex::Real>(
              rhoOld, amrex::max<amrex::Real>(0.0, rhoY(iv, n)));
            rhoNew += rhoY(iv, n);
          }
          rho(iv) = rhoNew;

          // Keep kinetic energy, recompute, rhoe, rhoE, and rhoU
          const amrex::Real kinNRG =
            0.5 * rhoOld_inv * rhoOld_inv *
            (AMREX_D_TERM(
              (rhoU(iv, 0) * rhoU(iv, 0)), +(rhoU(iv, 1) * rhoU(iv, 1)),
              +(rhoU(iv, 2) * rhoU(iv, 2))));
          const amrex::Real eOld = (rhoE(iv) * rhoOld_inv) - kinNRG;
          rhoe(iv) = rhoNew * eOld;
          rhoE(iv) = rhoNew * eOld + rhoNew * kinNRG;
          for (int n = 0; n < AMREX_SPACEDIM; n++) {
            rhoU(iv, n) *= rhoNew * rhoOld_inv;
          }
        }
      }
    });
}
