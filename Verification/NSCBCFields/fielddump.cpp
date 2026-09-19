// ============================================================================
//  fielddump -- flatten one variable of a single-level AMReX plotfile onto a
//  regular array, for the NSCBC multi-dimensional metrics.  Builds at whatever
//  AMREX_SPACEDIM the AMReX it links against was configured with.
//
//  (Named fielddump rather than pltdump: the repository .gitignore carries a
//  `plt*` rule for plotfiles, which silently swallows any source file whose
//  name begins with "plt".)
//
//  The 2-D boundary tests need the whole field, not a line-out, because the
//  quantity of interest is the SHAPE of the wavefront.  fextract gives 1-D
//  slices only, so this exists.
//
//    fielddump <plotfile> <varname> <outfile>
//
//  Writes an ASCII header followed by the values, deliberately trivial to
//  parse:
//
//    2-D:  # nx ny xlo ylo dx dy time            then nx*ny  values, j slowest
//    3-D:  # nx ny nz xlo ylo zlo dx dy dz time  then nx*ny*nz values, k
//    slowest
//
//  The 2-D header is byte-for-byte what it always was, so every existing
//  reader keeps working; a reader that wants to handle both can branch on the
//  number of tokens on the header line.
// ============================================================================

#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_MultiFab.H>

#include <cstdio>
#include <string>
#include <vector>

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv, false);
  int rc = 0;
  {
    if (argc < 4) {
      amrex::Print() << "usage: fielddump <plotfile> <varname> <outfile>\n";
      amrex::Finalize();
      return 1;
    }
    const std::string pf(argv[1]);
    const std::string var(argv[2]);
    const std::string out(argv[3]);

    amrex::PlotFileData plotfile(pf);
    const int lev = 0; // these tests are single level by construction
    const amrex::Box dom = plotfile.probDomain(lev);
    const auto plo = plotfile.probLo();
    const auto dx = plotfile.cellSize(lev);

    bool found = false;
    for (const auto& n : plotfile.varNames()) {
      found = found || (n == var);
    }
    if (!found) {
      amrex::Print() << "fielddump: variable '" << var << "' not in " << pf
                     << "\n  available:";
      for (const auto& n : plotfile.varNames()) {
        amrex::Print() << " " << n;
      }
      amrex::Print() << "\n";
      amrex::Finalize();
      return 1;
    }

    const amrex::MultiFab mf = plotfile.get(lev, var);

    const int nx = dom.length(0);
    const int ny = (AMREX_SPACEDIM > 1) ? dom.length(1) : 1;
    const int nz = (AMREX_SPACEDIM > 2) ? dom.length(2) : 1;
    const size_t ntot = static_cast<size_t>(nx) * ny * nz;
    std::vector<double> a(ntot, std::numeric_limits<double>::quiet_NaN());

    // Index of (i,j,k) in the flattened array, with the last dimension
    // slowest.  Written once so the write loop below cannot disagree with it.
    const auto idx = [=](const int i, const int j, const int k) {
      return (static_cast<size_t>(k) * ny + j) * nx + i;
    };

    for (amrex::MFIter mfi(mf); mfi.isValid(); ++mfi) {
      const amrex::Box& bx = mfi.validbox();
      const auto& arr = mf.const_array(mfi);
      const auto lo = amrex::lbound(bx);
      const auto hi = amrex::ubound(bx);
      for (int k = lo.z; k <= hi.z; ++k) {
        for (int j = lo.y; j <= hi.y; ++j) {
          for (int i = lo.x; i <= hi.x; ++i) {
            a[idx(
              i - dom.smallEnd(0),
              (AMREX_SPACEDIM > 1) ? j - dom.smallEnd(1) : 0,
              (AMREX_SPACEDIM > 2) ? k - dom.smallEnd(2) : 0)] = arr(i, j, k);
          }
        }
      }
    }

    FILE* f = std::fopen(out.c_str(), "w");
#if AMREX_SPACEDIM == 3
    std::fprintf(
      f,
      "# nx ny nz xlo ylo zlo dx dy dz time\n"
      "%d %d %d %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
      nx, ny, nz, plo[0], plo[1], plo[2], dx[0], dx[1], dx[2], plotfile.time());
#else
    std::fprintf(
      f, "# nx ny xlo ylo dx dy time\n%d %d %.17g %.17g %.17g %.17g %.17g\n",
      nx, ny, plo[0], plo[1], dx[0], dx[1], plotfile.time());
#endif
    for (int k = 0; k < nz; ++k) {
      for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
          std::fprintf(f, "%.17e ", a[idx(i, j, k)]);
        }
        std::fprintf(f, "\n");
      }
    }
    std::fclose(f);
    amrex::Print() << "fielddump: wrote " << nx << " x " << ny << " x " << nz
                   << " '" << var << "' at t = " << plotfile.time() << " to "
                   << out << "\n";
  }
  amrex::Finalize();
  return rc;
}
