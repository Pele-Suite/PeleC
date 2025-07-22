#ifndef DIFFUSION_H
#define DIFFUSION_H

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

#include <AMReX_FArrayBox.H>

#include "PeleC.H"
#include "PeleCUtilities.H"
#include "GradUtil.H"
#include "Diffterm.H"
#include "MOL.H"
#include "EB.H"
#include "AMReX_EB_Redistribution.H"

// This header file contains functions and declarations for the diffusion
// portion of PeleC GPU. As per the convention of AMReX, inlined functions are
// defined here. Where as non-inline functions are declared here.

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
pc_diffextrap(
  const int i,
  const int j,
  const int k,
  const amrex::Array4<amrex::Real>& dif,
  const int /*ng*/,
  const int first,
  const int last,
  AMREX_D_DECL(const int lx, const int ly, const int lz),
  AMREX_D_DECL(const int hx, const int hy, const int hz),
  const amrex::Dim3 dlo,
  const amrex::Dim3 dhi)
{
  for (int n = first; n < last; n++) {
    //============================= Sides
    //============================================= low x
    if (dlo.x < lx) {
      if (i == lx) {
        dif(i - 1, j, k, n) = dif(i, j, k, n);
      }
    }
    // high x
    if (dhi.x - 1 > hx) {
      if (i == hx) {
        dif(i + 1, j, k, n) = dif(i, j, k, n);
      }
    }
#if AMREX_SPACEDIM > 1
    // low y
    if (dlo.y < ly) {
      if (j == ly) {
        dif(i, j - 1, k, n) = dif(i, j, k, n);
      }
    }
    // high y
    if (dhi.y - 1 > hy) {
      if (j == hy) {
        dif(i, j + 1, k, n) = dif(i, j, k, n);
      }
    }
#if AMREX_SPACEDIM == 3
    // low z
    if (dlo.z < lz) {
      if (k == lz) {
        dif(i, j, k - 1, n) = dif(i, j, k, n);
      }
    }
    // high z
    if (dhi.z - 1 > hz) {
      if (k == hz) {
        dif(i, j, k + 1, n) = dif(i, j, k, n);
      }
    }
#endif
    //============================= Edges
    //============================================== low x low y
    if (dlo.x < lx && dlo.y < ly) {
      if (i == lx && j == ly) {
        dif(i - 1, j - 1, k, n) = dif(i, j, k, n);
      }
    }
    // low x high y
    if (dlo.x < lx && dhi.y - 1 > hy) {
      if (i == lx && j == hy) {
        dif(i - 1, j + 1, k, n) = dif(i, j, k, n);
      }
    }
    // high x low y
    if (dhi.x - 1 > hx && dlo.y < ly) {
      if (i == hx && j == ly) {
        dif(i + 1, j - 1, k, n) = dif(i, j, k, n);
      }
    }
    // high x high y
    if (dhi.x - 1 > hx && dhi.y - 1 > hy) {
      if (i == hx && j == hy) {
        dif(i + 1, j + 1, k, n) = dif(i, j, k, n);
      }
    }
#if AMREX_SPACEDIM == 3
    // low x low z
    if (dlo.x < lx && dlo.z < lz) {
      if (i == lx && k == lz) {
        dif(i - 1, j, k - 1, n) = dif(i, j, k, n);
      }
    }
    // low x high z
    if (dlo.x < lx && dhi.z - 1 > hz) {
      if (i == lx && k == hz) {
        dif(i - 1, j, k + 1, n) = dif(i, j, k, n);
      }
    }
    // high x low z
    if (dhi.x - 1 > hx && dlo.z < lz) {
      if (i == hx && k == lz) {
        dif(i + 1, j, k - 1, n) = dif(i, j, k, n);
      }
    }
    // high x high z
    if (dhi.x - 1 > hx && dhi.z - 1 > hz) {
      if (i == hx && k == hz) {
        dif(i + 1, j, k + 1, n) = dif(i, j, k, n);
      }
    }

    // low y low z
    if (dlo.y < ly && dlo.z < lz) {
      if (j == ly && k == lz) {
        dif(i, j - 1, k - 1, n) = dif(i, j, k, n);
      }
    }
    // low y hi z
    if (dlo.y < ly && dhi.z - 1 > hz) {
      if (j == ly && k == hz) {
        dif(i, j - 1, k + 1, n) = dif(i, j, k, n);
      }
    }
    // hi y low z
    if (dhi.y - 1 > hy && dlo.z < lz) {
      if (j == hy && k == lz) {
        dif(i, j + 1, k - 1, n) = dif(i, j, k, n);
      }
    }
    // low y low z
    if (dhi.y - 1 > hy && dhi.z - 1 > hz) {
      if (j == hy && k == hz) {
        dif(i, j + 1, k + 1, n) = dif(i, j, k, n);
      }
    }
    //============================= Corners
    //====================================== low x low y low z
    if (dlo.x < lx && dlo.y < ly && dlo.z < lz) {
      if (i == lx && j == ly && k == lz) {
        dif(i - 1, j - 1, k - 1, n) = dif(i, j, k, n);
      }
    }
    // hi x low y low z
    if (dhi.x - 1 > hx && dlo.y < ly && dlo.z < lz) {
      if (i == hx && j == ly && k == lz) {
        dif(i + 1, j - 1, k - 1, n) = dif(i, j, k, n);
      }
    }

    // low x hi y low z
    if (dlo.x < lx && dhi.y - 1 > ly && dlo.z < lz) {
      if (i == lx && j == hy && k == lz) {
        dif(i - 1, j + 1, k - 1, n) = dif(i, j, k, n);
      }
    }

    // low x low y hi z
    if (dlo.x < lx && dlo.y < ly && dhi.z - 1 > hz) {
      if (i == lx && j == ly && k == hz) {
        dif(i - 1, j - 1, k + 1, n) = dif(i, j, k, n);
      }
    }

    // hi x hi y low z
    if (dhi.x - 1 > hx && dhi.y - 1 > hy && dlo.z < lz) {
      if (i == hx && j == hy && k == lz) {
        dif(i + 1, j + 1, k - 1, n) = dif(i, j, k, n);
      }
    }

    // hi x low y hi z
    if (dhi.x - 1 > hx && dlo.y < ly && dhi.z - 1 > lz) {
      if (i == hx && j == ly && k == hz) {
        dif(i + 1, j - 1, k + 1, n) = dif(i, j, k, n);
      }
    }

    // low x hi y hi z
    if (dlo.x < lx && dhi.y - 1 > hy && dlo.z - 1 > hz) {
      if (i == lx && j == hy && k == hz) {
        dif(i - 1, j + 1, k + 1, n) = dif(i, j, k, n);
      }
    }

    // hi x hi y hi z
    if (dhi.x - 1 > hx && dhi.y - 1 > hy && dlo.z - 1 > hz) {
      if (i == hx && j == hy && k == hz) {
        dif(i + 1, j + 1, k + 1, n) = dif(i, j, k, n);
      }
    }
#endif // AMREX_SPACEDIM == 3
#endif // AMREX_SPACEDIM > 1
  }
}

#endif
