#ifndef INDEX_DEFINES_H
#define INDEX_DEFINES_H

#include "mechanism.H"
#include "TransportTypes.H"

#define URHO 0
#define UMX 1
#define UMY 2
#define UMZ 3
#define UEDEN 4
#define UEINT 5
#define UTEMP 6
#define UFA 7

#define QRHO 0
#define QU 1
#define QV 2
#define QW 3
#define QGAME 4
#define QPRES 5
#define QREINT 6
#define QTEMP 7
#define QFA 8

#define QGAMC 0
#define QC 1
#define QCSML 2
#define QDPDR 3
#define QDPDE 4
#define QRSPEC 5

#define GDRHO 0
#define GDU 1
#define GDV 2
#define GDW 3
#define GDPRES 4
#define GDGAME 5

#define dComp_rhoD 0
#define dComp_mu (dComp_rhoD + NUM_SPECIES)
#define dComp_xi (dComp_mu + 1)
#define dComp_lambda (dComp_xi + 1)

// Optional multicomponent diffusion (pelec.use_multicomponent_diffusion):
// available for the ideal gas (Fuego) EOS with Simple transport, without
// sprays or soot, and only when PelePhysics provides the extended transport
// API (feature macro PELE_PHYSICS_TRANSPORT_EXTENDED_API in TransportTypes.H).
// The Soret and Dufour terms additionally need the LCT conductivity models,
// which PelePhysics compiles in with USE_LCT_TRANSPORT = TRUE (GNU make) or
// PELE_PHYSICS_ENABLE_LCT_TRANSPORT (CMake), i.e. PELE_TRANSPORT_LCT.
#if defined(USE_FUEGO_EOS) && defined(USE_SIMPLE_TRANSPORT) && \
  !defined(PELE_USE_SPRAY) && !defined(PELE_USE_SOOT) &&       \
  defined(PELE_PHYSICS_TRANSPORT_EXTENDED_API)
#define PELEC_MULTICOMPONENT_DIFFUSION 1
#endif
// Additional transport coefficients used by the multicomponent path:
// the flux diffusion matrix rho*Y_k*D_kj (row-major, k*NUM_SPECIES + j), the
// rescaled thermal diffusion ratios and the Gershgorin diffusivity bound of
// the multicomponent operator used by the species-diffusion time-step limit
#define dComp_Dmat (dComp_lambda + 1)
#define dComp_chit (dComp_Dmat + NUM_SPECIES * NUM_SPECIES)
#define dComp_mcDgersh (dComp_chit + NUM_SPECIES)
#define nCompTrMC (dComp_mcDgersh + 1)

#define comp_Cs2 0
#define comp_CI (comp_Cs2 + 1)
#define comp_PrT (comp_CI + 1)
#define comp_Cs2ovPrT comp_PrT
#define nCompC (comp_PrT + 1)

#define NTHERM 7
#define QTHERM 8

#ifndef NUM_ADV
#define NUM_ADV 0
#endif
#ifndef NUM_AUX
#define NUM_AUX 0
#endif

#define UFS (UFA + NUM_ADV)
#define UFX (UFS + NUM_SPECIES)
#define ULIN (UFX + NUM_AUX)
#ifdef PELE_USE_SOOT
#define NUM_LIN (NUM_SOOT_MOMENTS + 1)
#define UFSOOT ULIN
#else
#define NUM_LIN 0
#endif
#define NVAR (NTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN)

#define QFS (QFA + NUM_ADV)
#define QFX (QFS + NUM_SPECIES)
#define QLIN (QFX + NUM_AUX)
#ifdef PELE_USE_SOOT
#define QFSOOT QLIN
#endif
#define QVAR (QTHERM + NUM_ADV + NUM_SPECIES + NUM_AUX + NUM_LIN)
#define NQAUX 6
#define NGDNV 6

#define R_RHO 0
#define R_UN 1
#define R_UT1 2
#define R_UT2 3
#define R_P 4
#define R_ADV 5
#define R_Y (R_ADV + NUM_ADV)
#define R_LIN (R_Y + NUM_SPECIES + NUM_AUX)
#define R_NUM (5 + NUM_SPECIES + NUM_ADV + NUM_LIN + NUM_AUX)

#if AMREX_SPACEDIM == 1
// hopefully we won't do EB in 1D
#define NLSPTS 1
#define NEL_TRIMAT 1

#elif AMREX_SPACEDIM == 2

#define NLSPTS 8
#define NEL_TRIMAT 3

#elif AMREX_SPACEDIM == 3

#define NLSPTS 26
#define NEL_TRIMAT 6

#endif

#endif
