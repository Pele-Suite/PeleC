# The same problem as NSCBC-FlameOutflow, with CH4/air chemistry instead of
# H2/air.  It exists because the H2 case runs at Le ~ 0.3, where the imposed
# wrinkle is thermodiffusively unstable, and a reader is entitled to ask
# whether that contaminates the boundary measurement.  CH4/air at phi = 0.75
# has Le(CH4) = 0.97 -- neutral -- so it settles the question, and its larger
# expansion ratio makes it an independent test of the du_out/dn scaling.
#
# Sources are the parent case's; only the mechanism differs.
set(PELE_PHYSICS_EOS_MODEL Fuego)
set(PELE_PHYSICS_CHEMISTRY_MODEL drm19)
set(PELE_PHYSICS_TRANSPORT_MODEL Simple)
set(PELE_PHYSICS_ENABLE_SOOT OFF)
set(PELE_PHYSICS_ENABLE_SPRAY OFF)
set(PELE_PHYSICS_SPRAY_FUEL_NUM 0)
include(BuildExeAndLib)
