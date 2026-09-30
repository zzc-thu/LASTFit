# Blunt cone at 1 degree angle of attack

This case evaluates the asymmetric three-dimensional steady base flow and the
polar-axis treatment for a Mach 6 spherical-nosed cone.

## Included configuration

- Grid declared by `Config.cfg`: 120 x 151 x 40
- Mach number: 6
- Reynolds number: 10000
- Cone half-angle: 5 degrees
- Angle of attack: 1 degree

## Release status

The solver configuration, eight-piece initial grid, and matching eight-rank
flow/shock restart for 120 x 151 x 40 are included. The manuscript currently
reports 30 x 101 x 40, while the legacy PVTS file used with the plotting data
has a different extent. These records must not be presented as one run until
their provenance is resolved.

`output/` contains one matching three-dimensional field. The windward/leeward
wall-pressure samples and their manuscript-grid provenance still need to be
resolved before the comparison is fully reproducible.
