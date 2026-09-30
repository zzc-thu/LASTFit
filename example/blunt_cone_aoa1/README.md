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

The solver configuration and an initial grid with these dimensions have been
identified. The manuscript currently reports 30 x 101 x 40, while the legacy
PVTS file used with the plotting data has a different extent. These records
must not be presented as one run until their provenance is resolved.

The final release output should include one three-dimensional pressure field
and the windward/leeward wall-pressure samples used in the comparison.
