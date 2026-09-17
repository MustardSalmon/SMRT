# Changelog

All notable changes to SMRT are documented here.

## [1.0.0] — 2026-09-17

Initial open-source release.

### Added

- Gmsh-based square pin-cell geometry generation.
- Gmsh-based hexagonal pin-cell geometry generation.
- Structured mesh generation with explicit radial and azimuthal subdivisions.
- Automatic mesh generation with material-specific target mesh sizes.
- Gmsh `.msh` parsing for first-order two-dimensional triangle and quadrilateral meshes.
- Characteristic-line generation for multiple azimuthal directions.
- Full-cell and one-sixth symmetric hexagonal tracking modes.
- ALPHA-compatible `.trk` text-file export.
- Four documented examples:
  - `HexPinStructured`;
  - `HexPinAutomatic`;
  - `RectPinStructured`;
  - `RectPinAutomatic`.

### Scope

This release provides preprocessing data for the ALPHA neutron lattice physics code. It does not include the neutron transport solver. More complex geometry-generation modules are planned for future releases.
