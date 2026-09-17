# Contributing to SMRT

Thank you for contributing to SMRT. The project provides Gmsh-based preprocessing for two-dimensional Method of Characteristics neutron transport calculations and interfaces with the ALPHA neutron lattice physics code through the `.trk` format.

## Scope

Contributions should preserve the current responsibilities of the project:

- geometry and mesh generation;
- Gmsh mesh parsing;
- characteristic-line generation;
- ALPHA-compatible track-file export;
- documented square and hexagonal pin-cell examples.

The neutron transport solver itself is outside the scope of this repository.

## Development setup

Use Python 3.12.3 or a compatible Python 3 environment. Install the tested runtime dependencies with:

```powershell
python -m pip install -r requirements.txt
```

SciPy is optional. If it is installed, selected spatial searches may use its accelerated implementation.

## Before submitting a change

Please verify all four public examples:

```powershell
python examples/HexPinStructured/model.py
python examples/HexPinStructured/tracking.py

python examples/HexPinAutomatic/model.py
python examples/HexPinAutomatic/tracking.py

python examples/RectPinStructured/model.py
python examples/RectPinStructured/tracking.py

python examples/RectPinAutomatic/model.py
python examples/RectPinAutomatic/tracking.py
```

Please check that the `.msh`, `.svg`, and `.trk` outputs are generated successfully and that no machine-specific absolute paths are introduced into the source or documentation.

## Input-card and format changes

Changes to material ordering, mesh IDs, `.trk` sections, field order, units, azimuth conventions, or hexagonal reflection identifiers may affect ALPHA compatibility. Such changes must include:

- an update to `docs/FileFormats.md` or `docs/ALPHA_Interface.md`;
- an entry in `CHANGELOG.md`;
- a regenerated example output or an explicit compatibility note.

## Documentation and style

- Keep public documentation in English.
- Use clear, descriptive names for geometry and material parameters.
- Keep example input cards executable from their own directories.
- Do not commit credentials, private paths, generated cache directories, or unrelated local IDE files.
