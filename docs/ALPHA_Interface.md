# SMRT–ALPHA Interface

SMRT is the geometry, mesh, and characteristic-line preprocessing component used with the ALPHA neutron lattice physics code developed by the ALPHA Group, Nuclear Reactor Physics, School of Physics, Zhejiang University.

SMRT does not solve the neutron transport equation. Its role is to generate the mesh and characteristic-line data required by ALPHA. ALPHA reads the SMRT `.trk` file, stores the data in arrays, and performs the two-dimensional Method of Characteristics calculation.

The `.trk` files included with the current public examples have been successfully read by ALPHA.

## Data-flow contract

```text
SMRT input card
      │
      ├── GmshGeometryModeling
      │       └── two-dimensional .msh mesh
      │
      ├── RayTracing
      │       └── characteristic-line segmentation
      │
      └── Driver.TrackFile(..., "alpha")
              └── ALPHA-compatible .trk file
                              │
                              └── ALPHA arrays and MOC calculation
```

## Geometry conventions

- The model is two-dimensional.
- The geometric center of every pin-cell model is `(0, 0)`.
- Coordinates in the `.msh` file use the same units as the input geometry parameters.
- The `.trk` file uses the same length units for FSR areas, characteristic-line spacing, segment lengths, and RTM dimensions.
- `GeomType 0` identifies square/rectangular geometry.
- `GeomType 1` identifies a complete hexagonal geometry.
- `GeomType 2` identifies a one-sixth symmetric hexagonal geometry.

## Material-ID contract

Material IDs in the `.trk` file are zero-based. They are assigned by the material order in the input card, not by the textual material name:

```python
uo2 = Material.Material("uo2")
zr = Material.Material("zr")
h2o = Material.Material("h2o")

materials = Material.Materials([uo2, zr, h2o])
```

The resulting IDs are:

| Material object position | Material ID |
| ---: | ---: |
| `uo2` | `0` |
| `zr` | `1` |
| `h2o` | `2` |

The mapping must remain unchanged between mesh creation and track-file generation. When adding or reordering materials, the ALPHA-side material data must be updated consistently.

## ALPHA `.trk` input

The current ALPHA text output contains four logical parts:

1. `PART_0_PIN_INFO`: pin/RTM information, FSR count, material IDs, FSR IDs, and mesh-face references.
2. `PART_1_PIN_LAYOUT_INFO`: pin-map and sublayer information.
3. `PART_2_RAY_PARAM_INFO`: azimuthal count, azimuthal angles, track-family counts, and effective characteristic-line spacing.
4. `PART_3_RTM_GEOM_INFO`: FSR areas, track counts, segment counts, segment lengths, swept FSR IDs, and hexagonal reflection/continuation IDs.

For every characteristic-line segment, ALPHA can obtain:

- the azimuthal direction;
- the characteristic-line ID;
- the segment ID;
- the segment length;
- the FSR/material region traversed by the segment.

For hexagonal geometry, the track record also stores the reflected/continuation azimuth ID and reflected/continuation track ID. The physical reflected angle is recovered through the azimuth table in `PART_2_RAY_PARAM_INFO`.

## Recommended handoff procedure

1. Run `model.py` to create the `.msh` mesh.
2. Confirm the material order in the input card.
3. Run `tracking.py` with the intended `spacing`, `azim`, and, for supported hexagonal cases, `scan_mode`.
4. Pass the resulting `.trk` file to ALPHA.
5. Use the same material-ID mapping and length units in the ALPHA calculation.

Example output names are:

```text
HexPin-full-0.03-24.trk
RectPin-0.03-24.trk
```

The name encodes the geometry, scan mode where applicable, nominal track spacing, and azimuthal count.

## Compatibility notes

The `.trk` format is an interface contract between SMRT and ALPHA. Changes to section names, field order, IDs, units, or reflection fields should be treated as format changes and documented in `CHANGELOG.md`. The SMRT and ALPHA versions used together should therefore be recorded in scientific calculations.
