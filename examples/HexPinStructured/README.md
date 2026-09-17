# HexPin — Structured Mesh

This is a structured-mesh test case for a hexagonal pin cell. The geometry is created by `GmshGeometryModeling`, followed by characteristic-line tracking to generate the mesh and characteristic-line segment information.

## Files

- `model.py`: Creates the materials, hexagonal pin-cell geometry, and structured unstructured mesh. It writes `HexPin.msh` and `HexPin.svg`.
- `tracking.py`: Reads `HexPin.msh`, performs characteristic-line tracking, and writes `HexPin-full-0.03-24.trk`.

## Running the example

Run the following commands from this directory:

```powershell
python model.py
python tracking.py
```

Run `model.py` before `tracking.py`.

## Geometry and material parameters

The model contains three materials in the order `[uo2, zr, co2]`: fuel, cladding, and coolant. The main geometry parameters are:

- `PinPitch = 1.08`: hexagonal pin-cell pitch.
- `radius1 = 0.4225`: fuel radius.
- `radius2 = 0.487571`: outer cladding radius.
- `edge_width = PinPitch / sqrt(3)`: parameter defining the outer hexagonal boundary.

## Structured mesh parameters

```python
radial_meshes=[3, 1, 2]
azimuthal_meshes=[24, 24, 24]
mesh_mode="structured"
```

The lists follow the material order `[uo2, zr, co2]`:

- `radial_meshes`: number of radial layers in each material region.
- `azimuthal_meshes`: number of azimuthal sectors in each material region.

In this case, the fuel, cladding, and coolant regions contain 3, 1, and 2 radial layers, respectively, and each region has 24 azimuthal sectors. This mode is useful when the mesh topology and the number of subdivisions need to be explicitly controlled.

## Characteristic-line tracking parameters

The parameters in `tracking.py` are:

- `spacing = 0.03`: transverse spacing between characteristic lines.
- `azim = 24`: number of azimuthal directions.
- `num_proc = 4`: number of parallel processes.
- `scan_mode = "full"`: tracks the complete hexagonal pin cell. To test one-sixth symmetry tracking, set it to `"one_sixth"`.

The tracking result is written to `HexPin-full-0.03-24.trk`.
