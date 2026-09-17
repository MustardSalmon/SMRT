# RectPin — Structured Mesh

This is a structured-mesh test case for a square pin cell. The geometry is created by `GmshGeometryModeling`, followed by characteristic-line tracking to generate the mesh and characteristic-line segment information.

## Files

- `model.py`: Creates the materials, square pin-cell geometry, and structured unstructured mesh. It writes `RectPin.msh` and `RectPin.svg`.
- `tracking.py`: Reads `RectPin.msh`, performs characteristic-line tracking, and writes `RectPin-0.03-24.trk`.

## Running the example

Run the following commands from this directory:

```powershell
python model.py
python tracking.py
```

Run `model.py` before `tracking.py`.

## Geometry and material parameters

The model contains three materials in the order `[uo2, zr, h2o]`: fuel, cladding, and water. The main geometry parameters are:

- `PinPitch = 1.26`: width and height of the square pin cell.
- `radius1 = 0.46`: fuel radius.
- `radius2 = 0.54`: outer cladding radius.

## Structured mesh parameters

```python
radial_meshes=[3, 1, 2]
azimuthal_meshes=[24, 24, 24]
mesh_mode="structured"
```

The lists follow the material order `[uo2, zr, h2o]`:

- `radial_meshes`: number of radial layers in each material region.
- `azimuthal_meshes`: number of azimuthal sectors in each material region.

In this case, the fuel, cladding, and water regions contain 3, 1, and 2 radial layers, respectively, and each region has 24 azimuthal sectors. The outer water region is finally mapped to the square boundary, while the mesh remains controlled by these structured subdivision parameters.

## Characteristic-line tracking parameters

The parameters in `tracking.py` are:

- `spacing = 0.03`: transverse spacing between characteristic lines.
- `azim = 24`: number of azimuthal directions.
- `num_proc = 4`: number of parallel processes.

The tracking result is written to `RectPin-0.03-24.trk`.
