# RectPin — Automatic Mesh

This is an automatic-mesh test case for a square pin cell. The geometry is created by `GmshGeometryModeling`, followed by characteristic-line tracking to generate the mesh and characteristic-line segment information.

## Files

- `model.py`: Creates the materials, square pin-cell geometry, and automatic unstructured mesh. It writes `RectPin.msh` and `RectPin.svg`.
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

## Automatic mesh parameters

```python
mesh_sizes=[0.1, 0.2, 0.3]
mesh_mode="automatic"
```

The list follows the material order `[uo2, zr, h2o]`. Each value is the target characteristic mesh size for the corresponding material region, using the same units as the geometry:

- `0.1`: target element size in the fuel region. The smaller value produces a denser mesh.
- `0.2`: target element size in the cladding region.
- `0.3`: target element size in the water region. The larger value produces a coarser mesh.

These values are Gmsh target-size controls, not strict guarantees of element area or element count. The final mesh is also affected by geometry boundaries, curvature, and material interfaces.

## Characteristic-line tracking parameters

The parameters in `tracking.py` are:

- `spacing = 0.03`: transverse spacing between characteristic lines.
- `azim = 24`: number of azimuthal directions.
- `num_proc = 4`: number of parallel processes.

The tracking result is written to `RectPin-0.03-24.trk`.
