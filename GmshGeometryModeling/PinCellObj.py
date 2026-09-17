"""Pin-cell objects backed by Gmsh."""

from __future__ import annotations

from GmshGeometryModeling.GmshInterface import HexPinBuilder, RectPinBuilder


class HexRodLikePinCell:
    """One hexagonal cell containing fuel, cladding, and coolant regions."""

    def __init__(
        self,
        edge_length: float,
        radii: list[float] | tuple[float, float],
        materials,
        radial_meshes=None,
        azimuthal_meshes=None,
        mesh_mode: str = "structured",
        mesh_sizes=None,
    ):
        if len(radii) != 2:
            raise ValueError("radii must contain fuel and clad radii")
        if len(materials) != 3:
            raise ValueError("materials must contain fuel, clad, and coolant")
        if mesh_mode not in {"structured", "automatic"}:
            raise ValueError("mesh_mode must be 'structured' or 'automatic'")
        if mesh_mode == "structured":
            if radial_meshes is None or azimuthal_meshes is None:
                raise ValueError("structured mode requires radial_meshes and azimuthal_meshes")
            if len(radial_meshes) != 3 or len(azimuthal_meshes) != 3:
                raise ValueError("mesh settings must contain three region values")
        else:
            radial_meshes = [1, 1, 1] if radial_meshes is None else radial_meshes
            azimuthal_meshes = [12, 12, 12] if azimuthal_meshes is None else azimuthal_meshes
        self.edge_length = float(edge_length)
        self.radii = tuple(float(radius) for radius in radii)
        self.materials = list(materials)
        self.radial_meshes = list(radial_meshes)
        self.azimuthal_meshes = list(azimuthal_meshes)
        self.mesh_mode = mesh_mode
        self.orientation = "y"
        if mesh_sizes is None and mesh_mode == "automatic":
            mesh_sizes = [material.mesh_size for material in materials]
        if mesh_mode == "automatic":
            if len(mesh_sizes) != 3 or any(size is None or size <= 0 for size in mesh_sizes):
                raise ValueError("automatic mode requires three positive mesh_sizes")
            self.mesh_sizes = [float(size) for size in mesh_sizes]
        else:
            self.mesh_sizes = None
        self._builder = None

    def create(self) -> None:
        self._builder = HexPinBuilder(
            edge_length=self.edge_length,
            radii=self.radii,
            material_names=[material.name for material in self.materials],
            orientation=self.orientation,
            mesh_mode=self.mesh_mode,
            radial_meshes=self.radial_meshes,
            azimuthal_meshes=self.azimuthal_meshes,
            mesh_sizes=self.mesh_sizes,
        )
        self._builder.create_geometry()


class RectRodLikePinCell:
    """A square pin cell with concentric circular inner regions."""

    def __init__(
        self,
        width_X: float,
        width_Y: float,
        radii,
        materials,
        radial_meshes=None,
        azimuthal_meshes=None,
        mesh_mode: str = "structured",      # automatic
        mesh_sizes=None,
    ):
        if len(materials) != len(radii) + 1:
            raise ValueError("the number of materials must equal len(radii) + 1")
        if mesh_mode not in {"structured", "automatic"}:
            raise ValueError("mesh_mode must be 'structured' or 'automatic'")
        if mesh_mode == "structured":
            if radial_meshes is None or azimuthal_meshes is None:
                raise ValueError("structured mode requires radial_meshes and azimuthal_meshes")
            if len(radial_meshes) != len(materials) or len(azimuthal_meshes) != len(materials):
                raise ValueError("mesh settings must match the number of materials")
        else:
            radial_meshes = [1] * len(materials) if radial_meshes is None else radial_meshes
            azimuthal_meshes = [24] * len(materials) if azimuthal_meshes is None else azimuthal_meshes
        if mesh_sizes is None and mesh_mode == "automatic":
            mesh_sizes = [material.mesh_size for material in materials]
        if mesh_mode == "automatic":
            if len(mesh_sizes) != len(materials) or any(size is None or size <= 0 for size in mesh_sizes):
                raise ValueError("automatic mode requires one positive mesh_size per material")

        self.width_X = float(width_X)
        self.width_Y = float(width_Y)
        self.radii = [float(radius) for radius in radii]
        self.materials = list(materials)
        self.radial_meshes = list(radial_meshes)
        self.azimuthal_meshes = list(azimuthal_meshes)
        self.mesh_mode = mesh_mode
        self.mesh_sizes = None if mesh_sizes is None else [float(size) for size in mesh_sizes]
        self.orientation = "square"
        self._builder = None

    def setMeshSize(self, mesh_size: float) -> None:
        # Kept for API symmetry with HexRodLikePinCell. Automatic mode uses
        # the per-material mesh_sizes list instead.
        if mesh_size <= 0:
            raise ValueError("mesh size must be positive")

    def create(self) -> None:
        self._builder = RectPinBuilder(
            width_x=self.width_X,
            width_y=self.width_Y,
            radii=self.radii,
            material_names=[material.name for material in self.materials],
            mesh_mode=self.mesh_mode,
            radial_meshes=self.radial_meshes,
            azimuthal_meshes=self.azimuthal_meshes,
            mesh_sizes=self.mesh_sizes,
        )
        self._builder.create_geometry()
