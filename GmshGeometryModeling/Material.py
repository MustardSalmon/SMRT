"""Material definitions for the Gmsh input-card interface."""

from __future__ import annotations

from dataclasses import dataclass

from GmshGeometryModeling import GmshInterface


@dataclass
class Material:
    """A named material region used by the Gmsh input-card interface."""

    name: str
    color: tuple[float, float, float] | None = None
    mesh_size: float | None = None

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("material name must not be empty")

    def setColor(self, red: float, green: float, blue: float) -> None:
        self.color = (float(red), float(green), float(blue))

    def setMeshSize(self, mesh_size: float) -> None:
        if mesh_size <= 0.0:
            raise ValueError("material mesh size must be positive")
        self.mesh_size = float(mesh_size)


class Materials:
    """Collection used to initialize the active Gmsh model."""

    def __init__(self, materials: list[Material] | tuple[Material, ...]):
        self.materials = list(materials)

    def create(self) -> None:
        names = [material.name for material in self.materials]
        if len(names) != len(set(names)):
            raise ValueError("material names must be unique")
        GmshInterface.initialize_model("GmshModel")

    def getNames(self) -> list[str]:
        return [material.name for material in self.materials]
