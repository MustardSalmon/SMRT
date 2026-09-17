"""Small helpers for Gmsh input cards."""

from __future__ import annotations

from pathlib import Path

from GmshGeometryModeling import GmshInterface


def save(file_path: str | Path, boundary_geom: str = "hex") -> Path:
    return GmshInterface.save_gmsh_file(file_path, boundary_geom)


def reset() -> None:
    GmshInterface.reset_model()
