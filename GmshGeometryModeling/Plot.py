"""Plot Gmsh meshes using colors supplied by the input-card materials."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from RayTracing.GmshParser import GmshFileTranslator


plt.rcParams.update(
    {
        "font.family": "Times New Roman",
        "font.serif": ["Times New Roman"],
        "mathtext.fontset": "stix",
        "svg.fonttype": "none",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)


def _material_colors(materials) -> dict[str, tuple[float, ...] | str]:
    """Return a material-name/color mapping without importing an input card."""

    if materials is None:
        return {}
    material_list = getattr(materials, "materials", materials)
    return {
        material.name: material.color
        for material in material_list
        if getattr(material, "color", None) is not None
    }


def plot_mesh(
    mesh_path: str | Path,
    materials=None,
    geom_type: str = "hex",
    output_path: str | Path | None = None,
    output_format: str = "svg",
    title: str | None = None,
) -> Path:
    """Render a mesh using colors from the supplied Material objects.

    ``materials`` may be either a list of ``Material`` objects or the
    ``Material.Materials`` collection. It is passed by the input card, so
    this function never imports or executes ``model.py``.
    """

    mesh_path = Path(mesh_path).resolve()
    output_format = output_format.lower().lstrip(".")
    if output_format not in {"svg", "pdf", "png"}:
        raise ValueError("output_format must be 'svg', 'pdf', or 'png'")
    if not mesh_path.exists():
        raise FileNotFoundError(f"Mesh file not found: {mesh_path}")

    if output_path is None:
        output_path = mesh_path.with_suffix(f".{output_format}")
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)

    mesh = GmshFileTranslator(mesh_path, geom_type=geom_type).extractDataFromGmshFile()
    colors = _material_colors(materials)
    fig, axis = plt.subplots(figsize=(7, 7))
    material_names = []

    for face_index in range(sum(mesh.getNumFaces())):
        face_id = mesh.getFaceID(face_index)
        coordinates = mesh.getFaceNodeCoordinates(face_id)
        material = mesh.getPhysicalName(face_id)
        if material not in material_names:
            material_names.append(material)
        axis.fill(
            coordinates[0],
            coordinates[1],
            facecolor=colors.get(material, "white"),
            edgecolor="#777777",
            linewidth=0.15,
        )

    axis.set_aspect("equal")
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_title(title or f"Gmsh {mesh_path.stem} mesh")
    axis.legend(
        handles=[
            Patch(facecolor=colors.get(name, "white"), edgecolor="0.4", label=name)
            for name in material_names
        ],
        loc="upper right",
    )
    fig.tight_layout()
    fig.savefig(output_path, format=output_format)
    plt.close(fig)
    return output_path
