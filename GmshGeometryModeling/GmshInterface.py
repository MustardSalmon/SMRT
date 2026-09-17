"""Gmsh backend with structured and automatic mesh modes."""

from __future__ import annotations

import math
from pathlib import Path

import gmsh


_ACTIVE_MODEL = False
_MESH_MODE = "structured"
_REGIONS: list[tuple[int, str, float | None]] = []
_BOUNDARY_LINES: list[int] = []
_STRUCTURED_MESH: dict | None = None


def initialize_model(model_name: str = "GmshModel") -> None:
    global _ACTIVE_MODEL, _REGIONS, _BOUNDARY_LINES, _STRUCTURED_MESH
    if _ACTIVE_MODEL:
        reset_model()
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.ElementOrder", 1)
    gmsh.model.add(model_name)
    _REGIONS = []
    _BOUNDARY_LINES = []
    _STRUCTURED_MESH = None
    _ACTIVE_MODEL = True


def reset_model() -> None:
    global _ACTIVE_MODEL, _REGIONS, _BOUNDARY_LINES, _STRUCTURED_MESH
    if _ACTIVE_MODEL:
        gmsh.finalize()
    _ACTIVE_MODEL = False
    _REGIONS = []
    _BOUNDARY_LINES = []
    _STRUCTURED_MESH = None


def _set_automatic_mesh_fields() -> None:
    fields = []
    for surface, _material, mesh_size in _REGIONS:
        if mesh_size is None:
            raise ValueError("automatic mode requires a mesh size for every material")
        constant = gmsh.model.mesh.field.add("Constant")
        gmsh.model.mesh.field.setNumber(constant, "VIn", mesh_size)
        gmsh.model.mesh.field.setNumber(constant, "VOut", 1.0e22)
        gmsh.model.mesh.field.setNumbers(constant, "SurfacesList", [surface])
        fields.append(constant)

    minimum = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(minimum, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(minimum)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)


def _write_structured_msh(file_path: str | Path) -> Path:
    if _STRUCTURED_MESH is None:
        raise RuntimeError("structured mesh data has not been created")
    output = Path(file_path).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    data = _STRUCTURED_MESH
    names = data["material_names"]
    physical_names = [(2, index + 1, name) for index, name in enumerate(names)]
    physical_names.append((1, 10, "_boundary"))

    elements = []
    for nodes in data["boundary_lines"]:
        elements.append((1, 10, 10, nodes))
    for material_id, nodes in data["quads"]:
        elements.append((3, material_id, material_id, nodes))
    for material_id, nodes in data["tris"]:
        elements.append((2, material_id, material_id, nodes))

    with output.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        stream.write("$PhysicalNames\n%d\n" % len(physical_names))
        for dimension, tag, name in physical_names:
            stream.write(f'{dimension} {tag} "{name}"\n')
        stream.write("$EndPhysicalNames\n")
        stream.write("$Nodes\n%d\n" % len(data["nodes"]))
        for node_id, (x, y) in enumerate(data["nodes"], start=1):
            stream.write(f"{node_id} {x:.16g} {y:.16g} 0\n")
        stream.write("$EndNodes\n")
        stream.write("$Elements\n%d\n" % len(elements))
        for element_id, (element_type, physical, entity, nodes) in enumerate(elements, start=1):
            node_text = " ".join(str(node) for node in nodes)
            stream.write(f"{element_id} {element_type} 2 {physical} {entity} {node_text}\n")
        stream.write("$EndElements\n")
    return output


def save_gmsh_file(file_path: str | Path, boundary_geom: str = "hex") -> Path:
    """Mesh and save the active model using its selected mesh mode."""

    if not _ACTIVE_MODEL:
        raise RuntimeError("no active Gmsh model; call Materials.create() first")
    if boundary_geom.lower() not in {"hex", "rect"}:
        raise ValueError("boundary_geom must be 'hex' or 'rect'")
    try:
        if _MESH_MODE == "structured":
            return _write_structured_msh(file_path)

        gmsh.model.geo.synchronize()
        _set_automatic_mesh_fields()
        for index, (surface, material_name, _mesh_size) in enumerate(_REGIONS, start=1):
            group = gmsh.model.addPhysicalGroup(2, [surface], index)
            gmsh.model.setPhysicalName(2, group, material_name)
        boundary_group = gmsh.model.addPhysicalGroup(1, _BOUNDARY_LINES, 10)
        gmsh.model.setPhysicalName(1, boundary_group, "_boundary")
        gmsh.model.mesh.generate(2)
        output = Path(file_path).resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        gmsh.write(str(output))
        return output
    finally:
        reset_model()


class HexPinBuilder:
    """Build one hexagonal pin cell in structured or automatic mode."""

    def __init__(
        self,
        edge_length: float,
        radii,
        material_names,
        orientation: str = "y",
        model_name: str = "HexPin",
        mesh_mode: str = "structured",
        radial_meshes=None,
        azimuthal_meshes=None,
        mesh_sizes=None,
    ):
        if len(radii) != 2 or not (0.0 < radii[0] < radii[1]):
            raise ValueError("radii must satisfy 0 < fuel_radius < clad_radius")
        if len(material_names) != 3 or len(set(material_names)) != 3:
            raise ValueError("three distinct material names are required")
        if edge_length <= 0:
            raise ValueError("edge_length must be positive")
        if orientation not in {"x", "y"}:
            raise ValueError("orientation must be 'x' or 'y'")
        if mesh_mode not in {"structured", "automatic"}:
            raise ValueError("mesh_mode must be 'structured' or 'automatic'")
        if mesh_mode == "structured" and (radial_meshes is None or azimuthal_meshes is None):
            raise ValueError("structured mode requires radial_meshes and azimuthal_meshes")
        if radial_meshes is None:
            radial_meshes = [1, 1, 1]
        if azimuthal_meshes is None:
            azimuthal_meshes = [12, 12, 12]
        if len(radial_meshes) != 3 or len(azimuthal_meshes) != 3:
            raise ValueError("mesh settings must contain three region values")
        if any(int(value) <= 0 for value in list(radial_meshes) + list(azimuthal_meshes)):
            raise ValueError("mesh counts must be positive")
        if mesh_mode == "structured":
            if len(set(azimuthal_meshes)) != 1:
                raise ValueError("concentric regions must use the same azimuthal mesh count")
            if azimuthal_meshes[2] % 6:
                raise ValueError("the outer azimuthal mesh count must be divisible by 6")
        if mesh_mode == "automatic":
            if mesh_sizes is None or len(mesh_sizes) != 3 or any(size is None or size <= 0 for size in mesh_sizes):
                raise ValueError("automatic mode requires three positive mesh_sizes")

        self.edge_length = float(edge_length)
        self.fuel_radius, self.clad_radius = map(float, radii)
        self.material_names = list(material_names)
        self.orientation = orientation
        self.model_name = model_name
        self.mesh_mode = mesh_mode
        self.radial_meshes = [int(value) for value in radial_meshes]
        self.azimuthal_meshes = [int(value) for value in azimuthal_meshes]
        self.mesh_sizes = None if mesh_sizes is None else [float(size) for size in mesh_sizes]
        inradius = math.sqrt(3.0) * self.edge_length / 2.0
        if self.clad_radius >= inradius:
            raise ValueError("clad radius must be smaller than the hexagon inradius")

    def _add_hexagon(self):
        edge = self.edge_length
        if self.orientation == "x":
            angles = [-2.0 * math.pi / 3.0, -math.pi / 3.0, 0.0,
                      math.pi / 3.0, 2.0 * math.pi / 3.0, math.pi]
        else:
            angles = [-5.0 * math.pi / 6.0, -math.pi / 2.0, -math.pi / 6.0,
                      math.pi / 6.0, math.pi / 2.0, 5.0 * math.pi / 6.0]
        vertices = [(edge * math.cos(angle), edge * math.sin(angle)) for angle in angles]
        points = [gmsh.model.geo.addPoint(x, y, 0.0) for x, y in vertices]
        lines = [gmsh.model.geo.addLine(points[index], points[(index + 1) % 6]) for index in range(6)]
        return gmsh.model.geo.addCurveLoop(lines), lines

    @staticmethod
    def _add_circle(radius: float):
        center = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
        points = [
            gmsh.model.geo.addPoint(radius * math.cos(angle), radius * math.sin(angle), 0.0)
            for angle in (0.0, math.pi / 2.0, math.pi, 3.0 * math.pi / 2.0)
        ]
        arcs = [
            gmsh.model.geo.addCircleArc(points[index], center, points[(index + 1) % 4])
            for index in range(4)
        ]
        return gmsh.model.geo.addCurveLoop(arcs)

    def _create_automatic_geometry(self) -> None:
        fuel_loop = self._add_circle(self.fuel_radius)
        clad_loop = self._add_circle(self.clad_radius)
        hex_loop, boundary_lines = self._add_hexagon()
        fuel_surface = gmsh.model.geo.addPlaneSurface([fuel_loop])
        clad_surface = gmsh.model.geo.addPlaneSurface([clad_loop, -fuel_loop])
        coolant_surface = gmsh.model.geo.addPlaneSurface([hex_loop, -clad_loop])
        for surface, material, size in zip(
            (fuel_surface, clad_surface, coolant_surface), self.material_names, self.mesh_sizes
        ):
            _REGIONS.append((surface, material, size))
        _BOUNDARY_LINES.extend(boundary_lines)

    @staticmethod
    def _corrected_radius(radius: float, sectors: int) -> float:
        return math.sqrt(math.pi * radius * radius * 2.0 / (sectors * math.sin(2.0 * math.pi / sectors)))

    @staticmethod
    def _radii_by_equal_area(inner: float, outer: float, count: int) -> list[float]:
        return [math.sqrt(inner * inner + (outer * outer - inner * inner) * index / count)
                for index in range(count + 1)]

    def _create_structured_mesh(self) -> None:
        sectors = self.azimuthal_meshes[0]
        start_angle = -2.0 * math.pi / 3.0 if self.orientation == "x" else -5.0 * math.pi / 6.0
        fuel_radius = self._corrected_radius(self.fuel_radius, sectors)
        clad_radius = self._corrected_radius(self.clad_radius, sectors)
        nodes: list[tuple[float, float]] = []
        node_map: dict[tuple[float, float], int] = {}
        tris: list[tuple[int, tuple[int, int, int]]] = []
        quads: list[tuple[int, tuple[int, int, int, int]]] = []
        boundary_lines: list[tuple[int, int]] = []

        def add_node(x: float, y: float) -> int:
            key = (round(x, 13), round(y, 13))
            if key not in node_map:
                node_map[key] = len(nodes) + 1
                nodes.append((x, y))
            return node_map[key]

        def circular_ring(radius: float) -> list[int]:
            return [add_node(radius * math.cos(start_angle + 2.0 * math.pi * index / sectors),
                             radius * math.sin(start_angle + 2.0 * math.pi * index / sectors))
                    for index in range(sectors)]

        center = add_node(0.0, 0.0)
        fuel_rings = []
        for radius in self._radii_by_equal_area(0.0, fuel_radius, self.radial_meshes[0])[1:]:
            fuel_rings.append(circular_ring(radius))
        for index, node in enumerate(fuel_rings[0]):
            tris.append((1, (center, node, fuel_rings[0][(index + 1) % sectors])))
        for ring_index in range(1, len(fuel_rings)):
            previous, current = fuel_rings[ring_index - 1], fuel_rings[ring_index]
            for index in range(sectors):
                quads.append((1, (previous[index], current[index], current[(index + 1) % sectors],
                                   previous[(index + 1) % sectors])))
        fuel_outer = fuel_rings[-1]

        clad_rings = [fuel_outer]
        clad_radii = self._radii_by_equal_area(fuel_radius, clad_radius, self.radial_meshes[1])
        for radius in clad_radii[1:]:
            clad_rings.append(circular_ring(radius))
        for ring_index in range(len(clad_rings) - 1):
            previous, current = clad_rings[ring_index], clad_rings[ring_index + 1]
            for index in range(sectors):
                quads.append((2, (previous[index], current[index], current[(index + 1) % sectors],
                                   previous[(index + 1) % sectors])))
        clad_outer = clad_rings[-1]

        outer_angles = [start_angle + 2.0 * math.pi * index / 6.0 for index in range(6)]
        outer_vertices = [(self.edge_length * math.cos(angle), self.edge_length * math.sin(angle))
                          for angle in outer_angles]
        outer_boundary = []
        per_edge = sectors // 6
        for index in range(sectors):
            edge_index, local_index = divmod(index, per_edge)
            v0 = outer_vertices[edge_index]
            v1 = outer_vertices[(edge_index + 1) % 6]
            fraction = local_index / per_edge
            outer_boundary.append(add_node(v0[0] + fraction * (v1[0] - v0[0]),
                                           v0[1] + fraction * (v1[1] - v0[1])))

        coolant_rings = [clad_outer]
        for layer in range(1, self.radial_meshes[2] + 1):
            fraction = layer / self.radial_meshes[2]
            ring = []
            for inner_node, outer_node in zip(clad_outer, outer_boundary):
                ix, iy = nodes[inner_node - 1]
                ox, oy = nodes[outer_node - 1]
                ring.append(add_node(ix + fraction * (ox - ix), iy + fraction * (oy - iy)))
            coolant_rings.append(ring)
        for ring_index in range(len(coolant_rings) - 1):
            previous, current = coolant_rings[ring_index], coolant_rings[ring_index + 1]
            for index in range(sectors):
                quads.append((3, (previous[index], current[index], current[(index + 1) % sectors],
                                   previous[(index + 1) % sectors])))
        for index, node in enumerate(outer_boundary):
            boundary_lines.append((node, outer_boundary[(index + 1) % sectors]))

        _STRUCTURED_MESH = {
            "nodes": nodes,
            "tris": tris,
            "quads": quads,
            "boundary_lines": boundary_lines,
            "material_names": self.material_names,
        }
        globals()["_STRUCTURED_MESH"] = _STRUCTURED_MESH

    def create_geometry(self) -> None:
        if not _ACTIVE_MODEL:
            raise RuntimeError("call Materials.create() before creating geometry")
        global _MESH_MODE
        _MESH_MODE = self.mesh_mode
        if self.mesh_mode == "structured":
            self._create_structured_mesh()
        else:
            self._create_automatic_geometry()

    def build(self, output_path: str | Path) -> Path:
        initialize_model(self.model_name)
        try:
            self.create_geometry()
            return save_gmsh_file(output_path)
        except Exception:
            reset_model()
            raise


class RectPinBuilder:
    """Build a square pin cell with circular inner material regions."""

    def __init__(
        self,
        width_x: float,
        width_y: float,
        radii,
        material_names,
        mesh_mode: str = "structured",
        radial_meshes=None,
        azimuthal_meshes=None,
        mesh_sizes=None,
    ):
        if width_x <= 0 or width_y <= 0:
            raise ValueError("rectangle widths must be positive")
        if len(material_names) != len(radii) + 1:
            raise ValueError("the number of materials must equal len(radii) + 1")
        if any(radius <= 0 for radius in radii) or any(a >= b for a, b in zip(radii, radii[1:])):
            raise ValueError("radii must be strictly increasing and positive")
        if mesh_mode not in {"structured", "automatic"}:
            raise ValueError("mesh_mode must be 'structured' or 'automatic'")
        if radial_meshes is None:
            radial_meshes = [1] * len(material_names)
        if azimuthal_meshes is None:
            azimuthal_meshes = [24] * len(material_names)
        if len(radial_meshes) != len(material_names) or len(azimuthal_meshes) != len(material_names):
            raise ValueError("mesh settings must match the number of materials")
        if any(int(value) <= 0 for value in list(radial_meshes) + list(azimuthal_meshes)):
            raise ValueError("mesh counts must be positive")
        if mesh_mode == "structured":
            if len(set(azimuthal_meshes)) != 1:
                raise ValueError("concentric regions must use the same azimuthal mesh count")
            if azimuthal_meshes[-1] % 4:
                raise ValueError("the outer azimuthal mesh count must be divisible by 4")
        elif mesh_sizes is None or len(mesh_sizes) != len(material_names) or any(
            size is None or size <= 0 for size in mesh_sizes
        ):
            raise ValueError("automatic mode requires one positive mesh_size per material")

        self.width_x = float(width_x)
        self.width_y = float(width_y)
        self.radii = [float(radius) for radius in radii]
        self.material_names = list(material_names)
        self.mesh_mode = mesh_mode
        self.radial_meshes = [int(value) for value in radial_meshes]
        self.azimuthal_meshes = [int(value) for value in azimuthal_meshes]
        self.mesh_sizes = None if mesh_sizes is None else [float(size) for size in mesh_sizes]
        self.mesh_size = max(self.mesh_sizes) if self.mesh_sizes else 0.03

    def _add_rectangle(self):
        points = [
            gmsh.model.geo.addPoint(-self.width_x / 2.0, -self.width_y / 2.0, 0.0, self.mesh_size),
            gmsh.model.geo.addPoint(self.width_x / 2.0, -self.width_y / 2.0, 0.0, self.mesh_size),
            gmsh.model.geo.addPoint(self.width_x / 2.0, self.width_y / 2.0, 0.0, self.mesh_size),
            gmsh.model.geo.addPoint(-self.width_x / 2.0, self.width_y / 2.0, 0.0, self.mesh_size),
        ]
        lines = [gmsh.model.geo.addLine(points[index], points[(index + 1) % 4]) for index in range(4)]
        return gmsh.model.geo.addCurveLoop(lines), lines

    @staticmethod
    def _add_circle(radius: float):
        center = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
        points = [
            gmsh.model.geo.addPoint(radius * math.cos(angle), radius * math.sin(angle), 0.0)
            for angle in (0.0, math.pi / 2.0, math.pi, 3.0 * math.pi / 2.0)
        ]
        arcs = [
            gmsh.model.geo.addCircleArc(points[index], center, points[(index + 1) % 4])
            for index in range(4)
        ]
        return gmsh.model.geo.addCurveLoop(arcs)

    def _create_automatic_geometry(self) -> None:
        loops = [self._add_circle(radius) for radius in self.radii]
        rectangle_loop, boundary_lines = self._add_rectangle()
        surfaces = [gmsh.model.geo.addPlaneSurface([loops[0]])]
        for index in range(1, len(loops)):
            surfaces.append(gmsh.model.geo.addPlaneSurface([loops[index], -loops[index - 1]]))
        surfaces.append(gmsh.model.geo.addPlaneSurface([rectangle_loop, -loops[-1]]))
        for surface, material, size in zip(surfaces, self.material_names, self.mesh_sizes):
            _REGIONS.append((surface, material, size))
        _BOUNDARY_LINES.extend(boundary_lines)

    @staticmethod
    def _corrected_radius(radius: float, sectors: int) -> float:
        return math.sqrt(math.pi * radius * radius * 2.0 / (sectors * math.sin(2.0 * math.pi / sectors)))

    @staticmethod
    def _radii_by_equal_area(inner: float, outer: float, count: int) -> list[float]:
        return [math.sqrt(inner * inner + (outer * outer - inner * inner) * index / count)
                for index in range(count + 1)]

    def _create_structured_mesh(self) -> None:
        sectors = self.azimuthal_meshes[0]
        start_angle = -3.0 * math.pi / 4.0
        corrected_radii = [self._corrected_radius(radius, sectors) for radius in self.radii]
        nodes: list[tuple[float, float]] = []
        node_map: dict[tuple[float, float], int] = {}
        tris: list[tuple[int, tuple[int, int, int]]] = []
        quads: list[tuple[int, tuple[int, int, int, int]]] = []
        boundary_lines: list[tuple[int, int]] = []

        def add_node(x: float, y: float) -> int:
            key = (round(x, 13), round(y, 13))
            if key not in node_map:
                node_map[key] = len(nodes) + 1
                nodes.append((x, y))
            return node_map[key]

        def circular_ring(radius: float) -> list[int]:
            return [add_node(radius * math.cos(start_angle + 2.0 * math.pi * index / sectors),
                             radius * math.sin(start_angle + 2.0 * math.pi * index / sectors))
                    for index in range(sectors)]

        center = add_node(0.0, 0.0)
        disk_rings = []
        for radius in self._radii_by_equal_area(0.0, corrected_radii[0], self.radial_meshes[0])[1:]:
            disk_rings.append(circular_ring(radius))
        for index, node in enumerate(disk_rings[0]):
            tris.append((1, (center, node, disk_rings[0][(index + 1) % sectors])))
        for ring_index in range(1, len(disk_rings)):
            previous, current = disk_rings[ring_index - 1], disk_rings[ring_index]
            for index in range(sectors):
                quads.append((1, (previous[index], current[index], current[(index + 1) % sectors],
                                   previous[(index + 1) % sectors])))
        previous_ring = disk_rings[-1]

        for region_index in range(1, len(self.material_names) - 1):
            region_rings = [previous_ring]
            radii = self._radii_by_equal_area(
                corrected_radii[region_index - 1], corrected_radii[region_index],
                self.radial_meshes[region_index]
            )
            for radius in radii[1:]:
                region_rings.append(circular_ring(radius))
            for ring_index in range(len(region_rings) - 1):
                previous, current = region_rings[ring_index], region_rings[ring_index + 1]
                for index in range(sectors):
                    quads.append((region_index + 1, (
                        previous[index], current[index], current[(index + 1) % sectors],
                        previous[(index + 1) % sectors]
                    )))
            previous_ring = region_rings[-1]

        vertices = [
            (-self.width_x / 2.0, -self.width_y / 2.0),
            (self.width_x / 2.0, -self.width_y / 2.0),
            (self.width_x / 2.0, self.width_y / 2.0),
            (-self.width_x / 2.0, self.width_y / 2.0),
        ]
        per_edge = sectors // 4
        outer_boundary = []
        for index in range(sectors):
            edge_index, local_index = divmod(index, per_edge)
            v0, v1 = vertices[edge_index], vertices[(edge_index + 1) % 4]
            fraction = local_index / per_edge
            outer_boundary.append(add_node(v0[0] + fraction * (v1[0] - v0[0]),
                                           v0[1] + fraction * (v1[1] - v0[1])))

        region_rings = [previous_ring]
        for layer in range(1, self.radial_meshes[-1] + 1):
            fraction = layer / self.radial_meshes[-1]
            ring = []
            for inner_node, outer_node in zip(previous_ring, outer_boundary):
                ix, iy = nodes[inner_node - 1]
                ox, oy = nodes[outer_node - 1]
                ring.append(add_node(ix + fraction * (ox - ix), iy + fraction * (oy - iy)))
            region_rings.append(ring)
        for ring_index in range(len(region_rings) - 1):
            previous, current = region_rings[ring_index], region_rings[ring_index + 1]
            for index in range(sectors):
                quads.append((len(self.material_names), (
                    previous[index], current[index], current[(index + 1) % sectors],
                    previous[(index + 1) % sectors]
                )))
        for index, node in enumerate(outer_boundary):
            boundary_lines.append((node, outer_boundary[(index + 1) % sectors]))

        globals()["_STRUCTURED_MESH"] = {
            "nodes": nodes,
            "tris": tris,
            "quads": quads,
            "boundary_lines": boundary_lines,
            "material_names": self.material_names,
        }

    def create_geometry(self) -> None:
        if not _ACTIVE_MODEL:
            raise RuntimeError("call Materials.create() before creating geometry")
        global _MESH_MODE
        _MESH_MODE = self.mesh_mode
        if self.mesh_mode == "structured":
            self._create_structured_mesh()
        else:
            self._create_automatic_geometry()

    def build(self, output_path: str | Path) -> Path:
        initialize_model("RectPin")
        try:
            self.create_geometry()
            return save_gmsh_file(output_path, "rect")
        except Exception:
            reset_model()
            raise
