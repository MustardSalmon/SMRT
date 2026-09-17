"""Convert a Gmsh 2-D mesh into the internal ray-tracing mesh objects."""

from __future__ import annotations

from collections import OrderedDict
from pathlib import Path

import gmsh
import numpy as np

from RayTracing.MeshBasicObj import Edges, Faces, Mesh, Nodes


class GmshFileTranslator:
    """Read first-order Gmsh triangles/quadrilaterals and boundary lines."""

    def __init__(self, file_path: str | Path, geom_type: str = "hex"):
        self.file_path = str(Path(file_path).resolve())
        self.geom_type = geom_type

    @staticmethod
    def _physical_element_map(dim: int, skip_boundary: bool = True):
        element_to_group = {}
        group_names = OrderedDict()
        for group_dim, group_tag in gmsh.model.getPhysicalGroups(dim):
            name = gmsh.model.getPhysicalName(group_dim, group_tag)
            if not name or (skip_boundary and name == "_boundary") or "TALLY" in name:
                continue
            group_names[group_tag] = name
            for entity_tag in gmsh.model.getEntitiesForPhysicalGroup(dim, group_tag):
                element_types, element_tags, _ = gmsh.model.mesh.getElements(dim, entity_tag)
                for tags in element_tags:
                    for element_tag in tags:
                        element_to_group[int(element_tag)] = name
        return element_to_group, list(group_names.values())

    @staticmethod
    def _collect_elements(dim: int, wanted_types: set[int]):
        records = []
        for entity_dim, entity_tag in gmsh.model.getEntities(dim):
            element_types, element_tags, node_tags = gmsh.model.mesh.getElements(dim, entity_tag)
            for element_type, tags, flat_nodes in zip(element_types, element_tags, node_tags):
                if int(element_type) not in wanted_types:
                    continue
                _, _, _, num_nodes, _, _ = gmsh.model.mesh.getElementProperties(int(element_type))
                connectivity = np.asarray(flat_nodes, dtype=np.int64).reshape((-1, num_nodes))
                for element_tag, nodes in zip(tags, connectivity):
                    records.append((int(element_tag), tuple(int(node) for node in nodes)))
        return records

    def extractDataFromGmshFile(self) -> Mesh:
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            gmsh.open(self.file_path)

            element_to_material, physical_groups = self._physical_element_map(2)
            quad_records = self._collect_elements(2, {3})
            tri_records = self._collect_elements(2, {2})
            edge_groups, _ = self._physical_element_map(1, skip_boundary=False)
            edge_records = [
                record for record in self._collect_elements(1, {1})
                if edge_groups.get(record[0]) == "_boundary"
            ]

            if not quad_records and not tri_records:
                raise ValueError("Gmsh file contains no first-order triangle or quadrilateral elements")

            used_nodes = set()
            for _, nodes in quad_records + tri_records + edge_records:
                used_nodes.update(nodes)

            node_tags, coordinates, _ = gmsh.model.mesh.getNodes()
            coordinate_map = {
                int(tag): (float(coordinates[3 * idx]), float(coordinates[3 * idx + 1]))
                for idx, tag in enumerate(node_tags)
            }
            ordered_nodes = sorted(used_nodes)
            node_id_map = {old_id: new_id for new_id, old_id in enumerate(ordered_nodes, start=1)}
            node_ids = np.arange(1, len(ordered_nodes) + 1, dtype=np.int32)
            x_coords = np.array([coordinate_map[tag][0] for tag in ordered_nodes], dtype=np.float64)
            y_coords = np.array([coordinate_map[tag][1] for tag in ordered_nodes], dtype=np.float64)
            nodes = Nodes(node_ids, x_coords, y_coords)

            quad_ids = np.arange(1, len(quad_records) + 1, dtype=np.int32)
            tri_ids = np.arange(len(quad_records) + 1, len(quad_records) + len(tri_records) + 1, dtype=np.int32)
            quad_connectivity = np.array(
                [[node_id_map[node] for node in record[1]] for record in quad_records], dtype=np.int32
            ).reshape((-1, 4))
            tri_connectivity = np.array(
                [[node_id_map[node] for node in record[1]] for record in tri_records], dtype=np.int32
            ).reshape((-1, 3))

            quads = Faces(quad_ids, quad_connectivity, "quad")
            tris = Faces(tri_ids, tri_connectivity, "tri")
            for face_id, (old_tag, _) in zip(quad_ids, quad_records):
                quads.setPhysicalName(int(face_id), element_to_material.get(old_tag, "unassigned"))
            for face_id, (old_tag, _) in zip(tri_ids, tri_records):
                tris.setPhysicalName(int(face_id), element_to_material.get(old_tag, "unassigned"))

            edge_ids = np.arange(1, len(edge_records) + 1, dtype=np.int32)
            edge_connectivity = np.array(
                [[node_id_map[node] for node in record[1]] for record in edge_records], dtype=np.int32
            ).reshape((-1, 2))
            edges = Edges(edge_ids, edge_connectivity)

            mesh = Mesh(nodes, edges, quads, tris, self.geom_type, physical_groups)

            # Prefer the explicit Gmsh boundary physical group over the
            # bounding-box heuristic. This is important for both hexagon
            # orientations and for structured meshes generated from input cards.
            if edge_records:
                boundary_node_ids = np.unique(edge_connectivity.reshape(-1))
                nodes.setBoundaryNodeIDs(boundary_node_ids)
                boundary_node_set = set(int(node) for node in boundary_node_ids)
                def boundary_face_ids(face_set):
                    if face_set is None or face_set.getNumFaces() == 0:
                        return np.array([], dtype=np.int32)
                    boundary_ids = []
                    for face_index, face_id in enumerate(face_set.getIDs(np.arange(face_set.getNumFaces()))):
                        face_nodes = face_set.getNodeConnectivity(int(face_id))
                        # Keep all faces incident to an explicit boundary node.
                        # This matches the boundary convention and gives the
                        # ray tracer a robust entry-face candidate near corners.
                        if set(int(node) for node in face_nodes) & boundary_node_set:
                            boundary_ids.append(int(face_id))
                    return np.asarray(boundary_ids, dtype=np.int32)

                tris.setBoundaryFaces(boundary_face_ids(tris))
                quads.setBoundaryFaces(boundary_face_ids(quads))

            return mesh
        finally:
            gmsh.finalize()
