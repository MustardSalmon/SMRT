from RayTracing import MeshBasicObj
from RayTracing import RunTime as rt
from RayTracing import MathFuncLib as mfl
from RayTracing.MacroConfig import RUNNING, INIT, WARNING, ERROR, CASE_NAME
import numpy as np


class MeshSplitter:
    def __init__(self, origin_mesh: MeshBasicObj.Mesh):
        self.__origin_mesh = origin_mesh
        self.__face_ids = None
        self.__node_ids = None
        self.__bdry_face_ids = None
        self.__bdry_node_ids = None
        # Boundary edges are stored explicitly so tracking can intersect rays with
        # the real mesh outline instead of an ideal analytic bounding shape.
        self.__bdry_edges = None
        self.split_type = None

        # (xmin, xmax, ymin, ymax)
        self.__target_bounding_box = origin_mesh.getBoundingBox()

    def __getFaceNodeConnectivity(self, face_id: int) -> tuple:
        quads, tris = self.__origin_mesh.getFaces()
        if self.__origin_mesh.getFaceType(face_id) == 'quad':
            return tuple(int(node_id) for node_id in quads.getNodeConnectivity(face_id))
        else:
            return tuple(int(node_id) for node_id in tris.getNodeConnectivity(face_id))

    @staticmethod
    def __iterFaceEdges(node_ids: tuple):
        num_nodes = len(node_ids)
        for node_idx in range(num_nodes):
            edge = (node_ids[node_idx], node_ids[(node_idx + 1) % num_nodes])
            yield tuple(sorted(edge))

    def __getBoundaryByTopology(self, face_ids: np.ndarray) -> tuple:
        # A boundary edge of the selected split domain is an edge owned by only
        # one selected face. This avoids relying on geometric tolerances for
        # hand-built meshes whose boundary nodes are close to, but not exactly
        # on, the analytic symmetry lines.
        edge_owner_faces = {}
        for face_id in face_ids:
            face_id = int(face_id)
            node_ids = self.__getFaceNodeConnectivity(face_id)
            for edge in self.__iterFaceEdges(node_ids):
                edge_owner_faces.setdefault(edge, []).append(face_id)

        boundary_face_ids = set()
        boundary_node_ids = set()
        boundary_edges = []
        for edge, owner_faces in edge_owner_faces.items():
            if len(owner_faces) == 1:
                boundary_face_ids.add(owner_faces[0])
                boundary_node_ids.update(edge)
                boundary_edges.append(edge)

        return (np.array(sorted(boundary_face_ids), dtype=face_ids.dtype),
                np.array(sorted(boundary_node_ids), dtype=np.int32),
                tuple(sorted(boundary_edges)))

    def setTargetBoundingBox(self, coords: tuple = None, node_ids: tuple = None):
        if coords:
            self.__target_bounding_box = coords
        elif node_ids:
            all_nodes = self.__origin_mesh.getNodes()
            xs = []
            ys = []
            for node_id in node_ids:
                x, y = all_nodes.getCoords(node_id)
                xs.append(x)
                ys.append(y)
            self.__target_bounding_box = (min(xs), max(xs), min(ys), max(ys))
        else:
            rt.printMessage(CASE_NAME, ERROR, 'coords or nodes must be given')
            raise RuntimeError

    def getTargetBoundingBox(self):
        return self.__target_bounding_box

    def getFaceIDs(self):
        return self.__face_ids

    def getFaceID(self, face_idx):
        return self.__face_ids[face_idx]

    def getNodeIDs(self):
        return self.__node_ids

    def getBoundaryFaceIDs(self):
        return self.__bdry_face_ids

    def getBoundaryNodeIDs(self):
        return self.__bdry_node_ids

    def getBoundaryEdges(self):
        return self.__bdry_edges

    def getNumFaces(self):
        return len(self.__face_ids)

    def getFaceIndex(self, face_id):
        return np.where(self.__face_ids == face_id)[0][0]

    def isFaceExisted(self, face_id: int) -> bool:
        if face_id in self.getFaceIDs():
            return True
        else:
            return False

    def split(self, split_type='1/6'):
        self.split_type = split_type
        if self.__origin_mesh.getBoundaryGeomType() == 'hex':
            rt.printMessage(CASE_NAME, INIT, 'Splitting mesh by %s ...' % (split_type))

            # find the nodes
            nodes = self.__origin_mesh.getNodes()
            x_coords, y_coords = nodes.getAllCoords()
            sect_node_idx_set, sect_boundary_node_idx_set = [], []
            if split_type == '1/6':
                sect_node_idx_set, sect_boundary_node_idx_set = mfl.getPointsIndexInTriSectOfHex(x_coords, y_coords,
                                                                                                 self.__target_bounding_box)
            elif split_type == 'full':
                sect_node_idx_set, sect_boundary_node_idx_set = mfl.getPointsIndexInHex(x_coords, y_coords,
                                                                                        self.__target_bounding_box)
            else:
                rt.printMessage(CASE_NAME, ERROR, 'the given split_type should be "1/6" or "full"')

            sect_node_id_set = nodes.getIDs(sect_node_idx_set)
            sect_boundary_node_id_set = nodes.getIDs(sect_boundary_node_idx_set)

            # find the faces
            quads, tris = self.__origin_mesh.getFaces()
            sect_face_id_set = []
            for node_id in sect_node_id_set:
                node_faces = nodes.getFaceConnectivity(node_id)
                for face_id in node_faces:
                    is_in = True
                    if self.__origin_mesh.getFaceType(face_id) == 'quad':
                        face_nodes = quads.getNodeConnectivity(face_id)
                    else:
                        face_nodes = tris.getNodeConnectivity(face_id)
                    for node_id1 in face_nodes:
                        if node_id1 not in sect_node_id_set:
                            is_in = False
                            break
                    if is_in:
                        sect_face_id_set.append(face_id)

            sect_face_id_set = np.array(sect_face_id_set)
            sect_face_id_set = np.unique(sect_face_id_set)

            # Build boundary faces/nodes from face-edge topology. The previous
            # node-on-analytic-boundary test can miss true exterior faces when
            # Imported mesh coordinates can have small boundary offsets.
            sect_boundary_face_id_set, sect_boundary_node_id_set, sect_boundary_edges = \
                self.__getBoundaryByTopology(sect_face_id_set)

            self.__face_ids = sect_face_id_set
            self.__bdry_face_ids = sect_boundary_face_id_set
            self.__node_ids = sect_node_id_set
            self.__bdry_node_ids = sect_boundary_node_id_set
            self.__bdry_edges = sect_boundary_edges

class MeshPartitioner:
    def __init__(self, mesh: MeshBasicObj.Mesh, num_part: int = 4, mesh_splitter: MeshSplitter = None):
        self.__mesh = mesh
        self.__mesh_splitter = mesh_splitter
        self.__num_part = num_part
        self.__part_face_ids = None
        self.__part_node_ids = None
        self.__part_bdry_face_ids = None
        self.__part_bdry_node_ids = None

    def part(self):
        import pymetis
        if self.__mesh_splitter:
            num_faces = self.__mesh_splitter.getNumFaces()
        else:
            num_quads, num_tris = self.__mesh.getNumFaces()
            num_faces = num_quads + num_tris

        # collect the adjacency, pop the self!
        adjacency = []
        for face_idx in range(num_faces):
            adj_faces = []
            if self.__mesh_splitter:
                face_id = self.__mesh_splitter.getFaceID(face_idx)
                _adj_faces = list(self.__mesh.getAdjacentFace(face_id))
                for adj_face_id in _adj_faces:
                    if adj_face_id != face_id and self.__mesh_splitter.isFaceExisted(adj_face_id):
                        adj_faces.append(self.__mesh_splitter.getFaceIndex(adj_face_id))
                    else:
                        continue
            else:
                face_id = self.__mesh.getFaceID(face_idx)
                _adj_faces = list(self.__mesh.getAdjacentFace(face_id))
                for adj_face_id in _adj_faces:
                    if adj_face_id == face_id:
                        continue
                    else:
                        adj_faces.append(self.__mesh.getFaceIndex(adj_face_id))

            adjacency.append(adj_faces)

        (edgecuts, parts) = pymetis.part_graph(nparts=self.__num_part, adjacency=adjacency)

        self.__part_face_ids = []
        for n in range(self.__num_part):
            coll = []
            for i in range(len(parts)):
                if parts[i] == n:
                    if self.__mesh_splitter:
                        face_id = self.__mesh_splitter.getFaceID(i)
                    else:
                        face_id = self.__mesh.getFaceID(i)
                    coll.append(face_id)
            self.__part_face_ids.append(coll)


        tf = open('test.txt', 'w')
        tf.write('delete block all\n')
        num_quads, num_tris = self.__mesh.getNumFaces()

        for n in range(self.__num_part):
            coll = []
            for i in range(len(parts)):
                if parts[i] == n:
                    if self.__mesh_splitter:
                        face_id = self.__mesh_splitter.getFaceID(i)
                    else:
                        face_id = self.__mesh.getFaceID(i)
                    coll.append(face_id)

            print(len(coll))

            tf.write('block %d add face ' % (n + 1))
            for f in coll:
                if f <= num_quads:
                    tf.write(' %d ' % f)
            tf.write('\n')

            tf.write('block %d add tri ' % (n + 1))
            for f in coll:
                if f > num_quads:
                    tf.write(' %d ' % f)
            tf.write('\n')

        tf.write('draw block all\n')
        tf.close()
