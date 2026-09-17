import numpy as np
from RayTracing import MathFuncLib as mfl
from RayTracing import RunTime as rt
from RayTracing.MacroConfig import CASE_NAME, INIT, WARNING, ERROR


class Nodes:
    def __init__(self, node_id_set: np.ndarray, X_coords: np.ndarray, Y_coords: np.ndarray):
        self.__ids = node_id_set
        # self.__num_nodes = len(node_id_set)
        self.__num_nodes = np.max(node_id_set)
        self.__X_coords = X_coords
        self.__Y_coords = Y_coords
        self.__face_connectivity = [[] for _ in range(self.__num_nodes)]

        self.__boundary_node_ids = []

    def getCoords(self, node_id: int) -> tuple:
        node_idx = node_id - self.__ids[0]
        return self.__X_coords[node_idx], self.__Y_coords[node_idx]

    def getAllCoords(self) -> tuple:
        return self.__X_coords, self.__Y_coords

    def setBoundaryNodeIDs(self, node_ids: np.ndarray):
        self.__boundary_node_ids = node_ids

    def getID(self, node_idx: int):
        return self.__ids[node_idx]

    def getIDs(self, node_idx_set: np.ndarray):
        return self.__ids[node_idx_set]

    def getBoundaryIDs(self):
        return self.__boundary_node_ids

    def setFaceConnectivity(self, node_id: int, face_id: int):
        node_idx = node_id - self.__ids[0]
        if face_id not in self.__face_connectivity[node_idx]:
            self.__face_connectivity[node_idx].append(face_id)
        else:
            rt.printMessage(CASE_NAME, WARNING,
                            f"the face {face_id} has been added in nodes.FaceConnectivity, jump the operation")

    def getFaceConnectivity(self, node_id: int):
        node_idx = node_id - self.__ids[0]
        return self.__face_connectivity[node_idx]


class Edges:
    def __init__(self, edge_id_set: np.ndarray, node_connectivity: np.ndarray):
        self.__ids = edge_id_set
        self.__node_connectivity = node_connectivity
        self.__face_connectivity = []

        self.__arc_ids = []
        self.__arc_center_X_coords = []
        self.__arc_center_Y_coords = []

    def getNodeConnectivity(self, edge_id: int):
        edge_idx = edge_id - self.__ids[0]
        return self.__node_connectivity[edge_idx]

    def getEdgeIDs(self):
        return self.__ids

    def getNodesConnectivity(self):
        return self.__node_connectivity


class Faces:
    def __init__(self, face_id_set: np.ndarray, node_connectivity: np.ndarray, face_type: str):
        self.__ids = face_id_set
        self.__node_connectivity = node_connectivity
        self.__face_type = face_type
        self.__num_faces = len(face_id_set)
        self.__physical_names = ['unassigned' for _ in range(self.__num_faces)]
        self.__edge_connectivity = None
        self.__areas = None

        self.__boundary_faces = []

    def getNumFaces(self):
        return self.__num_faces

    def getNodeConnectivity(self, face_id: int):
        face_idx = face_id - self.__ids[0]
        return self.__node_connectivity[face_idx]

    def getAllNodeConnectivity(self):
        return self.__node_connectivity

    def getID(self, face_idx: int):
        return self.__ids[face_idx]

    def getIDs(self, face_idx_set: np.ndarray):
            return self.__ids[face_idx_set]

    def getPhysicalName(self, face_id: int):
        face_idx = face_id - self.__ids[0]
        return self.__physical_names[face_idx]

    def getBoundaryFaces(self):
        return self.__boundary_faces

    def setPhysicalName(self, face_id: int, physical_name: str):
        face_idx = face_id - self.__ids[0]
        self.__physical_names[face_idx] = physical_name

    def setBoundaryFaces(self, boundary_face_id_set: np.ndarray):
        self.__boundary_faces = boundary_face_id_set


class Mesh:
    def __init__(self, nodes: Nodes, edges: Edges, quads: Faces or None, tris: Faces or None, geom_type: str,
                 physical_groups: list):
        self.__nodes = nodes
        self.__edges = edges
        self.__quads = quads
        self.__tris = tris
        self.__geom_type = geom_type

        self.__bounding_box = None

        self.__physical_groups = physical_groups

        # set boundary flag
        self.__assignBoundaryFlagToNodes()
        self.__assignBoundaryFlagToFaces()

        # set face connectivity for nodes
        self.__assignFaceConnectivityForNodes()

        # set physical name for faces

    def __assignBoundaryFlagToNodes(self):
        # get bounding box
        x_coords, y_coords = self.__nodes.getAllCoords()
        self.__bounding_box = (x_coords.min(), x_coords.max(), y_coords.min(), y_coords.max())

        # find boundary node
        if self.__geom_type == 'hex':
            boundary_node_idx_set = mfl.getBoundaryPointsIndexInHex(x_coords, y_coords, self.__bounding_box)
        elif self.__geom_type == 'rect':
            boundary_node_idx_set = mfl.getBoundaryPointsIndexInRect(x_coords, y_coords, self.__bounding_box)
        else:
            raise ValueError('The geometry type must be rect or hex, plz check your geometry model!')

        boundary_node_id_set = self.__nodes.getIDs(boundary_node_idx_set)
        self.__nodes.setBoundaryNodeIDs(boundary_node_id_set)

    def __assignBoundaryFlagToFaces(self):
        boundary_node_ids = self.__nodes.getBoundaryIDs()

        if self.__tris:
            all_tri_node_ids = self.__tris.getAllNodeConnectivity()
        else:
            all_tri_node_ids = []
        if self.__quads:
            all_quad_node_ids = self.__quads.getAllNodeConnectivity()
        else:
            all_quad_node_ids = []

        boundary_tri_id_set = np.array([], dtype=np.int32)
        boundary_quad_id_set = np.array([], dtype=np.int32)

        for boundary_node_id in boundary_node_ids:
            if self.__tris:
                boundary_tri_idx_set, node_idx_in_tri = np.where(all_tri_node_ids == boundary_node_id)
                if boundary_tri_idx_set.shape[0] != 0:
                    boundary_tri_id_set = np.concatenate(
                        (boundary_tri_id_set, self.__tris.getIDs(boundary_tri_idx_set))
                    )

            if self.__quads:
                boundary_quad_idx_set, node_idx_in_quad = np.where(all_quad_node_ids == boundary_node_id)
                if boundary_quad_idx_set.shape[0] != 0:
                    boundary_quad_id_set = np.concatenate(
                        (boundary_quad_id_set, self.__quads.getIDs(boundary_quad_idx_set))
                    )

        if self.__tris:
            self.__tris.setBoundaryFaces(boundary_tri_id_set)
        if self.__quads:
            self.__quads.setBoundaryFaces(boundary_quad_id_set)

    def __assignFaceConnectivityForNodes(self):
        # quads then tris
        quad_node_id_set = self.__quads.getAllNodeConnectivity() if self.__quads else []

        for quad_idx in range(len(quad_node_id_set)):
            quad_id = self.__quads.getID(quad_idx)

            # get the relative node id of current quad
            node_id_0, node_id_1, node_id_2, node_id_3 = quad_node_id_set[quad_idx]

            self.__nodes.setFaceConnectivity(node_id_0, quad_id)
            self.__nodes.setFaceConnectivity(node_id_1, quad_id)
            self.__nodes.setFaceConnectivity(node_id_2, quad_id)
            self.__nodes.setFaceConnectivity(node_id_3, quad_id)

        tri_node_id_set = self.__tris.getAllNodeConnectivity() if self.__tris else []

        for tri_idx in range(len(tri_node_id_set)):
            tri_id = self.__tris.getID(tri_idx)

            # get the relative node id of current quad
            node_id_0, node_id_1, node_id_2 = tri_node_id_set[tri_idx]

            self.__nodes.setFaceConnectivity(node_id_0, tri_id)
            self.__nodes.setFaceConnectivity(node_id_1, tri_id)
            self.__nodes.setFaceConnectivity(node_id_2, tri_id)

    def getBoundaryGeomType(self) -> str:
        return self.__geom_type

    def getBoundingBox(self) -> tuple:
        return self.__bounding_box

    def getBoundaryQuads(self):
        return self.__quads.getBoundaryFaces() if self.__quads else np.array([], dtype=np.int32)

    def getBoundaryTris(self):
        return self.__tris.getBoundaryFaces() if self.__tris else np.array([], dtype=np.int32)

    def getPointLocationFlag(self, face_id: int, x: float, y: float) -> bool:
        vert_x, vert_y = self.getFaceNodeCoordinates(face_id)
        return mfl.isPointInPoly(len(vert_x), vert_x, vert_y, x, y)

    def getFaceNodeCoordinates(self, face_id) -> tuple:
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if num_quads and face_id <= self.__quads.getID(num_quads - 1):
            # is quad
            node1, node2, node3, node4 = self.__quads.getNodeConnectivity(face_id)
            x1, y1 = self.__nodes.getCoords(node1)
            x2, y2 = self.__nodes.getCoords(node2)
            x3, y3 = self.__nodes.getCoords(node3)
            x4, y4 = self.__nodes.getCoords(node4)
            return (x1, x2, x3, x4), (y1, y2, y3, y4)
        else:
            # is tri
            node1, node2, node3 = self.__tris.getNodeConnectivity(face_id)
            x1, y1 = self.__nodes.getCoords(node1)
            x2, y2 = self.__nodes.getCoords(node2)
            x3, y3 = self.__nodes.getCoords(node3)
            return (x1, x2, x3), (y1, y2, y3)

    def getAdjacentFace(self, face_id: int) -> tuple:
        adj_faces = []
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if num_quads and face_id <= self.__quads.getID(num_quads - 1):
            # is quad
            nodes = self.__quads.getNodeConnectivity(face_id)
        else:
            # is tri
            nodes = self.__tris.getNodeConnectivity(face_id)

        for node in nodes:
            adj_faces += self.__nodes.getFaceConnectivity(node)

        return tuple(set(adj_faces))

    def getNumFaces(self):
        if self.__tris:
            num_tris = self.__tris.getNumFaces()
        else:
            num_tris = 0

        if self.__quads:
            num_quads = self.__quads.getNumFaces()
        else:
            num_quads = 0

        return num_quads, num_tris

    def getPhysicalName(self, face_id: int):
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if num_quads and face_id <= self.__quads.getID(num_quads - 1):
            return self.__quads.getPhysicalName(face_id)
        else:
            return self.__tris.getPhysicalName(face_id)

    def getPhysicalGroups(self):
        return self.__physical_groups

    def getFaceID(self, face_idx: int):
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if face_idx < num_quads:
            return self.__quads.getID(face_idx)
        else:
            return self.__tris.getID(face_idx - num_quads)

    def getFaceIndex(self, face_id: int):
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if num_quads and face_id <= self.__quads.getID(num_quads - 1):
            return face_id - self.__quads.getID(0)
        if self.__tris and self.__tris.getNumFaces():
            return num_quads + face_id - self.__tris.getID(0)
        raise ValueError('face_id does not exist in this mesh')

    def getArea(self, face_id: int):
        vertx, verty = self.getFaceNodeCoordinates(face_id)
        if len(vertx) == 3:
            x1, x2, x3 = vertx
            y1, y2, y3 = verty

            # calc tri area
            area = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))

            return abs(area)

        else:
            x1, x2, x3, x4 = vertx
            y1, y2, y3, y4 = verty

            area1 = 0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
            area2 = 0.5 * (x1 * (y4 - y3) + x4 * (y3 - y1) + x3 * (y1 - y4))

            return abs(area2) + abs(area1)

    def getFaceType(self, face_id: int):
        num_quads = self.__quads.getNumFaces() if self.__quads else 0
        if num_quads and face_id <= self.__quads.getID(num_quads - 1):
            return 'quad'
        else:
            return 'tri'

    def getNodes(self):
        return self.__nodes

    def getFaces(self):
        return self.__quads, self.__tris

    def getAllEdge(self):
        return self.__edges
