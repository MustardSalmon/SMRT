from RayTracing.MeshBasicObj import Mesh
from RayTracing.MeshTreat import MeshSplitter
from RayTracing.Track import TrackFactory, SegmentationTally


class RayTracingModule:
    def __init__(self, mesh: Mesh, track_factory: TrackFactory, mesh_splitter: MeshSplitter = None):
        self.mesh = mesh
        self.track_factory = track_factory
        self.physical_id = 0
        self.rtm_id = 0
        self.mesh_splitter = mesh_splitter
        self.seg_tally = None

    def setPhysicalID(self, id):
        self.physical_id = id

    def setRTMID(self, id):
        self.rtm_id = id

    def getPhysicalID(self):
        return self.physical_id

    def getRTMID(self):
        return self.rtm_id

    def getSegmentationTally(self):
        if self.seg_tally == None:
            self.seg_tally = SegmentationTally(self.mesh, self.track_factory, mesh_spliter=self.mesh_splitter)
        return self.seg_tally


class RTMMap:
    def __init__(self, layout: list, width_Z: float, axial_mesh: int = 1):
        self.layout = layout
        self.dim_Y = len(layout)
        self.dim_X = len(layout[0])
        self.width_Z = width_Z
        self.axial_mesh = axial_mesh
