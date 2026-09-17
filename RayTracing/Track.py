import os
import struct

import numpy as np
# import MathFuncLib as mfl
from RayTracing import MathFuncLib as mfl
from RayTracing import RunTime as rt
from RayTracing.MacroConfig import CASE_NAME, INIT, WARNING, ERROR, RUNNING, SQRT_3, SMALL_VALUE, BIG_VALUE
from RayTracing.MeshBasicObj import Mesh
from RayTracing.MeshTreat import MeshSplitter
from RayTracing.FunctionalFuncLib import Line, OutputPassBottomRayIDFileHandle
import math
import time
import multiprocessing as mp
from matplotlib import pyplot as plt

MAX_OFFSET_GROW = 1000  # adaptive offset iter


class AngularTrackWrap:
    def __init__(self, mesh: Mesh, spacing: float, azim_ang: float, num_tracks: int,
                 mesh_splitter: MeshSplitter = None, isOutputStatus: bool = False):
        self.__mesh = mesh
        self.__mesh_splitter = mesh_splitter
        self.__azim_ang = azim_ang
        self.__spacing = spacing
        self.__num_tracks = num_tracks
        self.__starts = [(0., 0.) for _ in range(num_tracks)]
        self.__ends = [(0., 0.) for _ in range(num_tracks)]
        self.__track_length = [0. for _ in range(num_tracks)]
        self.__seg_length = [[] for _ in range(num_tracks)]
        self.__swept_fsr = [[] for _ in range(num_tracks)]
        self.__num_seg = [0 for _ in range(num_tracks)]

        self.__direction = (math.cos(azim_ang), math.sin(azim_ang))
        self.__k = math.sin(azim_ang) / math.cos(azim_ang)
        if self.__k > BIG_VALUE:
            self.__k = math.inf

        self.__isOutputStatus = isOutputStatus
        if self.__mesh_splitter and self.__mesh_splitter.split_type == '1/6' and isOutputStatus:
            self.__exitRayID = {'incident': [], 'emergent': [], 'totalTracks': self.__num_tracks}

    def setStart(self, track_idx: int, start: tuple):
        if track_idx >= self.__num_tracks:
            rt.printMessage(CASE_NAME, ERROR, "the given track_idx larger than number of track")
            raise ValueError
        else:
            self.__starts[track_idx] = start

    def setEnd(self, track_idx: int, end: tuple):
        if track_idx >= self.__num_tracks:
            rt.printMessage(CASE_NAME, ERROR, "the given track_idx larger than number of track")
            raise ValueError
        else:
            self.__ends[track_idx] = end

    def setTrackLen(self, track_idx: int, track_len: float):
        if track_idx >= self.__num_tracks:
            rt.printMessage(CASE_NAME, ERROR, "the given track_idx larger than number of track")
            raise ValueError
        else:
            self.__track_length[track_idx] = track_len

    def getStart(self, track_idx: int):
        return self.__starts[track_idx]

    def getEnd(self, track_idx: int):
        return self.__ends[track_idx]

    def getTrackLen(self, track_idx: int):
        return self.__track_length[track_idx]

    def getNumSeg(self):
        return sum(self.__num_seg)

    def getSegLength(self, track_idx):
        return self.__seg_length[track_idx]

    def getSweptFace(self, track_idx):
        return self.__swept_fsr[track_idx]

    def getExitRay(self):
        """Get the ID of the line that passes through the bottom edge"""
        if len(self.__exitRayID['emergent']) != 0:  # only incident ray have value in [0, pi]
            raise ValueError
        return self.__exitRayID['totalTracks'], self.__exitRayID['incident']

    def tracking(self):
        for t in range(self.__num_tracks):
            self.__trackingForOneTrack(t)

    def __findEnrtyAndExitFaceID(self, track_idx: int) -> tuple:
        # calc the boundary hit point
        # update the track geometric param
        incident, emergent = self.__calcBoundaryHit(track_idx)

        # todo: calc 1/6 geometry needs and ray id through vacuum boundary
        if self.__mesh_splitter and self.__mesh_splitter.split_type == '1/6' and self.__isOutputStatus:
            # # get all ray id
            # self.__exitRayID['incident'].append(Line(track_idx, 0, incident, emergent))

            # get y_min for this geometry
            x_min, x_max, y_min, y_max = self.__mesh_splitter.getTargetBoundingBox()
            # Use a tolerance here because boundary hits are computed from mesh
            # edges and may not be bitwise equal to the analytic y_min value.
            if abs(incident[-1] - y_min) < SMALL_VALUE:
                self.__exitRayID['incident'].append(Line(track_idx, 0, incident, emergent))
            elif abs(emergent[-1] - y_min) < SMALL_VALUE:
                self.__exitRayID['emergent'].append(Line(track_idx, 0, incident, emergent))

        track_len = mfl.calc2DVectorMode(incident[0], incident[1], emergent[0], emergent[1])
        self.setTrackLen(track_idx, track_len)

        # for calc the intersection easier on the entry and exit mesh
        incident_offset_X = incident[0] - self.__direction[0] * SMALL_VALUE
        incident_offset_Y = incident[1] - self.__direction[1] * SMALL_VALUE
        emergent_offset_X = emergent[0] + self.__direction[0] * SMALL_VALUE
        emergent_offset_Y = emergent[1] + self.__direction[1] * SMALL_VALUE
        self.setStart(track_idx, (incident_offset_X, incident_offset_Y))
        self.setEnd(track_idx, (emergent_offset_X, emergent_offset_Y))

        # find the entry FSR
        incident_offset_X = incident[0]
        incident_offset_Y = incident[1]
        emergent_offset_X = emergent[0]
        emergent_offset_Y = emergent[1]

        entry_face = None
        exit_face = None

        if self.__mesh_splitter:
            boundary_faces = self.__mesh_splitter.getBoundaryFaceIDs()
        else:
            boundary_quads = self.__mesh.getBoundaryQuads()
            boundary_tris = self.__mesh.getBoundaryTris()
            boundary_faces = list(boundary_quads) + list(boundary_tris)

        for _i in range(MAX_OFFSET_GROW):
            incident_offset_X += self.__direction[0] * 5e-6
            incident_offset_Y += self.__direction[1] * 5e-6
            for face_id in boundary_faces:
                if self.__mesh.getPointLocationFlag(face_id, incident_offset_X, incident_offset_Y):
                    entry_face = face_id
                    break

            if entry_face:
                break

        for _i in range(MAX_OFFSET_GROW):
            emergent_offset_X -= self.__direction[0] * 5e-6
            emergent_offset_Y -= self.__direction[1] * 5e-6
            for face_id in boundary_faces:
                if self.__mesh.getPointLocationFlag(face_id, emergent_offset_X, emergent_offset_Y):
                    exit_face = face_id
                    break

            if exit_face:
                break

        if not entry_face:
            rt.printMessage(CASE_NAME, ERROR, "Can not find the entry FSR for track %d" % (track_idx))
            raise RuntimeError("Can not find the entry FSR for track %d" % (track_idx))

        if not exit_face:
            rt.printMessage(CASE_NAME, ERROR, "Can not find the exit FSR for track %d" % (track_idx))
            raise RuntimeError("Can not find the exit FSR for track %d" % (track_idx))

        return entry_face, exit_face

    def __calcSegment(self, track_idx: int, face_id: int) -> tuple:
        """
        根据特征线ID & 特征线所在平源区, 计算特征线段长
        :param track_idx: 特征线ID
        :param face_id: 特征线穿过平源区
        :return: 特征线长度 & 下一个平源区特征线入口点的坐标
        :bugReason: 1. 设置zero判定条件过于宽松, 当有线段有长度且长度小于1e-6时, 会自动置零, 因此会出现线与面只有一个交点
                    2. 当修改zero判定条件时(1e-6 -> 1e-8), 出现特征线经过角点(平源区两条边的角点)的情况, 会计算出线段长(角点的两条边), 此时会自动退出__calcSegment,
                       在 next fsr 查找时, fsr与特征线出口的交点有关, 由于返回了错误的出口交点, 所以会出现无法找到 next fsr 的问题
                    上述两个问题相互影响, 但是由于python float的精度只有 1e-6, 因此不能将zero判定条件修改, 采用其他修改方法
        :modify: 当计算的线段长度小于1e-6时, 计算特征线与平源区所有边的交点, 将线段长和交点分别存储至_length, intersectX, intersectY
        等变量, 通过找到两两交点之间线段长的最大值, 此最大值就是真正的线段长, 其中的两个点即为特征线的入口和出口点
        """
        p0 = self.getStart(track_idx)
        p1 = self.getEnd(track_idx)

        x0, y0, x1, y1 = 0., 0., 0., 0.

        # ZouHang Writing
        # vert_x, vert_y = self.__mesh.getFaceNodeCoordinates(face_id)
        # hit_cnt, num_vert,  = 0, len(vert_x)
        #
        # for vert in range(num_vert):
        #     # find the edge node for intersecting
        #     q0 = (vert_x[vert], vert_y[vert])
        #     if vert == num_vert - 1:
        #         q1 = (vert_x[0], vert_y[0])
        #     else:
        #         q1 = (vert_x[vert + 1], vert_y[vert + 1])
        #
        #     # calc the intersection
        #     status, xi, yi = mfl.calcLineIntersection(p0, p1, q0, q1)
        #     if status == 'intersect':
        #         if hit_cnt == 0:
        #             x0, y0 = xi, yi
        #             hit_cnt += 1
        #         elif math.sqrt((yi - y0) ** 2 + (xi - x0) ** 2) > 1e-6:
        #             hit_cnt += 1
        #             x1, y1 = xi, yi
        #             break
        #
        #         if hit_cnt > 1:
        #             rt.printMessage(CASE_NAME, ERROR,
        #                             'find more than 2 intersection for track %d and face %d, plz check!' % (
        #                                 track_idx, face_id))
        #             raise RuntimeError

        # CaoWei Modifying
        vert_x, vert_y = self.__mesh.getFaceNodeCoordinates(face_id)
        hit_cnt, num_vert, isDistance = 0, len(vert_x), 'normal'

        for vert in range(num_vert):
            # find the edge node for intersecting
            q0 = (vert_x[vert], vert_y[vert])
            if vert == num_vert - 1:
                q1 = (vert_x[0], vert_y[0])
            else:
                q1 = (vert_x[vert + 1], vert_y[vert + 1])

            # calc the intersection
            status, xi, yi = mfl.calcLineIntersection(p0, p1, q0, q1)
            if status == 'intersect':
                if hit_cnt == 0:
                    x0, y0 = xi, yi
                    hit_cnt += 1
                elif math.sqrt((yi - y0) ** 2 + (xi - x0) ** 2) > 1e-6:
                    hit_cnt += 1
                    x1, y1 = xi, yi
                    break
                else:
                    isDistance = 'too short'
                    break

        if isDistance == 'too short':
            _length, intersectX, intersectY = {0.0: ()}, [], []
            for vert in range(num_vert):
                q0 = (vert_x[vert], vert_y[vert])
                if vert == num_vert - 1:
                    q1 = (vert_x[0], vert_y[0])
                else:
                    q1 = (vert_x[vert + 1], vert_y[vert + 1])

                status, xi, yi = mfl.calcLineIntersection(p0, p1, q0, q1)
                if status == 'intersect':
                    intersectX.append(xi)
                    intersectY.append(yi)

            for _i in range(len(intersectX)):
                for _j in range(_i+1, len(intersectX)):
                    _length[mfl.calc2DVectorMode(intersectX[_i], intersectY[_i], intersectX[_j], intersectY[_j])] = (_i, _j)

            # Look for the segment with the largest length
            _lengthMax = sorted(_length.keys(), reverse=True)
            p1ID, p2ID = _length[_lengthMax[0]][0], _length[_lengthMax[0]][1]

            # Determine the point of incidence and emergence of the line segment
            if abs(intersectX[p1ID] - x0) < 1e-6 and abs(intersectY[p1ID] - y0) < 1e-6:
                x1, y1 = intersectX[p2ID], intersectY[p2ID]
            else:
                x1, y1 = intersectX[p1ID], intersectY[p1ID]
            hit_cnt += 1

        if hit_cnt == 2:
            if y1 > y0:
                return mfl.calc2DVectorMode(x0, y0, x1, y1), (x1, y1)
            else:
                return mfl.calc2DVectorMode(x0, y0, x1, y1), (x0, y0)
        else:
            rt.printMessage(CASE_NAME, ERROR,
                            'failed to calc intersection for track %d and face %d, plz check!' % (
                                track_idx, face_id))
            rt.printMessage(CASE_NAME, ERROR, 'AZIM VALUE : %f' % (self.__azim_ang * 180. / math.pi))
            rt.printMessage(CASE_NAME, ERROR, 'TRACK INDEX: %d' % track_idx)
            rt.printMessage(CASE_NAME, ERROR, 'SWEPT FACE : %d' % face_id)
            rt.printMessage(CASE_NAME, ERROR, 'TRACK START: (%f, %f)' % (p0[0], p0[1]))
            rt.printMessage(CASE_NAME, ERROR, 'TRACK END  : (%f, %f)' % (p1[0], p1[1]))
            rt.printMessage(CASE_NAME, ERROR, 'HIT COUNT  : %d' % hit_cnt)
            raise RuntimeError

    def __findNextFaceID(self, track_id: int, incident: tuple, src_face: int) -> int:
        incident_offset_X = incident[0]
        incident_offset_Y = incident[1]
        for _i in range(MAX_OFFSET_GROW):
            incident_offset_X += self.__direction[0] * 1e-6
            incident_offset_Y += self.__direction[1] * 1e-6

            adj_faces = self.__mesh.getAdjacentFace(src_face)

            for face in adj_faces:
                if face == src_face:
                    continue
                else:
                    if self.__mesh.getPointLocationFlag(face, incident_offset_X, incident_offset_Y):
                        return face

            for face_idx in range(sum(self.__mesh.getNumFaces())):
                face = self.__mesh.getFaceID(face_idx)
                if face == src_face:
                    continue
                else:
                    if self.__mesh.getPointLocationFlag(face, incident_offset_X, incident_offset_Y):
                        return face

        wrnmsg = '\nAZIM VALUE : %f \n' % (self.__azim_ang * 180. / math.pi)
        wrnmsg += 'TRACK INDEX: %d \n' % track_id
        wrnmsg += 'CURR FACE  : %d \n' % src_face
        wrnmsg += 'INCIDENT   : (%f, %f) \n' % (incident[0], incident[1])
        rt.printMessage(CASE_NAME, WARNING, wrnmsg)
        rt.printMessage(CASE_NAME, ERROR, 'Can not find next FSR for track %d, plz check!' % (track_id))

    def __trackingForOneTrack(self, track_idx: int):
        entry, _exit = self.__findEnrtyAndExitFaceID(track_idx)

        # start calculating segmentation
        curr_face = entry
        is_hit_exit = False

        while not is_hit_exit:
            # calc seg len of current FSR
            seg_length, (x_in, y_in) = self.__calcSegment(track_idx, curr_face)
            self.__seg_length[track_idx].append(seg_length)
            self.__swept_fsr[track_idx].append(curr_face)
            self.__num_seg[track_idx] += 1

            if curr_face == _exit:
                is_hit_exit = True
            else:
                # find next FSR
                curr_face = self.__findNextFaceID(track_idx, (x_in, y_in), curr_face)

        # length check
        tracking_length = sum(self.__seg_length[track_idx])
        reference_length = self.getTrackLen(track_idx)
        bias = tracking_length - reference_length
        if abs(bias) > 1e-5:
            msg = 'The sum of seg_length is not equal to track_length in track %d, bias %.8f' % (track_idx, bias)
            rt.printMessage(CASE_NAME, WARNING, msg)

        # recover the start and end
        start = self.getStart(track_idx)
        end = self.getEnd(track_idx)
        self.setStart(track_idx,
                      (start[0] + self.__direction[0] * SMALL_VALUE, start[1] + self.__direction[1] * SMALL_VALUE))
        self.setEnd(track_idx, (end[0] - self.__direction[0] * SMALL_VALUE, end[1] - self.__direction[1] * SMALL_VALUE))

    def __calcBoundaryHit(self, track_idx: int) -> tuple:

        def __calcSplitterBoundaryHit(track_idx: int) -> tuple:
            # When a mesh splitter is active, the transport domain is the
            # selected mesh outline. Intersecting with those real boundary edges
            # keeps track lengths consistent with hand-built or imported meshes.
            boundary_edges = self.__mesh_splitter.getBoundaryEdges()
            if not boundary_edges:
                rt.printMessage(CASE_NAME, ERROR, 'can not find boundary edges in mesh splitter, plz check!')
                raise RuntimeError

            def __sub(p0: tuple, p1: tuple):
                return p1[0] - p0[0], p1[1] - p0[1]

            def __appendUnique(points: list, point: tuple):
                for old_point in points:
                    if mfl.calc2DVectorMode(old_point[0], old_point[1], point[0], point[1]) < 1e-7:
                        return
                points.append(point)

            def __lineBoundaryEdgeIntersection(p0: tuple, direction: tuple, q0: tuple, q1: tuple) -> list:
                # Treat the ray as an infinite line and only constrain the
                # boundary-edge parameter. Very shallow tracks can start far
                # outside the mesh, so a finite artificial ray segment can miss
                # the exit boundary.
                r = direction
                s = __sub(q0, q1)
                q2p = __sub(p0, q0)
                rxs = mfl.vectorCross2D(r, s)
                if abs(rxs) < 1e-12:
                    # If the track lies on a boundary edge, keep both edge
                    # endpoints; sorting by projection below selects the domain
                    # span and unique filtering removes duplicated vertices.
                    if abs(mfl.vectorCross2D(q2p, r)) < 1e-10:
                        return [q0, q1]
                    return []

                t = mfl.vectorCross2D(q2p, s) / rxs
                u = mfl.vectorCross2D(q2p, r) / rxs
                if -1e-8 <= u <= 1.0 + 1e-8:
                    return [(p0[0] + t * r[0], p0[1] + t * r[1])]
                return []

            p0 = self.__starts[track_idx]

            nodes = self.__mesh.getNodes()
            hit_points = []
            for node_id_0, node_id_1 in boundary_edges:
                q0 = nodes.getCoords(node_id_0)
                q1 = nodes.getCoords(node_id_1)
                for hit_point in __lineBoundaryEdgeIntersection(p0, self.__direction, q0, q1):
                    __appendUnique(hit_points, hit_point)

            if len(hit_points) < 2:
                rt.printMessage(CASE_NAME, ERROR,
                                'number of boundary hit points less than 2 for track %d, plz check!' % track_idx)
                raise RuntimeError

            def __projection(point: tuple):
                return (point[0] - p0[0]) * self.__direction[0] + (point[1] - p0[1]) * self.__direction[1]

            hit_points = sorted(hit_points, key=__projection)
            return hit_points[0], hit_points[-1]

        def __calcRectHit(track_idx: int) -> tuple:
            #           p1
            #      p2 /
            # +------/-+
            # |     /  |
            # |    /   |
            # +---/----+
            #    / p0
            # p3
            x_min, x_max, y_min, y_max = self.__mesh.getBoundingBox()
            x0, y0 = self.__starts[track_idx]

            k = self.__k
            b = mfl.calcIntercept(x0, y0, self.__k)

            p0 = ((y_min - b) / k, y_min)
            p1 = (x_max, x_max * k + b)
            p2 = ((y_max - b) / k, y_max)
            p3 = (x_min, x_min * k + b)

            hit_points = []

            for p_coord in [p0, p1, p2, p3]:
                if p_coord[0] < x_min or p_coord[0] > x_max:
                    continue
                elif p_coord[1] < y_min or p_coord[1] > y_max:
                    continue
                else:
                    hit_points.append(p_coord)

            if len(hit_points) > 2:
                rt.printMessage(CASE_NAME, ERROR, 'number of boundary hit points more than 2, plz check')
                raise RuntimeError

            hit_x_coord = [hit_point[0] for hit_point in hit_points]
            if self.__direction[0] > 0:
                incident_idx = hit_x_coord.index(min(hit_x_coord))
                emergent_idx = hit_x_coord.index(max(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            elif self.__direction[0] < 0:
                incident_idx = hit_x_coord.index(max(hit_x_coord))
                emergent_idx = hit_x_coord.index(min(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            else:
                rt.printMessage(CASE_NAME, ERROR, 'direction of track ' + str(track_idx) + ' is error, plz check')
                raise ValueError

            return incident_point, emergent_point

        def __calcHexHit(track_idx: int) -> tuple:
            if self.__mesh_splitter and self.__mesh_splitter.split_type == 'full':
                x_min, x_max, y_min, y_max = self.__mesh_splitter.getTargetBoundingBox()
            else:
                x_min, x_max, y_min, y_max = self.__mesh.getBoundingBox()

            x0, y0 = self.__starts[track_idx]
            b = mfl.calcIntercept(x0, y0, self.__k)

            hex_center = [(x_max + x_min) / 2, (y_min + y_max) / 2]
            hex_edge_width = (x_max - x_min) / 2

            #
            #                /
            #           /---/---------\
            #          /   /           \
            #         /   /             \
            #        /   /               \
            #        \  /                /
            #         \/                /
            #         /\               /
            #           \-------------/
            #
            # v0 = [hex_center[0] + hex_edge_width / 2, y_min]
            # v1 = [x_max, hex_center[1]]
            # v2 = [hex_center[0] + hex_edge_width / 2, y_max]
            # v3 = [hex_center[0] - hex_edge_width / 2, y_max]
            # v4 = [x_min, hex_center[1]]
            # v5 = [hex_center[0] - hex_edge_width / 2, y_min]

            k = self.__k

            if math.isinf(self.__k):
                # todo: fix the bug of num_azim = 6 in hex ray tracing
                pass
            else:
                p0 = ((y_min - b) / k, y_min)

                x0 = (hex_center[1] - SQRT_3 * x_max - b) / (k - SQRT_3)
                y0 = k * x0 + b
                p1 = (x0, y0)

                x0 = (hex_center[1] + SQRT_3 * x_max - b) / (k + SQRT_3)
                y0 = k * x0 + b
                p2 = (x0, y0)

                p3 = ((y_max - b) / k, y_max)

                x0 = (hex_center[1] - SQRT_3 * x_min - b) / (k - SQRT_3)
                y0 = k * x0 + b
                p4 = (x0, y0)

                x0 = (hex_center[1] + SQRT_3 * x_min - b) / (k + SQRT_3)
                y0 = k * x0 + b
                p5 = (x0, y0)

            x_lf = hex_center[0] - hex_edge_width / 2
            x_lr = hex_center[0] + hex_edge_width / 2

            hit_points = []

            for p_coord in [p0, p1, p2, p3, p4, p5]:
                if p_coord[0] < x_min - SMALL_VALUE or p_coord[0] > x_max + SMALL_VALUE:
                    continue
                elif p_coord[1] < y_min - SMALL_VALUE or p_coord[1] > y_max + SMALL_VALUE:
                    continue
                elif abs(p_coord[1] - y_max) < SMALL_VALUE or abs(p_coord[1] - y_min) < SMALL_VALUE:
                    if p_coord[0] < x_lf - SMALL_VALUE or p_coord[0] > x_lr + SMALL_VALUE:
                        continue
                    else:
                        hit_points.append(p_coord)
                else:
                    hit_points.append(p_coord)

            if len(hit_points) > 2:
                rt.printMessage(CASE_NAME, ERROR, 'number of boundary hit points more than 2, plz check!')

            incident_point, emergent_point = 0, 0
            hit_x_coord = [hit_point[0] for hit_point in hit_points]
            if self.__direction[0] > 0:
                incident_idx = hit_x_coord.index(min(hit_x_coord))
                emergent_idx = hit_x_coord.index(max(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            elif self.__direction[0] < 0:
                incident_idx = hit_x_coord.index(max(hit_x_coord))
                emergent_idx = hit_x_coord.index(min(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            else:
                rt.printMessage(CASE_NAME, ERROR, 'direction of ray ' + str(track_idx) + ' is error, plz check')

            return incident_point, emergent_point

        def __calcTriHit(track_idx: int) -> tuple:
            # triangle on the bottom of hex
            if self.__mesh_splitter and self.__mesh_splitter.split_type == '1/6':
                x_min, x_max, y_min, y_max = self.__mesh_splitter.getTargetBoundingBox()
            else:
                x_min, x_max, y_min, y_max = self.__mesh.getBoundingBox()
            x0, y0 = self.__starts[track_idx]
            b = mfl.calcIntercept(x0, y0, self.__k)

            hex_center = [(x_max + x_min) / 2, (y_min + y_max) / 2]
            hex_edge_width = (x_max - x_min) / 2

            #              v2
            #              /\
            #             /  \/
            #            /   /\
            #           /   /  \
            #          /___/____\
            #         v0  /     v1
            #
            # v0 = [hex_center[0] - hex_edge_width / 2, y_min]
            # v1 = [hex_center[0] + hex_edge_width / 2, y_min]
            # v2 = [hex_center[0], hex_center[1]]

            k = self.__k

            if math.isinf(self.__k):
                # todo: fix the bug of num_azim = 6 in hex ray tracing
                pass
            else:
                p0 = ((y_min - b) / k, y_min)  # v0 - v1

                # v1 - v2
                x0 = (hex_center[1] + SQRT_3 * hex_center[0] - b) / (k + SQRT_3)
                y0 = k * x0 + b
                p1 = (x0, y0)

                # v2 - v0
                x0 = (hex_center[1] - SQRT_3 * hex_center[0] - b) / (k - SQRT_3)
                y0 = k * x0 + b
                p2 = (x0, y0)

            x_lf = hex_center[0] - hex_edge_width / 2
            x_lr = hex_center[0] + hex_edge_width / 2

            hit_points = []

            for p_coord in [p0, p1, p2]:
                if p_coord[0] < x_lf - SMALL_VALUE or p_coord[0] > x_lr + SMALL_VALUE:
                    continue
                elif p_coord[1] < y_min - SMALL_VALUE or p_coord[1] > hex_center[1] + SMALL_VALUE:
                    continue
                else:
                    hit_points.append(p_coord)

            if len(hit_points) > 2:
                rt.printMessage(CASE_NAME, ERROR, 'number of boundary hit points more than 2, plz check!')
                raise RuntimeError

            hit_x_coord = [hit_point[0] for hit_point in hit_points]
            if self.__direction[0] > 0:
                incident_idx = hit_x_coord.index(min(hit_x_coord))
                emergent_idx = hit_x_coord.index(max(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            elif self.__direction[0] < 0:
                incident_idx = hit_x_coord.index(max(hit_x_coord))
                emergent_idx = hit_x_coord.index(min(hit_x_coord))
                incident_point = hit_points[incident_idx]
                emergent_point = hit_points[emergent_idx]

            else:
                rt.printMessage(CASE_NAME, ERROR, 'direction of ray ' + str(track_idx) + ' is error, plz check')
                raise RuntimeError

            return incident_point, emergent_point

        if self.__mesh_splitter:
            # Split domains may have real mesh boundaries that deviate slightly
            # from the analytic rect/hex/tri equations.
            return __calcSplitterBoundaryHit(track_idx)
        elif self.__mesh.getBoundaryGeomType() == 'rect':
            return __calcRectHit(track_idx)
        else:
            if self.__mesh_splitter and self.__mesh_splitter.split_type == '1/6':
                return __calcTriHit(track_idx)
            else:
                return __calcHexHit(track_idx)


class TrackFactory:
    def __init__(self, mesh: Mesh, usr_spacing: float, num_azim: int, num_proc: int = 1,
                 mesh_spliter: MeshSplitter = None, output_vacuum_file: bool = False):
        self.__mesh = mesh
        self.__mesh_spliter = mesh_spliter
        self.__usr_spacing = usr_spacing
        self.__num_azim = num_azim
        self.__geom_type = mesh.getBoundaryGeomType()
        self.__output_vacuum_file = output_vacuum_file

        if self.__geom_type == 'rect':
            if self.__num_azim % 4 != 0:
                rt.printMessage(CASE_NAME, ERROR, "the num_azim must be the times of 4")
                raise ValueError
        else:
            if self.__num_azim % 6 != 0 or self.__num_azim == 6:
                rt.printMessage(CASE_NAME, ERROR, "the num_azim must be the times of 6")
                raise ValueError

        # geometry type is hex                 geometry type is rectangle
        #       num U                                num U
        #     /------\                          +-------------+
        #    /        \ num U + num V           |             |
        #   /          \                        |             | num V
        #   \          /                        |             |
        #    \        / num V                   |             |
        #     \------/                          +-------------+

        self.__num_tracks_U = []
        self.__num_tracks_V = []
        self.__num_tracks = []

        self.__spacings = []
        self.__azim_angs = []

        self.__angular_track_wraps = []

        self.__corrTrackParam()

        self.__num_proc = num_proc

        if mesh_spliter and mesh_spliter.split_type == '1/6' and (output_vacuum_file is True):
            self._mode = 'b'  # b means binary / w means decimal
            self._rayMessage = {i: [] for i in range(self.__num_azim)}
            self.OutputFileOption = OutputPassBottomRayIDFileHandle(self.__num_azim)

    def getNumAzim(self):
        return self.__num_azim

    def getAzimAngs(self):
        return self.__azim_angs

    def getSpacings(self):
        return self.__spacings

    def getNumTracks(self):
        return self.__num_tracks, self.__num_tracks_U, self.__num_tracks_V

    def getNumSegs(self):
        num_segs = [0 for _ in range(self.__num_azim)]
        for azim_idx in range(self.__num_azim // 2):
            num_segs[azim_idx] = self.__angular_track_wraps[azim_idx].getNumSeg()
        return num_segs

    def getGeomType(self):
        return self.__geom_type

    def __corrTrackParam(self):
        def __calcEffTrackParamInRect():
            x_min, x_max, y_min, y_max = self.__mesh.getBoundingBox()

            width_X = x_max - x_min
            width_Y = y_max - y_min

            num_azim_half_space = int(round(self.__num_azim / 2))

            for i in range(num_azim_half_space):
                phi = 2 * math.pi / self.__num_azim * float(i + 0.5)
                num_ray_X = int(math.ceil(width_X / self.__usr_spacing * math.sin(phi)))
                num_ray_Y = int(math.ceil(width_Y / self.__usr_spacing * abs(math.cos(phi))))
                num_ray = num_ray_Y + num_ray_X

                phi_eff = math.atan(width_Y / width_X * num_ray_X / num_ray_Y)
                eff_azim_ang = math.acos(math.cos(phi) / abs(math.cos(phi)) * math.cos(phi_eff))
                eff_spacing = width_X / num_ray_X * math.sin(eff_azim_ang)

                self.__azim_angs.append(eff_azim_ang)
                self.__spacings.append(eff_spacing)
                self.__num_tracks_U.append(num_ray_X)
                self.__num_tracks_V.append(num_ray_Y)
                self.__num_tracks.append(num_ray)

                angular_tracks = AngularTrackWrap(self.__mesh, eff_spacing, eff_azim_ang, num_ray)
                self.__angular_track_wraps.append(angular_tracks)

            # calc ray direction and location and create Ray
            ray_id = 0
            for i in range(num_azim_half_space):
                rays_azim = []  # rays with the direction of azimuthal angle i
                eff_spacing_Y = width_Y / self.__num_tracks_V[i]
                direction = [math.cos(self.__azim_angs[i]), math.sin(self.__azim_angs[i])]
                length_tmp = width_X / abs(math.cos(self.__azim_angs[i]))
                for j in range(self.__num_tracks[i]):
                    if self.__azim_angs[i] < math.pi / 2:
                        location = (
                            x_min, y_min + eff_spacing_Y * (j + 0.5) - length_tmp * math.sin(self.__azim_angs[i]))
                    else:
                        location = (
                            x_max, y_min + eff_spacing_Y * (j + 0.5) - length_tmp * math.sin(self.__azim_angs[i]))

                    self.__angular_track_wraps[i].setStart(j, location)

        def __calcEffTrackParamInTri():
            # The selected one-sixth domain is the bottom triangle of an
            # x-oriented regular hexagon.  Build tracks by projecting the
            # three triangle vertices onto the normal of each direction.  The
            # previous implementation placed every track on an extended
            # bottom-edge line, so some shallow tracks missed the triangle
            # completely and produced fewer than two boundary intersections.
            x_min, x_max, y_min, y_max = self.__mesh_spliter.getTargetBoundingBox()
            center = ((x_max + x_min) / 2, (y_min + y_max) / 2)
            half_bottom_edge = (x_max - x_min) / 4
            triangle_vertices = (
                (center[0] - half_bottom_edge, y_min),
                (center[0] + half_bottom_edge, y_min),
                center,
            )

            num_azim_half_space = int(round(self.__num_azim / 2))
            num_azim_one_sixth_space = int(round(num_azim_half_space / 3))

            for sector_index in range(3):
                for i in range(num_azim_one_sixth_space):
                    eff_azim_ang = 2 * np.pi / self.__num_azim * float(i + 0.5)
                    eff_azim_ang += sector_index * math.pi / 3
                    direction = (math.cos(eff_azim_ang), math.sin(eff_azim_ang))
                    normal = (-direction[1], direction[0])
                    projections = [
                        normal[0] * vertex[0] + normal[1] * vertex[1]
                        for vertex in triangle_vertices
                    ]
                    projection_min = min(projections)
                    projection_max = max(projections)
                    projection_width = projection_max - projection_min
                    num_tracks = max(1, int(math.ceil(projection_width / self.__usr_spacing)))
                    eff_spacing = projection_width / num_tracks

                    self.__num_tracks_U.append(num_tracks)
                    self.__num_tracks_V.append(0)
                    self.__num_tracks.append(num_tracks)
                    self.__azim_angs.append(eff_azim_ang)
                    self.__spacings.append(eff_spacing)

                    angular_tracks = AngularTrackWrap(self.__mesh, eff_spacing, eff_azim_ang, num_tracks,
                                                      self.__mesh_spliter, self.__output_vacuum_file)
                    self.__angular_track_wraps.append(angular_tracks)

                    normal_center_projection = normal[0] * center[0] + normal[1] * center[1]
                    for track_index in range(num_tracks):
                        projection = projection_min + eff_spacing * (track_index + 0.5)
                        offset = projection - normal_center_projection
                        location = (
                            center[0] + normal[0] * offset,
                            center[1] + normal[1] * offset,
                        )
                        angular_tracks.setStart(track_index, location)

        def __calcEffTrackParamInHex():
            # hex orientation is X
            if self.__mesh_spliter and self.__mesh_spliter.split_type == 'full':
                x_min, x_max, y_min, y_max = self.__mesh_spliter.getTargetBoundingBox()
            else:
                x_min, x_max, y_min, y_max = self.__mesh.getBoundingBox()

            hex_center = [(x_max + x_min) / 2, (y_min + y_max) / 2]
            hex_edge_width = (x_max - x_min) / 2

            num_azim_half_space = int(round(self.__num_azim / 2))
            num_azim_one_sixth_space = int(round(num_azim_half_space / 3))

            for j in range(3):
                for i in range(num_azim_one_sixth_space):
                    phi = 2 * np.pi / self.__num_azim * float(i + 0.5)
                    n1 = int(math.ceil(hex_edge_width / self.__usr_spacing * math.sin(phi)))
                    n2 = int(math.ceil(hex_edge_width / self.__usr_spacing * math.sin(math.pi / 3 - phi)))
                    n3 = (n1 + n2) * 2

                    phi_eff = math.atan(math.sqrt(3) / (2 * float(n2) / float(n1) + 1))
                    eff_azim_ang = phi_eff + j * math.pi / 3
                    eff_spacing = hex_edge_width * math.sin(phi_eff) / n1

                    self.__num_tracks_U.append(n1)
                    self.__num_tracks_V.append(n2)
                    self.__num_tracks.append(n3)
                    self.__azim_angs.append(eff_azim_ang)
                    self.__spacings.append(eff_spacing)

                    angular_tracks = AngularTrackWrap(self.__mesh, eff_spacing, eff_azim_ang, n3, self.__mesh_spliter)
                    self.__angular_track_wraps.append(angular_tracks)

            # To calculate the location and direction of each ray, it is easy to implement that rotate the ray in 0~30° to
            # target region
            x_hex_lr = hex_center[
                           0] + hex_edge_width / 2  # x coord of lower right point, which is the start point of 0~60
            for k in range(3):
                for i in range(num_azim_one_sixth_space):
                    ang_idx = k * num_azim_one_sixth_space + i
                    delta_x = self.__spacings[ang_idx] / math.sin(self.__azim_angs[ang_idx] - k * math.pi / 3)
                    for j in range(self.__num_tracks[ang_idx]):
                        loc_x = x_hex_lr - delta_x * (j + 0.5) - hex_center[0]
                        loc_y = y_min - hex_center[1]
                        # rotate the location
                        # -- --   --              --  --  --
                        # | x | = | cos(a) -sin(a) |  | x' |
                        # | y |   | sin(a)  cos(a) |  | y' |
                        # -- --   --              --  --  --
                        theta = k * math.pi / 3
                        loc_x1 = loc_x * math.cos(theta) - loc_y * math.sin(theta)
                        loc_y1 = loc_x * math.sin(theta) + loc_y * math.cos(theta)

                        location = (loc_x1 + hex_center[0], loc_y1 + hex_center[1])
                        self.__angular_track_wraps[ang_idx].setStart(j, location)

        rt.printMessage(CASE_NAME, INIT, "Generating track ...")

        if self.__geom_type == 'rect':
            __calcEffTrackParamInRect()
        else:
            if self.__mesh_spliter and self.__mesh_spliter.split_type == '1/6':
                __calcEffTrackParamInTri()
            else:
                __calcEffTrackParamInHex()

    def trackingForOneAzim(self, azim_idx: int):
        rt.printMessage(CASE_NAME, RUNNING,
                        "Tracking for segmentation in parallel [%d/%d] ..." % ((azim_idx + 1), self.__num_azim // 2))
        self.__angular_track_wraps[azim_idx].tracking()

        # todo: output ray id for pass bottom line in file(mult proc) -- CaoWei
        if self.__mesh_spliter and self.__mesh_spliter.split_type == '1/6' and self.__output_vacuum_file:
            _numTrack, _rays = self.__angular_track_wraps[azim_idx].getExitRay()

            _mode = 'wb' if 'b' in self._mode else 'w'

            self.OutputFileOption.genePassBottomLineFileFactory(str(azim_idx), azim_idx, _rays, _numTrack, _mode)

        return self.__angular_track_wraps[azim_idx]

    def tracking(self):
        num_quad, num_tri = self.__mesh.getNumFaces()
        num_trk = sum(self.__num_tracks)
        rt.printMessage(CASE_NAME, RUNNING, '-----------------------------------')
        rt.printMessage(CASE_NAME, RUNNING, '#MESH %d  #QUAD %d   #TRI %d' % (num_quad + num_tri, num_quad, num_tri))
        rt.printMessage(CASE_NAME, RUNNING,
                        '#AZIM %d  #TRACK %d  #PROC %d' % (self.__num_azim, num_trk, self.__num_proc))
        rt.printMessage(CASE_NAME, RUNNING, '-----------------------------------')

        if self.__mesh_spliter:
            rt.printMessage(CASE_NAME, RUNNING, 'Tracking in %s mode' % self.__mesh_spliter.split_type)

        # prepare segmentation in half space
        num_azim_hs = self.__num_azim // 2

        if self.__num_proc <= 0:
            self.__num_proc = 1

        if self.__num_proc > num_azim_hs:
            self.__num_proc = num_azim_hs

        if self.__num_proc == 1:
            for azim_idx in range(num_azim_hs):
                rt.printMessage(CASE_NAME, RUNNING,
                                "Tracking for segmentation [%d/%d] ..." % ((azim_idx + 1), num_azim_hs))
                self.__angular_track_wraps[azim_idx].tracking()

                # todo: output ray id for pass bottom line in file(one proc) -- CaoWei
                if self.__mesh_spliter and self.__mesh_spliter.split_type == '1/6' and self.__output_vacuum_file:
                    _numTrack, _rays = self.__angular_track_wraps[azim_idx].getExitRay()

                    _mode = 'wb' if 'b' in self._mode else 'w'

                    azim_id = azim_idx
                    self.OutputFileOption.genePassBottomLineFileFactory(str(azim_id), azim_id, _rays, _numTrack, _mode)
        else:
            pool = mp.Pool(processes=self.__num_proc)

            rets = []
            for azim_idx in range(num_azim_hs):
                rets.append(pool.apply_async(func=self.trackingForOneAzim, args=(azim_idx,)))

            pool.close()
            pool.join()

            for azim_idx in range(num_azim_hs):
                self.__angular_track_wraps[azim_idx] = rets[azim_idx].get()

        # todo: generator pass bottom line ray ID file -- CaoWei
        if self.__mesh_spliter and self.__mesh_spliter.split_type == '1/6' and self.__output_vacuum_file:
            _numTracksCurrent, self._rayMessage = 0, {i: [] for i in range(self.__num_azim)}

            _mode = 'rb' if 'b' in self._mode else 'r'

            for azim_idx in range(self.__num_azim // 2):
                self.OutputFileOption._rayMessage = self._rayMessage
                _numTracksAzim = self.OutputFileOption.readPassBottomLineFileFactory(azim_idx, _numTracksCurrent, _mode)
                _numTracksCurrent = _numTracksCurrent + _numTracksAzim

                self._rayMessage[azim_idx] = sorted(self._rayMessage[azim_idx], key=lambda line: line.id, reverse=True)

                if _mode == 'rb':
                    __mode = 'wb' if azim_idx == 0 else 'ab'
                else:
                    __mode = 'w' if azim_idx == 0 else 'a'

                _suffix = '.dat' if 'b' in _mode else '.txt'

                self.OutputFileOption.genePassBottomLineFileFactory('file', azim_idx + self.__num_azim // 2,
                                                                    self._rayMessage[azim_idx], _numTracksAzim, __mode)
                os.remove(f'./{azim_idx}' + _suffix)

        rt.printMessage(CASE_NAME, RUNNING, 'Tracking success !')

    def getAngularTrackWrap(self, azim_idx):
        if azim_idx < self.__num_azim // 2:
            return self.__angular_track_wraps[azim_idx]
        else:
            rt.printMessage(CASE_NAME, ERROR, 'the azim_idx is larger than the number of eff azim, plz check!')


class SegmentationTally:
    def __init__(self, mesh: Mesh, track_factory: TrackFactory, mesh_spliter: MeshSplitter = None):
        self.mesh = mesh
        self.track_factory = track_factory
        self.mesh_spliter = mesh_spliter
        self.num_azim = track_factory.getNumAzim()
        self.spacings = track_factory.getSpacings()
        self.azim_angs = track_factory.getAzimAngs()
        self.angular_num_tracks, self.angular_num_tracks_U, self.angular_num_tracks_V = track_factory.getNumTracks()

        self.num_tracks = 0
        self.num_segs = 0

        self.num_faces = sum(mesh.getNumFaces())

        # azim_track_idx means the index of track at current azimuthal angle
        self.segmentation = {'azim_track_idx': [],
                             'azim_idx': [],
                             'start': [],
                             'end': [],
                             'num_seg': [],
                             'seg_length': [],
                             'swept_face': [],
                             'mirror_azim_idx': [],
                             'mirror_track_idx': []}

        self.__completeSegmentation()

        if self.track_factory.getGeomType() == 'hex':
            self.__getMirrorTrackIndex()

        self.__corrSegLenByVol()

    def __completeSegmentation(self):
        self.spacings += self.spacings
        self.azim_angs += [azim_ang + math.pi for azim_ang in self.azim_angs]

        self.angular_num_tracks += self.angular_num_tracks
        self.angular_num_tracks_U += self.angular_num_tracks_U
        self.angular_num_tracks_V += self.angular_num_tracks_V

        self.num_tracks = sum(self.angular_num_tracks)

        for azim_idx in range(self.num_azim):
            if azim_idx < self.num_azim // 2:
                track_wrap = self.track_factory.getAngularTrackWrap(azim_idx)
            else:
                track_wrap = self.track_factory.getAngularTrackWrap(azim_idx - self.num_azim // 2)

            aizm_num_track = self.angular_num_tracks[azim_idx]
            for azim_track_idx in range(aizm_num_track):

                start = track_wrap.getStart(azim_track_idx)
                end = track_wrap.getEnd(azim_track_idx)

                seg_length = track_wrap.getSegLength(azim_track_idx)
                swept_face = track_wrap.getSweptFace(azim_track_idx)

                self.num_segs += len(seg_length)

                self.segmentation['azim_track_idx'].append(azim_track_idx)
                self.segmentation['azim_idx'].append(azim_idx)

                if azim_idx < self.num_azim // 2:
                    self.segmentation['start'].append(start)
                    self.segmentation['end'].append(end)
                    self.segmentation['seg_length'].append(seg_length)
                    self.segmentation['swept_face'].append(swept_face)
                    self.segmentation['num_seg'].append(len(seg_length))

                else:
                    # all the data is reversed in pi~2pi
                    isrt_idx = sum(self.angular_num_tracks[:azim_idx])
                    self.segmentation['start'].insert(isrt_idx, end)
                    self.segmentation['end'].insert(isrt_idx, start)
                    self.segmentation['seg_length'].insert(isrt_idx, seg_length[::-1])
                    self.segmentation['swept_face'].insert(isrt_idx, swept_face[::-1])
                    self.segmentation['num_seg'].insert(isrt_idx, len(seg_length))

    def __getMirrorTrackIndex(self):
        # Mirror matching was originally an exact coordinate lookup. With
        # real-edge boundary hits, symmetric endpoints can differ by a few
        # thousandths due to hand-built boundary offsets, so match the nearest
        # start point in the target azimuth while keeping the tolerance below
        # the spacing scale to avoid jumping to a neighboring ray.
        min_track_spacing = min(self.spacings) if self.spacings else SMALL_VALUE
        if self.mesh_spliter and self.mesh_spliter.split_type == '1/6':
            # Independent track families in the triangular symmetry domain
            # can have different projected grids on a slanted side.  The
            # nearest physical boundary point can therefore differ by more
            # than a quarter spacing, while still identifying the same mirror
            # ray.
            mirror_match_tol = max(SMALL_VALUE, 1.5 * min_track_spacing)
        else:
            mirror_match_tol = max(SMALL_VALUE, 0.25 * min_track_spacing)
        relaxed_match_cnt = 0
        max_match_dist = 0.0

        azim_start_offsets = [0]
        for num_track in self.angular_num_tracks:
            azim_start_offsets.append(azim_start_offsets[-1] + num_track)

        azim_start_coords = []
        for azim_idx in range(self.num_azim):
            start_idx = azim_start_offsets[azim_idx]
            end_idx = azim_start_offsets[azim_idx + 1]
            azim_start_coords.append(
                np.array(self.segmentation['start'][start_idx:end_idx], dtype=float)
            )

        try:
            from scipy.spatial import cKDTree
        except ImportError:
            cKDTree = None

        # Use a KDTree when SciPy is available; the numpy fallback preserves the
        # same nearest-point behavior for environments without SciPy.
        if cKDTree:
            azim_start_trees = [cKDTree(coords) if len(coords) else None for coords in azim_start_coords]
        else:
            azim_start_trees = [None for _ in azim_start_coords]

        def __findNearestMirrorTrack(mirror_azim_idx: int, hit_coord: tuple) -> tuple:
            coords = azim_start_coords[mirror_azim_idx]
            if len(coords) == 0:
                return -1, math.inf

            if azim_start_trees[mirror_azim_idx] is not None:
                match_dist, mirror_track_idx = azim_start_trees[mirror_azim_idx].query(hit_coord)
                return int(mirror_track_idx), float(match_dist)

            hit_coord_np = np.array(hit_coord, dtype=float)
            dist_sq = np.sum((coords - hit_coord_np) ** 2, axis=1)
            mirror_track_idx = int(np.argmin(dist_sq))
            return mirror_track_idx, float(math.sqrt(dist_sq[mirror_track_idx]))

        def __getHitHexBdrySideIndex(azim_track_idx: int, azim_idx: int) -> int:
            #                  2
            #        /-------------------\
            #       /                     \
            #   3  /                       \
            #     /                         \  1
            #    /                           \
            #   /                             \
            #   \                             /
            #    \                           /
            #     \                         /
            #   4  \                       /  0
            #       \                     /
            #        \-------------------/
            #                  5
            # params: n1: number of tracks in the third boundary of the hit boundary set
            #         n2: number of tracks in the first boundary of the hit boundary set
            #         n3: number of tracks in the second boundary of the hit boundary set
            if azim_idx < self.num_azim / 2:
                n1 = self.angular_num_tracks_U[azim_idx]
                n2 = self.angular_num_tracks_V[azim_idx]
            else:
                n1 = self.angular_num_tracks_U[azim_idx - self.num_azim // 2]
                n2 = self.angular_num_tracks_V[azim_idx - self.num_azim // 2]
            n3 = n1 + n2
            # get hit boundary set
            if azim_idx < self.num_azim / 6:
                hit_boundary_set = (0, 1, 2)
            elif self.num_azim / 6 <= azim_idx < self.num_azim / 3:
                hit_boundary_set = (1, 2, 3)
            elif self.num_azim / 3 <= azim_idx < self.num_azim / 2:
                hit_boundary_set = (2, 3, 4)
            elif self.num_azim / 2 <= azim_idx < self.num_azim * 2 / 3:
                hit_boundary_set = (3, 4, 5)
            elif self.num_azim * 2 / 3 <= azim_idx < self.num_azim * 5 / 6:
                hit_boundary_set = (4, 5, 0)
            elif self.num_azim * 5 / 6 <= azim_idx < self.num_azim:
                hit_boundary_set = (5, 0, 1)
            else:
                rt.printMessage(CASE_NAME, ERROR, 'azim_idx do not exist')
                raise RuntimeError
            # get hit boundary id
            if azim_track_idx < n2:
                hit_boundary_id = hit_boundary_set[0]
            elif n2 <= azim_track_idx < n3 + n2:
                hit_boundary_id = hit_boundary_set[1]
            elif n3 + n2 <= azim_track_idx < 2 * n3:
                hit_boundary_id = hit_boundary_set[2]
            else:
                rt.printMessage(CASE_NAME, ERROR, 'track_idx do not exit')
                raise RuntimeError

            return hit_boundary_id

        def __getHitTriBdrySideIndex(track_idx: int) -> int:
            #
            #       /\
            #  0 v /  \ 1
            #     /____\
            #        2 u
            if self.mesh_spliter:
                x_min, x_max, y_min, y_max = self.mesh_spliter.getTargetBoundingBox()
            else:
                x_min, x_max, y_min, y_max = self.mesh.getBoundingBox()

            xend, yend = self.segmentation['end'][track_idx]
            if abs(yend - y_min) < SMALL_VALUE:
                hit_boundary_id = 2
            elif xend > (x_max + x_min) / 2 + SMALL_VALUE:
                hit_boundary_id = 1
            elif xend < (x_max + x_min) / 2 + SMALL_VALUE:
                hit_boundary_id = 0
            else:
                rt.printMessage(CASE_NAME, ERROR, 'Can not found the hit boundary of track_idx=%d' % (track_idx))
                raise RuntimeError

            return hit_boundary_id

        def __getMirrorAzimAngIndex(track_idx: int, azim_track_idx: int, azim_idx: int) -> int:
            if self.mesh_spliter and self.mesh_spliter.split_type == '1/6':
                # Reflect the transport direction across the side reached by
                # the track.  The old index formula was tied to the original
                # full-hex angular ordering and did not match the regular
                # angular sampling used by the one-sixth triangle.
                hit_boundary_id = __getHitTriBdrySideIndex(track_idx)
                boundary_angles = {
                    0: math.pi / 3,     # left slanted side
                    1: -math.pi / 3,    # right slanted side
                    2: 0.0,             # bottom side
                }
                incident_angle = self.azim_angs[azim_idx]
                reflected_angle = (
                    2.0 * boundary_angles[hit_boundary_id] - incident_angle
                ) % (2.0 * math.pi)

                def angular_distance(angle_a: float, angle_b: float) -> float:
                    delta = abs(angle_a - angle_b) % (2.0 * math.pi)
                    return min(delta, 2.0 * math.pi - delta)

                return min(
                    range(self.num_azim),
                    key=lambda candidate: angular_distance(
                        reflected_angle, self.azim_angs[candidate]
                    ),
                )

            # get the azim_id upper limit of 120, 240
            azim_id_limit_120 = int(self.num_azim / 3 - 1)
            azim_id_limit_240 = int(self.num_azim * 2 / 3 - 1)

            if self.mesh_spliter and self.mesh_spliter.split_type == '1/6':
                hit_boundary_id = __getHitTriBdrySideIndex(track_idx)
            else:
                hit_boundary_id = __getHitHexBdrySideIndex(azim_track_idx, azim_idx)

            if hit_boundary_id == 0 or hit_boundary_id == 3:
                if azim_idx <= azim_id_limit_120:
                    next_azim_id = azim_id_limit_120 - azim_idx
                elif azim_idx > azim_id_limit_120:
                    next_azim_id = self.num_azim + azim_id_limit_120 - azim_idx
                else:
                    rt.printMessage(CASE_NAME, ERROR, 'azim_id do not exist')
                    raise RuntimeError

            elif hit_boundary_id == 1 or hit_boundary_id == 4:
                if azim_idx <= azim_id_limit_240:
                    next_azim_id = azim_id_limit_240 - azim_idx
                elif azim_id_limit_240 < azim_idx:
                    next_azim_id = self.num_azim + azim_id_limit_240 - azim_idx
                else:
                    rt.printMessage(CASE_NAME, ERROR, 'azim_id do not exist')
                    raise RuntimeError

            elif hit_boundary_id == 2 or hit_boundary_id == 5:
                next_azim_id = self.num_azim - azim_idx - 1

            else:
                return -1

            return next_azim_id

        for track_idx in range(self.num_tracks):
            azim_idx = self.segmentation['azim_idx'][track_idx]
            azim_track_idx = self.segmentation['azim_track_idx'][track_idx]

            mirror_azim_idx = __getMirrorAzimAngIndex(track_idx, azim_track_idx, azim_idx)

            self.segmentation['mirror_azim_idx'].append(mirror_azim_idx)

            hit_coord = self.segmentation['end'][track_idx]
            mirror_azim_track_idx, match_dist = __findNearestMirrorTrack(mirror_azim_idx, hit_coord)

            if match_dist < mirror_match_tol:
                self.segmentation['mirror_track_idx'].append(mirror_azim_track_idx)
                max_match_dist = max(max_match_dist, match_dist)
                if match_dist > SMALL_VALUE:
                    relaxed_match_cnt += 1
            else:
                errmsg = 'can not find the mirror track, plz check!\n'
                errmsg += 'GLOBAL_TRACK_INDEX: %d \n' % track_idx
                errmsg += 'AZIM_INDEX        : %d \n' % azim_idx
                errmsg += 'AZIM_TRACK_INDEX  : %d \n' % azim_track_idx
                errmsg += 'MIRROR_AZIM_INDEX : %d \n' % mirror_azim_idx
                errmsg += 'NEAREST_TRACK     : %d \n' % mirror_azim_track_idx
                errmsg += 'NEAREST_DISTANCE  : %.10e \n' % match_dist
                errmsg += 'MATCH_TOLERANCE   : %.10e \n' % mirror_match_tol
                errmsg += 'TRACK_END_X       : %.10f \n' % hit_coord[0]
                errmsg += 'TRACK_END_Y       : %.10f' % hit_coord[1]
                rt.printMessage(CASE_NAME, ERROR, errmsg)

        if relaxed_match_cnt > 0:
            msg = 'Mirror track matching used nearest-boundary matching for %d tracks, max distance %.10e' % (
                relaxed_match_cnt, max_match_dist)
            rt.printMessage(CASE_NAME, WARNING, msg)

    def __corrSegLenByVol(self):
        # fsr area = A
        # in azimuthal angle a:
        #   fsr area in tracing = A' = sum(seg_len) * spacing
        #   A = A'
        #   A = sum(seg_len) * spacing * f
        #   f = A / A'
        #   seg_len_corr = seg_len * f

        azim_seg_len_in_face = [[0. for _i in range(self.num_faces)] for _j in range(self.num_azim)]

        for track_idx in range(self.num_tracks):
            azim_idx = self.segmentation['azim_idx'][track_idx]

            for seg_idx in range(self.segmentation['num_seg'][track_idx]):
                seg_len = self.segmentation['seg_length'][track_idx][seg_idx]

                swept_face = self.segmentation['swept_face'][track_idx][seg_idx]
                face_idx = self.mesh.getFaceIndex(swept_face)

                azim_seg_len_in_face[azim_idx][face_idx] += seg_len

        vol_corr_factor = [[0. for _i in range(self.num_faces)] for _j in range(self.num_azim)]

        for azim_idx in range(self.num_azim):
            for face_idx in range(self.num_faces):
                face_id = self.mesh.getFaceID(face_idx)

                sum_seg_len = azim_seg_len_in_face[azim_idx][face_idx]
                if sum_seg_len > SMALL_VALUE:
                    f = self.mesh.getArea(face_id) / (sum_seg_len * self.spacings[azim_idx])
                else:
                    f = 1

                vol_corr_factor[azim_idx][face_idx] = f

        for track_idx in range(self.num_tracks):
            azim_idx = self.segmentation['azim_idx'][track_idx]

            for seg_idx in range(self.segmentation['num_seg'][track_idx]):
                swept_face = self.segmentation['swept_face'][track_idx][seg_idx]
                face_idx = self.mesh.getFaceIndex(swept_face)

                self.segmentation['seg_length'][track_idx][seg_idx] *= vol_corr_factor[azim_idx][face_idx]
