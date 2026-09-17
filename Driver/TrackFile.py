import RayTracing as smrt
from Driver.PhysicsObj import RayTracingModule, RTMMap
import math
import h5py as h5
import numpy as np


class AlphaTrackFileFormatter:
    def __init__(self, file_path: str):
        self.__output_path = file_path

    def pinInfoOutPut(self, rtms: list):
        trk_file = open(self.__output_path, 'w')

        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)

        geom_type = rtms[0].mesh.getBoundaryGeomType()

        if geom_type == 'rect':
            trk_file.write('GeomType 0\n')
        elif geom_type == 'hex':
            if rtms[0].mesh_splitter and rtms[0].mesh_splitter.split_type == '1/6':
                trk_file.write('GeomType 2\n')
            else:
                trk_file.write('GeomType 1\n')
        else:
            print('Geometry type in track file is error, plz check !')
            raise ValueError

        trk_file.write('PART_0_PIN_INFO\n\n')

        num_of_pins = len(rtms)
        trk_file.write('NumOfPins ' + str(num_of_pins) + '\n')
        trk_file.write('!pin_id  rtm_id  num_of_fsr  mat_id  fsr_id  //face_id\n')

        rtm_idx = -1
        for rtm in rtms:
            rtm_idx += 1
            assert isinstance(rtm, RayTracingModule)
            if rtm.mesh_splitter:
                num_faces = rtm.mesh_splitter.getNumFaces()
            else:
                num_quad, num_tri = rtm.mesh.getNumFaces()
                num_faces = num_quad + num_tri

            trk_file.write(str(rtm.physical_id).rjust(7, ' '))  # pin id
            trk_file.write(str(rtm.rtm_id).rjust(8, ' '))  # rtm id
            trk_file.write(str(num_faces).rjust(12, ' '))  # number of FSR

            mat_names = rtm.mesh.getPhysicalGroups()

            for face_idx in range(num_faces):

                if rtm.mesh_splitter:
                    face_id = rtm.mesh_splitter.getFaceIDs()[face_idx]
                    fsr_id = rtm.mesh_splitter.getFaceIndex(face_id)
                else:
                    face_id = rtm.mesh.getFaceID(face_idx)
                    fsr_id = face_idx

                mat_name = rtm.mesh.getPhysicalName(face_id)
                mat_idx = mat_names.index(mat_name)

                if face_idx == 0:
                    trk_file.write(str(mat_idx).rjust(8, ' '))
                    trk_file.write(str(fsr_id).rjust(8, ' '))
                    trk_file.write('\t!' + str(face_id).rjust(8, ' '))
                    trk_file.write('\n')
                else:
                    trk_file.write(str(mat_idx).rjust(35, ' '))
                    trk_file.write(str(fsr_id).rjust(8, ' '))
                    trk_file.write('\t!' + str(face_id).rjust(8, ' '))
                    trk_file.write('\n')

            trk_file.write('!'.ljust(60, '-') + '\n')
        trk_file.write('\n\n')
        trk_file.close()

    # todo: complete the lattice function
    def layoutInfoOutPut(self, rtm_maps: list):
        track_file = open(self.__output_path, 'a')

        for rtm_map in rtm_maps:
            assert isinstance(rtm_map, RTMMap)

        track_file.write('PART_1_PIN_LAYOUT_INFO\n\n')

        track_file.write('NumOfPinMaps')
        track_file.write(str(len(rtm_maps)).rjust(6, ' ') + '\n')

        track_file.write('NumOfSublayers')
        num_of_sublayers = sum([rtm_map.axial_mesh for rtm_map in rtm_maps])
        track_file.write(str(num_of_sublayers).rjust(4, ' ') + '\n')

        track_file.write('!id  dim_X  dim_Y  width_Z  mesh_Z\n')

        # find the largest gap of neighboring Pin ID
        str_len_old = 2
        # for k in range(len(rtm_maps)):
        #     str_len = len(str(max(map(max, rtm_maps[k].layout))))
        #     if str_len > str_len_old:
        #         str_len_old = str_len

        rtm_idx = 0
        for rtm in rtm_maps:
            assert isinstance(rtm, RTMMap)
            track_file.write(str(rtm_idx).rjust(3, ' '))
            dim_X = rtm.dim_X
            dim_Y = rtm.dim_Y
            track_file.write(str(dim_X).rjust(7, ' '))
            track_file.write(str(dim_Y).rjust(7, ' '))

            width_Z = rtm.width_Z
            mesh_Z = rtm.axial_mesh
            track_file.write(str(format(width_Z, '.2f')).rjust(9, ' '))
            track_file.write(str(mesh_Z).rjust(8, ' ') + '\n')

            for i in range(dim_Y):
                for j in range(dim_X):
                    track_file.write(' ')
                    track_file.write(str(0).rjust(str_len_old, ' '))
                track_file.write('\n')
            track_file.write('\n')

            rtm_idx += 1

        track_file.write('\n')

    def rayParamOutPut(self, track_factory: smrt.Track.TrackFactory):
        trk_file = open(self.__output_path, 'a')

        eff_azims = track_factory.getAzimAngs()
        eff_spacings = track_factory.getSpacings()
        num_tracks, num_tracks_U, num_tracks_V = track_factory.getNumTracks()

        trk_file.write('PART_2_RAY_PARAM_INFO\n\n')
        trk_file.write('NumOfAzimAng  ' + str(track_factory.getNumAzim()) + '\n')
        trk_file.write('!azim  numU  numV  effAzim               effSpacing\n')

        for azim_id in range(track_factory.getNumAzim()):
            str_tmp = str(azim_id).rjust(5, ' ')

            if azim_id < track_factory.getNumAzim() / 2:
                local_azim_id = azim_id
                azim_diff = 0.0
            else:
                local_azim_id = azim_id - int(track_factory.getNumAzim() / 2)
                azim_diff = math.pi

            str_tmp += str(num_tracks_U[local_azim_id]).rjust(6, ' ')
            str_tmp += str(num_tracks_V[local_azim_id]).rjust(6, ' ')
            str_tmp += str(format(eff_azims[local_azim_id] + azim_diff, '.18f')).rjust(22, ' ')
            str_tmp += str(format(eff_spacings[local_azim_id], '.18f')).rjust(22, ' ') + '\n'

            trk_file.write(str_tmp)
        trk_file.write('\n\n')
        trk_file.close()

    def rtmGeomInfoOutPut(self, rtms: list):
        trk_file = open(self.__output_path, 'a')
        trk_file.write('PART_3_RTM_GEOM_INFO\n\n')
        trk_file.write('TotalNumOfRTM  ' + str(len(rtms)) + '\n')

        tot_num_seg = 0
        tot_num_trk = 0
        tot_num_fsr = 0

        # todo: ignore the rtm with the same rtm_id
        rtm_num_fsrs = []
        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)
            seg_tally = rtm.getSegmentationTally()

            if rtm.mesh_splitter:
                num_faces = rtm.mesh_splitter.getNumFaces()
            else:
                num_quad, num_tri = rtm.mesh.getNumFaces()
                num_faces = num_quad + num_tri
            rtm_num_fsrs.append(num_faces)
            tot_num_fsr += num_faces

            tot_num_trk += seg_tally.num_tracks
            tot_num_seg += seg_tally.num_segs

        trk_file.write('TotalNumOfFSR  ' + str(tot_num_fsr) + '\n')
        trk_file.write('TotalNumOfTrk  ' + str(tot_num_trk) + '\n')
        trk_file.write('TotalNumOfSeg  ' + str(tot_num_seg) + '\n')
        xmin, xmax, ymin, ymax = rtms[0].mesh.getBoundingBox()
        rtm_width_X = xmax - xmin
        rtm_width_Y = ymax - ymin
        trk_file.write('RTMWidthX      ' + str(rtm_width_X) + '\n')
        trk_file.write('RTMWidthY      ' + str(rtm_width_Y) + '\n\n\n')

        rtm_id_collect = []

        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)
            seg_tally = rtm.getSegmentationTally()

            if rtm.rtm_id in rtm_id_collect:
                continue
            else:
                rtm_id_collect.append(rtm.rtm_id)

            trk_file.write('RTM_ID  ' + str(rtm.rtm_id) + '\n\n')
            trk_file.write('FSRAreaInfo\n')
            trk_file.write('NumOfFSR  ' + str(rtm_num_fsrs[rtms.index(rtm)]) + '\n')
            trk_file.write('!fsr_id  fsr_area\n')

            if rtm.mesh_splitter:
                num_faces = rtm.mesh_splitter.getNumFaces()
            else:
                num_quad, num_tri = rtm.mesh.getNumFaces()
                num_faces = num_quad + num_tri

            for face_idx in range(num_faces):
                if rtm.mesh_splitter:
                    face_id = rtm.mesh_splitter.getFaceIDs()[face_idx]
                else:
                    face_id = rtm.mesh.getFaceID(face_idx)
                trk_file.write(str(face_idx).rjust(7, ' '))
                trk_file.write(str(format(rtm.mesh.getArea(face_id), '.16f')).rjust(20, ' ') + '\n')
            trk_file.write('\n')

            trk_file.write('TrackInfo\n')

            if rtm.mesh.getBoundaryGeomType() != 'hex':
                trk_file.write('!azim_id  num_ray  ray_id  num_seg  seg_id  seg_length          fsr_id\n')
            else:
                trk_file.write(
                    '!azim_id  num_ray  ray_id  next_azim_id  next_ray_id  num_seg  seg_id  seg_length          fsr_id\n')

            for track_idx in range(seg_tally.num_tracbks):
                azim_idx = seg_tally.segmentation['azim_idx'][track_idx]
                azim_num_track = seg_tally.angular_num_tracks[azim_idx]
                azim_track_idx = seg_tally.segmentation['azim_track_idx'][track_idx]
                trk_num_segs = seg_tally.segmentation['num_seg'][track_idx]

                if track_idx == 0 or azim_idx > seg_tally.segmentation['azim_idx'][track_idx - 1]:
                    trk_file.write(str(azim_idx).rjust(8, ' '))
                    trk_file.write(str(azim_num_track).rjust(9, ' '))

                if azim_track_idx == 0:
                    trk_file.write(str(azim_track_idx).rjust(8, ' '))
                else:
                    trk_file.write(str(azim_track_idx).rjust(25, ' '))

                if rtm.mesh.getBoundaryGeomType() == 'hex':
                    mirror_azim_ang_idx = seg_tally.segmentation['mirror_azim_idx'][track_idx]
                    mirror_track_idx = seg_tally.segmentation['mirror_track_idx'][track_idx]
                    trk_file.write(str(mirror_azim_ang_idx).rjust(14, ' '))
                    trk_file.write(str(mirror_track_idx).rjust(13, ' '))

                trk_file.write(str(trk_num_segs).rjust(9, ' '))

                for i in range(trk_num_segs):
                    if i == 0:
                        trk_file.write(str(i).rjust(8, ' '))
                    else:
                        if rtm.mesh.getBoundaryGeomType() != 'hex':
                            size = 42
                        else:
                            size = 42 + 19 + 8
                        trk_file.write(str(i).rjust(size, ' '))

                    seg_len = seg_tally.segmentation['seg_length'][track_idx][i]
                    swept_face = seg_tally.segmentation['swept_face'][track_idx][i]

                    if rtm.mesh_splitter:
                        fsr_id = rtm.mesh_splitter.getFaceIndex(swept_face)
                    else:
                        fsr_id = rtm.mesh.getFaceIndex(swept_face)

                    trk_file.write(str(format(seg_len, '.16f')).rjust(20, ' '))
                    trk_file.write(str(fsr_id).rjust(8, ' ') + '\n')

                if rtm.mesh.getBoundaryGeomType() != 'hex':
                    size = 81
                else:
                    size = 81 + 28
                trk_file.write('!'.ljust(size, '-') + '\n')

            trk_file.write('\n\n')


class ZMOCTrackFileFormatter:
    def __init__(self, file_path: str):
        self.__output_path = file_path

    def pinInfoOutput(self, rtms: list):
        f = h5.File(self.__output_path, 'w')
        geom_type = rtms[0].mesh.getBoundaryGeomType()
        if geom_type == 'rect':
            f.create_dataset('GeomType', data=np.array([0]))
        elif geom_type == 'hex':
            if rtms[0].mesh_splitter:
                f.create_dataset('GeomType', data=np.array([2]))
            else:
                f.create_dataset('GeomType', data=np.array([1]))
        else:
            print('Geometry type in track file is error, plz check !')
            raise ValueError

        pin_group = f.create_group('PinInfo')

        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)
            group_name = 'Pin' + str(rtm.rtm_id)
            pin_subgroup = pin_group.create_group(group_name)
            pin_subgroup.create_dataset('PinID', data=np.array([rtm.rtm_id]))
            pin_subgroup.create_dataset('RTMID', data=np.array([rtm.rtm_id]))

            if rtm.mesh_splitter:
                num_faces = rtm.mesh_splitter.getNumFaces()
            else:
                num_quad, num_tri = rtm.mesh.getNumFaces()
                num_faces = num_quad + num_tri
            pin_subgroup.create_dataset('#FSRs', data=np.array([num_faces]))

            mat_names = rtm.mesh.getPhysicalGroups()
            mat_ids = []
            face_ids = []
            fsr_ids = []
            areas = []
            for face_idx in range(num_faces):
                if rtm.mesh_splitter:
                    face_id = rtm.mesh_splitter.getFaceIDs()[face_idx]
                    fsr_id = rtm.mesh_splitter.getFaceIndex(face_id)
                else:
                    face_id = rtm.mesh.getFaceID(face_idx)
                    fsr_id = face_idx
                mat_name = rtm.mesh.getPhysicalName(face_id)
                mat_idx = mat_names.index(mat_name)

                area = rtm.mesh.getArea(face_id)

                mat_ids.append(mat_idx)
                face_ids.append(face_id)
                fsr_ids.append(fsr_id)
                areas.append(area)

            pin_subgroup.create_dataset('MatIDs', data=np.array(mat_ids))
            pin_subgroup.create_dataset('FSRIDs', data=np.array(fsr_ids))
            pin_subgroup.create_dataset('FaceIDs', data=np.array(face_ids))
            pin_subgroup.create_dataset('Areas', data=np.array(areas))

        f.close()

    def rayParamOutput(self, track_factory: smrt.Track.TrackFactory):
        f = h5.File(self.__output_path, 'a')

        eff_azims = track_factory.getAzimAngs()
        eff_spacings = track_factory.getSpacings()
        num_tracks, num_tracks_U, num_tracks_V = track_factory.getNumTracks()

        ray_group = f.create_group('RayParamInfo')
        ray_group.create_dataset('#Azimuth', data=np.array([track_factory.getNumAzim()]))

        ray_group.create_dataset('#TracksU', data=np.array(num_tracks_U * 2))
        ray_group.create_dataset('#TracksV', data=np.array(num_tracks_V * 2))
        ray_group.create_dataset('Spacings', data=np.array(eff_spacings * 2))
        ray_group.create_dataset('Azimuth', data=np.array(eff_azims + [azim + math.pi for azim in eff_azims]))

        f.close()

    def rtmGeomInfoOutput(self, rtms: list):
        f = h5.File(self.__output_path, 'a')
        seg_group = f.create_group('Segmentation')

        tot_num_fsr = 0
        tot_num_seg = 0
        tot_num_trk = 0

        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)
            seg_tally = rtm.getSegmentationTally()

            if rtm.mesh_splitter:
                num_faces = rtm.mesh_splitter.getNumFaces()
            else:
                num_quad, num_tri = rtm.mesh.getNumFaces()
                num_faces = num_quad + num_tri
            tot_num_fsr += num_faces

            tot_num_trk += seg_tally.num_tracks
            tot_num_seg += seg_tally.num_segs

        seg_group.create_dataset('#TotalFSR', data=np.array([tot_num_fsr]))
        seg_group.create_dataset('#TotalTrack', data=np.array([tot_num_trk]))
        seg_group.create_dataset('#TotalSegments', data=np.array([tot_num_seg]))

        for rtm in rtms:
            assert isinstance(rtm, RayTracingModule)
            rtm_seg_group = seg_group.create_group('RTM'+str(rtm.rtm_id))
            seg_tally = rtm.getSegmentationTally()

            for track_idx in range(seg_tally.num_tracks):
                azim_idx = seg_tally.segmentation['azim_idx'][track_idx]
                azim_num_track = seg_tally.angular_num_tracks[azim_idx]
                azim_track_idx = seg_tally.segmentation['azim_track_idx'][track_idx]
                trk_num_segs = seg_tally.segmentation['num_seg'][track_idx]

                if (azim_idx == 0 or azim_idx > seg_tally.segmentation['azim_idx'][track_idx - 1]) and azim_track_idx == 0:
                    azim_group = rtm_seg_group.create_group('AZIM' + str(azim_idx))
                    azim_group.create_dataset('#NumTrackInAzim', data=np.array([azim_num_track]))

                ray_group = azim_group.create_group('RAY' + str(azim_track_idx))

                if rtm.mesh.getBoundaryGeomType() == 'hex':
                    mirror_azim_ang_idx = seg_tally.segmentation['mirror_azim_idx'][track_idx]
                    mirror_track_idx = seg_tally.segmentation['mirror_track_idx'][track_idx]
                    ray_group.create_dataset("#MirrorAzim", data=np.array([mirror_azim_ang_idx]))
                    ray_group.create_dataset("#MirrorTrack", data=np.array([mirror_track_idx]))

                ray_group.create_dataset("#NumSegInTrack", data=np.array([trk_num_segs]))

                seg_len, sweep_fsr, sweep_face = [], [], []
                for i in range(trk_num_segs):
                    seg_len.append(seg_tally.segmentation['seg_length'][track_idx][i])

                    faceID = seg_tally.segmentation['swept_face'][track_idx][i]
                    sweep_face.append(faceID)

                    if rtm.mesh_splitter:
                        sweep_fsr.append(rtm.mesh_splitter.getFaceIndex(faceID))
                    else:
                        sweep_fsr.append(rtm.mesh.getFaceIndex(faceID))

                ray_group.create_dataset("SegLength", data=np.array(seg_len))
                ray_group.create_dataset("FSR", data=np.array(sweep_fsr))
                # ray_group.create_dataset("Face", data=np.array(sweep_face))

        f.close()

        pass


class TrackFile:
    def __init__(self, file_path: str, fmt: str = 'zmoc'):
        self.file_path = file_path
        self.fmt = fmt

    def __alphaFmtFactory(self, rtms: list, rtm_maps: list):
        alpha_fmtter = AlphaTrackFileFormatter(self.file_path)

        alpha_fmtter.pinInfoOutPut(rtms)
        alpha_fmtter.layoutInfoOutPut(rtm_maps)
        alpha_fmtter.rayParamOutPut(rtms[0].track_factory)
        alpha_fmtter.rtmGeomInfoOutPut(rtms)

    def __zmocFmtFactory(self, rtms: list, rtm_maps: list):
        zmoc_fmtter = ZMOCTrackFileFormatter(self.file_path)
        zmoc_fmtter.pinInfoOutput(rtms)
        zmoc_fmtter.rayParamOutput(rtms[0].track_factory)
        zmoc_fmtter.rtmGeomInfoOutput(rtms)

    def export(self, rtms: list, rtm_maps: list):
        if self.fmt == 'alpha':
            self.__alphaFmtFactory(rtms, rtm_maps)
        elif self.fmt == 'zmoc':
            self.__zmocFmtFactory(rtms, rtm_maps)
