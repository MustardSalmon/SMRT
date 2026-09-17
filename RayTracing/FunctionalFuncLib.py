# coding=utf-8
# @Author   : Mr.Cao
# @Email    : mrcao@zju.edu.cn
# @Time     : 20/08/2024 16:24
import struct
from dataclasses import dataclass
from typing import List


@dataclass
class Line:
    id: int
    globalID: int
    incident: List
    emergent: List

class OutputPassBottomRayIDFileHandle:
    def __init__(self, _num_azim):
        self._rayMessage = {}
        self.__num_azim = _num_azim

    def genePassBottomLineFileForAzimBinary(self, _fileName, _azimIdx, _rays, _numTrack, _mode):
        """
        generator Binary pass bottom Line file in azim_id
        :param _fileName: file name
        :param _azimIdx: generator azim ID of this file
        :param _rays: number of track pass bottom line in this azim
        :param _numTrack: number of track in this azim
        :param _mode: generator file mode -- wb || ab
        :return: None
        """
        f = open('./' + _fileName + '.dat', _mode)
        _azimIdxBin = struct.pack('i', _azimIdx)
        f.write(_azimIdxBin)
        _raysBin = struct.pack('i', len(_rays))
        f.write(_raysBin)
        _numTrackBin = struct.pack('i', _numTrack)
        f.write(_numTrackBin)

        for i in range(len(_rays)):
            _id = _rays[i].id if _azimIdx < self.__num_azim // 2 else _numTrack - 1 - _rays[i].id
            _ray_id = struct.pack('i', _id)
            f.write(_ray_id)

        f.close()

    def genePassBottomLineFileForAzimDecimal(self, _fileName, _azimIdx, _rays, _numTrack, _mode):
        f = open('./' + _fileName + '.txt', _mode)
        f.write(f'! azimID numTracksPassBottom numTracks\n')
        f.write(f'{_azimIdx} {len(_rays)} {_numTrack} \n')

        f.write(f'! rayID incidentX incidentY emergentX emergentY\n')
        for i in range(len(_rays)):
            _id = _rays[i].id if _azimIdx < self.__num_azim // 2 else _numTrack - 1 - _rays[i].id
            _incident = _rays[i].incident if _azimIdx < self.__num_azim // 2 else _rays[i].emergent
            _emergent = _rays[i].emergent if _azimIdx < self.__num_azim // 2 else _rays[i].incident

            # only generator ray ID
            f.write(str(_id).rjust(8, ' ') + '\n')

            # # generator ray ID and incident, emergent coordinate
            # f.write(str(_id).rjust(8, ' '))
            # f.write(str(format(_incident[0], '.16f')).rjust(20, ' '))
            # f.write(str(format(_incident[1], '.16f')).rjust(20, ' '))
            # f.write(str(format(_emergent[0], '.16f')).rjust(20, ' '))
            # f.write(str(format(_emergent[1], '.16f')).rjust(20, ' ') + '\n')

        f.close()

    def genePassBottomLineFileFactory(self, _fileName, _azimIdx, _rays, _numTrack, _mode):
        if 'b' in _mode:
            self.genePassBottomLineFileForAzimBinary(_fileName, _azimIdx, _rays, _numTrack, _mode)
        else:
            self.genePassBottomLineFileForAzimDecimal(_fileName, _azimIdx, _rays, _numTrack, _mode)

    def readPassBottomLineFileForAzimBinary(self, azim_idx, _numTracksCurrent, _mode):
        _numTracksCurrent = -1

        f = open('./' + str(azim_idx) + '.dat', _mode)
        binaryData = f.read()
        unpackedData = struct.unpack('i' * (len(binaryData) // 4), binaryData)
        f.close()

        _numTracksAzim = unpackedData[2]
        for data in unpackedData[3::]:
            line = Line(data, data + _numTracksCurrent, [], [])
            self._rayMessage[azim_idx].append(line)

        return _numTracksAzim

    def readPassBottomLineFileForAzimDecimal(self, azim_idx, _numTracksCurrent, _mode):
        _numTracksAzim = -1

        f = open('./' + str(azim_idx) + '.txt', _mode)
        lines = f.readlines()

        for line in lines:
            if '!' in line:
                continue
            lineList = line.split()

            if len(lineList) == 3:
                _numTracksAzim = int(lineList[-1])
            else:
                # # read generator ray ID and incident, emergent coordinate file
                # line = Line(int(lineList[0]), int(lineList[0]) + _numTracksCurrent,
                #             [float(lineList[1]), float(lineList[2])], [float(lineList[3]), float(lineList[4])])

                # read only generator ray ID file
                line = Line(int(lineList[0]), int(lineList[0]) + _numTracksCurrent, [], [])

                self._rayMessage[azim_idx].append(line)

        f.close()

        return _numTracksAzim

    def readPassBottomLineFileFactory(self, azim_idx, _numTracksCurrent, _mode):
        if 'b' in _mode:
            _numTracksAzim = self.readPassBottomLineFileForAzimBinary(azim_idx, _numTracksCurrent, _mode)
        else:
            _numTracksAzim = self.readPassBottomLineFileForAzimDecimal(azim_idx, _numTracksCurrent, _mode)

        return _numTracksAzim



@dataclass
class Pin:
    id: int
    row_id: int
    col_id: int
    surface: List
    face: List


