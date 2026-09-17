import math

import numpy as np

from RayTracing.MacroConfig import DEFAULT_BOUNDARY_OVERLAPPED_CRIT, SQRT_3, SMALL_VALUE, CASE_NAME, WARNING
from RayTracing import RunTime as rt


def getBoundaryPointsIndexInRect(x_coords: np.ndarray, y_coords: np.ndarray, bounding_box: tuple):
    x_min, x_max, y_min, y_max = bounding_box

    x_coords_test = x_coords - x_min
    pts_idxes_x_min, = np.where((x_coords_test >= 0) & (x_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))

    x_coords_test = x_max - x_coords
    pts_idxes_x_max, = np.where((x_coords_test >= 0) & (x_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))

    y_coords_test = y_coords - y_min
    pts_idxes_y_min, = np.where((y_coords_test >= 0) & (y_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))

    y_coords_test = y_max - y_coords
    pts_idxes_y_max, = np.where((y_coords_test >= 0) & (y_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))

    bdry_pts_idxes = np.concatenate((pts_idxes_x_min, pts_idxes_x_max, pts_idxes_y_min, pts_idxes_y_max))
    bdry_pts_idxes = np.unique(bdry_pts_idxes)

    return bdry_pts_idxes


def getBoundaryPointsIndexInHex(x_coords: np.ndarray, y_coords: np.ndarray, bounding_box: tuple):
    x_min, x_max, y_min, y_max = bounding_box
    hex_x = (x_min + x_max) * 0.5
    hex_y = (y_min + y_max) * 0.5
    hex_edge_width = (x_max - x_min) * 0.5

    # classify the inscribed circle
    vec_modes_sq = np.square(x_coords - hex_x) + np.square(y_coords - hex_y)
    # minus the SMALL_VALUE to include the boundary node overlapped with the inscribed radius
    classified_pts_idxes, = np.where(vec_modes_sq >= ((hex_edge_width * SQRT_3 / 2) ** 2 - SMALL_VALUE))
    x_coords_classified = x_coords[classified_pts_idxes]
    y_coords_classified = y_coords[classified_pts_idxes]
    vec_modes_sq_classified = vec_modes_sq[classified_pts_idxes]

    # classify by azimuthal angle and the obtained index is relative with inscribed circle classified points
    with np.errstate(divide='ignore'):
        azim_angs = np.arctan((y_coords_classified - hex_y) / (x_coords_classified - hex_x))
    pts_idxes_0_60, = np.where((0 <= azim_angs) & (azim_angs < np.pi / 3))
    pts_idxes_60_0, = np.where((-np.pi / 3 <= azim_angs) & (azim_angs < 0))
    pts_idxes_60_90, = np.where((azim_angs > np.pi / 3) | (azim_angs < -np.pi / 3))

    # find boundary point indexes in upper and bottom
    # the obtained index is relative with the azimuthal angle classified points
    y_coords_60_90 = y_coords_classified[pts_idxes_60_90]
    y_coords_test = y_coords_60_90 - y_min
    pts_idxes_y_min, = np.where((y_coords_test >= 0) & (y_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))
    pts_idxes_y_min = pts_idxes_60_90[pts_idxes_y_min]
    pts_idxes_y_min = classified_pts_idxes[pts_idxes_y_min]

    y_coords_test = y_max - y_coords_60_90
    pts_idxes_y_max, = np.where((y_coords_test >= 0) & (y_coords_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT))
    pts_idxes_y_max = pts_idxes_60_90[pts_idxes_y_max]
    pts_idxes_y_max = classified_pts_idxes[pts_idxes_y_max]

    # find boundary point indexes between 0~60
    vec_modes_0_60 = np.sqrt(vec_modes_sq_classified[pts_idxes_0_60])
    azim_angs_0_60 = azim_angs[pts_idxes_0_60]
    vec_modes_0_60_project = vec_modes_0_60 * np.cos(azim_angs_0_60 - np.pi / 6)
    vec_modes_test = (hex_edge_width * SQRT_3 / 2) - vec_modes_0_60_project
    cond = (-DEFAULT_BOUNDARY_OVERLAPPED_CRIT < vec_modes_test) & (vec_modes_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT)
    pts_idxes_0_60_bdry, = np.where(cond)
    pts_idxes_0_60_bdry = pts_idxes_0_60[pts_idxes_0_60_bdry]
    pts_idxes_0_60_bdry = classified_pts_idxes[pts_idxes_0_60_bdry]

    # find boundary point indexes between -60~0
    vec_modes_60_0 = np.sqrt(vec_modes_sq_classified[pts_idxes_60_0])
    azim_angs_60_0 = azim_angs[pts_idxes_60_0]
    vec_modes_60_0_project = vec_modes_60_0 * np.cos(azim_angs_60_0 + np.pi / 6)
    vec_modes_test = (hex_edge_width * SQRT_3 / 2) - vec_modes_60_0_project
    cond = (-DEFAULT_BOUNDARY_OVERLAPPED_CRIT < vec_modes_test) & (vec_modes_test < DEFAULT_BOUNDARY_OVERLAPPED_CRIT)
    pts_idxes_60_0_bdry, = np.where(cond)
    pts_idxes_60_0_bdry = pts_idxes_60_0[pts_idxes_60_0_bdry]
    pts_idxes_60_0_bdry = classified_pts_idxes[pts_idxes_60_0_bdry]

    bdry_pts_idxes = np.concatenate((pts_idxes_y_min, pts_idxes_y_max,
                                     pts_idxes_0_60_bdry, pts_idxes_60_0_bdry))
    bdry_pts_idxes = np.unique(bdry_pts_idxes)

    return bdry_pts_idxes


def getPointsIndexInTriSectOfHex(x_coords: np.ndarray, y_coords: np.ndarray, bounding_box: tuple):
    x_min, x_max, y_min, y_max = bounding_box
    hex_x = (x_min + x_max) * 0.5
    hex_y = (y_min + y_max) * 0.5

    # control function
    # y < sqrt(3) * (x - hex_x) + hex_y
    # y < -sqrt(3) * (x - hex_x) + hex_y
    # y > y_min

    # y - sqrt(3) * x < hex_y - sqrt(3) * hex_x
    # y + sqrt(3) * x < hex_y + sqrt(3) * hex_x
    # y > y_min

    y_minus_sqrt3_x = y_coords - SQRT_3 * x_coords
    y_plus_sqrt3_x = y_coords + SQRT_3 * x_coords

    classified_pts_idxes_1, = np.where(y_minus_sqrt3_x < hex_y - SQRT_3 * hex_x + SMALL_VALUE)
    classified_pts_idxes_2, = np.where(y_plus_sqrt3_x < hex_y + SQRT_3 * hex_x + SMALL_VALUE)
    classified_pts_idxes_3, = np.where(y_coords > y_min - SMALL_VALUE)
    all_pts_idxes = np.intersect1d(classified_pts_idxes_1, classified_pts_idxes_2)
    all_pts_idxes = np.intersect1d(all_pts_idxes, classified_pts_idxes_3)

    cond1 = np.abs(y_minus_sqrt3_x - (hex_y - SQRT_3 * hex_x)) < SMALL_VALUE
    cond2 = np.abs(y_plus_sqrt3_x - (hex_y + SQRT_3 * hex_x)) < SMALL_VALUE
    cond3 = np.abs(y_coords - y_min) < SMALL_VALUE
    classified_pts_idxes_4, = np.where(cond1)
    classified_pts_idxes_5, = np.where(cond2)
    classified_pts_idxes_6, = np.where(cond3)
    bdry_pts_idxes_1 = np.intersect1d(all_pts_idxes, classified_pts_idxes_4)
    bdry_pts_idxes_2 = np.intersect1d(all_pts_idxes, classified_pts_idxes_5)
    bdry_pts_idxes_3 = np.intersect1d(all_pts_idxes, classified_pts_idxes_6)
    # bdry_pts_idxes_3, = np.where(cond3)
    bdry_pts_idxes = np.concatenate((bdry_pts_idxes_1, bdry_pts_idxes_2, bdry_pts_idxes_3))

    return all_pts_idxes, bdry_pts_idxes


def getPointsIndexInHex(x_coords: np.ndarray, y_coords: np.ndarray, bounding_box: tuple):
    x_min, x_max, y_min, y_max = bounding_box
    hex_x = (x_min + x_max) * 0.5
    hex_y = (y_min + y_max) * 0.5

    # control function
    # y - hex_y < -sqrt(3) * (x - x_max)
    # y - hex_y > -sqrt(3) * (x - x_min)
    # y - hex_y < sqrt(3) * (x - x_min)
    # y - hex_y > sqrt(3) * (x - x_max)
    # y < y_max
    # y > y_min

    # y + sqrt(3) * x < hex_y + sqrt(3) * x_max
    # y + sqrt(3) * x > hex_y + sqrt(3) * x_min
    # y - sqrt(3) * x < hex_y - sqrt(3) * x_min
    # y - sqrt(3) * x > hex_y - sqrt(3) * x_max
    # y < y_max
    # y > y_min

    y_minus_sqrt3_x = y_coords - SQRT_3 * x_coords
    y_plus_sqrt3_x = y_coords + SQRT_3 * x_coords

    classified_pts_idxes_1, = np.where(y_plus_sqrt3_x < hex_y + SQRT_3 * x_max + SMALL_VALUE)
    classified_pts_idxes_2, = np.where(y_plus_sqrt3_x > hex_y + SQRT_3 * x_min - SMALL_VALUE)
    classified_pts_idxes_3, = np.where(y_minus_sqrt3_x < hex_y - SQRT_3 * x_min + SMALL_VALUE)
    classified_pts_idxes_4, = np.where(y_minus_sqrt3_x > hex_y - SQRT_3 * x_max - SMALL_VALUE)
    classified_pts_idxes_5, = np.where(y_coords > y_min - SMALL_VALUE)
    classified_pts_idxes_6, = np.where(y_coords < y_max + SMALL_VALUE)

    all_pts_idxes = np.intersect1d(classified_pts_idxes_1, classified_pts_idxes_2)
    all_pts_idxes = np.intersect1d(all_pts_idxes, classified_pts_idxes_3)
    all_pts_idxes = np.intersect1d(all_pts_idxes, classified_pts_idxes_4)
    all_pts_idxes = np.intersect1d(all_pts_idxes, classified_pts_idxes_5)
    all_pts_idxes = np.intersect1d(all_pts_idxes, classified_pts_idxes_6)

    cond1 = np.abs(y_plus_sqrt3_x - (hex_y + SQRT_3 * x_max)) < SMALL_VALUE
    cond2 = np.abs(y_plus_sqrt3_x - (hex_y + SQRT_3 * x_min)) < SMALL_VALUE
    cond3 = np.abs(y_minus_sqrt3_x - (hex_y - SQRT_3 * x_min)) < SMALL_VALUE
    cond4 = np.abs(y_minus_sqrt3_x - (hex_y - SQRT_3 * x_max)) < SMALL_VALUE
    cond5 = np.abs(y_coords - y_min) < SMALL_VALUE
    cond6 = np.abs(y_coords - y_max) < SMALL_VALUE

    classified_pts_idxes_11, = np.where(cond1)
    classified_pts_idxes_12, = np.where(cond2)
    classified_pts_idxes_13, = np.where(cond3)
    classified_pts_idxes_14, = np.where(cond4)
    classified_pts_idxes_15, = np.where(cond5)
    classified_pts_idxes_16, = np.where(cond6)
    bdry_pts_idxes_1 = np.intersect1d(all_pts_idxes, classified_pts_idxes_11)
    bdry_pts_idxes_2 = np.intersect1d(all_pts_idxes, classified_pts_idxes_12)
    bdry_pts_idxes_3 = np.intersect1d(all_pts_idxes, classified_pts_idxes_13)
    bdry_pts_idxes_4 = np.intersect1d(all_pts_idxes, classified_pts_idxes_14)
    bdry_pts_idxes_5 = np.intersect1d(all_pts_idxes, classified_pts_idxes_15)
    bdry_pts_idxes_6 = np.intersect1d(all_pts_idxes, classified_pts_idxes_16)
    bdry_pts_idxes = np.concatenate(
        (bdry_pts_idxes_1, bdry_pts_idxes_2, bdry_pts_idxes_3, bdry_pts_idxes_4, bdry_pts_idxes_5, bdry_pts_idxes_6))

    return all_pts_idxes, bdry_pts_idxes


def calcIntercept(x0: float, y0: float, k: float) -> float:
    if k == np.inf or k == math.inf:
        return x0
    else:
        return y0 - k * x0


def calc2DVectorMode(x0: float, y0: float, x1: float, y1: float) -> float:
    return math.sqrt((y1 - y0) ** 2 + (x1 - x0) ** 2)


def isPointInPoly(num_vert: int, vert_x: list, vert_y: list, test_x: float, test_y: float) -> bool:
    '''
    PNPoly algorithm
    I DON'T UNDERSTAND
    '''
    if ((test_x < min(vert_x)) or (test_x > max(vert_x)) or (test_y < min(vert_y)) or (test_y > max(vert_y))):
        return False
    else:
        is_in = False
        i = 0
        j = num_vert - 1

        while i < num_vert:
            if (((vert_y[i] > test_y) != (vert_y[j] > test_y)) and (
                    test_x < (vert_x[j] - vert_x[i]) * (test_y - vert_y[i]) / (vert_y[j] - vert_y[i]) + vert_x[i])):
                is_in = not is_in

            j = i
            i += 1

        return is_in


def vectorDot2D(v0: tuple, v1: tuple) -> float:
    return v0[0] * v1[0] + v0[1] * v1[1]


def vectorCross2D(v0: tuple, v1: tuple) -> float:
    # A x B = |x0  y0| = x0y1 - x1y0
    #         |x1  y1|
    return v0[0] * v1[1] - v1[0] * v0[1]


def calcLineIntersection(p0: tuple, p1: tuple, q0: tuple, q1: tuple) -> tuple:
    '''
    calc the intersection
            / q+s
           /
    p ----*---- p+r
         /
        /
       q

    the expr of intersection can be written as
    p + tr = q + us
    while t and u are scalar parameters

    solve t:
        (p + tr) x s = (q + us) x s
                   t = (q - p) x s / (r x s)

    solve u:
        (p + tr) x r = (q + us) x r
                   u = (p - q) x r / (s x r)
                   u = (q - p) x r / (r x s)

    case1:
        if r x s = 0 and (q - p) x r = 0, the two lines are collinear

    case2:
        if r x s = 0 and (q - p) x r != 0, the two lines are parallel and non-intersecting

    case3:
        if r x s != 0 and 0 <= t <=1 and 0 <= u <= 1, the two lines segments meet at the point p + tr  = q + us
        return p + tr

    case4:
        otherwise, the two lines are not parallel but do not intersect
    '''

    def pointSubstract(c0: tuple, c1: tuple):
        return (c1[0] - c0[0], c1[1] - c0[1])

    p = p0
    r = pointSubstract(p0, p1)
    q = q0
    s = pointSubstract(q0, q1)

    q2p = pointSubstract(p, q)
    rxs = vectorCross2D(r, s)

    q2pxs = vectorCross2D(q2p, s)
    q2pxr = vectorCross2D(q2p, r)

    t = q2pxs / rxs
    u = q2pxr / rxs

    if abs(rxs) < 1e-8:
        if abs(q2pxr) < 1e-8:
            rt.printMessage(CASE_NAME, WARNING, "meet the collinear in handling intersection, ignore the result")
            return ('collinear', -9999., -9999.)
        else:
            rt.printMessage(CASE_NAME, WARNING, "meet the parallel in handling intersection, ignore the result")
            return ('parallel', -9999., -9999.)

    elif abs(q2pxr) > 1e-8 and (-1e-7 <= t <= 1. + 1e-7) and (-1e-7 <= u <= 1. + 1e-7):
        return ('intersect', p[0] + t * r[0], p[1] + t * r[1])

    else:
        return ('nointersect', -9999., -9999.)
