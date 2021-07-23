import numpy as np


def _planeBoxOverlap(normal, vert, maxbox):
    '''
    vmin = np.zeros(3)
    vmax = np.zeros(3)

    for q in range(3):
        v = vert[q]
        if normal[q] > 0:
            vmin[q] = -maxbox[q] - v
            vmax[q] = maxbox[q] - v
        else:
            vmin[q] = maxbox[q] - v
            vmax[q] = -maxbox[q] - v
    '''
    a = maxbox - vert
    b = -maxbox - vert

    vmin = np.where(normal > 0, b, a)
    vmax = np.where(normal > 0, a, b)

    if normal.dot(vmin) > 0:
        return False

    if normal.dot(vmax) >= 0:
        return True

    return False

# X-tests
def _axistest_X01(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p0 = a * v0[1] - b * v0[2]
    p2 = a * v2[1] - b * v2[2]
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2]

    if p0 > p2:
        p2, p0 = p0, p2

    return (p0 > rad) or (p2 < -rad)


def _axistest_X2(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p0 = a * v0[1] - b * v0[2]
    p1 = a * v1[1] - b * v1[2]
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2]

    if p0 > p1:
        p1, p0 = p0, p1

    return (p0 > rad) or (p1 < -rad)

# Y-tests
def _axistest_Y02(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p0 = -a * v0[0] + b * v0[2]
    p2 = -a * v2[0] + b * v2[2]
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2]
    if p0 > p2:
        p2, p0 = p0, p2

    return (p0 > rad) or (p2 < -rad)


def _axistest_Y1(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p0 = -a * v0[0] + b * v0[2]
    p1 = -a * v1[0] + b * v1[2]
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[2]
    if p0 > p1:
        p1, p0 = p0, p1

    return (p0 > rad) or (p1 < -rad)

# Z-tests
def _axistest_Z12(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p1 = a * v1[0] - b * v1[1]
    p2 = a * v2[0] - b * v2[1]
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1]
    if p2 > p1:
        p1, p2 = p2, p1

    return (p2 > rad) or (p1 < -rad)


def _axistest_Z0(a, b, fa, fb, v0, v1, v2, boxhalfsize):
    p0 = a * v0[0] - b * v0[1]
    p1 = a * v1[0] - b * v1[1]
    rad = fa * boxhalfsize[0] + fb * boxhalfsize[1]
    if p0 > p1:
        p1, p0 = p0, p1

    return (p0 > rad) or (p1 < -rad)


def intersects(aabb, tri):
    '''
    AABB-Triangle intersection code taken from
    http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
    '''
    '''
    # use separating axis theorem to test overlap between triangle and box
    # need to test for overlap in these directions:
    # 1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
    #    we do not even need to test these)
    # 2) normal of the triangle
    # 3) crossproduct(edge from tri, {x,y,z}-directin)
    #    this gives 3x3=9 more tests
    '''

    boxhalfsize = aabb.widths * 0.5
    boxcenter = aabb.origin + boxhalfsize

    # move everything so that the boxcenter is in (0,0,0)
    v0 = tri[0] - boxcenter
    v1 = tri[1] - boxcenter
    v2 = tri[2] - boxcenter

    # compute triangle edges
    e0 = v1 - v0
    e1 = v2 - v1
    e2 = v0 - v2

    # Bullet 3:
    # test the 9 tests first (this was faster)
    fex = abs(e0[0])
    fey = abs(e0[1])
    fez = abs(e0[2])
    if _axistest_X01(e0[2], e0[1], fez, fey, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Y02(e0[2], e0[0], fez, fex, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Z12(e0[1], e0[0], fey, fex, v0, v1, v2, boxhalfsize):
        return False

    fex = abs(e1[0])
    fey = abs(e1[1])
    fez = abs(e1[2])
    if _axistest_X01(e1[2], e1[1], fez, fey, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Y02(e1[2], e1[0], fez, fex, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Z0(e1[1], e1[0], fey, fex, v0, v1, v2, boxhalfsize):
        return False

    fex = abs(e2[0])
    fey = abs(e2[1])
    fez = abs(e2[2])
    if _axistest_X2(e2[2], e2[1], fez, fey, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Y1(e2[2], e2[0], fez, fex, v0, v1, v2, boxhalfsize):
        return False
    if _axistest_Z12(e2[1], e2[0], fey, fex, v0, v1, v2, boxhalfsize):
        return False

    # Bullet 1:
    #  first test overlap in the {x,y,z}-directions
    #  find min, max of the triangle each direction, and test for overlap in
    #  that direction -- this is equivalent to testing a minimal AABB around
    #  the triangle against the AABB

    # test in X-direction
    min = np.min((v0[0], v1[0], v2[0]))
    max = np.max((v0[0], v1[0], v2[0]))
    if (min > boxhalfsize[0]) or (max < -boxhalfsize[0]):
        return False

    # test in Y-direction
    min = np.min((v0[1], v1[1], v2[1]))
    max = np.max((v0[1], v1[1], v2[1]))
    if (min > boxhalfsize[1]) or (max < -boxhalfsize[1]):
        return False

    # test in Z-direction
    min = np.min((v0[2], v1[2], v2[2]))
    max = np.max((v0[2], v1[2], v2[2]))
    if (min > boxhalfsize[2]) or (max < -boxhalfsize[2]):
        return False

    # Bullet 2:
    #  test if the box intersects the plane of the triangle
    #  compute plane equation of triangle: normal*x+d=0
    normal = np.cross(e0, e1)

    if not _planeBoxOverlap(normal, v0, boxhalfsize):
        return False

    return True  # box and triangle overlaps
