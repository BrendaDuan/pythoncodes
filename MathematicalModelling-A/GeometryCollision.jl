module GeometryCollision

using GeometryBasics: Rect, Ngon
using LinearAlgebra: dot, cross
using StaticArrays: MVector

function _planeBoxOverlap(normal, vert, maxbox)
    vmin = zero(MVector{3,eltype(maxbox)})
    vmax = zero(MVector{3,eltype(maxbox)})
    for q in Base.OneTo(3)
        v = vert[q]
        if normal[q] > 0
            vmin[q] = -maxbox[q] - v
            vmax[q] = maxbox[q] - v
        else
            vmin[q] = maxbox[q] - v
            vmax[q] = -maxbox[q] - v
        end
    end
    if dot(normal, vmin) > 0
        return false
    end
    if dot(normal, vmax) >= 0
        return true
    end
    return false
end

# X-tests

function _axistest_X01(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p0 = a * v0[2] - b * v0[3]
    p2 = a * v2[2] - b * v2[3]
    if p0 < p2
        min = p0
        max = p2
    else
        min = p2
        max = p0
    end
    rad = fa * boxhalfsize[2] + fb * boxhalfsize[3]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

function _axistest_X2(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p0 = a * v0[2] - b * v0[3]
    p1 = a * v1[2] - b * v1[3]
    if p0 < p1
        min = p0
        max = p1
    else
        min = p1
        max = p0
    end
    rad = fa * boxhalfsize[2] + fb * boxhalfsize[3]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

# Y-tests

function _axistest_Y02(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p0 = -a * v0[1] + b * v0[3]
    p2 = -a * v2[1] + b * v2[3]
    if p0 < p2
        min = p0
        max = p2
    else
        min = p2
        max = p0
    end
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[3]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

function _axistest_Y1(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p0 = -a * v0[1] + b * v0[3]
    p1 = -a * v1[1] + b * v1[3]
    if p0 < p1
        min = p0
        max = p1
    else
        min = p1
        max = p0
    end
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[3]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

# Z-tests

function _axistest_Z12(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p1 = a * v1[1] - b * v1[2]
    p2 = a * v2[1] - b * v2[2]
    if p2 < p1
        min = p2
        max = p1
    else
        min = p1
        max = p2
    end
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

function _axistest_Z0(a, b, fa, fb, v0, v1, v2, boxhalfsize)
    p0 = a * v0[1] - b * v0[2]
    p1 = a * v1[1] - b * v1[2]
    if p0 < p1
        min = p0
        max = p1
    else
        min = p1
        max = p0
    end
    rad = fa * boxhalfsize[1] + fb * boxhalfsize[2]
    if (min > rad) || (max < -rad)
        return true
    else
        return false
    end
end

"""
    intersects(aabb::Rect{3, T}, tri::Ngon{3, T, 3}) where T

AABB-Triangle intersection code taken from
http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/
"""
function intersects(aabb::Rect{3,T}, tri::Ngon{3,T,3}) where {T}

    # use separating axis theorem to test overlap between triangle and box
    # need to test for overlap in these directions:
    # 1) the {x,y,z}-directions (actually, since we use the AABB of the triangle
    #    we do not even need to test these)
    # 2) normal of the triangle
    # 3) crossproduct(edge from tri, {x,y,z}-directin)
    #    this gives 3x3=9 more tests

    boxhalfsize = aabb.widths * T(0.5)
    boxcenter = aabb.origin + boxhalfsize

    # move everything so that the boxcenter is in (0,0,0)
    v0 = tri[1] - boxcenter
    v1 = tri[2] - boxcenter
    v2 = tri[3] - boxcenter

    # compute triangle edges
    e0 = v1 - v0
    e1 = v2 - v1
    e2 = v0 - v2

    # Bullet 3:
    # test the 9 tests first (this was faster)
    fex = abs(e0[1])
    fey = abs(e0[2])
    fez = abs(e0[3])
    _axistest_X01(e0[3], e0[2], fez, fey, v0, v1, v2, boxhalfsize) && return false
    _axistest_Y02(e0[3], e0[1], fez, fex, v0, v1, v2, boxhalfsize) && return false
    _axistest_Z12(e0[2], e0[1], fey, fex, v0, v1, v2, boxhalfsize) && return false

    fex = abs(e1[1])
    fey = abs(e1[2])
    fez = abs(e1[3])
    _axistest_X01(e1[3], e1[2], fez, fey, v0, v1, v2, boxhalfsize) && return false
    _axistest_Y02(e1[3], e1[1], fez, fex, v0, v1, v2, boxhalfsize) && return false
    _axistest_Z0(e1[2], e1[1], fey, fex, v0, v1, v2, boxhalfsize) && return false

    fex = abs(e2[1])
    fey = abs(e2[2])
    fez = abs(e2[3])
    _axistest_X2(e2[3], e2[2], fez, fey, v0, v1, v2, boxhalfsize) && return false
    _axistest_Y1(e2[3], e2[1], fez, fex, v0, v1, v2, boxhalfsize) && return false
    _axistest_Z12(e2[2], e2[1], fey, fex, v0, v1, v2, boxhalfsize) && return false

    # Bullet 1:
    #  first test overlap in the {x,y,z}-directions
    #  find min, max of the triangle each direction, and test for overlap in
    #  that direction -- this is equivalent to testing a minimal AABB around
    #  the triangle against the AABB

    # test in X-direction
    min = minimum((v0[1], v1[1], v2[1]))
    max = maximum((v0[1], v1[1], v2[1]))
    if (min > boxhalfsize[1]) || (max < -boxhalfsize[1])
        return false
    end

    # test in Y-direction
    min = minimum((v0[2], v1[2], v2[2]))
    max = maximum((v0[2], v1[2], v2[2]))
    if (min > boxhalfsize[2]) || (max < -boxhalfsize[2])
        return false
    end

    # test in Z-direction
    min = minimum((v0[3], v1[3], v2[3]))
    max = maximum((v0[3], v1[3], v2[3]))
    if (min > boxhalfsize[3]) || (max < -boxhalfsize[3])
        return false
    end

    # Bullet 2:
    #  test if the box intersects the plane of the triangle
    #  compute plane equation of triangle: normal*x+d=0
    normal = cross(e0, e1)

    if !_planeBoxOverlap(normal, v0, boxhalfsize)
        return false
    end

    return true  # box and triangle overlaps
end

end  # module
