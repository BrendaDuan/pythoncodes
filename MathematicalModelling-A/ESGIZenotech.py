import io
import struct
import math
import numpy
import numpy as np
import gc

import GeometryCollision


# https://pypi.org/project/numpy-stl/
def run_numpy_stl_demo(filename):
    from stl import mesh
    from mpl_toolkits import mplot3d
    from matplotlib import pyplot

    # Using an existing stl file:
    stl_mesh = mesh.Mesh.from_file(filename)

    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    axes.add_collection3d(mplot3d.art3d.Poly3DCollection(stl_mesh.vectors))

    # Auto scale to the mesh size
    scale = stl_mesh.points.flatten()
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def readUInt32(file):
    return struct.unpack('I', file.read(4))[0]


def readFloat32(file):
    return struct.unpack('f', file.read(4))[0]

def readNextFloat32(file):
    float_data = file.read(4)
    if not float_data:
        return None, False
    return struct.unpack('f', float_data)[0], True


# TODO: import from opengl?
class GLTriangleFace(object):
    def __init__(self, *items):
        self.items = items

    def __iter__(self):
        yield from self.items


class FRect(object):
    def __init__(self, origin, widths):
        self.origin = origin
        self.widths = widths


# https://github.com/JuliaGeometry/GeometryBasics.jl/blob/master/src/basic_types.jl
class Mesh(object):
    '''
    Base.getindex(mesh::Mesh, i::Integer) = getfield(mesh, :simplices)[i]
    '''
    def __init__(self, vertices, normals):
        self.vertices = vertices
        self.normals = normals
        self.count = len(vertices) // 3
        # self.faces = faces
        # import pdb
        # pdb.set_trace()

    def __iter__(self):
        for i in range(self.count):
            tri = numpy.array([
                self.vertices[i * 3 + 0],
                self.vertices[i * 3 + 1],
                self.vertices[i * 3 + 2]])
            yield tri


def load_binary_stl(stlfile):
    #Binary STL
    #https://en.wikipedia.org/wiki/STL_%28file_format%29#Binary_STL
    stlfile.read(80) # throw out header

    triangle_count = readUInt32(stlfile)
    # faces = []     # triangle_count
    vertices = []  # triangle_count * 3
    normals = []   # triangle_count * 3

    i = 0
    while True:
        # check eof
        first_float, ok = readNextFloat32(stlfile)
        if not ok:
            break

        # faces.append(GLTriangleFace(i * 3 + 0, i * 3 + 1, i * 3 + 2))
        normal = (first_float, readFloat32(stlfile), readFloat32(stlfile))
        normals.append(normal)
        normals.append(normal) # hurts, but we need per vertex normals
        normals.append(normal)

        # import pdb
        # pdb.set_trace()
        #
        point = numpy.array([readFloat32(stlfile), readFloat32(stlfile), readFloat32(stlfile)])
        vertices.append(point)

        point = numpy.array([readFloat32(stlfile), readFloat32(stlfile), readFloat32(stlfile)])
        vertices.append(point)

        point = numpy.array([readFloat32(stlfile), readFloat32(stlfile), readFloat32(stlfile)])
        vertices.append(point)

        stlfile.read(2) # throwout 16bit attribute
        i += 1

    return Mesh(vertices=vertices, normals=normals)


class AABBModel(object):

    def __init__(self, centre, size, min, max):
        # "centre of the cube"
        self.centre = centre
        # "size of the cube (half the width/length/depth)"
        self.size = size
        # "minimum of bounding box"
        self.min = min
        # "maximum of bounding box"
        self.max = max


def AABB(centre=None, size=None, min=None, max=None):
    # TODO: use dot. and numpy
    if centre is None:
        # vector
        # import pdb
        # pdb.set_trace()
        _centre = 0.5*(min + max)
        _size = 0.5*(max - min)
        return AABBModel(_centre, _size, min, max)
    else:
        _min = centre - size
        _max = centre + size
        return AABBModel(centre, size, _min, _max)


def intersects(c1, c2):
    return np.all(np.abs(c1.centre - c2.centre) < (c1.size + c2.size))
    # import pdb
    # pdb.set_trace()
    # flags = abs(c1.centre - c2.centre) < (c1.size + c2.size)
    # return all(flags)


# Extents [-1500, -1700, 2.125] to [1500, 1600, 292.75]
def bounding_box(tri):
    tri_min = np.array([math.inf, math.inf, math.inf])
    tri_max = np.array([-math.inf, -math.inf, -math.inf])

    # import pdb
    # pdb.set_trace()
    for pt in tri:
        for idx in range(3):
            tri_min[idx] = min(tri_min[idx], pt[idx])
            tri_max[idx] = max(tri_max[idx], pt[idx])

    return AABB(min=tri_min, max=tri_max)


def clamp3(x3, lo, hi):
    nx3 = np.copy(x3)
    for i in range(3):
        if x3[i] < lo:
            nx3[i] = lo
        elif x3[i] > hi:
            nx3[i] = hi
        else:
            nx3[i] = x3[i]

    return nx3

def CartesianIndices3(x, y, z):
    # import pdb
    # pdb.set_trace()
    for i in range(x[0], x[1]+1):
        for j in range(y[0], y[1]+1):
            for k in range(z[0], z[1]+1):
                yield (i, j, k)


def CartesianIndices32(a, b, c):
    yield from CartesianIndices3([0, a], [0, b], [0, c])


def voxelize(stl, cube, depth=8):
    # See http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/tribox_tam.pdf
    n_cubes = 2**depth
    occupied = np.zeros((n_cubes, n_cubes, n_cubes))
    smalldims = 2*cube.size / n_cubes

    for tri in stl:
        # gc.collect()
        bbox = bounding_box(tri)
        if intersects(bbox, cube):
            '''
            min = clamp.(floor.(Int32, 1 .+ 0.5n_cubes .* (bbox.min .- cube.min) ./ cube.size), 1, n_cubes)
            max = clamp.(ceil.(Int32, 1 .+ 0.5n_cubes .* (bbox.max .- cube.min) ./ cube.size), 1, n_cubes)
            for xyz in CartesianIndices((min[1]:max[1], min[2]:max[2], min[3]:max[3]))
                origin = cube.min .+ 2 .* (Tuple(xyz) .- 1) ./ n_cubes .* cube.size
                smallcube = FRect(origin, smalldims)
                occupied[xyz] = GeometryCollision.intersects(smallcube, tri)
            '''
            # 0, 2 --> 1, 3 --> 1, 2, 3 ---> 0, 1, 2
            # import pdb
            # pdb.set_trace()
            a1 = 0 + 0.5*n_cubes * (bbox.min - cube.min) / cube.size
            _min = clamp3(np.floor(a1), 0, n_cubes - 1)
            a2 = 0 + 0.5*n_cubes * (bbox.max - cube.min) / cube.size
            _max = clamp3(np.ceil(a2), 0, n_cubes - 1)
            X = (int(_min[0]), int(_max[0]))
            Y = (int(_min[1]), int(_max[1]))
            Z = (int(_min[2]), int(_max[2]))
            # print(X, Y, Z, _min, _max)
            for x, y, z in CartesianIndices3(X, Y, Z):
                # import pdb
                # pdb.set_trace()
                origin = cube.min + 2 * (np.array([x, y, z]) - 1) / n_cubes * cube.size
                smallcube = FRect(origin, smalldims)
                # print(origin, smalldims)
                occupied[x][y][z] = GeometryCollision.intersects(smallcube, tri)

    return occupied


def vtk_hexahedron(cube):
    x = [cube.widths[0], 0, 0]
    y = [0, cube.widths[1], 0]
    z = [0, 0, cube.widths[2]]

    vertex1 = cube.origin
    vertex2 = vertex1 + x
    vertex3 = vertex2 + y
    vertex4 = vertex1 + y

    vertex5 = vertex1 + z
    vertex6 = vertex5 + x
    vertex7 = vertex6 + y
    vertex8 = vertex5 + y

    vtk = np.zeros((3, 8))
    vtk[:, 0] = vertex1
    vtk[:, 1] = vertex2
    vtk[:, 2] = vertex3
    vtk[:, 3] = vertex4
    vtk[:, 4] = vertex5
    vtk[:, 5] = vertex6
    vtk[:, 6] = vertex7
    vtk[:, 7] = vertex8
    return vtk


def findall3(occupied):
    X, Y, Z = occupied.shape
    for i in range(X):
        for j in range(Y):
            for k in range(Z):
                if occupied[i][j][k]:
                    yield (i, j, k)


def triangles_in_bbox(stl, cube):
    count = 0
    for tri in stl:
        if intersects(bounding_box(tri), cube):
            count += 1
    return count


def enumerate3(indexA, indexB):
    pass

def _bitstream(stream, levels, level, idx_prev, terminator):
    if level > len(levels):
        terminator(stream, idx_prev)
        return stream

    curr = levels[level]
    idx_curr = ((2*idx_prev - [1, 1, 1]), 2*idx_prev)

    value = 0
    for (i, idx) in enumerate(idx_curr):
        if curr[idx]:
            value =  value | (1 << (i-1))

    stream.append(value)
    for (i, idx) in enumerate(idx_curr):
        if curr[idx]:
            _bitstream(stream, levels, level + 1, idx, terminator)

    return stream


def nothing(*args, **kwargs):
    pass


def bitstream(occupied, terminator=nothing):
    occupied_shape = occupied.shape
    assert occupied_shape[0] == occupied_shape[1] == occupied_shape[2]

    depth = int(np.log2(occupied_shape[0]))

    levels = [
        np.zeros((2**level, 2**level, 2**level))
        for level in range(1, depth)
    ]
    levels.append(occupied)

    # import pdb
    # pdb.set_trace()

    '''
    for level in reverse(1:depth-1)
        curr = levels[level]
        prev = levels[level+1]
        for idx_curr in CartesianIndices(size(curr))
            idx_prev = (2idx_curr - CartesianIndex(1, 1, 1)):2idx_curr
            curr[idx_curr] = any(prev[idx_prev])
        end
    '''

    for level in range(depth-2, 0, -1):
        curr = levels[level]
        prev = levels[level+1]
        a, b, c = curr.shape
        # 3, 3, 3
        # 1,2,3 --> 0,1,2-->
        for x, y, z in CartesianIndices32(a-1, b-1, c-1):
            # (1,2,3) --> (2,4,6) --> (1,3,5):(2,4,6)
            # (0, 1, 2) --> (0, 2, 4):(1,3,5),(2,4,6)
            sub_part = prev[2*x:2*x+2, 2*y:2*y+2, 2*z:2*z+2]
            curr[x][y][z] = np.any(sub_part)

    stream = _bitstream([], levels, 1, [0, 0, 0], terminator)
    return stream


def bitstream2(stl, cube, max_length=1, voxelize_depth=10):
    # Calculate number of levels required
    depth = int(np.ceil(np.log2(2*np.max(cube.shape) / max_length)))

    if depth > voxelize_depth:
        rem_depth = depth - voxelize_depth
        n_cubes = 2**rem_depth
        bitstreams = np.zeros((n_cubes, n_cubes, n_cubes))
        # TODO: finish this reusing the mechanisms within _bitstream (use terminator)
    else:
        return bitstream(voxelize(stl, cube, depth))


def _bitstream_tri(stream, levels, full_levels, idx_prev, terminator):
    curr = levels[level]
    full_curr = full_levels[level]
    x, y, z = idx_prev

    # idx_curr = ((2*idx_prev - [1, 1, 1]), 2*idx_prev)
    X = (2*x, 2*x + 1)
    Y = (2*y, 2*y + 1)
    Z = (2*z, 2*z + 1)

    # (1, 1, 1) --> (1,1,1):(2,2,2)
    # (0, 0, 0) --> (0,0,0):(1,1,1)
    # (2, 2, 2) --> (3,3,3):(4,4,4)

    value = 0
    for (i, idx) in enumerate(CartesianIndices32(X, Y, Z)):
        xi, yi, zi = idx
        if curr[xi][yi][zi]:
            value = value | (1 << (i - 1))

    stream.append(value)

    if level < len(levels):
        for (i, idx) in enumerate(CartesianIndices32(X, Y, Z)):
            xi, yi, zi = idx
            if full_curr[xi][yi][zi]:
                stream.append(0)
            elif curr[xi][yi][zi]:
                _bitstream_tri(stream, levels, full_levels, level + 1, idx, terminator)
    else:
        terminator(stream, idx_prev)

    return stream


# 原始 bitstream
def bitstream_tri(occupied, terminator=nothing):
    occupied_shape = occupied.shape
    assert occupied_shape[0] == occupied_shape[1] == occupied_shape[2]

    depth = int(np.log2(occupied_shape[0]))
    levels = [
        np.zeros((2**level, 2**level, 2**level))
        for level in range(1, depth)
    ]

    levels.append(occupied)

    full_levels = [
        np.zeros((2**level, 2**level, 2**level))
        for level in range(1, depth)
    ]
    full_levels.append(occupied)


    for level in range(depth-1, 0, -1):
        curr = levels[level]
        prev = levels[level+1]
        full_curr = full_levels[level]
        full_prev = full_levels[level+1]
        for idx_curr in CartesianIndices3(curr.shape):
            idx_prev = ((2*idx_curr - [1, 1, 1]), 2*idx_curr)
            curr[idx_curr] = any(prev[idx_prev])
            full_curr[idx_curr] = all(full_prev[idx_prev])

    stream = _bitstream_tri([], levels, full_levels, 1, [1, 1, 1], terminator)
    return stream


v1_indexes = [
    (0, 0, 0),
    (0, 0, 1),
    (0, 1, 0),
    (0, 1, 1),
    (1, 0, 0),
    (1, 0, 1),
    (1, 1, 0),
    (1, 1, 1),
]

def bitstream_v1(occupied):
    '''
    cube --> 8*child
    8bit
    '''
    length = occupied.shape[0]
    # 优化手段
    # 当cube 尺寸较小时,先 sum 下，是不是全是0，
    # 这样全都是空白点，就不需要后续流程了
    if length < 16:
        if not np.sum(occupied):
            return [], False

    # 递归退出条件
    # 8个点组成一个 2*2*2 cube --->
    '''
    2*2*2---> return one byte 8bit, flag
    '''
    if length <= 2:
        if not np.sum(occupied):
            return [], False

        value = 0
        # 1. 将 2*2*2 occupied reshape为 8 size 的 array
        # 2. enumerate 遍历，同时拿到索引和值
        for i, bit_value in enumerate(occupied.reshape(length*length*length)):
            # 设置 value 的 i 位 为 1（看起来确实是从右边开始算）
            if bit_value:
                value = value | (1 << i)
        # 如果全为 0, 同样返回空
        if not value:
            return [], False
        return [value], True

    # 分割为 8 个 child
    # split to 8 child
    child_length = length // 2

    # 左右分界点
    # 比如 length 为 32, 则 child_length 为 16
    # 左右分界点为 [0, 16], [16,32]
    left_right_indexes = [
        (0, child_length),
        (child_length, length)
    ]

    # value 为控制位
    value = 0
    # 子 cube 的 stream 列表
    sub_streams = []
    # 从 0 数到 7
    for i in range(8):
        # 根据编号拿到 (0,0,0),(0,0,1),... 这样的索引
        ii, ij, ik = v1_indexes[i]
        # 如果 ii 为 0 说明, X 轴需要取左半边
        # 如果 ij 为 0 说明, Y 轴需要取前半边
        # 如果 ik 为 0 说明, Z 轴需要取下半边
        # left right
        iil, iir = left_right_indexes[ii]
        ijl, ijr = left_right_indexes[ij]
        ikl, ikr = left_right_indexes[ik]
        #import pdb
        #pdb.set_trace()
        # 取得 子 child
        # 一个例子: child = occupied[0: 16, 16: 32, 0: 16]
        # one child
        child = occupied[iil: iir, ijl: ijr, ikl: ikr]
        # 递归调用,得到子 cube 的 stream
        sub_stream, flag = bitstream_v1(child)
        # 如果 子 cube 为空,则跳过
        if not flag:
            continue
        # 同上
        # child has detail
        value = value | (1 << i)
        # add child detail
        sub_streams.append(sub_stream)

    if value:
        # join child stream
        stream = [value]
        last = len(sub_streams) - 1
        print(length, len(sub_streams))
        for i, sub_stream in enumerate(sub_streams):
            if length <= 4:
                if i == 0:
                    stream.append(0)
            # 使用 extend 拼接子 stream
            stream.extend(sub_stream)

            if length > 4:
                # 在中间拼接 0, 如果不是最后一个,则不需要
                if i < last:
                    stream.append(0)
        return stream, True
    return [], False


def bitstream_v2(occupied):
    # TODO
    pass
    '''
    stream = []
    # return (N, 3) array
    points_index =  np.argwhere(occupied > 0)
    points_x = points_index[:, 0]
    points_y = points_index[:, 1]
    points_z = points_index[:, 2]

    min_x = np.min(points_x)
    max_x = np.max(points_x)
    min_y = np.min(points_y)
    max_y = np.max(points_y)
    min_z = np.min(points_z)
    max_z = np.max(points_z)
    '''

'''
def _reconstruct(occupied, idx_prev, depth, max_depth, stream, stream_pos):
    # TODO: make this work for tri state encodings
    value = stream[stream_pos]
    stream_pos += 1
    idx_curr = ((2*idx_prev - [1, 1, 1]), 2*idx_prev)
    if depth == max_depth:
        for (i, idx) in enumerate3(idx_curr):
            if (value & (1 << (i-1))) != 0:
                occupied[idx] = True

    else:
        for (i, idx) in enumerate3(idx_curr):
            if (value & (1 << (i-1))) != 0:
                stream_pos = _reconstruct(occupied, idx, depth+1, max_depth, stream, stream_pos)

    return stream_pos
'''

'''
def reconstruct(stream, max_depth):
    occupied = np.zeros((2**max_depth, 2**max_depth, 2**max_depth))
    _reconstruct(occupied, CartesianIndex(1, 1, 1), 1, max_depth, stream, 0)
    return occupied
'''

def run_demo(files=["canary.stl"], depths=[7, 8, 9, 10], plot=None):
    root = AABB(centre=[0,512,0], size=[512,512,512])

    for file in files:
        stl = load_data(file)
        tri_count = triangles_in_bbox(stl, root)
        for depth in depths:
            occupied = voxelize(stl, root, depth)
            stream = bitstream_tri(occupied)

            if plot:
                plot(f'{file}_{depth}', root, stl, occupied)

            print(f'Done {file} at depth {depth}')


def spy3(occupied):
    idx = np.array(list(findall3(occupied)))

    x = idx[:, 0]
    y = idx[:, 1]
    z = idx[:, 2]
    return (x, y, z)


def dir_vec(i):
    if i == 1:
        return [1, 0, 0]
    if i == 2:
        return [0, 1, 0]

    return [0, 0, 1]


def determine_dir(occupied):
    sz = occupied.shape

    sumx = 0
    for k in range(sz[2]):
        for j in range(sz[1]):
            for i in range(0, sz[0], 2):
                sumx += all(occupied[i, j, k], occupied[i+1, j, k])

    sumy = 0
    for k in range(sz[2]):
        for j in range(0, sz[1], 2):
            for i in range(sz[0]):
                sumy += all(occupied[i, j, k], occupied[i, j+1, k])

    sumz = 0
    for k in range(0, sz[2], 2):
        for j in range(sz[1]):
            for i in range(sz[0]):
                sumz += all(occupied[i, j, k], occupied[i, j, k+1])

    print(sumx, sumy, sumz)
    if (sumx == 0) and (sumy == 0) and(sumz == 0):
        return np.argmax(sz)
    else:
        return np.argmax((sumx, sumy, sumz))


def condense_dir(occupied, dir):
    return condense_dir3(occupied, dir_vec(dir))

def condense_dir3(occupied, dir):
    sz = occupied[0].shape
    new_sz = sz / (1 + dir)
    new_occupied = np.zeros(new_sz)
    new_occupied_all = np.zeros(new_sz)
    for i in CartesianIndices(new_sz):
        ii = CartesianIndex(Tuple(i) * (1 + dir))
        new_occupied[i] = any((occupied[1][ii - CartesianIndex(dir)], occupied[1][ii]))
        new_occupied_all[i] = all((occupied[2][ii - CartesianIndex(dir)], occupied[2][ii]))

    return (new_occupied, new_occupied_all)

def condense(occupied):
    sz = occupied.shape

    depth = np.log2(sz)
    depth_sum = int(sum(depth))

    level = occupied
    level_all = occupied

    levels = []
    dirs = []
    levels.append((level, level))

    for i in range(2, depth_sum):
        dir = determine_dir(level)
        dirs.append(dir)
        (level, level_all) = condense_dir((level, level_all), dir)
        levels.append((level, level_all))

    dirs.append(np.argmax(level.shape))
    return (levels, dirs)


def tree_size(levels, dirs, level, idx_prev):
    if level == 0:
        return 0

    dir = dir_vec(dirs[level])
    curr = levels[level]
    counter = 1
    idx_curr = idx_prev * (1 + dir) - dir
    if curr[2][idx_curr]:
        counter += 1
    elif curr[1][idx_curr]:
        counter += tree_size(levels, dirs, level-1, idx_curr)

    idx_curr += dir
    if curr[2][idx_curr]:
        counter += 1
    elif curr[1][idx_curr]:
        counter += tree_size(levels, dirs, level-1, idx_curr)

    return counter


def tree_size2(levels, dirs):
    return tree_size(levels, dirs, length(levels), [1, 1, 1])
