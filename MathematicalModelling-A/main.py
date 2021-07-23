import sys
import io
import numpy as np
import time

from ESGIZenotech import load_binary_stl
from ESGIZenotech import AABB
from ESGIZenotech import voxelize
from ESGIZenotech import bitstream
from ESGIZenotech import run_demo
from ESGIZenotech import run_numpy_stl_demo
from ESGIZenotech import bitstream_v1
from ESGIZenotech import bitstream_v2


def run_local(filename):
    s = time.time()
    with io.open(filename, 'rb') as stl_f:
        stl = load_binary_stl(stl_f)

    print(stl)
    print(stl.count)
    # todo plot stl
    print("load time:", time.time() - s)
    s = time.time()

    root = AABB(centre=np.array([0, 0, 0]), size=np.array([512, 512, 512]))
    occupied = voxelize(stl, root, 10)

    print("occupied time:", time.time() - s)
    s = time.time()

    # import pdb
    # pdb.set_trace()
    st = bitstream(occupied)
    # stream size
    print(len(st))
    print("bitstream time:", time.time() - s)


# 优化后的 bitstream
def run_bitstream_v1(filename, level=10):
    s = time.time()
    with io.open(filename, 'rb') as stl_f:
        stl = load_binary_stl(stl_f)

    print("triangle count", stl.count)
    # todo plot stl
    print("load time: ", time.time() - s, "s")
    s = time.time()

    root = AABB(centre=np.array([0, 0, 0]), size=np.array([512, 512, 512]))
    occupied = voxelize(stl, root, level)
    length = occupied.shape[0]

    print("cube level: ", level)
    print("cube size: ", length*length*length)
    print("cube point count: ", sum(sum(sum(occupied))))
    print("calc occupied time: ", time.time() - s, "s")
    s = time.time()

    #import pdb
    #pdb.set_trace()
    stream, flag = bitstream_v1(occupied)
    # stream size
    print("after compress bitstream size: ", len(stream))
    print("calc bitstream_v1 time:", time.time() - s, "s")

    # show cube
    show_cube(occupied)


# 优化后的 bitstream
def run_bitstream_v2(filename, level=10):
    s = time.time()
    with io.open(filename, 'rb') as stl_f:
        stl = load_binary_stl(stl_f)

    print("triangle count", stl.count)
    # todo plot stl
    print("load time: ", time.time() - s, "s")
    s = time.time()

    root = AABB(centre=np.array([0, 0, 0]), size=np.array([512, 512, 512]))
    occupied = voxelize(stl, root, level)
    length = occupied.shape[0]

    print("cube level: ", level)
    print("cube size: ", length*length*length)
    print("cube point count: ", sum(sum(sum(occupied))))
    print("calc occupied time: ", time.time() - s, "s")
    s = time.time()

    # import pdb
    # pdb.set_trace()
    stream, flag = bitstream_v2(occupied)
    # stream size
    print("after compress bitstream size: ", len(stream))
    print("calc bitstream_v1 time:", time.time() - s, "s")

    show_cube(occupied)

    import pdb
    pdb.set_trace()
    print('ok')


def show_cube(occupied):
    # show cube
    from matplotlib import pyplot
    from mpl_toolkits.mplot3d import Axes3D

    length = occupied.shape[0]
    points_index =  np.argwhere(occupied > 0)

    X = np.arange(length)
    Y = np.arange(length)
    Z = np.zeros((length, length))
    for pi, pj, pk in points_index:
        Z[pi][pj] = max(Z[pi][pj], pk)

    xx, yy = np.meshgrid(X, Y)
    X, Y = xx.ravel(), yy.ravel()
    dz = Z.ravel()
    # Construct arrays with the dimensions for the 16 bars.
    dx = dy = 1

    #
    dz_index = np.argwhere(dz>0)
    X1 = X[dz_index][:, 0]
    Y1 = Y[dz_index][:, 0]
    dz1 = dz[dz_index][:, 0]

    dz2 = (dz1 - length // 2)

    figure = pyplot.figure()
    axes = figure.gca(projection='3d')
    axes.bar3d(X1, Y1, 0, dx, dy, dz2, zsort='min')
    axes.set_xlabel('X')
    axes.set_ylabel('Y')
    axes.set_zlabel('Z(value)')

    pyplot.show()


def run_show_demo():
    from stl import mesh
    import math
    import numpy

    # Create 3 faces of a cube
    data = numpy.zeros(6, dtype=mesh.Mesh.dtype)

    # Top of the cube
    data['vectors'][0] = numpy.array([[0, 1, 1],
                                      [1, 0, 1],
                                      [0, 0, 1]])
    data['vectors'][1] = numpy.array([[1, 0, 1],
                                      [0, 1, 1],
                                      [1, 1, 1]])
    # Front face
    data['vectors'][2] = numpy.array([[1, 0, 0],
                                      [1, 0, 1],
                                      [1, 1, 0]])
    data['vectors'][3] = numpy.array([[1, 1, 1],
                                      [1, 0, 1],
                                      [1, 1, 0]])
    # Left face
    data['vectors'][4] = numpy.array([[0, 0, 0],
                                      [1, 0, 0],
                                      [1, 0, 1]])
    data['vectors'][5] = numpy.array([[0, 0, 0],
                                      [0, 0, 1],
                                      [1, 0, 1]])

    # Since the cube faces are from 0 to 1 we can move it to the middle by
    # substracting .5
    data['vectors'] -= .5

    # Generate 4 different meshes so we can rotate them later
    meshes = [mesh.Mesh(data.copy()) for _ in range(4)]

    # Rotate 90 degrees over the Y axis
    meshes[0].rotate([0.0, 0.5, 0.0], math.radians(90))

    # Translate 2 points over the X axis
    meshes[1].x += 2

    # Rotate 90 degrees over the X axis
    meshes[2].rotate([0.5, 0.0, 0.0], math.radians(90))
    # Translate 2 points over the X and Y points
    meshes[2].x += 2
    meshes[2].y += 2

    # Rotate 90 degrees over the X and Y axis
    meshes[3].rotate([0.5, 0.0, 0.0], math.radians(90))
    meshes[3].rotate([0.0, 0.5, 0.0], math.radians(90))
    # Translate 2 points over the Y axis
    meshes[3].y += 2


    # Optionally render the rotated cube faces
    from matplotlib import pyplot
    from mpl_toolkits import mplot3d

    # Create a new plot
    figure = pyplot.figure()
    axes = mplot3d.Axes3D(figure)

    # Render the cube faces
    for m in meshes:
        # import pdb
        # pdb.set_trace()
        axes.add_collection3d(mplot3d.art3d.Poly3DCollection(m.vectors))

    # Auto scale to the mesh size
    scale = numpy.concatenate([m.points for m in meshes]).flatten()
    axes.auto_scale_xyz(scale, scale, scale)

    # Show the plot to the screen
    pyplot.show()


def main():
    if sys.argv[1] == 'showdemo':
        run_show_demo()
        return
    elif sys.argv[1] == 'numpystl':
        filename = sys.argv[2]
        # 使用 numpy-stl 库
        run_numpy_stl_demo(filename)
    elif sys.argv[1] == 'local':
        filename = sys.argv[2]
        run_local(filename)
    elif sys.argv[1] == 'demo':
        filename = sys.argv[2]
        run_demo(files=[filename])
    elif sys.argv[1] == 'bitv1':
        filename = sys.argv[2]
        # size (1/2/4/8)
        length = int(sys.argv[3])
        # default is 1
        level = 10 - int(np.log2(length))
        run_bitstream_v1(filename, level)
    elif sys.argv[1] == 'bitv2':
        filename = sys.argv[2]
        # size (1/2/4/8)
        length = int(sys.argv[3])
        # default is 1
        level = 10 - int(np.log2(length))
        run_bitstream_v2(filename, level)


if __name__ == '__main__':
    main()
