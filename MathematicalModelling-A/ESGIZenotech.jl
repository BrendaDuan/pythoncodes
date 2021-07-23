module ESGIZenotech

using GeometryBasics: Point3f0, Vec3f0, GLTriangleFace, Rect, Ngon, Mesh, meta, unit, FRect
using StaticArrays: MVector
using WriteVTK

include("GeometryCollision.jl")

using .GeometryCollision

export load_binary_stl, load_data, triangles_in_bbox, run_demo
export AABB, voxelize, bitstream, bitstream_tri, reconstruct, spy3

# stl = load_data(filename)
# root = AABB(centre=Vec3f0(0,0,0), size=Vec3f0(512,512,512))
# occupied = voxelize(stl, root, 10);
# st = bitstream(occupied)

# Taken from MeshIO.jl with modifications
function load_binary_stl(io::IO, facetype=GLTriangleFace, pointtype=Point3f0,
              normaltype=Vec3f0)
    #Binary STL
    #https://en.wikipedia.org/wiki/STL_%28file_format%29#Binary_STL
    read(io, 80) # throw out header
    triangle_count = read(io, UInt32)

    faces = Array{facetype}(undef, triangle_count)
    vertices = Array{pointtype}(undef, triangle_count * 3)
    normals = Array{normaltype}(undef, triangle_count * 3)

    i = 0
    while !eof(io)
        if length(faces) <= i
            # in case the triangle count was reported incorrectly
            resize!(faces, length(faces)*2)
            resize!(vertices, length(vertices)*2)
            resize!(normals, length(normals)*2)
        end
        faces[i+1] = GLTriangleFace(i * 3 + 1, i * 3 + 2, i * 3 + 3)
        normal = (read(io, Float32), read(io, Float32), read(io, Float32))

        normals[i*3+1] = normaltype(normal...)
        normals[i*3+2] = normals[i*3+1] # hurts, but we need per vertex normals
        normals[i*3+3] = normals[i*3+1]

        vertices[i*3+1] = pointtype(read(io, Float32), read(io, Float32), read(io, Float32))
        vertices[i*3+2] = pointtype(read(io, Float32), read(io, Float32), read(io, Float32))
        vertices[i*3+3] = pointtype(read(io, Float32), read(io, Float32), read(io, Float32))

        skip(io, 2) # throwout 16bit attribute
        i += 1
    end
    resize!(faces, i)
    resize!(vertices, i*3)
    resize!(normals, i*3)

    return Mesh(meta(vertices; normals=normals), faces)
end

function load_data(file)
    load_binary_stl(open(joinpath(@__DIR__, "..", "data", "stl", file)))
end

"""
    AABB

An axis-aligned bounding box
"""
struct AABB
    "centre of the cube"
    centre::Vec3f0
    "size of the cube (half the width/length/depth)"
    size::Vec3f0
    "minimum of bounding box"
    min::Vec3f0
    "maximum of bounding box"
    max::Vec3f0
end

function AABB(; centre=nothing, size=nothing, min=nothing, max=nothing)
    if centre === nothing
        _centre = 0.5.*(min .+ max)
        _size = 0.5.*(max .- min)
        return AABB(_centre, _size, min, max)
    else
        _min = centre - size
        _max = centre + size
        return AABB(centre, size, _min, _max)
    end
end

intersects(c1::AABB, c2::AABB) = all(abs.(c1.centre .- c2.centre) .< (c1.size .+ c2.size))

# Extents [-1500, -1700, 2.125] to [1500, 1600, 292.75]

function bounding_box(tri::Ngon{3})
    tri_min = MVector{3, Float32}(Inf32, Inf32, Inf32)
    tri_max = MVector{3, Float32}(-Inf32, -Inf32, -Inf32)
    for pt in tri
        for idx in 1:3
            tri_min[idx] = tri_min[idx] < pt[idx] ? tri_min[idx] : pt[idx]
            tri_max[idx] = tri_max[idx] > pt[idx] ? tri_max[idx] : pt[idx]
        end
    end
    return AABB(min=Vec3f0(tri_min), max=Vec3f0(tri_max))
end

function voxelize(stl::Mesh, cube::AABB, depth::Integer=8)
    # See http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/tribox_tam.pdf
    n_cubes = 2^depth
    occupied = falses(n_cubes, n_cubes, n_cubes)
    smalldims = 2cube.size ./ n_cubes
    for tri in stl
        bbox = bounding_box(tri)
        if intersects(bbox, cube)
            min = clamp.(floor.(Int32, 1 .+ 0.5n_cubes .* (bbox.min .- cube.min) ./ cube.size), 1, n_cubes)
            max = clamp.(ceil.(Int32, 1 .+ 0.5n_cubes .* (bbox.max .- cube.min) ./ cube.size), 1, n_cubes)
            for xyz in CartesianIndices((min[1]:max[1], min[2]:max[2], min[3]:max[3]))
                origin = cube.min .+ 2 .* (Tuple(xyz) .- 1) ./ n_cubes .* cube.size
                smallcube = FRect(origin, smalldims)
                occupied[xyz] = GeometryCollision.intersects(smallcube, tri)
            end
        end
    end
    return occupied
end

function vtk_hexahedron(cube::Rect{3, T}) where T
    x = Vec{3, T}(cube.widths[1], 0, 0)
    y = Vec{3, T}(0, cube.widths[2], 0)
    z = Vec{3, T}(0, 0, cube.widths[3])

    vertex1 = cube.origin
    vertex2 = vertex1 + x
    vertex3 = vertex2 + y
    vertex4 = vertex1 + y

    vertex5 = vertex1 + z
    vertex6 = vertex5 + x
    vertex7 = vertex6 + y
    vertex8 = vertex5 + y

    return reduce(hcat, (vertex1, vertex2, vertex3, vertex4, vertex5, vertex6, vertex7, vertex8))
end

function write_voxelized(filename, occupied, cube::AABB)
    n_pts = sum(occupied)
    points = zeros(Float32, 3, 8n_pts)
    for xyz in findall(occupied)
        origin = cube.min .+ 2 .* (Tuple(xyz) .- 1) ./ n_cubes .* cube.size
        smallcube = FRect(origin, smalldims)
    end
    cells = [MeshCell(VTKCellTypes.VTK_HEXAHEDRON, 8*(i-1)+1:8i) for i in 1:n_pts]
    vtkfile = vtk_grid(filename, points, cells)
end

function voxelize_old(stl::Mesh, cube::AABB, depth::Integer=8)
    # See http://fileadmin.cs.lth.se/cs/personal/tomas_akenine-moller/code/tribox_tam.pdf
    n_cubes = 2^depth
    occupied = falses(n_cubes, n_cubes, n_cubes)
    smalldims = 2cube.size ./ n_cubes
    for tri in stl
        bbox = bounding_box(tri)
        if intersects(bbox, cube)
            min = clamp.(floor.(Int32, 1 .+ 0.5n_cubes .* (bbox.min .- cube.min) ./ cube.size), 1, n_cubes)
            max = clamp.(ceil.(Int32, 1 .+ 0.5n_cubes .* (bbox.max .- cube.min) ./ cube.size), 1, n_cubes)
            occupied[min[1]:max[1], min[2]:max[2], min[3]:max[3]] .= true
        end
    end
    return occupied
end

function triangles_in_bbox(stl::Mesh, cube::AABB)
    sum(intersects(bounding_box(tri), cube) for tri in stl)
end

function _bitstream!(stream::Vector{UInt8}, levels::Vector{BitArray{N}}, level, idx_prev, terminator) where N
    if level > length(levels)
        terminator(stream, idx_prev)
        return stream
    end
    curr = levels[level]
    idx_curr = (2idx_prev - CartesianIndex(1, 1, 1)):2idx_prev
    value = zero(UInt8)
    for (i, idx) in enumerate(idx_curr)
        if curr[idx]
            value |= one(UInt8) << (i-1)
        end
    end
    push!(stream, value)
    for (i, idx) in enumerate(idx_curr)
        if curr[idx]
            _bitstream!(stream, levels, level+1, idx, terminator)
        end
    end
    return stream
end

function bitstream(occupied::BitArray, terminator=(stream, idx)->nothing)
    @assert (size(occupied, 1) == size(occupied, 2)) && (size(occupied, 2) == size(occupied, 3))
    depth = Int32(log2(size(occupied, 1)))
    levels = [BitArray(undef, 2^level, 2^level, 2^level) for level in 1:depth-1]
    push!(levels, occupied)
    for level in reverse(1:depth-1)
        curr = levels[level]
        prev = levels[level+1]
        for idx_curr in CartesianIndices(size(curr))
            idx_prev = (2idx_curr - CartesianIndex(1, 1, 1)):2idx_curr
            curr[idx_curr] = any(prev[idx_prev])
        end
    end
    return _bitstream!(UInt8[], levels, 1, CartesianIndex(1, 1, 1), terminator)
end

function bitstream(stl::Mesh, cube::AABB; max_length=1, voxelize_depth=10)
    # Calculate number of levels required
    depth = ceil(Int32, log2(2maximum(cube.size)/max_length))
    if depth > voxelize_depth
        rem_depth = depth - voxelize_depth
        n_cubes = 2^rem_depth
        bitstreams = Array{3, UInt8}(undef, n_cubes, n_cubes, n_cubes)
        # TODO: finish this reusing the mechanisms within _bitstream (use terminator)
    else
        return bitstream(voxelize(stl, cube, depth))
    end
end

function _bitstream_tri!(stream::Vector{UInt8},
    levels::Vector{BitArray{N}},
    full_levels::Vector{BitArray{N}},
    level,
    idx_prev,
    terminator
) where N
    curr = levels[level]
    full_curr = full_levels[level]
    idx_curr = (2idx_prev - CartesianIndex(1, 1, 1)):2idx_prev
    value = zero(UInt8)
    for (i, idx) in enumerate(idx_curr)
        if curr[idx]
            value |= one(UInt8) << (i-1)
        end
    end
    push!(stream, value)
    if level < length(levels)
        for (i, idx) in enumerate(idx_curr)
            if full_curr[idx]
                push!(stream, zero(UInt8))
            elseif curr[idx]
                _bitstream_tri!(stream, levels, full_levels, level+1, idx, terminator)
            end
        end
    else
        terminator(stream, idx_prev)
    end
    return stream
end

function bitstream_tri(occupied::BitArray, terminator=(stream, idx)->nothing)
    @assert (size(occupied, 1) == size(occupied, 2)) && (size(occupied, 2) == size(occupied, 3))
    depth = Int32(log2(size(occupied, 1)))
    levels = [BitArray(undef, 2^level, 2^level, 2^level) for level in 1:depth-1]
    push!(levels, occupied)
    full_levels = [BitArray(undef, 2^level, 2^level, 2^level) for level in 1:depth-1]
    push!(full_levels, occupied)
    for level in reverse(1:depth-1)
        curr = levels[level]
        prev = levels[level+1]
        full_curr = full_levels[level]
        full_prev = full_levels[level+1]
        for idx_curr in CartesianIndices(size(curr))
            idx_prev = (2idx_curr - CartesianIndex(1, 1, 1)):2idx_curr
            curr[idx_curr] = any(prev[idx_prev])
            full_curr[idx_curr] = all(full_prev[idx_prev])
        end
    end
    return _bitstream_tri!(sizehint!(UInt8[], 2*1048576), levels, full_levels, 1, CartesianIndex(1, 1, 1), terminator)
end

function _reconstruct!(occupied::BitArray, idx_prev::CartesianIndex, depth, max_depth, stream::Vector{UInt8}, stream_pos)
    # TODO: make this work for tri state encodings
    value = stream[stream_pos]
    stream_pos += 1
    idx_curr = (2idx_prev - CartesianIndex(1, 1, 1)):2idx_prev
    if depth == max_depth
        for (i, idx) in enumerate(idx_curr)
            if (value & (1 << (i-1))) != 0
                occupied[idx] = true
            end
        end
    else
        for (i, idx) in enumerate(idx_curr)
            if (value & (1 << (i-1))) != 0
                stream_pos = _reconstruct!(occupied, idx, depth+1, max_depth, stream, stream_pos)
            end
        end
    end
    return stream_pos
end

function reconstruct(stream::Vector{UInt8}, max_depth::Integer)
    occupied = falses(2^max_depth, 2^max_depth, 2^max_depth)
    _reconstruct!(occupied, CartesianIndex(1, 1, 1), 1, max_depth, stream, 1)
    return occupied
end

# bitstream(stl, AABB(centre=Vec3f0(0, 0, 150), size=Vec3f0(1500, 1700, 150)))

function run_demo(; files=["shear_0p001", "shear_0p010", "shear_0p100"], depths=[7, 8, 9, 10], plot=nothing)
    root = AABB(centre=Vec3f0(0,512,0), size=Vec3f0(512,512,512))
    results = []
    for file in files
        stl = load_data("$file.stl")
        tri_count = triangles_in_bbox(stl, root)
        for depth in depths
            occupied = voxelize(stl, root, depth)
            stream = bitstream_tri(occupied)
            push!(results, (file=file, scale=2root.size[1]/2^depth, depth=depth, tri_count=tri_count, tri_size=tri_count*36, cube_count=sum(occupied), cube_size=length(stream), factor=tri_count*36/length(stream)))
            if plot !== nothing
                plot("$(file)_$depth", root, stl, occupied)
            end
            println("Done $file at depth $depth")
        end
    end
    return results
end

function spy3(occupied)
    idx = findall(occupied)
    x = [i.I[1] for i in idx]
    y = [i.I[2] for i in idx]
    z = [i.I[3] for i in idx]
    return (x, y, z)
end

dir_vec(i::Integer) = i == 1 ? (1, 0, 0) : (i == 2 ? (0, 1, 0) : (0, 0, 1))

function determine_dir(occupied)
    sz = size(occupied)
    sumx = 0
    for k in 1:sz[3], j in 1:sz[2], i in 1:2:sz[1]-1
        sumx += all((occupied[i, j, k], occupied[i+1, j, k]))
    end
    sumy = 0
    for k in 1:sz[3], j in 1:2:sz[2]-1, i in 1:sz[1]
        sumy += all((occupied[i, j, k], occupied[i, j+1, k]))
    end
    sumz = 0
    for k in 1:2:sz[3]-1, j in 1:sz[2], i in 1:sz[1]
        sumz += all((occupied[i, j, k], occupied[i, j, k+1]))
    end
    @show (sumx, sumy, sumz)
    if (sumx == 0) && (sumy == 0) && (sumz == 0)
        return argmax(sz)
    else
        return argmax((sumx, sumy, sumz))
    end
end

condense_dir(occupied, dir::Integer) = condense_dir(occupied, dir_vec(dir))

function condense_dir(occupied, dir::Tuple)
    sz = size(occupied[1])
    new_sz = sz .÷ (1 .+ dir)
    new_occupied = BitArray(undef, new_sz)
    new_occupied_all = BitArray(undef, new_sz)
    for i in CartesianIndices(new_sz)
        ii = CartesianIndex(Tuple(i) .* (1 .+ dir))
        new_occupied[i] = any((occupied[1][ii - CartesianIndex(dir)], occupied[1][ii]))
        new_occupied_all[i] = all((occupied[2][ii - CartesianIndex(dir)], occupied[2][ii]))
    end
    return (new_occupied, new_occupied_all)
end

function condense(occupied)
    sz = size(occupied)
    depth = Int.(log2.(sz))
    level = occupied
    level_all = occupied
    levels = Vector{Tuple{typeof(level), typeof(level)}}()
    dirs = Vector{UInt8}()
    push!(levels, (level, level))
    for i in 2:sum(depth)
        dir = UInt8(determine_dir(level))
        push!(dirs, dir)
        (level, level_all) = condense_dir((level, level_all), dir)
        push!(levels, (level, level_all))
    end
    push!(dirs, argmax(size(level)))
    return (levels, dirs)
end

function tree_size(levels, dirs, level, idx_prev::CartesianIndex)
    (level == 0) && return 0
    dir = dir_vec(dirs[level])
    curr = levels[level]
    counter = 1
    idx_curr = CartesianIndex(Tuple(idx_prev) .* (1 .+ dir)) - CartesianIndex(dir)
    if curr[2][idx_curr]
        counter += 1
    elseif curr[1][idx_curr]
        counter += tree_size(levels, dirs, level-1, idx_curr)
    end
    idx_curr += CartesianIndex(dir)
    if curr[2][idx_curr]
        counter += 1
    elseif curr[1][idx_curr]
        counter += tree_size(levels, dirs, level-1, idx_curr)
    end
    return counter
end

function tree_size(levels, dirs)
    return tree_size(levels, dirs, length(levels), CartesianIndex(1, 1, 1))
end

end # module
