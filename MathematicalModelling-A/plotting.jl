using Makie, GLMakie

GLMakie.set_window_config!(vsync = false)

function plot_data(name, cube, stl, occupied)
    len = size(occupied, 1)
    sc = Scene(resolution = (1920, 1080))

    mesh!(sc, stl, show_axis=false)

    lookat = Vec3f0(389.37955, 63.58121, 215.1531)
    upvector = Vec3f0(-0.21262191, 0.2826073, 0.93537426)
    eyeposition = Vec3f0(804.97406, -488.80692, 476.5174)
    update_cam!(sc, eyeposition, lookat, upvector)
    sc.center = false # prevent to recenter on display (IMPORTANT!!)

    save("C:\\Users\\db9052\\Julia\\ESGIZenotech\\outputs\\$name-base.png", sc)

    xyz = (1:len) .- 0.5len .- 0.5
    cube_size = 2cube.size[1] / len
    pts = [cube.centre + cube_size .* Vec3f0(xyz[pt[1]], xyz[pt[2]], xyz[pt[3]]) for pt in findall(occupied)]
    meshscatter!(sc, pts, marker = Rect3D(Vec3f0(0), Vec3f0(1)), markersize = cube_size, color = :skyblue2, alpha=0.1)

    save("C:\\Users\\db9052\\Julia\\ESGIZenotech\\outputs\\$name-cubes.png", sc)
    return sc
end

function octree_visual()
    sc = mesh(Rect(Vec3f0(0), Vec3f0(1,1,1)), color=(:red,0.4), transparency=true, show_axis=false)
    mesh!(sc, Rect(Vec3f0(0), Vec3f0(0.5,0.5,0.5)), color=(:blue,0.4), transparency=true, show_axis=false)
    mesh!(sc, Rect(Vec3f0(0), Vec3f0(0.25,0.5,0.25)), color=(:green,0.4), transparency=true, show_axis=false)
    mesh!(sc, Rect(Vec3f0(0,0,0.25), Vec3f0(0.25,0.25,0.25)), color=(:green,0.4), transparency=true, show_axis=false)
    return sc
end

"""
    digitsep(value::Integer; separator=",", per_separator=3)
Convert an integer to a string, separating each `per_separator` digits by
`separator`.
    digitsep(12345678)  # "12,345,678"
    digitsep(12345678, seperator= "'")  # "12'345'678"
    digitsep(12345678, seperator= "-", per_separator=4)  # "1234-5678"
"""
function digitsep(value::Integer; seperator=",", per_separator=3)
    isnegative = value < zero(value)
    value = string(abs(value))  # Stringify, no seperators.
    # Figure out last character index of each group of digits.
    group_ends = reverse(collect(length(value):-per_separator:1))
    groups = [value[max(end_index - per_separator + 1, 1):end_index]
              for end_index in group_ends]
    return (isnegative ? "-" : "") * join(groups, seperator)
end

function print_tables()
    # bbox-based voxelisation
    # results = Any[(file = "shear_0p001", scale = 8.0f0, depth = 7, tri_count = 1392302, tri_size = 50122872, cube_count = 90133, cube_size = 19953, factor = 2512.0469102390616), (file = "shear_0p001", scale = 4.0f0, depth = 8, tri_count = 1392302, tri_size = 50122872, cube_count = 399968, cube_size = 83038, factor = 603.6136708494906), (file = "shear_0p001", scale = 2.0f0, depth = 9, tri_count = 1392302, tri_size = 50122872, cube_count = 1870981, cube_size = 369363, factor = 135.7008471341201), (file = "shear_0p001", scale = 1.0f0, depth = 10, tri_count = 1392302, tri_size = 50122872, cube_count = 9229389, cube_size = 1653652, factor = 30.310411138498306), (file = "shear_0p010", scale = 8.0f0, depth = 7, tri_count = 371475, tri_size = 13373100, cube_count = 31604, cube_size = 9175, factor = 1457.558583106267), (file = "shear_0p010", scale = 4.0f0, depth = 8, tri_count = 371475, tri_size = 13373100, cube_count = 103220, cube_size = 29002, factor = 461.1095786497483), (file = "shear_0p010", scale = 2.0f0, depth = 9, tri_count = 371475, tri_size = 13373100, cube_count = 369749, cube_size = 94642, factor = 141.30195896113776), (file = "shear_0p010", scale = 1.0f0, depth = 10, tri_count = 371475, tri_size = 13373100, cube_count = 1512647, cube_size = 339730, factor = 39.36390663173697), (file = "shear_0p100", scale = 8.0f0, depth = 7, tri_count = 13197, tri_size = 475092, cube_count = 3275, cube_size = 1525, factor = 311.535737704918), (file = "shear_0p100", scale = 4.0f0, depth = 8, tri_count = 13197, tri_size = 475092, cube_count = 6689, cube_size = 3140, factor = 151.3031847133758), (file = "shear_0p100", scale = 2.0f0, depth = 9, tri_count = 13197, tri_size = 475092, cube_count = 15883, cube_size = 6384, factor = 74.41917293233082), (file = "shear_0p100", scale = 1.0f0, depth = 10, tri_count = 13197, tri_size = 475092, cube_count = 47884, cube_size = 14960, factor = 31.757486631016043)]
    # proper voxelisation
    results = Any[(file = "shear_0p001", scale = 8.0f0, depth = 7, tri_count = 1392302, tri_size = 50122872, cube_count = 20755, cube_size = 11441, factor = 4380.986976662879), (file = "shear_0p001", scale = 4.0f0, depth = 8, tri_count = 1392302, tri_size = 50122872, cube_count = 95417, cube_size = 50560, factor = 991.3542721518987), (file = "shear_0p001", scale = 2.0f0, depth = 9, tri_count = 1392302, tri_size = 50122872, cube_count = 432547, cube_size = 222793, factor = 224.97507551853064), (file = "shear_0p001", scale = 1.0f0, depth = 10, tri_count = 1392302, tri_size = 50122872, cube_count = 1772687, cube_size = 902715, factor = 55.52458084777587), (file = "shear_0p010", scale = 8.0f0, depth = 7, tri_count = 371475, tri_size = 13373100, cube_count = 6841, cube_size = 4828, factor = 2769.904722452361), (file = "shear_0p010", scale = 4.0f0, depth = 8, tri_count = 371475, tri_size = 13373100, cube_count = 23459, cube_size = 15456, factor = 865.236801242236), (file = "shear_0p010", scale = 2.0f0, depth = 9, tri_count = 371475, tri_size = 13373100, cube_count = 86301, cube_size = 52496, factor = 254.74512343797622), (file = "shear_0p010", scale = 1.0f0, depth = 10, tri_count = 371475, tri_size = 13373100, cube_count = 330197, cube_size = 187236, factor = 71.42376466064218), (file = "shear_0p100", scale = 8.0f0, depth = 7, tri_count = 13197, tri_size = 475092, cube_count = 590, cube_size = 758, factor = 626.7704485488126), (file = "shear_0p100", scale = 4.0f0, depth = 8, tri_count = 13197, tri_size = 475092, cube_count = 1291, cube_size = 1499, factor = 316.9392928619079), (file = "shear_0p100", scale = 2.0f0, depth = 9, tri_count = 13197, tri_size = 475092, cube_count = 3455, cube_size = 3124, factor = 152.07810499359795), (file = "shear_0p100", scale = 1.0f0, depth = 10, tri_count = 13197, tri_size = 475092, cube_count = 10726, cube_size = 7776, factor = 61.09722222222222)]
    res = Dict{Any, Any}()
    for result in results
        res[(result.file, Int32(result.scale))] = result
    end

    files = ["shear_0p001", "shear_0p010", "shear_0p100"]
    scales = Int32[1, 2, 4, 8]

    print("| Counts |")
    for file in files
        print(" $file |")
    end
    println()
    print("| - |")
    for file in files
        print(" - |")
    end
    println()
    print("| Triangles |")
    for file in files
        print(" $(digitsep(res[(file, scales[1])].tri_count)) |")
    end
    println()
    for scale in scales
        print("| $scale metre cubes |")
        for file in files
            print(" $(digitsep(res[(file, scale)].cube_count)) |")
        end
        println()
    end
    println()

    print("| Storage size (bytes) |")
    for file in files
        print(" $file |")
    end
    println()
    print("| - |")
    for file in files
        print(" - |")
    end
    println()
    print("| Triangles |")
    for file in files
        print(" $(digitsep(res[(file, scales[1])].tri_size)) |")
    end
    println()
    for scale in scales
        print("| $scale metre cubes |")
        for file in files
            print(" $(digitsep(res[(file, scale)].cube_size)) |")
        end
        println()
    end
    println()

    print("| Compression factor |")
    for file in files
        print(" $file |")
    end
    println()
    print("| - |")
    for file in files
        print(" - |")
    end
    println()
    for scale in scales
        print("| $scale metre cubes |")
        for file in files
            print(" $(digitsep(round(Int, res[(file, scale)].factor))) |")
        end
        println()
    end
    println()
end
