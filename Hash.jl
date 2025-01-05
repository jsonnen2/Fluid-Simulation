
# O(n) hashing operation which places particles into a cell.
# A particle only needs to lookup particles in adjacent cells

module Hash

using ..GfxBase
using ..Mesh
using StaticArrays
using LinearAlgebra

export adjacent_cells_int, adjacent_cells_tuple
export coordinate_to_cell_int, coordinate_to_cell_tuple
export cell_tuple_to_int, cell_int_to_tuple
export find_neighbors


function generate_hash(position::Vector{Vec3})

    # Generate a dictionary which stores each particle that is within a cell. Cells are unit cubes in R^3.
    # key: cell number => value: list of particle indices
    # Runtime = O(num_particles)

    
    # Ensure boundary conditions are met.
    @assert all(x -> all(x .>= bounding_box.min), position) "Coordinate is out of bounds (less than min)"
    @assert all(x -> all(x .<= bounding_box.max), position) "Coordinate is out of bounds (greater than max)"

    max_cell = Hash.coordinate_to_cell_int(bounding_box.max)
    storage = Dict(i => Int[] for i in 0:max_cell)

    # TODO: broadcast to find cell from position
    num_cells = Int.(ceil.((bounding_box.max .- bounding_box.min) ./ smoothing_radius))
    range = bounding_box.max .- bounding_box.min
    discretized_position = (position .- Scalar(bounding_box.min)) ./ Scalar(range) .* Scalar(num_cells)
    discretized_position = [Int.(floor.(pos)) for pos in discretized_position]
    # TODO: convert position to scale uniformly from 1 to num_cells
    cells = getindex.(discretized_position, 1) .* Scalar(num_cells[2] * num_cells[3]) .+
            getindex.(discretized_position, 2) .* Scalar(num_cells[3]) .+
            getindex.(discretized_position, 3)

    for (index, cell_id) in enumerate(cells)
        if cell_id == 2197 ## DEBUG
            println(position[index])
        elseif cell_id < 0
            println(position[index])
        end
        push!(storage[cell_id], index)
    end

    return storage
end

function find_neighbors(cell::Int, hash::Dict{Int, Vector{Int}}) :: Vector{Int}
    # Given cell::Int and hash::Dict{Int, Vector{Int}}
    # Return Vector{Int} which are all indices of all neighbors AND itself
    # TODO: slow. ~0.1 seconds for 10k

    nearby_indices = adjacent_cells_int(cell)
    storage = []
    for idx in nearby_indices
        # if idx == 2198 1099 157 7 2
        #     println(cell)
        #     println(nearby_indices)
        # end
        append!(storage, hash[idx])
    end
    return storage
end

function adjacent_cells_tuple(cell::Tuple{Int, Int, Int}) :: Vector{Tuple{Int, Int, Int}}
    x, y, z = cell
    adjacent_cells = [

        (x-1, y-1, z-1), (x-1, y-1, z), (x-1, y-1, z+1),
        (x-1, y, z-1), (x-1, y, z), (x-1, y, z+1),
        (x-1, y+1, z-1), (x-1, y+1, z), (x-1, y+1, z+1),
        
        (x, y-1, z-1), (x, y-1, z), (x, y-1, z+1),
        (x, y, z-1), (x, y, z), (x, y, z+1),
        (x, y+1, z-1), (x, y+1, z), (x, y+1, z+1),
        
        (x+1, y-1, z-1), (x+1, y-1, z), (x+1, y-1, z+1),
        (x+1, y, z-1), (x+1, y, z), (x+1, y, z+1),
        (x+1, y+1, z-1), (x+1, y+1, z), (x+1, y+1, z+1),
    ]
    min_boundary = coordinate_to_cell_tuple(bounding_box.min)
    max_boundary = coordinate_to_cell_tuple(bounding_box.max)
    mask = findall(x -> all(min_boundary .<= x .<= max_boundary), adjacent_cells)
    return adjacent_cells[mask]
end

function adjacent_cells_int(i::Int) :: Vector{Int}
    # TODO: there is a boundary problem. Basically, I don't remove cells which "wrap-around" the box.
    cell = cell_int_to_tuple(i)
    adjacent_cells = adjacent_cells_tuple(cell)
    adjacent_ints = [cell_tuple_to_int(x) for x in adjacent_cells]
    return adjacent_ints
end

####################
# HELPER FUNCTIONS #
####################

function coordinate_to_cell_tuple(coord::Vec3)::Tuple{Int, Int, Int}
    # 0 indexed.
    # Returns tuples on the range: bounding_box.min to num_cells-1
    # Returns a tuple, which allows the bounding box to be rectangular

    # Ensure boundary conditions are met.
    @assert all(coord .>= bounding_box.min) "Coordinate is out of bounds (less than min)"
    @assert all(coord .<= bounding_box.max) "Coordinate is out of bounds (greater than max)"
    
    cell = floor.((coord .- bounding_box.min) ./ smoothing_radius)
    return Tuple(Int.(cell))
end

function coordinate_to_cell_int(coord::Vec3)::Int
    cell = coordinate_to_cell_tuple(coord)
    return cell_tuple_to_int(cell)
end

function cell_int_to_tuple(i::Int)::Tuple{Int, Int, Int}
    num_cells = Int.(ceil.((bounding_box.max .- bounding_box.min) ./ smoothing_radius))
    x = div(i, num_cells[2] * num_cells[3])
    y = div(i % (num_cells[2] * num_cells[3]), num_cells[3])
    z = i % num_cells[3]
    return (x, y, z)
end

function cell_tuple_to_int(cell::Tuple{Int, Int, Int})::Int
    num_cells = Int.(ceil.((bounding_box.max .- bounding_box.min) ./ smoothing_radius))
    x, y, z = cell
    i = x * (num_cells[2] * num_cells[3]) + y * num_cells[3] + z
    return i
end

function init_position(volume::Box, type::String)
    if type == "uniform"
        n = round(Int, num_particles^(1/3))
        if n^3 < num_particles
            n += 1
        end
        spacing::Vec3 = (volume.max .- volume.min) ./ n
        x = collect(0:n-1) .* spacing[1] .+ spacing[1] / 2
        x = repeat(x, n^2)[1:num_particles]
        y = collect(0:n-1) .* spacing[2] .+ spacing[2] / 2
        y = repeat(repeat(y, inner=n), n)[1:num_particles]
        z = collect(0:n-1) .* spacing[3] .+ spacing[3] / 2
        z = repeat(z, inner=n^2)[1:num_particles]

        position = [@SVector [a,b,c] for (a,b,c) in zip(x,y,z)]
    elseif type == "random"
        # initialize position to a random location within my volume
        position = [@SVector rand(3) for _ in 1:num_particles]
        position = (.*).(position, Scalar((volume.max .- volume.min)))
        position .+= Scalar(volume.min)
    end
    return position
end

function convert_to_camera_coords(position::Vector{Vec3})

    # Convert position from world coordinates (x,y,z) into camera coords (i,j,k)
    j_direction = normalize(camera.up)
    k_direction = normalize(camera.view)
    i_direction = cross(j_direction, k_direction)
    M = SMatrix{3,3,Float64}([i_direction..., j_direction..., k_direction...])
    camera_coords = Scalar(M) .* (position .- Scalar(camera.origin))

    # camera.origin in camera coordinates
    origin = M * camera.origin
    # draw a ray from camera origin to position
    ray = camera_coords .- Scalar(origin)

    # find intersection with the screen- 
    # a plane defined by point camera.origin + camera.view*distance and normal -camera.view
    # t = dot(n, a - p) / dot(n, d)
    # TODO: can optimize for FLOPS
    point = origin + camera.view * camera.focal # TODO Some function of focal
    normal = -camera.view
    t = Scalar(dot(normal, point - origin)) ./ dot.(Scalar(normal), ray)

    # All points with t > 1 & t < 0 don't get hashed
    ray = ray[0 .< t .< 1]
    t = t[0 .< t .< 1]
    
    # position on the screen for each particle
    # Domain: (0,0) -> (1,1)
    screen_coords = Scalar(origin) .+ t .* ray

    # Store positions of particles at each pixel.
    # TODO: this init is VERY slow-- container of sorted vectors
    # 

    # hash = Dict((i+width*j) => [] for i in 0:width for j in 0:height)
    # z_buffer = Dict((i+width*j) => [] for i in 0:width for j in 0:height)

    # for (i, pos) in enumerate(screen_coords)
    #     # if pixel is on the screen
    #     if (0 < pos[1] < 1) && (0 < pos[2] < 1)
    #         pixel_coord = (Int(round(pos[1]*width)), Int(round(pos[2]*height)))
    #         push!(hash[pixel_coord], position[i])
    #         push!(z_buffer[pixel_coord], ray[i][3])
    #     end
    # end

    # TODO: use a sorting alg with benefits most from nearly sorted data.
    # # Sort hash according to z_buffer
    # for key in keys(hash)
    #     sorted_indices = sortperm(z_buffer[key])
    #     sorted_data = hash[key][sorted_indices]
    #     hash[key] = sorted_data
    # end

    return ray
end

function init_file(filepath::String)
    # create or clear a file located at filepath.
    if isfile(filepath)
        rm(filepath)
    end
    touch(filepath)

    folderpath = first(split(filepath, ".txt"))
    if ispath(folderpath)
        rm(folderpath, recursive=true)
    end
    mkdir(folderpath)
end

function save_to_file(position::Vector{Vec3}, velocity::Vector{Vec3}, filepath::String)
    # Use standard out to convert position::Vector{Vec3} to a string 

    mag = sqrt.(getindex.(velocity, 1).^2 .+ 
                getindex.(velocity, 2).^2 .+
                getindex.(velocity, 3).^2)
    color_scalar = mag ./ Scalar(maximum(mag))
    four_vector= [@SVector [pos[1], pos[2], pos[3], color] for (pos, color) in zip(position, color_scalar)]
    open(filepath, "a") do io
        println(io, four_vector)
    end
end

function save_scene_objects(objects::Vector{<:Shape}, local_folder::String)
    # create local_folder
    if !isdir(local_folder)
        mkdir(local_folder)
    end

    # Iterate through objects in my vector and save them
    for object in objects
        if object isa OBJ
            if isfile(object.filename)
                # load .obj file
                mesh::OBJMesh = Mesh.read_obj(object.filename)
                
                if object.rotate !== nothing
                    Mesh.rotate!(mesh, object.rotate)
                end
                if object.scale !== nothing
                    Mesh.scale!(mesh, object.scale)
                end
                if object.translate !== nothing
                    Mesh.translate!(mesh, object.translate)
                end
                object_name = last(split(object.filename, '/'))
                object_name = first(split(object_name, ".obj"))
                filepath = "$local_folder/$object_name-$(object.translate)-$(object.scale)-$(object.rotate).obj"
                object.filename = filepath
                Mesh.write_obj(filepath, mesh)
            else
                error("Filepath $(object.filename) does not exist.")
            end
        elseif object isa Box
            scale_factor = (object.max .- object.min) ./ Scalar(2)
            translate_factor = object.min
            filepath = "$local_folder/cube-$scale_factor-$translate_factor.obj"
            object.filename = filepath

            cube::OBJMesh = Mesh.cube_mesh()
            Mesh.scale!(cube, scale_factor)
            Mesh.translate!(cube, translate_factor)
            Mesh.write_obj(filepath, cube)

        elseif object isa Sphere
            # TODO: save some data into a file that could be read by Unity for their sphere GameObject
        end
    end
end

end #module