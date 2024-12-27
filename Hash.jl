
# O(n) hashing operation which places particles into a cell.
# A particle only needs to lookup particles in adjacent cells

module Hash

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
using .GfxBase
# using ..GfxBase
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

    # PROBLEMS
    # wrap-around in Hash.adjacent_cells_int. Basically, I don't check if the cell exists on the 
    #   other side of the bounding box. 
    # overkill in the size of my dictionary. A particle could exist on the boundary line and it would
    #   be placed into its own cell that only boundary line paritcles could live on. 

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

    [push!(storage[c], i) for (i, c) in enumerate(cells)]

    return storage
end

function find_neighbors(cell::Int, hash::Dict{Int, Vector{Int}}) :: Vector{Int}
    # Given cell::Int and hash::Dict{Int, Vector{Int}}
    # Return Vector{Int} which are all indices of all neighbors AND itself

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
end

function save_to_file(position::Vector{Vec3}, filepath::String)
    # Use standard out to convert position::Vector{Vec3} to a string 
    open(filepath, "a") do io
        println(io, position)
    end
end

end #module