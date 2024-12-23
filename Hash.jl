
# O(n) hashing operation which places particles into a cell.
# A particle only needs to lookup particles in adjacent cells

module Hash

push!(LOAD_PATH, pwd())
include("GfxBase.jl")
using .GfxBase
# using ..GfxBase
using StaticArrays

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

end