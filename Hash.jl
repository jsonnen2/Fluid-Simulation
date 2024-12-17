
# O(n) hashing operation which places particles into a cell.
# A particle only needs to lookup particles in adjacent cells

module Hash

using ..GfxBase

export adjacent_cells_int, adjacent_cells_tuple
export coordinate_to_cell_int, coordinate_to_cell_tuple
export cell_int_to_int, cell_int_to_tuple
export find_neighbors

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
    num_cells = Int.(ceil.((bounding_box.max .- bounding_box.min) ./ smoothing_radius))
    # Filter for cells outside [min, max]
    filtered_adj_cells = []
    for box in adjacent_cells
        if all(collect(box) .>= [0,0,0] .&& collect(box) .< num_cells) && box !== (0,0,0) #TODO does projection work?
            push!(filtered_adj_cells, box)
        end
    end
    return filtered_adj_cells
end

function adjacent_cells_int(i::Int) :: Vector{Int}
    # TODO: there is a boundary problem. Basically, I don't remove cells which "wrap-around" the box.
    cell = cell_int_to_tuple(i)
    adjacent_cells = adjacent_cells_tuple(cell)
    adjacent_ints = [cell_tuple_to_int(x) for x in adjacent_cells]
    return adjacent_ints
end

function coordinate_to_cell_tuple(coord::Vec3)::Tuple{Int, Int, Int}
    # Ensure boundary conditions are met.
    for (i, c) in enumerate(coord)
        if any(c .< bounding_box.min)
            println(coord[i])
        elseif any(c .> bounding_box.max)
            println(coord[i])
        elseif any(isnan, c)
            println(coord[i])
        end
    end

    @assert all(coord .>= bounding_box.min) "Coordinate is out of bounds (less than min)"
    @assert all(coord .<= bounding_box.max) "Coordinate is out of bounds (greater than max)"
    
    cell = floor.((coord .- bounding_box.min) ./ smoothing_radius)
    return Tuple(Int.(cell))
end

function coordinate_to_cell_int(coord::Vec3)::Int
    cell = coordinate_to_cell_tuple(coord)
    return cell_tuple_to_int(cell) + 1
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


function find_neighbors(cell::Int, hash::Dict{Int, Vector{Int}}) :: Vector{Int}
    # Given cell::Int and hash::Dict{Int, Vector{Int}}
    # Return Vector{Int} which are all indices of all neighbors AND itself

    nearby_indices = adjacent_cells_int(cell)
    storage = []
    for idx in nearby_indices
        append!(storage, hash[idx])
    end
    return storage
end

end