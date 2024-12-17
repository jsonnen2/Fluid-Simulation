
module fluid_sim

push!(LOAD_PATH, pwd())

include("GfxBase.jl")
include("Mesh.jl")
include("Hash.jl")
include("Collision.jl")
include("Math.jl")

using .GfxBase
using .Math
using .Hash
using StaticArrays


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
    storage = Dict(i => Int[] for i in 1:max_cell)

    for (idx, pos) in enumerate(position)
        cell = Hash.coordinate_to_cell_int(pos)
        push!(storage[cell], idx)
    end
    return storage
end

function render_pipeline()

    # initialize position to a random location within my bounding_box
    position = [@SVector rand(3) for _ in 1:num_particles]
    position = (.*).(position, Scalar((bounding_box.max .- bounding_box.min)))
    position .+= Scalar(bounding_box.min)
    
    position_prev = copy(position)
    velocity = [@SVector zeros(3) for _ in 1:num_particles] # fluid begins at rest
    mass = fill(1000.0, num_particles) # fluid has unit volume

    # Collision.handle_collisions(position, position_prev, objects) # TODO: soln for comiling time
    for i in 1:simulation_steps
        println(i)

        # Store particles into cells using an indexing technique
        @time hashed_particles::Dict{Int, Vector{Int}} = generate_hash(position)
        
        # Calculate density and pressure
        @time density, pressure = Math.calc_density_and_pressure(position, mass, hashed_particles)
        
        # Calculate acceleration using the Navier Stokes equation for incompressible fluids
        @time acceleration = Math.apply_forces(position, velocity, density, pressure, mass, hashed_particles)
        
        # Update position and velocity of particles in the system. Also handles collisions with objects. 
        @time position, position_prev, velocity = Math.update_position_and_velocity(position, position_prev, velocity, acceleration, delta_time)
        
        @time position = Collision.handle_collisions(position, position_prev, objects)
        @time position = Collision.handle_collisions(position, position_prev, objects)
        # display particles-- wrap function with user input to pause/play, & step forward by 1 timestep

    end
end
@time render_pipeline()
end # module