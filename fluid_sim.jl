
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


function simulate_fluids()
    # initialize position to a random location within my bounding_box
    position = [@SVector rand(3) for _ in 1:num_particles]
    position = (.*).(position, Scalar((bounding_box.max .- bounding_box.min)))
    position .+= Scalar(bounding_box.min)

    position_prev = copy(position)
    velocity = [@SVector zeros(3) for _ in 1:num_particles] # fluid begins at rest
    mass = fill(100.0, num_particles) # fluid has unit volume

    Hash.init_file(global_filepath) # create or clear the file located at filepath

    # warm up. Load all files using Julia's just in time Compiler
    print("Warm up.................... ")
    @time Collision.handle_collisions(position, position_prev, objects)

    for i in 1:simulation_steps
        println("======")
        println("Step $i")

        # Store particles into cells using an indexing technique
        print("Hash particles.............. ")
        @time hashed_particles = Hash.generate_hash(position)
        
        # Calculate density and pressure
        print("Calc density & pressure..... ")
        @time density, pressure = Math.calc_density_and_pressure(position, mass, hashed_particles)
        
        # Calculate acceleration using the Navier Stokes equation for incompressible fluids
        print("Apply Navier-Stokes......... ")
        @time acceleration = Math.apply_forces(position, velocity, density, pressure, mass, hashed_particles)
        
        # Update position and velocity of particles in the system. 
        print("Update position............. ")
        @time position, position_prev, velocity = Math.update_position_and_velocity(position, position_prev, velocity, acceleration, delta_time)

        # A ray is traced between position and position_prev. Any ray which collides with
        # an object in the scene is reflected according to the surface normal. 
        print("Detect collisions........... ")
        @time position = Collision.handle_collisions(position, position_prev, objects)
        
        # display particles-- wrap function with user input to pause/play, & step forward by 1 timestep
        # print("Display particles........... ")
        # @time Hash.convert_to_camera_coords(position)

        # Save position to an output file
        print("Saving particles............ ")
        @time Hash.save_to_file(position, global_filepath)
        
    end
end
@time simulate_fluids()

end # module