
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

    print("Initialize Containers....... ")
@time begin
    position = Hash.init_position(spawn_cube, spawn_type)
    position_prev = copy(position)
    velocity = [@SVector zeros(3) for _ in 1:num_particles] # fluid begins at rest
    mass = fill(10.0, num_particles) # fluid has unit volume

    Hash.init_file(global_filepath) # create or clear the file located at global_filepath
        if save_obj_files
            Hash.save_scene_objects(objects, obj_save_folder)
        end
end
    # warm up. Load all files using Julia's Just In Time compiler
    print("Warm up..................... ")
    @time for _ in 1:1 # group together for timing
        hashed_particles = Hash.generate_hash(position)
        acceleration = Math.apply_forces(position, velocity, mass, hashed_particles)
        Math.update_position_and_velocity(position, position_prev, velocity, acceleration, delta_time)
        Collision.handle_collisions(position, position_prev, objects)
    end

    for i in 1:simulation_steps
        if !console_log # Turn off printing to console.
            if i % 100 == 0
                println(i)
            end

            hashed_particles = Hash.generate_hash(position)
            acceleration = Math.apply_forces(position, velocity, mass, hashed_particles)
            position, position_prev, velocity = Math.update_position_and_velocity(position, position_prev, velocity, acceleration, delta_time)
            position = Collision.handle_collisions(position, position_prev, objects)
            Hash.save_to_file(position, velocity, global_filepath)
        else 

        println("=======")
        println("Step $i")

        # TODO: graph distribution for each time trial on 10,000 particles
        # Store particles into cells using an indexing technique
        print("Hash particles.............. ")
        @time hashed_particles = Hash.generate_hash(position)

        # Calculate density and pressure
        # print("Calc density & pressure..... ")
        # @time density, pressure = Math.calc_density_and_pressure(position, mass, hashed_particles)
        
        # Calculate acceleration using the Navier Stokes equation for incompressible fluids
        print("Apply Navier-Stokes......... ")
        @time acceleration = Math.apply_forces(position, velocity, mass, hashed_particles)
        save_pos = copy(position)
        save_pos_prev = copy(position_prev)

        # Update position and velocity of particles in the system. 
        print("Update position............. ")
        @time position, position_prev, velocity = Math.update_position_and_velocity(position, position_prev, velocity, acceleration, delta_time)
        save_pos_future = copy(position)

        # A ray is traced between position and position_prev. Any ray which collides with
        # an object in the scene is reflected according to the surface normal. 
        print("Detect collisions........... ")
        @time position = Collision.handle_collisions(position, position_prev, objects)
        
        where_negative = findall(x -> any(y -> y < 0, x), position)
        if length(where_negative) > 0
            println(where_negative)
            println(position[where_negative])
            println()
            println(save_pos_future[where_negative])
            println()
            println(save_pos[where_negative])
            println(length(where_negative))
        end
# ???????? seconds
        # Create a datatype that is useful for rendering particles on the scene.
        # print("Camera hashing........... ")
        # @time Hash.convert_to_camera_coords(position)

# 0.042220 seconds
        # Save position to an output file
        print("Saving particles............ ")
        @time Hash.save_to_file(position, velocity, global_filepath)
        end
    end
end
@time simulate_fluids()

end # module