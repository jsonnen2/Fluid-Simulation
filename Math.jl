
# TODO:
# Better formula for velocity. Current equation has ϵ^1 (first degree error)
# Pressure equations might be sus.
# better solution for empty cells (density = 0 which causes divide by 0 error)

module Math

using ..GfxBase
using ..Hash
using ..Collision
using StaticArrays


function calc_density_and_pressure(position::Vector{Vec3}, mass::Vector{Float64}, hash::Dict{Int, Vector{Int}})
    # calculates density and pressure. These values are stored in the Particle object
    # O(n * E(# particles per cell)) time complexity

    sound = 343 # speed of sound in air. Units m/s^2
    rest_density = 1000 # density of water with units: gram per cubic centimeter 
    adiabatic_index = 7.15 # water is nearly incompressible
    atmospheric_pressure = 101_325 # Pascals ≡ Newton per square meter

    #TODO: there has got to be a better way of dealing with 0 density in a cell
    # I initialize to epsilon bc in later calculations I divide by density (avoiding divide by 0 errors)
    density = fill(1e-10, num_particles) # grams per cubic centimeter
    pressure = zeros(num_particles) # Pascal ≡ Newton per square meter

    hash_max = Hash.coordinate_to_cell_int(bounding_box.max)
    # iterate over cells
    for cell in 1:hash_max

        nearby_indices::Vector{Int} = Hash.find_neighbors(cell, hash)
        cell_indices::Vector{Int} = hash[cell]

        for i in cell_indices
            density[i] += sum(mass[nearby_indices] .* kernel_funct(position[nearby_indices], position[i]))
            # Cole equation 
            coeff = sound^2 / adiabatic_index
            pressure[i] = coeff * ((density[i] / rest_density)^adiabatic_index - 1) + atmospheric_pressure
        end
    end
    return density, pressure
end

function apply_forces(position::Vector{Vec3}, velocity::Vector{Vec3}, density::Vector{Float64}, 
                        pressure::Vector{Float64}, mass::Vector{Float64}, hash::Dict{Int, Vector{Int}})

    hash_max = Hash.coordinate_to_cell_int(bounding_box.max)
    acceleration = [@SVector zeros(3) for _ in 1:num_particles] #TODO: could this be slow?

    # iterate over cells
    for cell in 1:hash_max
        
        nearby::Vector{Int} = Hash.find_neighbors(cell, hash)
        cell_indices::Vector{Int} = hash[cell]

        for i in cell_indices # TODO: parallelize (thread) this loop

            # Use Navier Stokes equation for incompressible fluids.
            # Calculate the force applied to a particle from: gravity, pressure, & viscosity.
            visc_coeff = Scalar(1) # 1e-3 = viscosity coefficient of water 
            gravitational_force = @SVector [0, -9.81 * mass[i], 0] # gravitational force on Earth

            # (pressure[i] / density[i]^2 + pressure / density^2) * mass^2 * kernel_gradient() * density[i]
            summation = (pressure[nearby] ./ density[nearby].^2) .+ Scalar(pressure[i] / density[i]^2)
            kernel = kernel_gradient(position[nearby], position[i])
            pressure_force = summation .* mass[nearby] .* kernel .* density[i]

            # visc_coeff * mass^2 * (velocity - velocity[i]) / density * kernel()
            diff = velocity[nearby] .- Scalar(velocity[i])
            kernel = kernel_funct(position[nearby], position[i])
            viscosity_force = Scalar(visc_coeff) .* mass[nearby].^2 .* diff .* kernel ./ density[nearby]
            
            force = gravitational_force .- sum(pressure_force) .+ sum(viscosity_force)
            acceleration[i] = force ./ density[i] # TODO: should I use density or mass?
        end
    end
    return acceleration
end

function update_position_and_velocity(position::Vector{Vec3}, position_prev::Vector{Vec3}, velocity::Vector{Vec3}, 
                                        acceleration::Vector{Vec3}, delta_time::Float64)
    # 4th degree error term
    # 2*x(t) - x(t-1) + a(t)*dt^2 + ϵ^4
    acc_term = acceleration .* Scalar(delta_time^2)
    new_position = position .+ position .- position_prev .+ acc_term
    new_velocity = (new_position .- position) ./ Scalar(delta_time)
    
    # if new_position lies outside of my bounding box, perform reflections to move it back inside. 
    # @time new_position = Collision.handle_collisions(new_position, position, objects)
    # @time new_position = Collision.handle_collisions(new_position, position, objects)
    
    return new_position, position, new_velocity
end

function euc_distance(nearby::Vector{Vec3}, center::Vec3)
    diff = (nearby .- Scalar(center))
    dist = (.^).(diff, 2)
    add = (+).(getindex.(dist, 1), getindex.(dist, 2), getindex.(dist, 3))
    return sqrt.(add)
end

# TODO: kernel problem = When area under the kernel equals 1, then each point sampled on 
#       the kernel will be quite small. Scalar multiplier? Is there theory on how much to scale? 
function kernel_funct(nearby::Vector{Vec3}, center::Vec3)
    # Kernel Function from YouTuber Sebestian Lague
    # "Coding Adventure: Simulating Fluids". Timestamp= 20:21
    distance::Vector{Float64} = euc_distance(nearby, center)
    clamp!(distance, 0.0, smoothing_radius)
    distance .= ifelse.(distance .== 0.0, smoothing_radius, distance) #TODO: is this O(n)?
    
    coeff = (15/pi) / smoothing_radius^4
    return coeff .* (smoothing_radius .- distance).^3
end

function kernel_gradient(nearby::Vector{Vec3}, center::Vec3)
    direction::Vector{Vec3} = nearby .- Scalar(center)
    return kernel_funct(nearby, center) .* direction
end


end