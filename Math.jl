
# TODO:
# Better formula for velocity. Current equation has ϵ^1 (first degree error)
# Pressure equations might be sus.
# better solution for empty cells (density = 0 which causes divide by 0 error)

module Math

using ..GfxBase
using ..Hash
using ..Collision
using StaticArrays
using LinearAlgebra


function apply_forces(position::Vector{Vec3}, velocity::Vector{Vec3}, 
                        mass::Vector{Float64}, hash::Dict{Int, Vector{Int}})

    hash_max = Hash.coordinate_to_cell_int(bounding_box.max)
    acceleration = [@SVector zeros(3) for _ in 1:num_particles] 
    
    # iterate over cells
    for cell in 0:hash_max
        
        nearby::Vector{Int} = Hash.find_neighbors(cell, hash) # TODO: slow. ~0.1 seconds for 10k
        cell_indices::Vector{Int} = hash[cell]

        # for i in cell_indices 
        Threads.@threads for i in cell_indices

            gravitational_force = @SVector [0, -9.81 * mass[i], 0] # gravitational force on Earth

            kernel = kernel_gradient(position[nearby], position[i])
            pressure_force = mass[i] .* -kernel
            
            acceleration[i] = gravitational_force - sum(pressure_force)
        end
    end
    contains_nan = any(x -> any(isnan, x), acceleration)
    @assert !contains_nan "The vector contains NaN values!"

    # Reduce accelerations whose magnitudes greater than the cap
    mask = (cap_acceleration .< norm.(acceleration))
    acceleration[mask] .= cap_acceleration .* normalize.(acceleration[mask])
    return acceleration
end

function update_position_and_velocity(position::Vector{Vec3}, position_prev::Vector{Vec3}, velocity::Vector{Vec3}, 
                                        acceleration::Vector{Vec3}, delta_time::Float64)
    # 4th degree error term
    # 2*x(t) - x(t-1) + a(t)*dt^2 + ϵ^4
    acc_term = acceleration .* Scalar(delta_time^2)
    new_position = position .+ position .- position_prev .+ acc_term
    new_velocity = (new_position .- position) ./ Scalar(delta_time)
    
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
    # distance .= ifelse.(distance .== 0.0, smoothing_radius, distance) #TODO: is this O(n)?
    
    coeff = (15/pi) / smoothing_radius^4
    return coeff .* (smoothing_radius .- distance).^3
end

function kernel_gradient(nearby::Vector{Vec3}, center::Vec3)
    direction::Vector{Vec3} = nearby .- Scalar(center)
    # direction ./= norm.(direction) # normalize direction vector to have a length of 1
    return kernel_funct(nearby, center) .* direction
end


end