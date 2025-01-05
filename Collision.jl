
module Collision

using StaticArrays
using LinearAlgebra

using ..GfxBase
using ..Mesh

""" Collision
Objects which are used while checking for collisions. I model each particle as a Ray shooting 
    from position to new_position. If this ray intersects an object in my scene, I reflect the ray.
    HitRecord is used as a container of all useful intersection data. """
struct Ray 
    direction::Vector{Vec3}
    origin::Vector{Vec3}
end
mutable struct HitRecord
    t::Vector{Float64}
    intersect::Vector{Vec3}
    normal::Vector{Vec3}
end

function ray_intersects(object::OBJ, rays::Ray)

    if isfile(object.filename)
        # load .obj file
        mesh = Mesh.read_obj(object.filename)
        if !save_obj_files
            if object.rotate !== nothing
                Mesh.rotate!(mesh, object.rotate)
            end
            if object.scale !== nothing
                Mesh.scale!(mesh, object.scale)
            end
            if object.translate !== nothing
                Mesh.translate!(mesh, object.translate)
            end
        end
    else
        #TODO: run gen_mesh() with the string passed in.
        error("Filepath $(object.filename) does not exist.")
    end

    closest_hitrec = HitRecord(
            fill(Inf, length(rays.direction)), 
            [@SVector zeros(3) for _ in 1:length(rays.direction)], 
            [@SVector zeros(3) for _ in 1:length(rays.direction)])
    
    for (i, tri) in enumerate(mesh.triangles) #TODO: thread
        # Create triangle
        vertex = mesh.positions[tri.positions]
        uv = mesh.uvs[tri.uvs]
        normal = mesh.normals[tri.normals]
        scene_triangle = Triangle(vertex, uv, normal)
        
        # check ray intersects triangle
        mask, t, intersect, normal = ray_intersects(scene_triangle, rays)

        # Update HitRecord
        if !isempty(mask)
            compare = (t .< closest_hitrec.t[mask])
            closest_hitrec.t[mask[compare]] .= t[compare]
            closest_hitrec.intersect[mask[compare]] .= intersect[compare] #TODO: doesn't update.
            closest_hitrec.normal[mask[compare]] .= normal[compare]
        end
    end
    return closest_hitrec
end

function ray_intersects(tri::Triangle, rays::Ray)
    # insane speed. A little slower than Math.update_position_and_velocity()

    E1::Vec3 = tri.vertex[1] - tri.vertex[2]
    E2::Vec3 = tri.vertex[1] - tri.vertex[3]
    d::Vector{Vec3} = rays.direction
    rhs::Vector{Vec3} = Scalar(tri.vertex[1]) .- rays.origin

    # Calculate the determinant of a Vector of 3x3 matrix created by [E1 E2 d]
    # The 3x3 matrix is stacked and the determinant is calculated using broadcasting operations
    # det([E1 E2 d])
    M = E1[1] .* (E2[2] .* getindex.(d, 3) .- E2[3] .* getindex.(d, 2)) .-
        E2[1] .* (E1[2] .* getindex.(d, 3) .- E1[3] .* getindex.(d, 2)) .+
        getindex.(d, 1) .* (E1[2] * E2[3] - E1[3] * E2[2])
    mask = (M .!== 0.0) # ray parallel to plane

    # Calculate barycentric coordinates to ensure the intersection is within bounds of the triangle
    # Safely broadcast and return empty SVector{3, Float64}[] when the subset is empty
    rhs = something(rhs[mask], SVector{3, Float64}[])
    d = something(d[mask], SVector{3, Float64}[])
    M = something(M[mask], SVector{3, Float64}[])
    # det([rhs E2 d]) / M
    beta = (getindex.(rhs, 1) .* (E2[2] .* getindex.(d, 3) .- E2[3] .* getindex.(d, 2)) .-
            E2[1] .* (getindex.(rhs, 2) .* getindex.(d, 3) .- getindex.(rhs, 3) .* getindex.(d, 2)) .+
            getindex.(d, 1) .* (getindex.(rhs, 2) * E2[3] - getindex.(rhs, 3) * E2[2])) ./ M
    Bmask = (0.0 .<= beta .<= 1.0) # triangle intersection only occurs with beta ∃ (0,1)

    # Safely broadcast and return empty SVector{3, Float64}[] when the subset is empty
    rhs = something(rhs[Bmask], SVector{3, Float64}[])
    d = something(d[Bmask], SVector{3, Float64}[])
    M = something(M[Bmask], SVector{3, Float64}[])
    # det([E1 rhs d]) / M
    gamma = (E1[1] .* (getindex.(rhs, 2) .* getindex.(d, 3) .- getindex.(rhs, 3) .* getindex.(d, 2)) .-
            getindex.(rhs, 1) .* (E1[2] .* getindex.(d, 3) .- E1[3] .* getindex.(d, 2)) .+
            getindex.(d, 1) .* (E1[2] * getindex.(rhs, 3) - E1[3] * getindex.(rhs, 2))) ./ M
    Gmask = (0.0 .<= gamma .<= 1.0) # triangle intersection only occurs with gamma ∃ (0,1)
    
    # alpha = 1 - beta - gamma
    alpha = 1 .- beta[Bmask][Gmask] .- gamma[Gmask]
    Amask = (0.0 .<= alpha .<= 1.0)

    # Safely broadcast and return empty SVector{3, Float64}[] when rhs[Gmask][Amask] is empty
    rhs = something(rhs[Gmask][Amask], SVector{3, Float64}[])
    M = something(M[Gmask][Amask], SVector{3, Float64}[])

    # Calculate distance t along the ray and find the intersection point
    # det([E1 E2 rhs]) / M
    t = (E1[1] .* (E2[2] .* getindex.(rhs, 3) .- E2[3] .* getindex.(rhs, 2)) .-
        E2[1] .* (E1[2] .* getindex.(rhs, 3) .- E1[3] .* getindex.(rhs, 2)) .+
        getindex.(rhs, 1) .* (E1[2] * E2[3] - E1[3] * E2[2])) ./ M
    Tmask = (0.0 .< t .< 1.0) # TODO: tmin and tmax ??

    # If no intersection occurs, return empty Vectors
    if length(t[Tmask]) == 0
        return [], [], [], []
    end
    intersection = rays.origin[mask][Bmask][Gmask][Amask][Tmask] + t[Tmask] .* rays.direction[mask][Bmask][Gmask][Amask][Tmask]
    if interpolate_normals
        # interpolated normal
        normal = Ref(tri.normal[1]) .* alpha[Amask][Tmask] .+
                Ref(tri.normal[2]) .* gamma[Gmask][Amask][Tmask] .+
                Ref(tri.normal[3]) .* beta[Bmask][Gmask][Amask][Tmask]
    else
        # triangle normal
        v1 = tri.vertex[1] - tri.vertex[2]
        v2 = tri.vertex[1] - tri.vertex[3]
        normal = normalize(cross(v1, v2))
    end

    total_mask = collect(1:length(rays.direction))[mask][Bmask][Gmask][Amask][Tmask]
    return total_mask, t[Tmask], intersection, normal
end

function ray_intersects(sphere::Sphere, rays::Ray)

    hitrec = HitRecord(
            fill(Inf, num_particles), 
            [@SVector zeros(3) for _ in 1:num_particles], 
            [@SVector zeros(3) for _ in 1:num_particles])

    p::Vector{Vec3} = rays.origin .- Scalar(sphere.center)
    d::Vector{Vec3} = rays.direction

    discrim::Vector{Float64} = dot.(p, d).^2 .- dot.(d,d) .* (dot.(p, p) .- Scalar(sphere.radius^2))
    mask = findall(discrim .>= 0.0)

    t1 = (-dot.(p[mask], d[mask]) .+ sqrt.(discrim[mask])) ./ dot.(d[mask], d[mask])
    t2 = (-dot.(p[mask], d[mask]) .- sqrt.(discrim[mask])) ./ dot.(d[mask], d[mask])
    t = min.(t1, t2)

    # Make a new mask. A particle only collides with an object if t < 1.0
    mask = mask[1e-10 .< t .< 1.0] # TODO: tmin and tmax ??
    hitrec.t[mask] = t[1e-10 .< t .< 1.0]
    hitrec.intersect[mask] = hitrec.t[mask] .* d[mask] .+ p[mask] .+ Scalar(sphere.center)
    hitrec.normal[mask] = normalize.(hitrec.intersect[mask] .- Scalar(sphere.center))

    return hitrec
end

function handle_collisions(new_position::Vector{Vec3}, position::Vector{Vec3}, objects::Vector{<:Shape})
    save_new_pos = copy(new_position)
    # Initialize Objects
    rays = Ray(new_position .- position, position)
    global_mask = collect(1:num_particles)
    
    # Trace the path of each particle to determine if it intersects an object in the scene.
    while !isempty(rays.direction)
        # Initialize empty storage for ray intersection.
        closest_hitrec = HitRecord(
            fill(Inf, length(rays.direction)), 
            [@SVector zeros(3) for _ in 1:length(rays.direction)], 
            [@SVector zeros(3) for _ in 1:length(rays.direction)])
        # Iterate over every object
        for object in objects
            # Find intersection data
            hitrec = ray_intersects(object, rays) 
            # create a boolean mask of intersections which are closer than the current closest_hitrec
            mask = (hitrec.t .< closest_hitrec.t)
            # Update closest_hitrec for all positions which collided with an object sooner
            closest_hitrec.t[mask] .= hitrec.t[mask]
            closest_hitrec.intersect[mask] .= hitrec.intersect[mask]
            closest_hitrec.normal[mask] .= hitrec.normal[mask]
        end
        # Perform reflection for all rays which intersect an object
        collision_mask = findall(x -> x != Inf && x > 1e-6, closest_hitrec.t)
        # println(closest_hitrec.t[collision_mask])
        
        d = rays.direction[collision_mask]
        t = closest_hitrec.t[collision_mask]
        i = closest_hitrec.intersect[collision_mask]
        i = [round.(elem, digits=10) for elem in i]
        n = closest_hitrec.normal[collision_mask]
        r = (1 .- t) .* (d .- 2 .* n .* dot.(d, n)) # reflected rays

        # Update new_position and rays
        new_position[global_mask[collision_mask]] .= i .+ r
        global_mask = global_mask[collision_mask]

        rays = Ray(r, i)
    end 
    
    return new_position
end

end #module