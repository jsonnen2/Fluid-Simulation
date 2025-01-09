module GfxBase

# TODO: I wish I could just set all global, const, struct in this file to be accessible globally without this export nonsense.
export Vec3, Vec2
export OBJ, OBJTriangle, OBJMesh, Shape, Box, BoundingBox, Sphere, Triangle, Camera
export objects, num_particles, delta_time, simulation_steps, smoothing_radius, bounding_box, console_log
export interpolate_normals, cap_acceleration, camera, width, height, global_filepath, obj_save_folder, spawn_cube, spawn_type
export save_obj_files

using StaticArrays

# some type aliases:
const Vec3 = SVector{3, Float64} # 3-vector of floats
const Vec2 = SVector{2, Float64} # 2-vector of floats

""" OBJTriangle
A struct that represents a single triangle in a mesh. """
mutable struct OBJTriangle
    positions::Array{Int,1} # vertex position indices
    uvs::Array{Int,1} # vertex texture coordinate indices
    normals::Array{Int,1} # normal vector indices
end

""" OBJMesh
A struct that represents an indexed triangle mesh for reading from
or writing to OBJ format. """
mutable struct OBJMesh
    positions::Array{Vec3,1} # all vertex positions
    uvs::Array{Vec2,1} # all texture coordinates
    normals::Array{Vec3,1} # all vertex normals
    triangles::Array{OBJTriangle,1} # the OBJTriangles belonging to the mesh
end

struct Triangle # lack of inheritance from Shape is purposeful.
    vertex::Vector{Vec3}
    uv::Vector{Vec2}
    normal::Vector{Vec3}
end

""" Shape
Every type of object which can be placed in my scene. """
abstract type Shape end
mutable struct OBJ <: Shape
    filename::String
    translate::Union{Vec3, Nothing}
    scale::Union{Vec3, Nothing}
    rotate::Union{Tuple{Float64, Vec3, Vec3}, Nothing} # degrees, rotation axis, rotation center

    # Constructor for full OBJ
    OBJ(filename::String, translate::Union{Vec3, Nothing}, scale::Union{Vec3, Nothing}, rotate::Union{Tuple{Float64, Vec3, Vec3}, Nothing}) = 
        new(filename, translate, scale, rotate)

    # Constructor for simple OBJ
    OBJ(filename::String) = new(filename, nothing, nothing, nothing)
end
struct Box <: Shape # Surface normals point outside to keep particles out
    min::Vec3
    max::Vec3
end
struct BoundingBox <: Shape # Surface normals point inside the box to keep particles in
    min::Vec3
    max::Vec3
end
struct Sphere <: Shape
    center::Vec3
    radius::Float64
end


global spawn_cube = Box(Vec3(15,55,15), Vec3(35,70,35))
global spawn_type = "uniform"

global bounding_box = BoundingBox(Vec3(0,0,0), Vec3(150,70,150)) # bounding box for the fluid
global smoothing_radius::Float64 = 4.0 # hashing cell size
bounding_box_center::Vec3 = Scalar(0.5) .* (bounding_box.max .- bounding_box.min)

# Objects in the scene
global objects = [
    # OBJ("Mesh/inside_box.obj", bounding_box_center, bounding_box_center, nothing)
    bounding_box
    OBJ("Mesh/bunny.obj", Vec3(25,20,25), Vec3(30,30,30), nothing)
    # OBJ("Mesh/bunny.obj", Vec3(0,0,0), Vec3(10,10,10), nothing)
]
# toggle to use triangle surface normals, or to interpolate normals using the barycentric coordinates
global interpolate_normals = true

# Camera to view the scene
mutable struct Camera
    origin::Vec3
    view::Vec3
    up::Vec3
    focal::Float64
end
global camera = Camera(Vec3(25,25,-10), Vec3(0,0,1), Vec3(0,1,0), 1.0)
# width, height of image in pixels
global width = 400
global height = 400

# maximum magnitude of acceleration allowed by the particles
global cap_acceleration = 100

# Real Time Simulation
global delta_time::Float64 = 1/20
global num_particles::Int = 10000
global simulation_steps::Int = 5000

# File stuff
local_file = "$spawn_type-$num_particles"
local_file = "splish_splash_f"
global global_filepath::String = "C:/J/Unity/3D project/Assets/Saved_Fluid_Sims/$local_file.txt"
global obj_save_folder::String = "C:/J/Unity/3D project/Assets/Saved_Fluid_Sims/$local_file"
global save_obj_files = true
global console_log = true

end # module GfxBase
