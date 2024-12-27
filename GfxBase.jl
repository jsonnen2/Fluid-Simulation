module GfxBase

export Vec3, Vec2
export OBJ, OBJTriangle, OBJMesh, Shape, Box, Sphere, Triangle, Camera
export objects, num_particles, delta_time, simulation_steps, smoothing_radius, bounding_box
export interpolate_normals, cap_acceleration, camera, width, height, global_filepath

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
struct OBJ <: Shape
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
struct Box <: Shape
    min::Vec3
    max::Vec3
end
struct Sphere <: Shape
    center::Vec3
    radius::Float64
end

# also cell_size
global smoothing_radius::Float64 = 400.0

# bounding box for the fluid
global bounding_box = Box(Vec3(0,0,0), Vec3(50,50,50))
# TODO: default not nothing. 50 cell width : 4 radius
bounding_box_center = Scalar(0.5) .* (bounding_box.max .- bounding_box.min) # TODO: OBJ default values are from this equation

# Objects in the scene
global objects = [
    # TODO: User shouldn't computer bb_center. Set this to default
    OBJ("Mesh/inside_box.obj", bounding_box_center, bounding_box_center, nothing)
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
global camera = Camera(Vec3(0,0,10), Vec3(0,0,1), Vec3(0,1,0), 1.0)
# width, height of image in pixels
global width = 400
global height = 400

# maximum magnitude of acceleration allowed by the particles
global cap_acceleration = 100

# Real Time Simulation
global delta_time::Float64 = 1/60
global num_particles::Int = 10
global simulation_steps::Int = 10
global global_filepath::String = "C:/J/Unity/Assets/Saved_Fluid_Sims/particle_10.txt"


end # module GfxBase
