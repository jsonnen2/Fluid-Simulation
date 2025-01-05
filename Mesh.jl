
""" OBJ file read/write code was written by Scott Wehrwein of Western Washington University. """

#COMPLETE: translate, scale, rotate triangle meshes.

module Mesh

export read_obj, write_obj, gen_mesh
export translate, scale, rotate

# push!(LOAD_PATH, pwd())
# include("GfxBase.jl")
# using .GfxBase
using ..GfxBase

using StaticArrays
using LinearAlgebra

""" read_obj(obj_filename)
Read a mesh in OBJ format from file obj_filename."""
function read_obj(obj_filename)
    m = OBJMesh([], [], [], []) # create a mesh
    open(obj_filename) do f
        for (line_number, line) in enumerate(eachline(f))
            if line == "" || line[1] == "#"
                continue # skip comments
            end
            # Read the line and add its contents to the correct field of m:
            tokens = split(strip(line))
            if tokens[1] == "v" # vertex
                push!(m.positions, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vt" # vertex texture
                push!(m.uvs, Vec2([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "vn" # vertex normal
                push!(m.normals, Vec3([parse(Float64, x) for x in tokens[2:end]]...))
            elseif tokens[1] == "f"
                # create a OBJTriangle face:
                points = []
                uvs = []
                normals = []
                # handle faces with no texture and/or normals
                for corner in tokens[2:end]
                    indices = split(corner, '/')
                    if length(indices) == 3 # all 3 present, third is a normal
                        push!(normals, parse(Int, indices[3]))
                    end
                    if length(indices) >= 2 && indices[2] != ""
                        # if there are 2 or more and the second isn't blank, it's a texture
                        push!(uvs, parse(Int, indices[2]))
                    end
                    if length(indices) >= 1 # first value is the position
                        push!(points, parse(Int, indices[1]))
                    else # unless it has none, in which case it's not valid
                        error("in line $line_number: face vertex $corner could not be parsed")
                    end
                end
                # create the triangle and add it to the triangles array
                push!(m.triangles, OBJTriangle(points, uvs, normals))
            end
        end
    end
    return m
end

""" write_obj(obj_filename)
Write the given mesh in OBJ format to file obj_filename."""
function write_obj(obj_filename, mesh::OBJMesh)
    open(obj_filename, "w") do f
        # write all positions:
        for v in mesh.positions
            write(f, "v $(v[1]) $(v[2]) $(v[3])\n")
        end
        # write all texture coords:
        for v in mesh.uvs
            write(f, "vt $(v[1]) $(v[2])\n")
        end
        # write all normals:
        for v in mesh.normals
            write(f, "vn $(v[1]) $(v[2]) $(v[3])\n")
        end
        # write all triangles:
        for tri in mesh.triangles
            write(f, "f $(tri_vertex_str(tri))\n")
        end
    end
end

""" tri_vertex_str(triangle)
Helper function for write_obj()
Return a string with the indices of applicable positions, texture coordinates,
and normals for a given triangle according to the OBJ specification.
In particular, if p, u, and n are position, vertex and normal, each corner
of the triangle is represented as one of the following:
    p       (position only)
    p/u     (position and texture)
    p//n    (position and normal)
    p/u/n   (position, texture, and normal) """
function tri_vertex_str(triangle::OBJTriangle)
    # determine whether textures and normals are present:
    write_uv = length(triangle.uvs) == length(triangle.positions)
    write_normals = length(triangle.normals) == length(triangle.positions)
    corners = []
    for i = 1:3
        output = "$(triangle.positions[i])"
        if write_uv && !write_normals
            output = output * "/$(triangle.uvs[i])" # * does concatenation(!)
        elseif !write_uv && write_normals
            output = output * "//$(triangle.normals[i])"
        elseif write_uv && write_normals
            output = output * "/$(triangle.uvs[i])/$(triangle.normals[i])"
        end
        push!(corners, output)
    end
    join(corners, " ")
end

function translate!(mesh::OBJMesh, translation::Vec3)
    mesh.positions .+= Scalar(translation)
end

function scale!(mesh::OBJMesh, multiplier::Vec3)
    mesh.positions .= (.*).(mesh.positions, Scalar(multiplier))
    mesh.normals = normalize.((.*).(mesh.normals, Scalar(multiplier)))
end

function rotate!(mesh::OBJMesh, rotation::Tuple{Float64, Vec3, Vec3})
    # Rotate object θ degrees about an axis of rotation at a specified center.
    # Rodrigues' Formula
    # v_rot = v * cosθ + cross(axis, v) * sinθ + axis * dot(axis, v) * (1 - cosθ)

    θ, axis, center = rotation
    axis = normalize(axis)

    v = mesh.positions
    v .-= Scalar(center)
    v_rot = v .* Scalar(cos(θ)) .+ cross.(Scalar(axis), v) .* Scalar(sin(θ)) .+ Scalar(axis) .* dot.(Scalar(axis), v) .* Scalar(1 - cos(θ))
    v_rot .+= Scalar(center)
    mesh.positions .= v_rot

    v = mesh.normals
    v .-= Scalar(center)
    v_rot = v .* Scalar(cos(θ)) .+ cross.(Scalar(axis), v) .* Scalar(sin(θ)) .+ Scalar(axis) .* dot.(Scalar(axis), v) .* Scalar(1 - cos(θ))
    v_rot .+= Scalar(center)
    mesh.normals .= v_rot
end

""" gen_mesh(outfile, geom, divisionsU, divisionsV)
Generate a mesh and save the result in a file with name outfile.
geom may be "cube", "cylinder", or "sphere".
Cylinder requires divisionsU; sphere requires divisionsU and divisionsV. """
function gen_mesh(outfile, geom, divisionsU=0, divisionsV=0)
    if geom == "cube"
        mesh = inside_cube_mesh()
    elseif geom == "cylinder"
        mesh = cylinder_mesh(divisionsU)
    elseif geom == "sphere"
        mesh = sphere_mesh(divisionsU, divisionsV)
    elseif geom == "torus"
        mesh = torus_mesh(divisionsU, divisionsV)
    end
    write_obj(outfile, mesh)
end

""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. """
function cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3(1, -1, -1)) # 1 RBC
    push!(positions, Vec3(1, -1, 1)) # 2 RBF
    push!(positions, Vec3(-1, -1, 1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3(1, 1, -1)) # 5 RTC
    push!(positions, Vec3(1, 1, 1)) # 6 RTF
    push!(positions, Vec3(-1, 1, 1)) # 7 LTF
    push!(positions, Vec3(-1, 1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3(1, 0, 0)) # R
    push!(normals, Vec3(-1, 0, 0)) # L
    push!(normals, Vec3(0, 1, 0)) # U
    push!(normals, Vec3(0, -1, 0)) # D
    push!(normals, Vec3(0, 0, 1)) # C
    push!(normals, Vec3(0, 0, -1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1, 2, 3], [1, 2, 3], [4, 4, 4])) # bottom face 1
    push!(triangles, OBJTriangle([1, 3, 4], [1, 3, 4], [4, 4, 4])) # bottom face 2
    push!(triangles, OBJTriangle([1, 5, 6], [4, 1, 2], [1, 1, 1])) # right face 1
    push!(triangles, OBJTriangle([1, 6, 2], [4, 2, 3], [1, 1, 1])) # right face 2
    push!(triangles, OBJTriangle([2, 6, 7], [4, 1, 2], [5, 5, 5])) # far face 1
    push!(triangles, OBJTriangle([2, 7, 3], [4, 2, 3], [5, 5, 5])) # far face 2
    push!(triangles, OBJTriangle([3, 7, 8], [2, 3, 4], [2, 2, 2])) # left face 1
    push!(triangles, OBJTriangle([3, 8, 4], [2, 4, 1], [2, 2, 2])) # left face 2
    push!(triangles, OBJTriangle([4, 8, 5], [2, 3, 4], [6, 6, 6])) # far face 1
    push!(triangles, OBJTriangle([4, 5, 1], [2, 4, 1], [6, 6, 6])) # far face 2
    push!(triangles, OBJTriangle([5, 8, 7], [1, 2, 3], [3, 3, 3])) # top face 1
    push!(triangles, OBJTriangle([5, 7, 6], [1, 3, 4], [3, 3, 3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end

""" cube_mesh()
Return a new OBJMesh representing a 2x2x2 cube centered at the origin and
axis-aligned. Surface normals point inside the cube because it is a boundary box. """
function inside_cube_mesh()
    positions = []
    uvs = []
    normals = []
    triangles = []
    # key to comments:
    # L/R = x = right/left
    # B/T = y = top/bottom
    # C/F = z = close/far
    push!(positions, Vec3(1, -1, -1)) # 1 RBC
    push!(positions, Vec3(1, -1, 1)) # 2 RBF
    push!(positions, Vec3(-1, -1, 1)) # 3 LBF
    push!(positions, Vec3(-1, -1, -1)) # 4 LBC
    push!(positions, Vec3(1, 1, -1)) # 5 RTC
    push!(positions, Vec3(1, 1, 1)) # 6 RTF
    push!(positions, Vec3(-1, 1, 1)) # 7 LTF
    push!(positions, Vec3(-1, 1, -1)) # 8 LTC

    # texture coordinates:
    push!(uvs, Vec2(1, 1)) # TR
    push!(uvs, Vec2(0, 1)) # TL
    push!(uvs, Vec2(0, 0)) # BL
    push!(uvs, Vec2(1, 0)) # BR

    # normals:
    push!(normals, Vec3(-1, 0, 0)) # R
    push!(normals, Vec3(1, 0, 0)) # L
    push!(normals, Vec3(0, -1, 0)) # U
    push!(normals, Vec3(0, 1, 0)) # D
    push!(normals, Vec3(0, 0, -1)) # C
    push!(normals, Vec3(0, 0, 1)) # F

    # 8 faces, 2 triangles each
    push!(triangles, OBJTriangle([1, 2, 3], [1, 2, 3], [4, 4, 4])) # bottom face 1
    push!(triangles, OBJTriangle([1, 3, 4], [1, 3, 4], [4, 4, 4])) # bottom face 2
    push!(triangles, OBJTriangle([1, 5, 6], [4, 1, 2], [1, 1, 1])) # right face 1
    push!(triangles, OBJTriangle([1, 6, 2], [4, 2, 3], [1, 1, 1])) # right face 2
    push!(triangles, OBJTriangle([2, 6, 7], [4, 1, 2], [5, 5, 5])) # far face 1
    push!(triangles, OBJTriangle([2, 7, 3], [4, 2, 3], [5, 5, 5])) # far face 2
    push!(triangles, OBJTriangle([3, 7, 8], [2, 3, 4], [2, 2, 2])) # left face 1
    push!(triangles, OBJTriangle([3, 8, 4], [2, 4, 1], [2, 2, 2])) # left face 2
    push!(triangles, OBJTriangle([4, 8, 5], [2, 3, 4], [6, 6, 6])) # far face 1
    push!(triangles, OBJTriangle([4, 5, 1], [2, 4, 1], [6, 6, 6])) # far face 2
    push!(triangles, OBJTriangle([5, 8, 7], [1, 2, 3], [3, 3, 3])) # top face 1
    push!(triangles, OBJTriangle([5, 7, 6], [1, 3, 4], [3, 3, 3])) # top face 2

    # julia automatically returns the last value in the function:
    OBJMesh(positions, uvs, normals, triangles)

end


""" cylinder_mesh(n)
Return a new OBJMesh object approximation of a cylinder with radius 1 and
height 2, centered at the origin. The logitudinal axis is aligned with y, and
it is tesselated with n divisions arranged radially around the outer surface.
The ends of the cylinder are disc-shaped caps parallel to the xz plane. See the
assignment writeup for a diagram and details.
"""
function cylinder_mesh(divisionsU)
    # TODO - feel free to drop in your A1 solution code here
end


""" sphere_mesh(n, m)
Create a Latitude-Longitude-tesselated approximation of a sphere with radius 1
centered at the origin. There are n divisions around the equator and m
divisions from pole to pole along each line of longitude. The North pole is at
(0,1,0), the South pole at (0,-1,0), and points on the Greenwich meridian are
in the x = 0 plane with z > 0. The u texture coordinate depends on longitude,
with u=0 at 180 degrees West and u=1 at 180 degrees East. The v texture
coordinate varies with latitude with v=0 at the South pole and v=1 at the North
pole. Normals should be normal to the ideal sphere surface. See the assignment
for a diagram and further details. """
function sphere_mesh(N, M)
   # TODO
end

""" torus_mesh(n, m)
A torus is a doughnut-shaped surface defined by a major radius, affecting the size 
of the hole, and a minor radius, affecting the thickness of the ring. This code creates 
a torus with major radius 1 and minor radius r (controlled by an additional -r flag 
with a default of 0.25). Its u coordinates are like the sphere, and the v coordinate 
runs from 0 to 1 around the inside of the torus, with the direction arranged so that 
the texture is right-reading from the outside (i.e., the texture is not flipped when 
mapped to the surface). Like the sphere, it has a seam on the −z half of the yz-plane, 
and it has a similar seam around the inner surface of the doughnut hole; vertices along 
each seam share a pair of texture coordinates, and single vertex, at the position (0,0,r−1) 
where the seams meet, shares 4 texture coordinates. """
function torus_mesh(N, M)
    # N = MAJOR axis # points around the torus
    # M = MINOR axis # points on the tube
    # begin at Vec3(1.0, 0.0, 0.0) 

    inner_radius = 2.0
    outer_radius = 0.5
    
    normals = Matrix{Vec3}(undef, N, M)
    positions = Matrix{Vec3}(undef, N, M)
    
    domain = collect(1:N)
    center_line = Scalar(inner_radius) .* Vec3.(cos.(2*pi .* domain ./ N), 0, sin.(2*pi .* domain ./ N))
    # v is the center_line scaled to the radius of the tube
    v = center_line .* Scalar(outer_radius / inner_radius)
    # k is the axis of rotation in Rodrigues' formula
    k = normalize.(cross.(v, Scalar(Vec3(0.0, 1.0, 0.0))))
    for i in 1:M
        θ = 2*pi/M * i
        # Rodrigues' rotation formula
        # v_rot = v * cosθ + cross(k, v)*sinθ + k*dot(k, v)*(1 - cosθ)
        v_rot = v .* Scalar(cos(θ)) .+ cross.(k, v) .* Scalar(sin(θ)) .+ k .* dot.(k, v) .* Scalar(1 - cos(θ))
        positions[:, i] .=  v_rot .+ center_line
        normals[:, i] .= v_rot
    end
    domain = collect(1:M*N)
    x_plus_one = domain .+ Scalar(1-M) .+ Scalar(M) .* min.(domain .% M, 1)
    x_plus_M = (domain .+ M .- 1) .% (M*N) .+ 1
    x_plus_one_and_M = (x_plus_one .+ M .- 1) .% (M*N) .+ 1
    upper_tri = Vec3.(domain, x_plus_one, x_plus_M)
    lower_tri = Vec3.(x_plus_one, x_plus_one_and_M, x_plus_M)
    indices = vcat(upper_tri, lower_tri)

    tris = [OBJTriangle(i, Vec3(0.0,0.0,0.0), i) for i in indices]
    OBJMesh(vec(positions'), [Vec2(0,0)], vec(normals'), tris)
end

end #module