# Results

*Collision with a triangle mesh:*

![](Images/bunny_collision.gif)   


# Unity Instructions

I cloned my Assets folder. To view my particles inside Unity there are some steps to prepare your editor. First, you must create an empty GameObject and connect the `spawn_particles.cs` script to it. Next, ensure the variables on this GameObject are correct. The prefab must point to PointPrefab from the Prefabs folder. And filepath must point to the saved .txt file generated from `fluid_sim.jl`. Finally, navigate to Edit -> Project Settings -> Tags and Layers then add a "Particle" tag.