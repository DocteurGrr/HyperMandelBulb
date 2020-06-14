# HyperMandelBulb

![alt text](https://github.com/DocteurGrr/HyperMandelBulb/blob/master/MandelBulb4D.png?raw=true)

This very simple programs performs **iterative Fractal computations in high dimensions**.

The formula used is just an extension of the Mandelbulb fractal :
(https://en.wikipedia.org/wiki/Mandelbulb)
by using N-shperical coordinates.
(https://en.wikipedia.org/wiki/N-sphere#Spherical_coordinates)

For simplicity all the world axis are uniformly discretized and the fractal is estimated on discrete hyperpoints.
In fact, in order to visualize sets of 3D volumes, this program addresses dimensions 4, 5 and 6 by fixing 1, 2 or 3 coordinates. 
For instance, in dimension 4, with a grid of N^4 hyperpoints (i.e. hypervoxels) computed, N volumes of N^3 voxels can be generated and visualized. 
They are "cuts" in 3D dimensions along one of the 4 axis of the hyper-fractal.

In dimension 5, the hyperfractal can be cut by (5 2)= 5! / (3! 2!) = 10 planes, which represents 10 x N^2 possible 3D volumes.
In dimension 6, it can be cut by (6 3) = 6! / (3! 3!) = 20 volumes, which represents 20 x N^3 possible 3D volumes.

Of course, each volume is a fractal and could be zoomed in by adapting the voxel grid size and step in order to visualize more details. 
A GPU implemtentation might lead to interactively explore such a tremendous amount of possibilities. 

This startup code is written in C++.

It has three steps :
- it computes the hyperfractal by iterations on a 3D voxel grid. The result is stored in Truncated Signed Distance Function  (TSDF);
- it gets a triangle mesh by marching cube algorithm;
- it visualizes the mesh by a ray tracing computation.
 
HyperMandelBuld is a command line tool that takes in input a single argument : a json file describing the different parameters.
To read this file nlohmann-json header only library is used . (https://nlohmann.github.io/json/)

 Its other dependencies are :
- Open3D for uniform tsdf container and marching cube algorithm. (http://www.open3d.org/)
- Ospray (and ospcommon) for mesh visualization either by Scivis or Path Tracing renderers. (https://www.ospray.org/)
- TBB for parallelized fractal iteration computations. Note : This part could benefit from simd instructions usage.








