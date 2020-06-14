To launch HyperMandelBuld program, simply call :
HypeMandelBulb.exe ConfigMandelB.json

Config file is rather self explanatory.

You should probably modify it in order to change dirOutput !

ImageFileName is the output image file name (without .png extension).

width and height are output image width and height.

RenderingQuality can be either 0 (low) or 1 (high)

dim = Dimension stands between 3 and 6

fixedCoords describes with a triplet wich coordinates are fixed. Correct values are between 0 and dim-1, -1 stands for "ignore me", dim-3 values should not be -1.
For example, in 4D, to fix the second coordinate (y) use: [1, -1, -1]. Then a 3D volume depending on x, z and w will be computed.
In 5D, fixedCoords could be for example |1, 4, -1] in order to fix the second and fith coordinate.

valueCoords are the value of the fixed coordinates. That's also a triplet. Only 0 (dim = 3) up to 3 (dim = 6) of these values are used. 
The others should have value -1. Correct values are between 0 and N-1.

power is the fractal power (see Mandelbuld formulae). 

maxIter is the maximum iteration number (should be between 1 and 255)

bailout is the second iteration stopping criterion (see Mandelbuld formulae)

min and extent provides the Axis Aligned BoundingBox parameters of the 3d regular grid 

N is the number of grid point per axis (the computed grid has N^3 voxels)

 
