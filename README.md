# WireSculpt


## Directions

## Configuration Notes
- exedebug = way to debug our code as it builds an exe instead of an mll (just need a main and a if EXEDEBUG line encompassing the main)
- xcopy -> custom post build events where we are copying the libs that we need to the maya binary and to our debug folder so everything builds/runs properly
- admin mode -> need to open VS in this so xcopy has permission to copy to Maya app files
- %MAYA_LOCATION% -> needs this system variable to be set to device's Maya path for xcopy to work

## Resources

1. Generating a unit sphere (SphereMesh class): https://github.com/Erkaman/cute-deferred-shading/blob/master/src/main.cpp#L573

2. Extreme Points Resources:
code from here: https://github.com/FlorianTeich/concavity-aware-fields/tree/main

3. Uses this library (original linux version): https://github.com/DrTimothyAldenDavis/SuiteSparse

4. Used the windows version in this project: https://github.com/jlblancoc/suitesparse-metis-for-windows

5. Uses this library: https://libigl.github.io/tutorial/
6. Nearest Neighbor 2-Opt TSP Approximation: https://github.com/jimousse/tsp2D/blob/master/main.cpp
7. A* Implementation Reference: https://github.com/JDSherbert/A-Star-Pathfinding/blob/main/Pathfinder.cpp
8. Wikipedia A* Implementation Reference: https://en.wikipedia.org/wiki/A*_search_algorithm
