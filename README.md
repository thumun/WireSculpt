# WireSculpt

## Details
This plugin is based on the following SIGGRAPH paper published in 2021: WireRoom: Model-guided Explorative Design of Abstract Wire Art by Zhijin Yang, Pengfei Xu, Hongbo Fu, and Hui Huang.

This plugin was created by Claire Lu and Neha Thumu as a part of CIS660 at the University of Pennsylvania. 

## How to install and use our Plugin 
- First, make sure you have Maya 2022 or Maya 2024 installed (we have not tested this plug in for other versions). We also recommend Visual Studio 2022 to build the plugin.
- Clone our repo, then open it in Visual Studio in **admin mode** (without this, the project will **not** build properly). Once you build the project, the .mll file will be added to the debug folder of the project. Now, open your version of Maya and go to the plug in window (..)
- You can launch our GUI through the x window () or by opening the script editor and manually opening our wireUI.mel file and launching it from there.

INSERT IMAGE OF GUI WITH COLOR CODED BOXES!!

## Breakdown of Creating this Plugin 
### Maya Node Setup and General Code Setup 
- general Maya stuff
- file reading
- classes 
### Extreme Points Integration 
- used florian's repo as basis
- integrated necessary src code from libigl library
### A* Pathfinding 
- //
### Heat Map Distance 
- originally used CHoudrouge4's repo as basis 
- src code from libigl
### Contour Feature lines 
- //
### Integration 
- // 

## Resources
1. Generating a unit sphere (SphereMesh class): https://github.com/Erkaman/cute-deferred-shading/blob/master/src/main.cpp#L573
2. Extreme Points Resources:
code from here: https://github.com/FlorianTeich/concavity-aware-fields/tree/main
3. Uses this library (original linux version): https://github.com/DrTimothyAldenDavis/SuiteSparse
4. Used the windows version in this project: https://github.com/jlblancoc/suitesparse-metis-for-windows
5. Uses this library for heat map distance method and extreme points matrix computations: https://libigl.github.io/tutorial/
7. Nearest Neighbor 2-Opt TSP Approximation: https://github.com/jimousse/tsp2D/blob/master/main.cpp
9. A* Implementation Reference: https://github.com/JDSherbert/A-Star-Pathfinding/blob/main/Pathfinder.cpp
11. Wikipedia A* Implementation Reference: https://en.wikipedia.org/wiki/A*_search_algorithm

## Additional Notes 
### Configuration Settings
- exedebug = way to debug our code as it builds an exe instead of an mll (just need a main and a if EXEDEBUG line encompassing the main)
- xcopy -> custom post build events where we are copying the libs that we need to the maya binary and to our debug folder so everything builds/runs properly
- admin mode -> need to open VS in this so xcopy has permission to copy to Maya app files
- %MAYA_LOCATION% -> needs this system variable to be set to device's Maya path for xcopy to work
