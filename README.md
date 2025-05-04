# WireSculpt
<img src ="https://github.com/user-attachments/assets/f5162aa6-5210-4ac5-bf3d-00c368c34179" width="400" height="350">

## Details
This plugin is based on the following SIGGRAPH paper published in 2021: WireRoom: Model-guided Explorative Design of Abstract Wire Art by Zhijin Yang, Pengfei Xu, Hongbo Fu, and Hui Huang.

This plugin was created by Claire Lu and Neha Thumu as a part of CIS660 at the University of Pennsylvania. 

Our process can be found in this [design doc](https://github.com/thumun/WireSculpt/blob/neha_code_cleanup/WireSculpt%20Design%20Doc.pdf).

## How to install and use Wiresculpt
- First, make sure you have Maya 2022 or Maya 2024 installed (we have not tested this plug in for other versions). We also recommend Visual Studio 2022 to build the plugin.
- Clone our repo, then open it in Visual Studio in **admin mode** (without this, the project will **not** build properly). Once you build the project, the .mll file will be added to the debug folder of the project. To load the plugin, open your version of Maya and go to Windows->Settings/Preferences->Plug-in Manager and browse to open the file. 
<img src="https://github.com/user-attachments/assets/7f9e323c-2ce6-48ed-91f4-5b060ca0d886" width="300" height="300">

- You can launch our GUI through the menu at the top of the Maya UI or by typing 'CreateWireSculpt;' into the script editor.
![image](https://github.com/user-attachments/assets/8a2a553a-3e11-473e-aedc-5255254671c5)

To load a model into the plugin, click the load from file button. The 3D mesh has to be **triangulated**, **manifold**, and **obj**, or the plugin will **not** be able to process it.

<img src="https://github.com/user-attachments/assets/6af2361b-f5d7-4db1-b6ee-ad1244ec42a6" width="400" height="500">

The Abstraction option enables the wireframe sculpture to have a more abstract aesthetic by enabling more extreme points to be chosen as landmark vertices that get visited by the wire path algorithm. The extreme points parameters can be adjusted to change the verticies that are selected as extreme points. 
- Proximity: Controls how close new extreme points can be to existing ones
- MaxVal: threshold for determining whether a vertex qualifies as an extreme point based on its Laplacian value relative to its neighbors
- Filter threshold: Controls how many areas are considered concave (concave verticies are not considered extreme)

The feature attraction weights determine how closely the path follows the features.
 - Range (a): controls the strength of the attraction
 - Steepness (b): controls the range/locality of the attraction

The path repulsion weights aid in avoiding repetitions of verticies in the path.
- Range (a*): controls the strength of the repulsion
- Steepness (b*): controls the locality of repulsion

NURBs Curve Thickness parameter adjusts the width of the output line.

Contours are view dependent so it's shaped by the perspective.
- Front: captures contours from the front of the model 
- Side: captures contours from the side of the model

FOV (field of vision) is another parameter for affecting the look of the contours. It controls the field of view for which the occluding and suggestive contours get captured. By default, it is set at 0.7 which has worked well for most models.

Suggestive Contours can be switched on/off depending on wire output preference. By default, occluding contours and suggestive contours are used. Occluding contours capture the view-dependent silhouette of the model and suggestive contours capture finer details. 
- Trim SC: the amount that the suggestive contours are trimmed.

The Feature Vertices option allows control over how much of the models contours are captured in the wire path. A value of 1 adds all vertices from the contours to the path, making the wire sculpture match more closely with the model and appear more intricate, while 0 will cause the wire path to be solely constructed by the chosen extreme points, making the wire sculpture more minimalistic.

## Breakdown of Creating this Plugin 
### Maya Node Setup and General Code Setup 
- We setup the node following the examples provided by Autodesk. Our parameters are read in from the GUI and used as appropriate within our code. When the generate button is clicked and the parameters have been altered, the mesh will be recomputed based on our algorithms. Then the wire will be output to Maya. 
- We organized our code into plugin related classes (node information, helper files for outputting the wire) and classes to hold our data (verticies, edges, and faces). 
- We wrote a file reading method to process OBJ files and create our verticies, edges, and faces arrays (to be used in our algorithms). 
### Extreme Points Integration 
- We used the following repository in order to detect extreme points of meshes that are inputted by the user: https://github.com/FlorianTeich/concavity-aware-fields/tree/main
- We altered the code to parameterize the extreme points as the values used in the original implementation did not generalize to all inputs. 
- We also integrated necessary source code from the libigl library. (This can be found in our libigl and suitesparse folders) 
### A* / 2-Opt Pathfinding 
- We referenced the wikipedia page: https://en.wikipedia.org/wiki/A*_search_algorithm
- We used the following repository as a basis for our A* algorithm: https://github.com/JDSherbert/A-Star-Pathfinding/blob/main/Pathfinder.cpp
- We referenced the following repository for our 2-Opt Traveling Salesman Approximation: https://github.com/jimousse/tsp2D/blob/master/main.cpp
- The A* pathfinding algorithm is used to find the shortest path from one landmark vertex to another.
- The 2-Opt algorithm is used to find the optimal cycle of landmark vertices, i.e. the order in which landmark vertices are visited.
### Heat Map Distance 
- We originally translated this method: https://github.com/CHoudrouge4/HeatMethodForGeodesicDistance into C++ however this did not yield results that we could use. 
- We instead used the heat map method from the libigl library. 
### Contour Feature lines 
- We used the following source code to compute the contours of the model: https://gfx.cs.princeton.edu/gfx/proj/sugcon/
### Integration 
- We first read the file in using our process logic then computed the extreme points. These are set as our landmark verticies. 
- Then we get the feature lines of the mesh (the contours). 
- We then compute the heat map distance while using the feature lines as the sources.
- Next, we compute the feature attraction weights (so the wire path follows our features lines)
- The next step requires multiple passes of a pathfinding algorithm in order to properly get our output wire. After the first path, we calculate path repulsion weights in order to avoid backtracking of verticies. We then run A* once more to get our output.
- Our output is then shown in Maya!
## Resources
1. Generating a unit sphere (SphereMesh class): https://github.com/Erkaman/cute-deferred-shading/blob/master/src/main.cpp#L573
2. Extreme Points Resources:
code from here: https://github.com/FlorianTeich/concavity-aware-fields/tree/main
3. Uses this library (original linux version): https://github.com/DrTimothyAldenDavis/SuiteSparse
4. Used the windows version in this project: https://github.com/jlblancoc/suitesparse-metis-for-windows
5. Uses this library for heat map distance method and extreme points matrix computations: https://libigl.github.io/tutorial/
6. Nearest Neighbor 2-Opt TSP Approximation: https://github.com/jimousse/tsp2D/blob/master/main.cpp
7. A* Implementation Reference: https://github.com/JDSherbert/A-Star-Pathfinding/blob/main/Pathfinder.cpp
8. Wikipedia A* Implementation Reference: https://en.wikipedia.org/wiki/A*_search_algorithm
9. Suggestive Contours: https://gfx.cs.princeton.edu/gfx/proj/sugcon/
## Additional Notes 
### Models Used 
- Fox: https://skfb.ly/6GIrP
- Seagull: https://poly.pizza/m/6Tpj_vcWP3f
- Giraffe: https://sketchfab.com/3d-models/low-poly-giraffe-e9715aee26bd4e03b2740547b5b259d0
- Elephant: https://sketchfab.com/3d-models/low-poly-elephant-84fd98c561464b1ba5d3dd48ab161b9c
- Sailboat (sail alteratioons by Neha): https://poly.pizza/m/BgSZXwmm7k
- Saul Goodman Figure: https://sketchfab.com/3d-models/saul-goodman-531a84899eb44401a1ff5d8f735aa6ad
- Saul Goodman Head: https://www.artstation.com/artwork/blKOkn - saul head dineth
- Blender Add-on to make meshes manifold: https://extensions.blender.org/add-ons/print3d-toolbox/
- Metal shader: https://help.autodesk.com/view/MAYAUL/2025/ENU/?guid=GUID-426580A4-6079-46A2-B651-E3CE38A5DEB8 

### Configuration Settings
- exedebug = way to debug our code as it builds an exe instead of an mll (just need a main and a if EXEDEBUG line encompassing the main)
- xcopy -> custom post build events where we are copying the libs that we need to the maya binary and to our debug folder so everything builds/runs properly
- admin mode -> need to open VS in this so xcopy has permission to copy to Maya app files
- %MAYA_LOCATION% -> needs this system variable to be set to device's Maya path for xcopy to work
