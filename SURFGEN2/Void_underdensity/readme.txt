======================Readme for SURFGEN2 (for studying underdensity) ===========================================


<<<<<<<<<<<< input file >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
same as used in the percolation case.

<<<<<<<<<<<< usage >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
1. Compile using 'make'
2. usage: ./main <inputfilename> rho_th count_th NN_check <triangle_path> <mf_path>
./main = executable
<inputfilename>: give the name (with the path if needed) of the input file
rho_th (double): the isodensity (threshold value). SURFGEN2 will construct the surface using rho=rho_th (any scaler field in general). This code will study the individual underdensity regions (separated surfaces).
count_th (intiger): an isolated region required to have this minimum count of grid points/cells to be considered for triangulation and MFs/SFs calculation by SURFGEN2. The cluster statistics will not be affected by this though. 
NN_check (positive intiger): Cluster number to be saved (grid points as well as triangles)
<triangle_path>: path (and prefix) of the trianle file to be saved (Note NN_check has to entered properly)
<mf_path>: path (and prefix) of the MFs/SFs file to be saved. This is the main main result.
