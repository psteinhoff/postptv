This "package" is a small tool to convert your ptv.is-files into .vtk-files (unstructured grid)
The .vtk-files can be loaded into paraview, a powerfull post-processing tool.

add velo-trans.par, 3d_point_offset.par and ptv2vtk.config to your parameter folder.

To run the tool you have to options, a and b:

(a) Change the working folder in line 33 of ptv_is_unstructured_vtk.c to match youur working folder and compile it.
(b) Start the executable, press <return>, press <return> again and enter the path to your ptv2vtk.config file. (working folder has to be adjusted in ptv2vtk.config as well.

The tool offers different export options:

 'a': Enables export of 3d-points.
 'b': Enables export of tracked points with Ux,Uy,Uz. (like instantaneous vectors with changing positions)
 'c': Enables export of trajectories. (exports all trajectories)
 'd': Enables export of stacks. (exports trajectories with a specified number of time steps into one scene 
		-->  of course you can load different scenes at once into paraview)

 
The conversion is not optimized yet ( some useless calculation are performed). 
It takes some time to finish as each trajectory is tracked individually.