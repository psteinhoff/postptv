This "package" is a small tool to convert your ptv.is-files into .vtk-files (unstructured grid)
The .vtk-files can be loaded into paraview, a powerfull post-processing tool.

add velo-trans.par, 3d_point_offset.par to your parameter folder.


The tool offers different export options:

 'a': Enables export of 3d-points. (rt_is-files )
 'b': Enables export of tracked points with Ux,Uy,Uz. (like instantaneous vectors with changing positions)
 'c': Enables export of trajectories. (exports all trajectories)
 'd': Enables export of stacks. (exports trajectories with a specified number of time steps into one scene 
		-->  of course you can load different scenes at once into paraview)
 'e': Enables export of trajectories. (sorted and filtered by length - time consuming!!!)

