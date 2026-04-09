# DEM-ECP
here I share the code wrote in Python, which is able to generate the topology of the electrically conductive polymer, and calculate the corresponding electrical conductivity.

You can generate the linear particle using CNT_generator.py for LAMMPS input.

Run LAMMPS script.

Read the output of the LAMMPS with reader.py. here you can also call the segs_dis.py to get the distance matrix between nearest segments.

analyse the distance matrix and calculate the equivalent resistance with laplace_solver.py
