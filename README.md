## NPT-SLLOD 
This code implements NPT-SLLOD using Nose-Hoover thermostat for a system of soft-spheres. 
To compile the code use 
	g++ main_npt.cpp soft_sphere_nvt.cpp -o npt.o 
To run the code 
	./a.out {input_file} {BOXh} {PRESSURE} 
where BOXh is the half box length. 
