## NPT-SLLOD 
This code implements SLLOD in isobaric-isothermal ensemble for a system of soft-spheres. 

To compile the code use 

`g++ main_npt.cpp soft_sphere_nvt.cpp -o npt.o # NEEDS GCC >= 7 ` 

To run the code 

`./npt.o [input_file] [BOX] [PRESSURE]`

The format for the INPUT_FILE is as follows 

`particle_ID POS_X POS_Y RADIUS`

The simulation box is defined from 

``` 
	-[BOX]/2 [BOX]/2
	-[BOX]/2 [BOX]/2
```

Temperature and other related quantities are specified in the begining of 
soft_sphere_nvt.cpp. 

![image](/PRESSURE_TEMP_NPT_SLLOD.png)
