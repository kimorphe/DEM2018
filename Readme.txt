Source codes developed/extended originally from Mixutre/Src

(Working Directory and Test Data to add following functions)

New Features:
 v Nonuniform hydrated water distribution
 v Anisotropic Van der Waals force field
 v Monte Carlo moisture transport simulation
 v Simulation may be restarted
 

How to run: 

0. Minimal input dataset:
 gen_sheet.inp
 wall.dat (prepare dummy when wall particles is not used)
 dem.inp

1. gen_sheet <-- gen_sheet.inp
 --> Out_Dir/bbox.dat
 --> Out_Dir/gen_sheet.out
 --> Out_Dir/ptc.dat
 --> Out_Dir/sheet.dat

2. dem <-- dem.inp, ptc.dat, sheet.dat, bbox.dat, wall.dat (restart data file, cell_const.dat)
 --> ptc_nums.dat
 --> energy.out
 --> stress.out
 --> ptcl.out
 (--> cell_const.dat)
 --> x***.dat 

3. Test Data
 Move to DEM2018/TestData/Input

 Run following programs successively.

 ../../gen_sheet
 ../../dem

 Note1: 
   To able/disable water transport; 
   Set "mvw" parameter to 1/0 in "dem.inp"
 Note2: 
   To restart your simulation: 
    Set "rstat" parameter to 1 in "dem.inp".
    Also, save "x***.dat" as "restart0.dat" in your "TestData/Input" folder. 
    The DEM code imports the initial condition from "restart0.dat". 
    You would probably need to edit "dem.inp" to impose a new set of simulation 
	parameters, say boundary condition. 

 Simulation results (x***.dat) can shown in the following ways. 

(a) Use Pycodes/pplot.py 
    --> clay sheet distribution and deformations are drawn as a series of snapshots.  

(b) Create pixel image and then drawn by Python
    
  Run following programs succesively in DEM2018/TestData. 
  ../paint 
  ../Pycodes/imgs.py 
  
  To run "paint", input data file named "paint.inp" is required  
