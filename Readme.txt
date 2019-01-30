Source codes developed/extended originally from Mixutre/Src

(Working Directory and Test Data to add following functions)

New Features:
 v Nonuniform hydrated water distribution
 v Anisotropic Van der Waals force field
 v Monte Carlo moisture transport simulation
 

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

 Simulation results (x***.dat) can shown in the following ways. 

(a) Use Pycodes/pplot.py 
    --> clay sheet distribution and deformations are drawn as a series of snapshots.  

(b) Create pixel image and then drawn by Python
    
  Run following programs succesively in DEM2018/TestData. 
  ../paint 
  ../Pycodes/imgs.py 
  
  To run "paint", input data file named "paint.inp" is required  
