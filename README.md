To run,

   1. Open RUNME.m
   2. change experiment = 1, 2, 3, 4, 5, or 6
   3. run it
   4. run the corresponding FIGURE script (see thesis)

To fit,

   1. run `solve_exp1_and_exp2.m` or `solve_exp3_and_exp4.m`
   2. wait an eternity
   3. substitute the resulting `bestpar` in RUNME.m for the corresponding experiment, right where it says `startpar([...]) = [...]`
   4. run it

To use pre-generated data,

   1. open some FIGURE script 
   2. uncomment the first line and make it `load` the corresponding .mat file
   3. run the figure
