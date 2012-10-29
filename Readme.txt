PFI_GUI:	the start m-file.
PFI:		the GUI/User Input
delta_t:	Stability Condition Finder
PFI_CUDA:	Computations on GPU
PFI_CPU:	Computations on CPU
OMEGACALC:	GPU arrayfun to calculate Omega
PSICALC:	GPU arrayfun to calculate Psi  
VELCALC:	GPU arrayfun to calculate velocity components


Procedure:

In general the physical limits (Xmin/max Ymin/max) stay the same, but you'll note that changing
Nx or My will cause the value of dt to change.  dt then provides an upper limit 
to the value of dt to meet stability criteria.  We usually then truncate this value
to 1 or 2 sig digits to ensure that A) We are under the criteria and B) The value
is easy to remember.  You will also note that changing Re also changes the value of
dt.  The value in Re is a Reynolds number based on the radius of curvature of the leading
edge, which in our case is ~0.01578 of the chord length.  Thus to calculate the Reynolds number
based on chord length, we use  RE_chord = Re_le/.01578.

We need to capture about 10 grid lines in the y direction within the boundary layer.  As the
Reynolds number increases, the boundary layer thickness decreases.  Thus as the reynolds number
increases so to must our mesh size.  We have found that for Re = 700, we need about 400 grid spaces
in the y (My) direction, and that Re = 2100 requires about 800.

A Boeing 747 cruising at 35,000 feet and Mach 0.84 has a root chord Reynolds number of ~1e8 and 
has a tip chord Reynolds number of ~2.7e7.  A leading edge Reynolds number of 700 converts to ~4.5e4.
Hopefully this will help provide an idea of how large our domain, and thus a sparse matrix, needs to be for us
to get results.


Filewrite:  Enter negative values to only write out at the end of simulation.  
Report: Provides screen output showing:
Iteration number, iterations to convergence, residual in streamfunc, residual in vorticity, time spent on that iteration

Use inputs from the GUI to provide an idea for inputs to the Fortran code.