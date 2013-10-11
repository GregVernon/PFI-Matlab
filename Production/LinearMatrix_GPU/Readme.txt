This directory contains the MATLAB Production version of PFI_Linearized (I need to find a better name)

It uses the GPU for calculations of OMEGA, VELOCITY, and the b (or RHS) vector in the problem Ax=b

It solves for Psi with a sparse solver - TFQMR - on the CPU.

Additionally, I'd like to explore vtk as an output format.  I have generated a test file and deposited it here... try opening it with Paraview 4.0

Finally, this is still a work in progress... BoundaryCalc is still done on CPU (negligble time impact?)

PFI_FileWrite not yet to specification... want to include multiple different optional outputs:
Matlab Binary file
New ASCII Spec  (Clears up clutter, no?)
Old ASCII Spec  (What we've been doing)
VTK				(Paraview Compatible -- then we don't have to write/maintain our own post-processing code)


Anyways, give it a whirl!