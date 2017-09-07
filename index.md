# Welcome to Parabola Flow Interactive!
Parabola Flow Interactive (PFI) was initially developed by Wallace J. Morris II to investigate boundary layer flow phenomena, which he utilized in developing his theory for a universal prediction of airfoil stall onset.
![Stall Angle vs. Reynolds Number](http://i.imgur.com/PS6XSDi.png)

## Purpose 
This repository serves a slightly different purpose than the original which was to substantiate Dr. Morris' developing theory.  While this code does intend to continue this purpose by exploring higher Reynolds numbers, the primary purposes of this repository are to:
* Serve as an algorithm testbed for later conversion into high-performance languages
* Serve as a non-trivial, but simple, example of a CFD code -- even if it is limited to 2D incompressible flow
* Serve as a benchmark of sorts for testing Matlab performance and language development (string arrays!)
* Serve as an inspiration for those interested in using Matlab for their own simulation codes, even if they're not CFD

## Vision
* Improved spatial accuracy
   * Larger meshes
   * Higher-order finite differences
   * Compact finite difference schemes
* Improved temporal accuracy
   * Higher-order explicit integration
* Improved solution speeds 
   * Implicit integration
   * Implicit<->Explicit switching
   * Improved Preconditioners
   * Utilize parallel-compute
      * GPU Compute
      * MIC
      * Parallel CPU
      * CPU Vectorization
      * Parallel-in-time solvers
   * Input-file execution

   
## Code Style
* The code should emphasize clarity over performance.
   * Minimal nested functions
   * Composite functions must be clear _(should be "human-readable")_
       * `max(max(A))` is fine
       * `max(max(cell2mat(min(bsxfun(@hypot,sin(A),A)))))` is not
* The code should utilize latest language features as long as other goals are met
   * String Arrays
   * Implicit Expansion

## Background Publications
Morris, W. (2009). A universal prediction of stall onset for airfoils at a wide range of Reynolds number flows (Ph.D.). Rensselaer Polytechnic Institute. [http://digitool.rpi.edu:8881/R/FUG5ICDYFYEPMJF1MHBQIRVL5YMP1HUDBTII2T7XTP85UFY5LG-00157](http://digitool.rpi.edu:8881/R/FUG5ICDYFYEPMJF1MHBQIRVL5YMP1HUDBTII2T7XTP85UFY5LG-00157)

Morris, W. J., & Rusak, Z. (2013). Stall onset on aerofoils at low to moderately high Reynolds number flows. Journal of Fluid Mechanics, 733, 439–472. [doi:10.1017/jfm.2013.440](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/stall-onset-on-aerofoils-at-low-to-moderately-high-reynolds-number-flows/648F9A27BAEEBE84CF381225519749BC)

Rusak Z, Morris WJ, II. Stall Onset on Airfoils at Moderately High Reynolds Number Flows. ASME. J. Fluids Eng. 2011;133(11):111104-111104-12. [doi:10.1115/1.4005101.](http://fluidsengineering.asmedigitalcollection.asme.org/article.aspx?articleid=1439413)

Rusak Z, Morris WJ, II, Peles Y. Prediction of Leading-Edge Sheet Cavitation Inception on Hydrofoils at Low to Moderate Reynolds Number Flows. ASME. J. Fluids Eng. 2007;129(12):1540-1546. [doi:10.1115/1.2801350.](http://fluidsengineering.asmedigitalcollection.asme.org/article.aspx?articleID=1432841)


## Authors and Contributors
The original code can be found in Dr. Morris' Ph.D. thesis.  The original code was written in FORTRAN77.
Several versions of the code have since been developed, in various languages, but vary in their current development status.
This version of the code was translated from Dr. Morris' original code into Matlab by @GregVernon, who is also the main developer of this implementation.

![PFI](http://i.imgur.com/lR6CVo5.png)
