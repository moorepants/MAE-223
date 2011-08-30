Description
===========

These are all of the source code files that I wrote while taking Mont Hubbard's
multibody dynamics class in the Winter of 2006 at UC Davis. It includes
examples I used to learn [Autolev](http://www.autolev.com/) and my bicycle
dynamics class project. This bicycle dynamics project was my first attempt to
derive the equations of motion for the Whipple bicycle model.  These files may
be useful to someone but beware that there is no guarantee that these files
don't have bugs or are giving the correct answers.  In fact I don't think that
any of the bicycle files correctly derive the exact linear equations of motion
of the whipple bicycle as presented in Meijaard, et al. 2007, although they are
very close. My project paper on the bicycle can be downloaded [from my
website][paper].


[paper]: http://mae.ucdavis.edu/~biosport/jkm/bh_mypapers/Low%20Speed%20Bicycle%20Stability%20-%20Effects%20of%20Geometric%20Parameters%202006.pdf.

File Descriptions
=================


BicycleDynamicsProject
----------------------
These files were used to generate the (mostly correct) results in my [project
paper][paper].

BICYCLE[1-14].AL were all of my first attempts at deriving the equations of
motion of the non-linear and linear Whipple bicycle using Autolev. BICYCLE14.AL
and FINDQ8.AL were the final versions used to generate the data for the project
paper.

`bike_inertia.m`, `arm.m`, `leg.m`, `evalprimes.m` are Matlab source files that
use the equations of motion generated with BICYCLE14.AL and FINDQ8.AL to
generate the results in the paper in the project paper. It generates estimates
of the parameters for a bicycle and rider just based off of geometric
measurments and the rider's mass and then uses the equations of motion to plot
the eigenvalues versus speed for the bicycle.

thomasbike.al was written by another student, Thomas Engelhardt, and may
generate the equations of motion correctly. He was able to complete the
derivation of the Whipple model before I did and gave me this file to help me
debug my program.
