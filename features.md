PIMC Feature List
=================

### Local Permutation Cycle

Update Constant Singleton to use Associative Containers
-------------------------------------------------------
Goal: Clean up and re-factor the constants() singleton to use a map taking the
parameter name.  This is more consistent with the latest version of the code
and should simplify extensions.

### Concerns
 * map lookups are O(log(n)) should we be concerned by this?
 * would extra freedom of insertion into map lead to bugs?

DONE: Local Superfluid Density Estimators
-----------------------------------------

### Local Winding Number Estimator
### Local Area Estimator

 * Details provided in E. W. Draeger and D. M. Ceperley, Phys. Rev. Lett. 90, 065301 (2003).

MOSTLY DONE: Re-Factoring of Communicator
-----------------------------------------

Goal: Creation of a new class which contains all functionality needed for a
input/ouptut file.  Communicator could then just contain a map of filenames
to file objects

### Requirements
 * Data
  * File pointer
  * File label (string)
  * State (open/close)
 * Methods
  * open/close
  * reset

### Reset

Reset is used for those files that are not continually kept open but that we
wish to write over.  To be safe, we should write to a new temporary file, then
rename after successful completion.

### Possible Issues
 * fixed file functionality has not been fully tested

DONE: Generic Density Tensor Functionality
------------------------------------------

### Outline

We want a generic way to map a d-vec particle position to a bin in a 1d array
that is independent of dimension or container type

Define:
 * Discretization in each dimension: dx,dy,dz
 * Should we use a fixed number of bins?

### Prism
 * 1d: 

    index = (r[0] + 0.5*side[0])/dr[0]
    size = L/dr

 * 2d: 

    i[0] = (r[0] + 0.5*side[0])/dr[0]
    i[1] = (r[1] + 0.5*side[1])/dr[1]
    size = (Lx/dx)*(Ly/dy)
    index = i[0]*N[1] + i[1]

 * 3d:

    i[0] = (r[0] + 0.5*side[0])/dr[0]
    i[1] = (r[1] + 0.5*side[1])/dr[1]
    i[2] = (r[2] + 0.5*side[2])/dr[2]
    size = N[0]*N[1]*N[2]

    index = i[0]*N[1]*N[2] + i[1]*N[2] + i[2]

    dV = dx*dy*dz

### Cylinder:

3d: We will use cylindrical polar coordinates (r,theta,z)

	N[0] = R/dr
	N[1] = 2.0*pi/dtheta
	N[2] = side[2]/dz

	i[0] = (r[0]*r[0] + r[1]*r[1])/dr
	theta = atan2(r[1],r[0])
	if theta < 0: theta += 2.0*pi
	i[1] = theta/dtheta
	i[2] = (r[2] + 0.5*side[2])/dr[2]

	index = i[0]*N[1]*N[2] + i[1]*N[2] + i[2]

    dV(r) = r*dr*dtheta*dz


