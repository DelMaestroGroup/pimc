Description of Generic Density Functionality
	Outline
		- We want a generic way to map a d-vec particle position to a bin in a
		  1d array that is independent of dimension or container type
	Define:
		- Discretization in each dimension: dx,dy,dz
		- Should we use a fixed number of bins?

	Prism:
		1d: 
		index = (r[0] + 0.5*side[0])/dr[0]
		size = L/dr

		2d: 
		i[0] = (r[0] + 0.5*side[0])/dr[0]
		i[1] = (r[1] + 0.5*side[1])/dr[1]
		size = (Lx/dx)*(Ly/dy)
		index = i[0]*N[1] + i[1]

		3d:
		i[0] = (r[0] + 0.5*side[0])/dr[0]
		i[1] = (r[1] + 0.5*side[1])/dr[1]
		i[2] = (r[2] + 0.5*side[2])/dr[2]
		size = N[0]*N[1]*N[2]

		index = i[0]*N[1]*N[2] + i[1]*N[2] + i[2]

	Cylinder:
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


