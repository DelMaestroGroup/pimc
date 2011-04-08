#!/usr/bin/python
# hescrew.py
# Adrian Del Maestro
# 05.03.2010
# 
# Create a list of atomic positions representing a screw dislocation in solid helium-4

import os,sys,glob
from optparse import OptionParser
import math

# -----------------------------------------------------------------------------
def getXYLength(vec):
	'''Return the length of a 3-vector projected into the XY plane.'''
	r = 0.0
	for i in range(2):
		r += vec[i]*vec[i]
	return math.sqrt(r)

# -----------------------------------------------------------------------------
# Begin Main Program 
# -----------------------------------------------------------------------------
def main(): 

	# setup the command line parser options 
	parser = OptionParser() 
	parser.add_option("-a", "--lattice-constant", dest="a", type="float", \
			help="lattice constant [angstroms]") 
	parser.add_option("-n", "--number-basal-planes", dest="n", type="int", \
			help="number of basal planes in the c-axis (z) direction") 

	# parse the command line options and get the file name
	(options, args) = parser.parse_args() 
	if (not options.n) or (not options.a):
		parser.error("Need to specify both the lattice constant (a) and number of basal planes (n)!")

	# Create local copies of the options
	a = options.a
	c = math.sqrt(8.0/3.0)*a

	# Determine the system size in the z-direciton in terms of the number of basal planes
	L = (options.n)*c - 1.0E-6

	# The core radius is defined to contain 1/2 the volume of the simulation cell
	# i.e. pi R^2 L = L^3/2
	coreRadius = L/(math.sqrt(2.0*math.pi))
	
	# We determine how many lattice vectors we will need in each direction
	# We pad it a little bit, then reduce it later
	N = [0,0,0]
	N[0] = int(math.ceil(L/a)) + 4
	N[1] = N[0]
	N[2] = int(math.ceil(L/c)) + 4

	atomList = []
	# an hcp crystal is composed of two interpenetrating perfect triangular lattices
	# separated by c/2 in the z-direction and shifted from each other by (a/2,a/sqrt(12))
	# in the xy-plane.
	for i1 in range(-N[0],N[0]):
		for i2 in range(-N[1],N[1]):
			for i3 in range(-N[2],N[2]):

				i2even = not (i2 % 2)
				i3even = not (i3 % 2)

				if (i2even and i3even):
					shiftX = 0.0
					shiftY = 0.0
				elif (not i2even and i3even):
					shiftX = a/2.0
					shiftY = 0.0;
				elif (i2even and not i3even):
					shiftX = a/2.0
					shiftY = a/math.sqrt(12.0)
				else:
					shiftX = 0.0
					shiftY = a/math.sqrt(12.0)

				p = [0.0,0.0,0.0]

				p[0] = a*i1 + shiftX
				p[1] = a*math.sqrt(3.0)*i2/2.0 + shiftY
				p[2] = c*i3/2.0

				# Because we want the screw dislocation to occur at the origin, 
				# we shift all atoms in the y-direction
				p[1] -= 0.5*a*(0.5*math.sqrt(3.0)+1.0/math.sqrt(12.0))

				# Here we perform the shift that inserts the screw dislocation
				# See J.P. Hirth and J. Lothe, "Theory of Dislocations", page 60
				theta = math.atan2(p[1],p[0])
				if (theta < 0.0):
					theta += 2.0*math.pi
				p[2] += (c/(2.0*math.pi))*theta

				# We only keep positions that are inside our simulation box 
				# which has PBC in all directions
				inBox = 0
				for i in range(3):
					if (p[i] >= -0.5*L and p[i]<0.5*L):
						inBox += 1

				if inBox == 3:
					atomList.append(p)

	# Find out how many fixed and updateable particles there are
	numFixed = 0
	numUpdate = 0
	for atom in atomList:
		if getXYLength(atom) < coreRadius:
			numUpdate += 1
		else:
			numFixed += 1

	# Open the output file and write the header
	outFile = open('screw-%07.3f-%4.3f.dat' % (L,a),'w');
	outFile.write('#%15s\t%16s\t%16s\t%16s\n'%('Update or Fix','r_x','r_y','r_z'))

	# Go through the atom list.  If we are inside the core, than we have an __U__pdateable
	# particle, otherwise we have a __F__ixed particle.
	for atom in atomList:
		if getXYLength(atom) < coreRadius:
			atomState = 'U'
		else:
			atomState = 'F'

		outFile.write('%-16s\t%+16.8E\t%+16.8E\t%+16.8E\n' % (atomState,atom[0],atom[1],atom[2]))

	print 'numUpdate = %d' % numUpdate
	print 'numFixed  = %d' % numFixed

#	# shift the y positions
#	shiftAtomList = []
#	for atom in atomList:
#		shift = 0.0*a*(0.5*math.sqrt(3.0)+1.0/math.sqrt(12.0))
#		p = [atom[0],atom[1]-shift,atom[2]]
#		shiftAtomList.append(p)
#		if math.sqrt(p[0]*p[0]+p[1]*p[1]) < coreRadius:
#			print '%+12.6E\t%+12.6E\t%+12.6E\t%+12.6E' % (p[0],p[1],p[2],p[2]+shift)

	# Define the basis vector
#	a1 = [a,0.0,0.0]
#	a2 = [0.5*a,0.5*math.sqrt(3)*a,0.0]
#	a3 = [0.0,0.0,c]
#
#	shift = [0,0,0]
#	for i in range(3):
#		shift[i] = a1[i]/3.0 + a2[i]/3.0 + a3[i]/2.0
#
#	atomList = []
#	for n3 in range(0,N[2]/2):
#		for n1 in range(0,N[0]):
#			for n2 in range(0,N[1]):
#				p1 = [0,0,0]
#				p2 = [0,0,0]
#				for i in range(3):
#					p1[i] = n1*a1[i] + n2*a2[i] + n3*a3[i]
#					p2[i] = p1[i] + shift[i]
#				inBox = 0
#				for i in range(3):
#					if (p[i] >= -0.5*L) and (p[i] < 0.5*L):
##						inBox += 1;
#				if inBox==3:
#				atomList.append(p1)
#				atomList.append(p2)

#	for atom in atomList:
		#test = -2.690885
#		if abs(atom[2]-0.5*c) < 0.1:
#		shift = 0.5*a*(0.5*math.sqrt(3.0)+1.0/math.sqrt(12.0))
#		print 'He\t%12.6E\t%12.6E\t%12.6E' % (atom[0],atom[1]-shift,atom[2])


# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__": 
	main()
