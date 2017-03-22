#!/usr/bin/env python

"""
Script to calculate the area of a triangle defined by three atoms in an XDATCAR file.
Requires the structure_tools module to function.

Proper usage:
triangles.py [xdatcar] [atom ids for vertices of triangle]
"""


import sys
import structure_tools as st
import numpy as np


if len(sys.argv) != 5:
	print "Not sure about your arguments."
	print "Try:"
	print "triangles.py [xdatcar] [atom ids for vertices of triangle]"



xdatcar = sys.argv[1]

a1 = sys.argv[2]
a2 = sys.argv[3]
a3 = sys.argv[4]
multiple_triangles = False

print "reading "+xdatcar
imagelist = st.read_xdatcar(xdatcar)
print "now looking for atom ids {0}, {1}, {2}".format(a1,a2,a3)
arealist1 = []

for i, image in enumerate(imagelist[:]):
	# Make sure things are in cartesian coordinates
	image.d2c()

	# Calculate area of triangle formed by a1, a2, and a3.
	coordA = [ a.pos for a in image.atomlist if a.id == a1 ][0]
	coordB = [ a.pos for a in image.atomlist if a.id == a2 ][0]
	coordC = [ a.pos for a in image.atomlist if a.id == a3 ][0]
	print a1+"  "+str(coordA)
	print a2+"  "+str(coordB)
	print a3+"  "+str(coordC)

	AB = st.vecfromA2B( image.basis, coordA, coordB )
	AC = st.vecfromA2B( image.basis, coordA, coordC )
	print "length of AB: "+str( np.linalg.norm(AB) )
	print "length of AC: "+str( np.linalg.norm(AC) )

	# calculate area of the triangle
	area = 0.5 * np.linalg.norm( np.cross(AB,AC) )
	arealist1.append(area)
	print "triangle area = "+str(area)
	print ""


outfile = "area.dat"
out = open(outfile,'w')
for i, frame in enumerate(imagelist):
	out.write( "{0:6}  {1}\n".format(frame.framenum, arealist1[i]) )
print "Output written to "+outfile
	
