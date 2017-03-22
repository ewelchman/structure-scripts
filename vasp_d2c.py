#!/usr/bin/env python

"""
Translates POSCAR from direct to cartesian coordinates in Angstrom.
Writes a file "cartesianPOSCAR"
Usage:
vasp_d2c [poscar]
"""


import sys
import structure_tools as st


if ( len(sys.argv) == 2 ):
	readfile = sys.argv[1]
else:
	print '\n incorrect number of arguments \n'
	print 'You want to use: '
	print 'vasp_d2c [POSCAR]'
	sys.exit()

struct = st.Structure()
st.read_poscar(struct, readfile)
#print struct
#struct.d2c()
st.write_poscar(struct, "cartesianPOSCAR")
