#!/usr/bin/env python

"""
Script to extract a vasp POSCAR file from a QE input file

Correct usage:
QEtoVASP.py [file to read] [file to write]
"""

import sys
import structure_tools as st

###############  MAIN  #######################
if len(sys.argv) != 3:
	print "Wrong number of arguments"
	print "Correct usage:"
	print "QEtoVASP.py [file to read] [file to write]"
	sys.exit()

qefile = sys.argv[1]
print "looking in {0}".format(qefile)
vaspfile = sys.argv[2]

structure = st.Structure(qefile)
st.read_qe( structure, qefile )
st.write_poscar( structure, vaspfile )
