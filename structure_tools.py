#!/usr/bin/env python

"""
This is a module of methods to use for manipulating and analyzing
structures. Some of the methods are similar to those used in aselite,
but done in a way I understand.
"""

import os
import math
import numpy as np
import copy

chemical_symbols = ['X',  'H',  'He', 'Li', 'Be',
                    'B',  'C',  'N',  'O',  'F',
                    'Ne', 'Na', 'Mg', 'Al', 'Si',
                    'P',  'S',  'Cl', 'Ar', 'K',
                    'Ca', 'Sc', 'Ti', 'V',  'Cr',
                    'Mn', 'Fe', 'Co', 'Ni', 'Cu',
                    'Zn', 'Ga', 'Ge', 'As', 'Se',
                    'Br', 'Kr', 'Rb', 'Sr', 'Y',
                    'Zr', 'Nb', 'Mo', 'Tc', 'Ru',
                    'Rh', 'Pd', 'Ag', 'Cd', 'In',
                    'Sn', 'Sb', 'Te', 'I',  'Xe',
                    'Cs', 'Ba', 'La', 'Ce', 'Pr',
                    'Nd', 'Pm', 'Sm', 'Eu', 'Gd',
                    'Tb', 'Dy', 'Ho', 'Er', 'Tm',
                    'Yb', 'Lu', 'Hf', 'Ta', 'W',
                    'Re', 'Os', 'Ir', 'Pt', 'Au',
                    'Hg', 'Tl', 'Pb', 'Bi', 'Po',
                    'At', 'Rn', 'Fr', 'Ra', 'Ac',
                    'Th', 'Pa', 'U',  'Np', 'Pu',
                    'Am', 'Cm', 'Bk', 'Cf', 'Es',
                    'Fm', 'Md', 'No', 'Lr']

atomic_numbers = {}
for Z, symbol in enumerate(chemical_symbols):
	atomic_numbers[symbol] = Z

atomic_names = [
    '', 'Hydrogen', 'Helium', 'Lithium', 'Beryllium', 'Boron',
    'Carbon', 'Nitrogen', 'Oxygen', 'Fluorine', 'Neon', 'Sodium',
    'Magnesium', 'Aluminium', 'Silicon', 'Phosphorus', 'Sulfur',
    'Chlorine', 'Argon', 'Potassium', 'Calcium', 'Scandium',
    'Titanium', 'Vanadium', 'Chromium', 'Manganese', 'Iron',
    'Cobalt', 'Nickel', 'Copper', 'Zinc', 'Gallium', 'Germanium',
    'Arsenic', 'Selenium', 'Bromine', 'Krypton', 'Rubidium',
    'Strontium', 'Yttrium', 'Zirconium', 'Niobium', 'Molybdenum',
    'Technetium', 'Ruthenium', 'Rhodium', 'Palladium', 'Silver',
    'Cadmium', 'Indium', 'Tin', 'Antimony', 'Tellurium',
    'Iodine', 'Xenon', 'Caesium', 'Barium', 'Lanthanum',
    'Cerium', 'Praseodymium', 'Neodymium', 'Promethium',
    'Samarium', 'Europium', 'Gadolinium', 'Terbium',
    'Dysprosium', 'Holmium', 'Erbium', 'Thulium', 'Ytterbium',
    'Lutetium', 'Hafnium', 'Tantalum', 'Tungsten', 'Rhenium',
    'Osmium', 'Iridium', 'Platinum', 'Gold', 'Mercury',
    'Thallium', 'Lead', 'Bismuth', 'Polonium', 'Astatine',
    'Radon', 'Francium', 'Radium', 'Actinium', 'Thorium',
    'Protactinium', 'Uranium', 'Neptunium', 'Plutonium',
    'Americium', 'Curium', 'Berkelium', 'Californium',
    'Einsteinium', 'Fermium', 'Mendelevium', 'Nobelium',
    'Lawrencium', 'Unnilquadium', 'Unnilpentium', 'Unnilhexium']

atomic_masses = np.array([
   0.00000, # X
   1.00794, # H
   4.00260, # He
   6.94100, # Li
   9.01218, # Be
  10.81100, # B
  12.01100, # C
  14.00670, # N
  15.99940, # O
  18.99840, # F
  20.17970, # Ne
  22.98977, # Na
  24.30500, # Mg
  26.98154, # Al
  28.08550, # Si
  30.97376, # P
  32.06600, # S
  35.45270, # Cl
  39.94800, # Ar
  39.09830, # K
  40.07800, # Ca
  44.95590, # Sc
  47.88000, # Ti
  50.94150, # V
  51.99600, # Cr
  54.93800, # Mn
  55.84700, # Fe
  58.93320, # Co
  58.69340, # Ni
  63.54600, # Cu
  65.39000, # Zn
  69.72300, # Ga
  72.61000, # Ge
  74.92160, # As
  78.96000, # Se
  79.90400, # Br
  83.80000, # Kr
  85.46780, # Rb
  87.62000, # Sr
  88.90590, # Y
  91.22400, # Zr
  92.90640, # Nb
  95.94000, # Mo
    np.nan, # Tc
 101.07000, # Ru
 102.90550, # Rh
 106.42000, # Pd
 107.86800, # Ag
 112.41000, # Cd
 114.82000, # In
 118.71000, # Sn
 121.75700, # Sb
 127.60000, # Te
 126.90450, # I
 131.29000, # Xe
 132.90540, # Cs
 137.33000, # Ba
 138.90550, # La
 140.12000, # Ce
 140.90770, # Pr
 144.24000, # Nd
    np.nan, # Pm
 150.36000, # Sm
 151.96500, # Eu
 157.25000, # Gd
 158.92530, # Tb
 162.50000, # Dy
 164.93030, # Ho
 167.26000, # Er
 168.93420, # Tm
 173.04000, # Yb
 174.96700, # Lu
 178.49000, # Hf
 180.94790, # Ta
 183.85000, # W
 186.20700, # Re
 190.20000, # Os
 192.22000, # Ir
 195.08000, # Pt
 196.96650, # Au
 200.59000, # Hg
 204.38300, # Tl
 207.20000, # Pb
 208.98040, # Bi
    np.nan, # Po
    np.nan, # At
    np.nan, # Rn
    np.nan, # Fr
 226.02540, # Ra
    np.nan, # Ac
 232.03810, # Th
 231.03590, # Pa
 238.02900, # U
 237.04820, # Np
    np.nan, # Pu
    np.nan, # Am
    np.nan, # Cm
    np.nan, # Bk
    np.nan, # Cf
    np.nan, # Es
    np.nan, # Fm
    np.nan, # Md
    np.nan, # No
    np.nan])# Lw

# Covalent radii from:
#
#  Covalent radii revisited,
#  Beatriz Cordero, Veronica Gomez, Ana E. Platero-Prats, Marc Reves,
#  Jorge Echeverria, Eduard Cremades, Flavia Barragan and Santiago Alvarez,
#  Dalton Trans., 2008, 2832-2838 DOI:10.1039/B801115J 
missing = 0.2
covalent_radii = np.array([
    missing,  # X
    0.31,  # H
    0.28,  # He
    1.28,  # Li
    0.96,  # Be
    0.84,  # B
    0.76,  # C
    0.71,  # N
    0.66,  # O
    0.57,  # F
    0.58,  # Ne
    1.66,  # Na
    1.41,  # Mg
    1.21,  # Al
    1.11,  # Si
    1.07,  # P
    1.05,  # S
    1.02,  # Cl
    1.06,  # Ar
    2.03,  # K
    1.76,  # Ca
    1.70,  # Sc
    1.60,  # Ti
    1.53,  # V
    1.39,  # Cr
    1.39,  # Mn
    1.32,  # Fe
    1.26,  # Co
    1.24,  # Ni
    1.32,  # Cu
    1.22,  # Zn
    1.22,  # Ga
    1.20,  # Ge
    1.19,  # As
    1.20,  # Se
    1.20,  # Br
    1.16,  # Kr
    2.20,  # Rb
    1.95,  # Sr
    1.90,  # Y
    1.75,  # Zr
    1.64,  # Nb
    1.54,  # Mo
    1.47,  # Tc
    1.46,  # Ru
    1.42,  # Rh
    1.39,  # Pd
    1.45,  # Ag
    1.44,  # Cd
    1.42,  # In
    1.39,  # Sn
    1.39,  # Sb
    1.38,  # Te
    1.39,  # I
    1.40,  # Xe
    2.44,  # Cs
    2.15,  # Ba
    2.07,  # La
    2.04,  # Ce
    2.03,  # Pr
    2.01,  # Nd
    1.99,  # Pm
    1.98,  # Sm
    1.98,  # Eu
    1.96,  # Gd
    1.94,  # Tb
    1.92,  # Dy
    1.92,  # Ho
    1.89,  # Er
    1.90,  # Tm
    1.87,  # Yb
    1.87,  # Lu
    1.75,  # Hf
    1.70,  # Ta
    1.62,  # W
    1.51,  # Re
    1.44,  # Os
    1.41,  # Ir
    1.36,  # Pt
    1.36,  # Au
    1.32,  # Hg
    1.45,  # Tl
    1.46,  # Pb
    1.48,  # Bi
    1.40,  # Po
    1.50,  # At
    1.50,  # Rn
    2.60,  # Fr
    2.21,  # Ra
    2.15,  # Ac
    2.06,  # Th
    2.00,  # Pa
    1.96,  # U
    1.90,  # Np
    1.87,  # Pu
    1.80,  # Am
    1.69,  # Cm
    missing,  # Bk
    missing,  # Cf
    missing,  # Es
    missing,  # Fm
    missing,  # Md
    missing,  # No
    missing,  # Lr
    ])


electronegativity = np.array([
 np.nan, # X
   2.20, # H
 np.nan, # He
   0.98, # Li
   1.57, # Be
   2.04, # B
   2.55, # C
   3.04, # N
   3.44, # O
   3.98, # F
 np.nan, # Ne
   0.93, # Na
   1.31, # Mg
   1.61, # Al
   1.90, # Si
   2.19, # P
   2.58, # S
   3.16, # Cl
 np.nan, # Ar
   0.82, # K
   1.00, # Ca
   1.36, # Sc
   1.54, # Ti
   1.63, # V
   1.66, # Cr
   1.55, # Mn
   1.83, # Fe
   1.88, # Co
   1.91, # Ni
   1.90, # Cu
   1.65, # Zn
   1.81, # Ga
   2.01, # Ge
   2.18, # As
   2.55, # Se
   2.96, # Br
   3.00, # Kr
   0.82, # Rb
   0.95, # Sr
   1.22, # Y
   1.33, # Zr
   1.60, # Nb
   2.16, # Mo
   1.90, # Tc
   2.20, # Ru
   2.28, # Rh
   2.20, # Pd
   1.93, # Ag
   1.69, # Cd
   1.78, # In
   1.96, # Sn
   2.05, # Sb
   2.10, # Te
   2.66, # I
    2.6, # Xe
   0.79, # Cs
   0.89, # Ba
   1.10, # La
   1.12, # Ce
   1.13, # Pr
   1.14, # Nd
 np.nan, # Pm
   1.17, # Sm
 np.nan, # Eu
   1.20, # Gd
 np.nan, # Tb
   1.22, # Dy
   1.23, # Ho
   1.24, # Er
   1.25, # Tm
 np.nan, # Yb
   1.27, # Lu
   1.3 , # Hf
   1.5 , # Ta
   2.36, # W
   1.9 , # Re
   2.2 , # Os
   2.20, # Ir
   2.28, # Pt
   2.54, # Au
   2.0 , # Hg
   1.62, # Tl
   2.33, # Pb
   2.02, # Bi
   2.0 , # Po
   2.2 , # At
 np.nan, # Rn
   0.7 , # Fr
   0.9 , # Ra
   1.1 , # Ac
   1.3 , # Th
   1.5 , # Pa
   1.38, # U
   1.36, # Np
   1.28, # Pu
   1.3 , # Am
   1.3 , # Cm
   1.3 , # Bk
   1.3 , # Cf
   1.3 , # Es
   1.3 , # Fm
   1.3 , # Md
   1.3 , # No
 np.nan])# Lr


def rot_matrix(axis, theta):
	"""Get rotation matrix for rotating about an axis by theta"""
	axis = axis/np.sqrt( np.dot(axis,axis) )
	cos = math.cos(theta)
	sin = math.sin(theta)
	ux, uy, uz = axis
	return np.array([[cos+ux*ux*(1-cos), ux*uy*(1-cos)-uz*sin, ux*uz*(1-cos)+uy*sin],\
	                 [uy*ux*(1-cos)+uz*sin, cos+uy*uy*(1-cos), uy*uz*(1-cos)-ux*sin],\
	                 [uz*ux*(1-cos)-uy*sin, uz*uy*(1-cos)+ux*sin, cos+uz*uz*(1-cos)]])


def alignment_params(vecA,vecB):
	"""Get information needed to rotate vecA so that it aligns with vecB.
	vecA and vecB must be given as numpy arrays
	
	Returns angle in [unit] and a numpy vector of the rotation axis"""
	# Make sure vecA and vecB are unit vectors
	vecA = vecA / np.linalg.norm( vecA )
	vecB = vecB / np.linalg.norm( vecB )

	# Get angle between vectors and axis of rotation
	theta = math.acos( np.dot( vecA, vecB ) )
	rot_axis = np.cross( vecA, vecB )
	
	# arccos has two possible angles, so check if this is the one we want
	test_vec = np.dot( rot_matrix(rot_axis,theta), vecA )
	dot_prod = np.dot( test_vec, vecB )
	if math.fabs(dot_prod) < 0.9:
		# rotated vector is probably pointing the wrong way
		theta = -theta

	return rot_axis, theta


def rotate_coords(coord_list, axis, theta):
	""" Rotate an atomlist about a particular axis by theta degrees"""
#	print "rotating about "+str(axis)+" by "+str(theta)
	# Get a rotation matrix and coordinates that will be rotated
	rot_mat = rot_matrix( axis, theta )
	newcoords = copy.deepcopy(coord_list)
	for i, coord in enumerate(coord_list):
		newcoords[i] = np.dot( rot_mat, coord )
	return newcoords


def shift_coords(coord_list, vector):
	""" Shift coordinates of each atom in an atomlist by a particular vector"""
	newcoords = copy.deepcopy(coord_list)
	for i, coord in enumerate(coord_list):
		newcoords[i] = coord + vector
	return newcoords


def symbols_from_potcar(potcar_file):
	symbols = []
	if os.path.isfile(potcar_file):
		with open(potcar_file) as f:
			for line in f:
				if len(line.strip()) > 0 and (line.split())[0] == "TITEL":
					symbol = line.split()[3]
					symbols.append( symbol )

	else:
		symbols.append("X")
	
	return symbols


def read_xsf(system,filepath):
	"""Reads a .xsf crystal structure file"""
	data = open(filepath,'r')
	lines = data.readlines()
	data.close()

	# Start by looking for the lattice vectors
	basis = np.zeros( (3,3) )
	for i, line in enumerate(lines):
		if "PRIMVEC" in line.split()[0]:
			for j in range(3):
				for k in range(3):
					basis[j,k] = float( lines[i+1+j].split()[k] )
			system.scalefactor = 1
			system.basis = basis
			break
		
	# Now look for atomic coordinates with a similar procedure
	system.atomlist = []
	for l, line in enumerate(lines):
		if "PRIMCOORD" in line.split()[0]:
			system.numatoms = int(lines[l+1].split()[0])
			for i in range(system.numatoms):
				c = lines[l+2+i].split()
				system.atomlist.append( Atom( c[1], c[2], c[3] ) )
				system.atomlist[i].symbol = c[0]

	# Get system atomtypes and the number of each of them
	system.atomtypes = []
	system.atompop = []
	for a in system.atomlist:
		if a.symbol not in system.atomtypes:
			system.atomtypes.append(a.symbol)

	system.atompop = [ 0 for s in system.atomtypes ]
	for a in system.atomlist:
		for i, sym in enumerate(system.atomtypes):
			if sym == a.symbol:
				system.atompop[i] += 1

	
	# Get a symbol and unique label (id) for each atom
	atomtype = 0
	for i in range(system.numatoms):
		tot = sum( system.atompop[j] for j in range(atomtype+1) )
		if i == tot:
			atomtype += 1
		sum1 = sum( system.atompop[j] for j in range(atomtype))
		system.atomlist[i].symbol = system.atomtypes[atomtype]
		system.atomlist[i].id = "{0}{1}".format(system.atomtypes[atomtype],str(1+i-sum1))

	# At this point, assume we're in cartesian coords
	system.isdirect = False
	system.selective = False


def read_xdatcar(filepath):
	"""Reads a new (with atomic symbols) Vasp-style XDATCAR file.
	
	Returns a list of structures that can be accessed individually."""

	data = open(filepath,'r')
	lines = data.readlines()
	data.close()

	# Make a list for many systems
	imagelist = []

	# Start by reading the header and lattice vectors
	# They should be the same for each image in the file
	title = lines[0]
	scalefactor = float(lines[1].split()[0])
	basis = np.zeros( (3,3) )
	for i in range(3):
		for j in range(3):
			basis[i,j] = float( lines[i+2].split()[j] )
	# Assuming we're dealing with new vasp
	symbols = lines[5].split()
	populations = map(int, lines[6].split())
	numatoms = sum(populations)
	print("numatoms = {0}".format(str(numatoms)))

	# Check for direct or cartesian coords
	isdirect = True
	if lines[7].startswith( ("C","K") ):
		isdirect = False

	# Set up how to read sets of coordinates
	startline = 8
	finishline = 8 + numatoms

	while finishline <= len(lines):

		# Now read successive sets of coordinates
		image = Structure()
		image.title = title
		image.scalefactor = scalefactor
		image.basis = basis
		image.isdirect = isdirect
		image.atomtypes = symbols
		image.atompop = populations
		image.numatoms = numatoms
		image.selective = False
		image.atomlist = []

		# Read which frame in simulation this image corresponds to
		config = lines[startline-1].replace("=","").split()
#		print config
		image.framenum = int(config[2])

		# Actually read coordinates	
		for i, line in enumerate(lines[startline:finishline]):	
#			print "line {0} is {1}".format(i, line)
			c = line.split()
			image.atomlist.append( Atom( c[0], c[1], c[2] ))

		# Give unique id to each atom
		atomtype = 0
		for i in range(numatoms):
			tot = sum( image.atompop[j] for j in range(atomtype+1) )
			if i == tot:
				atomtype += 1
			sum1 = sum( image.atompop[j] for j in range(atomtype))
			image.atomlist[i].symbol = image.atomtypes[atomtype]
			image.atomlist[i].id = "{0}{1}".format(image.atomtypes[atomtype],str(1+i-sum1))
		
		# Add this structure to the list of images
		imagelist.append(image)

		# iterate things for the big loop
		startline += numatoms+1
		finishline += numatoms +1

	return imagelist


def read_poscar(system,filepath):
	"""Reads a Vasp-style structure file (POSCAR/CONTCAR)"""
	data = open(filepath,'r')
	lines = data.readlines()
	data.close()

	# Start by reading the header
	system.title = lines[0]
	system.scalefactor = float(lines[1].split()[0])
	# Now read basis vectors
	basis = np.zeros( (3,3) )
	for i in range(3):
		for j in range(3):
			basis[i,j] = float( lines[i+2].split()[j] )
	system.basis = system.scalefactor*basis
	
	# Check whether this is old or new VASP format
	# Old VASP doesn't list atomic symbols
	isoldVasp = False
	if len(lines[5].split()) != len(lines[6].split()):
		isoldVasp = True
		print "found old VASP file"
	else:
		print "found new VASP file"
	
	
	if not isoldVasp:
		symbols = lines[5].split()
		populations = map(int, lines[6].split())
	else:
		populations = map(int, lines[5].split())

		# Try getting the symbols from a POTCAR file
		# First we need to find the POTCAR
		symbols = symbols_from_potcar("POTCAR")
		if symbols[0] == "X":
			symbols = symbols_from_potcar("../CalcFold1/POTCAR")
		if symbols[0] == "X":
			print "No potcar found. Continusing with all symbols as X"

	system.atomtypes = symbols
	print "symbols = "+str(symbols)
	print "populations = "+str(populations)
	system.atompop = populations
	system.numatoms = sum(populations)

	# Check for selective dynamics
	if (lines[7].split()[0])[0] == 'S':
		system.selective = True
		if isoldVasp:	coordstartline = 8
		else:			coordstartline = 9
	else:
		system.selective = False
		if isoldVasp:	coordstartline = 7
		else:			coordstartline = 8
	
	# Now make a list of atoms
	system.atomlist = []
	for i in range(system.numatoms):
		c = lines[coordstartline+i].split()
		system.atomlist.append( Atom( c[0], c[1], c[2] ))
		if system.selective:
			system.atomlist[i].selective = ( c[3],c[4],c[5] )
	
	# Get a symbol and unique label (id) for each atom
	atomtype = 0
	for i in range(system.numatoms):
		tot = sum( system.atompop[j] for j in range(atomtype+1) )
		if i == tot:
			atomtype += 1
		sum1 = sum( system.atompop[j] for j in range(atomtype))
		system.atomlist[i].symbol = system.atomtypes[atomtype]
		system.atomlist[i].id = "{0}{1}".format(system.atomtypes[atomtype],str(1+i-sum1))
	
	# Check whether coodrinates are Fractional or Cartesian
	system.isdirect = True
	cartline = lines[coordstartline-1]
	if cartline.startswith( ("C","K") ):
		system.isdirect = False
	
	# Convert to Cartesian if needed
	if system.isdirect:
		system.d2c()
	
#	for atom in system.atomlist:
#		print atom.id+"  "+str(atom.pos)

	try:
		return system
	except:
		pass


def make_supercell(size):
	cells = []
	for i in range(-size,size+1):
		for j in range(-size,size+1):
			for k in range(-size,size+1):
				cells.append([i,j,k])
	return cells


def vecfromA2B(basis,coordA,coordB):
	"""Return the shortest vector from point A to point B in a periodic system"""
	image_cells = make_supercell(1)
	vectors = []
	for im in image_cells:
		Bpos = coordB + np.dot(im, basis)
		vectors.append( Bpos - coordA )
	vectors.sort( key = lambda vec: np.linalg.norm(vec) )
	return vectors[0]


def distfromA2B(basis,coordA,coordB):
	"""Return the distance between two points in a periodic system"""
	vec = vecfromA2B( basis, coordA, coordB )
	d = np.linalg.norm( vec )
	return d
			

def write_poscar(system,outfile):
	"""Write a vasp POSCAR file at outfile"""
	out = open(outfile,'w')
	title = " ".join(system.atomtypes)
	out.write(title+"\n")
	out.write("1.000\n")
	lat = system.basis
	for i in range(3):
		out.write('  {0:13.9f} {1:13.9f} {2:13.9f}\n'.format(lat[i,0],lat[i,1],lat[i,2]))
	out.write(' {0}\n'.format( "  ".join(system.atomtypes) ))
	out.write(' {0}\n'.format( "  ".join(map(str,system.atompop)) ))

	# Check for selective dynamics
	if system.selective:
		out.write('Selective Dynamics\n')

	# Check for Fractional/Cartesian coordinates
	if system.isdirect:
		out.write("Direct\n")
	else:
		out.write("Cartesian\n")

	# Write coordinates for each atom
	# Going through each type of atom one at a time
	txt2 = ""
	for sym in system.atomtypes:
		for atom in [a for a in system.atomlist if a.symbol == sym]:
			txt1 = '  {0:13.9f} {1:13.9f} {2:13.9f}'.format(atom.x,atom.y,atom.z)
			txt3 = "  !"+atom.id
			if system.selective:
				txt2 = " ".join( atom.selective )
			out.write( txt1 +"  "+ txt2 + txt3 +"\n" )

	out.close()
	print "Finished writing {0}".format(outfile)


def write_axsf(systemlist,outfile):
	"""Write a fixed-cell animated XSF file.

	Takes as arguments a list of structures and output file name"""

	out = open(outfile,'w')

	# Write header for the file with lattice from first listed structure
	num_steps = len(systemlist)
	out.write("ANIMSTEPS {0}\n".format(num_steps))
	out.write("CRYSTAL\n")
	out.write("PRIMVEC\n")
	lat = systemlist[0].basis
	for i in range(3):
		out.write('  {0:13.8f}  {1:13.8f}  {2:13.8f}\n'.format(lat[i,0],lat[i,1],lat[i,2]))
	
	# Now loop through listed structures and write coords for each atom
	for i, struct in enumerate(systemlist):

		# Make sure we're dealing with cartesian coords
		struct.d2c()
			
		out.write("PRIMCOORD {0}\n".format(i+1))
		out.write("{0} 1\n".format(struct.numatoms))
		for a in struct.atomlist:
			num = str(atomic_numbers[a.symbol])
			pos = '{0:13.9f} {1:13.9f} {2:13.9f}'.format(a.x, a.y, a.z)
			out.write("  {0:2}   {1}\n".format(num,pos))



def read_qe(system,filepath):
	"""Tries to parse a structure from a quantum espresso input file"""
	# Open the file
	data = open(filepath,'r')
	lines = data.readlines()
	data.close()

	ibrav = 0
	celldm = np.zeros(6)
	num_atoms = 0
	num_types = 0


	def parse_keyword(line, keyword):
		# Parse keywords from the lines
		chunks = line.replace("="," ").split()
		for j, chunk in enumerate(chunks):
			if keyword in chunk:  # Then we want the key from the next chunk
				var = chunks[j+1].replace(",","")
				return var


	for i, line in enumerate(lines):

		# Get type of Bravais lattice and lattice parameters
		if "ibrav" in line:
			ibrav = int(parse_keyword(line, "ibrav"))
			print "found ibrav {0}".format(ibrav)
			
		elif "celldm(1)" in line:
			celldm[0] = float(parse_keyword(line, "celldm(1)"))
		elif "celldm(2)" in line:
			celldm[1] = float(parse_keyword(line, "celldm(2)"))
		elif "celldm(3)" in line:
			celldm[2] = float(parse_keyword(line, "celldm(3)"))
		elif "celldm(4)" in line:
			celldm[3] = float(parse_keyword(line, "celldm(4)"))
		elif "celldm(5)" in line:
			celldm[4] = float(parse_keyword(line, "celldm(5)"))
		elif "celldm(6)" in line:
			celldm[5] = float(parse_keyword(line, "celldm(6)"))

		# Get number of atoms and atom types
		elif "nat" in line:
			num_atoms = int(parse_keyword(line, "nat"))
		elif "ntyp" in line:
			num_types = int(parse_keyword(line, "ntyp"))
			
	# Now figure out what the lattice vectors are
	basis = np.zeros( (3,3) )	
	# Convert from bohr to angstrom
	celldm[0] = celldm[0]*0.529177249
	a = celldm[0]

	if ibrav == 1:		# Cubic P (sc)
		basis[0][0] = a
		basis[1][1] = a
		basis[2][2] = a
	
	elif ibrav == 2:		# Cubic F (fcc)
		basis = np.array([[ -a/2, 0, a/2 ], \
						  [  0, a/2, a/2 ], \
						  [ -a/2, a/2, 0 ]])

	elif ibrav == 3:		# Cubic I (bcc)
		basis = np.array([[  a/2,  a/2,  a/2],\
						  [ -a/2,  a/2,  a/2],\
						  [ -a/2, -a/2,  a/2]])
	
	elif ibrav == 4:		# Hexagonal and Trigonal P
		c = a*celldm[2]
		basis = np.array([[ a, 0, 0], \
						  [ -a/2, a*math.sqrt(3)/2, celldm[2] ], \
						  [ 0, 0, celldm[2] ]])

	elif ibrav == 6:		# Tetragonal P (st)
		basis = np.array([[ a, 0, 0], \
						  [ 0, a, 0], \
						  [ 0, 0, celldm[2]]])
	
	elif ibrav == 7:		# Tetragonal I (bct)
		c = a*celldm[2]
		basis = np.array([[ a/2, -a/2, c/2 ],\
						  [ a/2,  a/2, c/2 ],\
						  [-a/2, -a/2, c/2]])

	elif ibrav == 8:		# Orthorhombic P
		b = a*celldm[1]
		c = a*celldm[2]
		basis = np.array([[ a, 0, 0],\
						  [ 0, b, 0],\
						  [ 0, 0, c]])

	elif ibrav == 9:		# Orthorhombic base-centered (bcc)
		b = a*celldm[1]
		c = a*celldm[2]
		basis = np.array([[ a/2, b/2, 0],\
						  [-a/2, b/2, 0],\
						  [   0,   0, c]])

	elif ibrav == 10:		# Orthorhombic face-centered
		b = a*celldm[1]
		c = a*celldm[2]
		basis = np.array([[ a/2, 0, c/2],\
						  [ a/2, b/2, 0],\
						  [ 0, b/2, c/2]])

	elif ibrav == 11:		# Orthorhombic body-centered
		b = a*celldm[1]
		c = a*celldm[2]
		cos_ac = celldm[4]; beta = math.acos(cos_ac)
		basis = np.array([[ a, 0, 0 ],\
						  [ 0, b, 0 ],\
						  [ c*cos_ac, 0, c*math.sin(beta) ]])

	elif ibrav == 12:		# Monoclinic P, unique axis c
		b = a*celldm[1]
		c = a*celldm[2]
		cos_ab = celldm[3]
		gamma = math.acos(cos_ab)
		basis = np.array([[                 a,                 0, 0 ],\
		                  [ b*math.cos(gamma), b*math.sin(gamma), 0 ],\
						  [                 0,                 0, c ]])

	elif ibrav == -12:      # Monoclinic P, unique axis b
		b = a*celldm[1]
		c = a*celldm[2]
		cos_ac = celldm[4]
		beta = math.acos(cos_ac)
		basis = np.array([[                a, 0,                0 ],\
		                  [                0, b,                0 ],\
						  [ c*math.cos(beta), 0, c*math.sin(beta) ]])

	elif ibrav == 13:		# Monoclinic base-centered
		b = a*celldm[1]
		c = a*celldm[2]
		cos_ab = celldm[3]; gamma = math.acos(cos_ab)
		basis = np.array([[    a/2,       0,           -c/2],\
						  [ b*cos_ab, b*math.sin(gamma),  0],\
						  [    a/2,       0,            c/2]])

	elif ibrav == 14:       # Triclinic
		a = celldm[0];      b = a*celldm[1];    c = a*celldm[2]
		cos_bc = celldm[3]; cos_ac = celldm[4]; cos_ab = celldm[5]
		alpha = math.acos(cos_bc)
		beta  = math.acos(cos_ac)
		gamma = math.acos(cos_ab)

		basis[0][0] = a
		basis[1][0] = b*cos_ab
		basis[1][1] = b*math.sin(gamma)
		basis[2][0] = c*cos_ac
		basis[2][1] = c * ( cos_bc - cos_ac*cos_ab ) / math.sin(gamma)
		basis[2][2] = c/math.sin(gamma) * math.sqrt( 1 + 2*cos_bc*cos_ac*cos_ab \
					                                - cos_bc**2 - cos_ac**2 - cos_ab**2 )	

	system.basis = basis
#	print basis

	# Now need to read atomic coordinates
	system.atomlist = []
	symbols = ["x" for x in xrange(num_atoms)]
	startline = 0
	# Make another pass through lines to figure out where they start
	for i, line in enumerate(lines):
		if "ATOMIC_POSITIONS" in line:
			startline = i+1
			# Read if coordinates are cartesian or fractional
			if "crystal" in line or "direct" in line:
				system.isdirect = True
			elif "cartesian" in line:
				system.isdirect = False
			else:
				print "Could not determine whether coordinates are direct or cartesian"

	# Now pull atomic coordinates from lines
	for j in range(num_atoms):
		coords = lines[ startline + j ].split()
#		print coords
		symbols[j] = coords[0]
		pos = coords[1:4]
		for k in range(3):
			pos[k] = float( pos[k] )
		system.atomlist.append( Atom(pos[0],pos[1],pos[2]) )

	# Assign each atom a symbol and unique id
	atomtypes = []
	type_pop = [0 for x in xrange(num_types)]
	for i, atom in enumerate(system.atomlist):
		atom.symbol = symbols[i]
		if atom.symbol not in atomtypes:
			atomtypes.append( atom.symbol )
		s = atomtypes.index( atom.symbol )
		type_pop[s] += 1
		atom.id = "{0}{1}".format(atom.symbol,type_pop[s])
#		print atom.id+"  "+str(atom.pos)

	# Make new atom list so that atoms of specific types are together
	sortedlist = []
	for s in atomtypes:
		for a in system.atomlist:
			if a.symbol == s:
				sortedlist.append(a)
	system.atomlist = sortedlist

	# Finalize some things that are necessary
	system.atomtypes = atomtypes
	system.atompop = type_pop
	system.numatoms = num_atoms
	system.selective = False
	print "Finished reading QE file"


class Atom:
	"""Class for representing a single atom"""
	def __init__(self,x,y,z):
		self.symbol = 'X'
		self.x = float(x)
		self.y = float(y)
		self.z = float(z)
		self.pos = np.array([self.x, self.y, self.z])
	
	def __repr__(self):
		txt = '{0:2}  {1:13.8f}  {2:13.8f} {3:13.8f}'.format(self.symbol,self.x,self.y,self.z)
		return txt

	def newpos(self,newpos):
		self.pos = newpos
		self.x, self.y, self.z = self.pos

	def find_nearest(self,species,system,num_results=1):
		"""Find the nearest atom(s) of a particular type in a particular system"""
		possible_nearest = []
		for a in (a for a in system.atomlist if (a.symbol == species and a.id != self.id)):
			dist = distfromA2B( system.basis, self.pos, a.pos )
			possible_nearest.append(( a, dist ))
		# Sorting this list may not scale well, I might come up with something else later
		possible_nearest.sort(key=lambda entry: entry[1])
		# Return as many results as are requested
		return possible_nearest[:num_results]


class Molecule:
	"""Class for representing a group of atoms"""
	def __init__(self):
		self.atomlist = []
	
	def get_boxface_centers(self):
		"""Define a box in cartesian space to enclose the molecule"""
		xmin = np.array([a.x for a in self.atomlist]).min()
		xmax = np.array([a.x for a in self.atomlist]).max()
		ymin = np.array([a.y for a in self.atomlist]).min()
		ymax = np.array([a.y for a in self.atomlist]).max()
		zmin = np.array([a.z for a in self.atomlist]).min()
		zmax = np.array([a.z for a in self.atomlist]).max()

		bcface1 = np.array([ xmin, (ymax+ymin)/2, (zmax+zmin)/2 ])
		bcface2 = np.array([ xmax, (ymax+ymin)/2, (zmax+zmin)/2 ])
		acface1 = np.array([ (xmin+xmax)/2, ymin, (zmin+zmax)/2 ])
		acface2 = np.array([ (xmin+xmax)/2, ymax, (zmin+zmax)/2 ])
		abface1 = np.array([ (xmin+xmax)/2, (ymin+ymax)/2, zmin ])
		abface2 = np.array([ (xmin+xmax)/2, (ymin+ymax)/2, zmax ])
		face_centers = [ bcface1, bcface2, acface1, acface2, abface1, abface2 ]
		self.face_names = ["bc1", "bc2", "ac1", "ac2", "ab1", "ab2"]
		return face_centers


def find_dihydrogen_bonds(system,dist_cutoff=2.5):
	pass	


class Dihydrogen_bond:
	"""Class for representing a dihydrogen bond"""
	def __init__(self,neighbor1,atom1,atom2,neighbor2):
		self.h1 = atom1
		self.h2 = atom2
		self.neighbor1 = neighbor1
		self.neighbor2 = neighbor2
		


class Structure:
	"""Class for representing a crystal structure"""
	def __init__(self,structurefile="POSCAR"):
		pass


	def update_atomtypes(self):
		"""update the atomtypes list"""
		atomtypes = []
		# Could get unique entries, but order matters
		for atom in self.atomlist:
			if atom.symbol not in atomtypes:
				atomtypes.append(atom.symbol)
		return atomtypes


	def update_atompop(self):
		"""Re-count atomic species and the number of each"""
		atompop = []
		# Count how many of each species of atom are in the system
		for sym in self.atomtypes:
			pop = sum([1 for a in self.atomlist if a.symbol == sym])
			atompop.append(pop)

		return atompop

	
	def d2c(self):
		"""Convert from fractional to cartesian coordinates if
		necessary"""
		if self.isdirect:
			basis = self.basis
			for atom in self.atomlist:
				atom.newpos( np.dot( atom.pos, basis ) )
			self.isdirect = False
		else:
			pass


	def get_bader_charges(self,acf_file):
		"""If an ACF.dat file is present, grab Bader charges from it"""
		print "getting bader charges from {0}".format(acf_file)
		chg_file = open(acf_file,'r')
		lines = chg_file.readlines()
		for i, atom in enumerate(self.atomlist):
			atom.charge = float( lines[i+2].split()[4] )
