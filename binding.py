#!/usr/bin/env python

"""
Try binding small molecules to the interior surface of a porous
material.
For a MOF, consider one linker and metal site. These must be specified
in the preamble. Try binding a small molecule (also specified in the
preamble) to particular binding sites (you guessed it, to be specified
in the preamble).

If the small molecule is considered to be within a rectangular prism in
cartesian space, make a structure where each face of that prism faces
the binding site.
Write these structures into a "docking" subdirectory

Proper usage:
python binding.py [empty mof poscar]
"""


import sys
import os
import numpy as np
import structure_tools as st
import copy

poscarfile = sys.argv[1]

### Define list of which atom ids are in linker of interest
linker_atom_ids = [
	"N3", "N4", "N9", "N10", "N15", "N16",
	"O3", "O4",
	"C3", "C4", "C9", "C10", "C15", "C16", "C21", "C22", "C27", "C28", "C33", "C34",
	"H3", "H4", "H9", "H10"
	]

# Which metal are you interested in adsorbing stuff to?
metal_ids = ["Mn1"]

# Is anything already adsorbed? We should allow it to relax
preadsorbed = False
preadsorbed_atom_ids = [
#	"O7", "H13"
	]

# Which atoms should we try adsorbing additional stuff to?
binding_site_ids = [
#	"N3", "N10", "N15",
#	"C4", "C10", "C15", "C22", "C27", "C34",
#	"O3", "O4"
	"N3", "N15", "C10", "C34"
	]

### Define the thing we're trying to bind to stuff
binding_atoms = []
binding_dist = 2.0

binding_atoms.append( ["Cl", "_A1",   3.2119,    2.7950,   -0.4902] )  
binding_atoms.append( ["Cl", "_A2",   3.5913,   -2.5846,    0.5305] )  
binding_atoms.append( ["Cl", "_A3",  -4.0248,    2.1408,    0.7695] )  
binding_atoms.append( ["Cl", "_A4",   5.0359,    0.2272,    0.0568] )  
binding_atoms.append( ["Cl", "_A5",  -5.4192,   -0.6169,   -0.1901] )  
binding_atoms.append( ["C" , "_A1",   0.5217,   -0.1094,   -0.0409] )  
binding_atoms.append( ["C" , "_A2",  -0.9209,   -0.2167,   -0.0716] )  
binding_atoms.append( ["C" , "_A3",   1.3002,   -1.2416,    0.1985] )  
binding_atoms.append( ["C" , "_A4",   1.1336,    1.1264,   -0.2507] )  
binding_atoms.append( ["C" , "_A5",  -1.7101,    0.8674,    0.3124] )  
binding_atoms.append( ["C" , "_A6",  -1.5224,   -1.4053,   -0.4854] )  
binding_atoms.append( ["C" , "_A7",   2.5243,    1.2300,   -0.2211] )  
binding_atoms.append( ["C" , "_A8",   2.6910,   -1.1380,    0.2281] )  
binding_atoms.append( ["C" , "_A9",  -3.1007,    0.7631,    0.2825] )  
binding_atoms.append( ["C" , "_A10",  -2.9130,   -1.5097,   -0.5152] )  
binding_atoms.append( ["C" , "_A11",   3.3030,    0.0978,    0.0182] )  
binding_atoms.append( ["C" , "_A12",  -3.7021,   -0.4254,   -0.1313] )  
binding_atoms.append( ["H" , "_A1",   0.8494,   -2.2110,    0.3971] )  
binding_atoms.append( ["H" , "_A2",   0.5534,    2.0191,   -0.4723] )  
binding_atoms.append( ["H" , "_A3",  -1.2696,    1.7948,    0.6708] )  
binding_atoms.append( ["H" , "_A4",  -0.9503,   -2.2648,   -0.8240] )  
binding_atoms.append( ["H" , "_A5",  -3.3655,   -2.4410,   -0.8469] )  

# Anything else that should be allowed to relax?
other_unfrozen_ids = [
	]


#########################################################################
#                            MAIN                                       #
#########################################################################

# Read MOF structure from file
mof = st.Structure()
mof = st.read_poscar(mof, poscarfile)

# Define linker as a molecule
linker = st.Molecule()

# Populate linker molecule with atoms listed in header
for atom in mof.atomlist:
	if atom.id in linker_atom_ids:
		linker.atomlist.append(atom)


print "Atoms in linker:"
for atom in linker.atomlist:
	print atom


# Get direction from metal to preadsorbed molecule
if preadsorbed:
	M = [atom for atom in mof.atomlist if atom.id == metal_ids[0]][0]
	adsorbant = [atom for atom in mof.atomlist if atom.id == preadsorbed_atom_ids[0]][0]
	vec_M_to_adsorbant = st.vecfromA2B( mof.basis, M.pos, adsorbant.pos )
	

linker_backbone = [atom for atom in linker.atomlist if atom.symbol != "H"]
for atom in linker_backbone:
	
	# Calculate plane through this atom and neighbors
	# Get list of vectors to everything else in backbone
	veclist = [
		st.vecfromA2B( mof.basis, atom.pos, n.pos ) for n in linker_backbone
		if atom.id != n.id
		]
	
	# Sort the list of vectors and pick the shortest two
	veclist.sort( key = lambda vec: np.linalg.norm(vec) )
	AB = veclist[0]
	AC = veclist[1]

	# Define normal of aforementioned plane from vectors to neighbors
	n = np.cross(AB,AC)
	atom.n_hat = n / np.linalg.norm(n)

	# Make sure n_hat points to same side of linker as any preadsorbed molecules
	if preadsorbed:
		dotproduct = np.dot( atom.n_hat, vec_M_to_adsorbant )
		if dotproduct < 0:
			atom.n_hat = -1*atom.n_hat
	
	print atom.id+"  "+str(atom.n_hat)


# Make a vector from M towards center of cavity
for m in [a for a in mof.atomlist if a.id in metal_ids]:
	# Find all of the atoms within sum of ionic radii
	vec_from_neighbors = []
	for a in [a for a in mof.atomlist if a.id != m.id]:
		vec = st.vecfromA2B( mof.basis, a.pos, m.pos )
		radius_a = st.covalent_radii[st.atomic_numbers[a.symbol]]
		radius_m = st.covalent_radii[st.atomic_numbers[m.symbol]]
		if np.linalg.norm(vec) <= radius_a + radius_m:
			vec_from_neighbors.append( vec/np.linalg.norm(vec) )

	sum_vec = np.zeros(3)
	for v in vec_from_neighbors:
		sum_vec = sum_vec + v
	# Metal has no neighbors on pore side,
	# so make the sum of vectors from neighbors to M a unit vector
	# this unit vector should point towards the pore
	sum_vec = sum_vec / np.linalg.norm(sum_vec)
	m.n_hat = copy.deepcopy(sum_vec)
	print m.id+" "+str(m.n_hat)


# Freeze everything except the linker
mof.selective = True
for atom in mof.atomlist:
	if atom in linker.atomlist:
	 	atom.selective = ( "T", "T", "T" )
	elif atom.id in other_unfrozen_ids:
		atom.selective = ( "T", "T", "T" )
	elif atom.id in metal_ids:
		atom.selective = ( "T", "T", "T" )
	elif atom.id in preadsorbed_atom_ids:
		atom.selective = ( "T", "T", "T" )
	else:
		atom.selective = ( "F", "F", "F" )


# Make sure output path is available
if not os.path.exists("docking"):
	os.system("mkdir docking")


# Make a molecule object to manipulate
binding_mol = st.Molecule()
for a in binding_atoms:
	binding_mol.atomlist.append( st.Atom( a[2], a[3], a[4] ) )
	added_atom = binding_mol.atomlist[-1:][0]
	added_atom.symbol = a[0]
	added_atom.id = a[0]+str(a[1])
for i, a in enumerate(binding_mol.atomlist):
	a.symbol = binding_atoms[i][0]
	a.selective = ( "T", "T", "T" )


# Figure out size of molecule box
binding_mol.face_centers = binding_mol.get_boxface_centers()

# Define axes of molecule box for desired placements
box_xAxis = np.array([1.0,0.0,0.0])
box_yAxis = np.array([0.0,1.0,0.0])
box_zAxis = np.array([0.0,0.0,1.0])

# For each of the faces of the box, rotate box so that the face's normal lines 
#  up with desired axis
axes = np.array([box_xAxis, -1*box_xAxis, box_yAxis, -1*box_yAxis, box_zAxis, -1*box_zAxis])
for site in [a for a in mof.atomlist if a.id in binding_site_ids+metal_ids]:
	for i, ax in enumerate(axes):
		# Get parameters for rotating binding molecule to line up with n_hat
		rot_axis, theta = st.alignment_params( ax, site.n_hat )
		# Rotate binding molecule
		mol_coords = np.array([a.pos for a in binding_mol.atomlist])
		print "rotating molecule coords"
		rotated_molecule_coords = st.rotate_coords( mol_coords, rot_axis, theta )
#		print rotated_molecule_coords
#		print "rotating molecule face centers"
		rotated_face_centers = st.rotate_coords( binding_mol.face_centers, rot_axis, theta )
		# Shift coordinates so that desired face center lines up with origin
		print "shifting molecule coords by "+str(-1*rotated_face_centers[i])
		shifted_molecule_coords = st.shift_coords( rotated_molecule_coords, -1*rotated_face_centers[i] )
#		shifted_face_centers = st.shift_coords( rotated_face_centers, -1*rotated_face_centers[i] )
		# Shift coordinates and to where they would be in structure
		print "shifting molecule to sit in mof by "+str(site.pos + binding_dist*site.n_hat)
		in_mof_coords = st.shift_coords( shifted_molecule_coords, site.pos + binding_dist*site.n_hat )

		# Make a new molecule with these coordinates
		new_mol = copy.deepcopy(binding_mol)
		for j, atom in enumerate(new_mol.atomlist):
			atom.newpos( in_mof_coords[j] )


		# Check for distance to other atoms in the mof
		# If too close to the mof, rotate binding molecule about site.n_hat
		# until an appropriate arrangement can be found
		molecule_ok = False
		rotation_angle = 0
		while not (molecule_ok or rotation_angle >= 2*np.pi):
			breaker = False
			for a in new_mol.atomlist:
				for b in mof.atomlist:
					dist = st.distfromA2B( mof.basis, a.pos, b.pos )
					radius_a = st.covalent_radii[st.atomic_numbers[a.symbol]]
					radius_b = st.covalent_radii[st.atomic_numbers[b.symbol]]
					# If we find a short enough distance, break out of both for loops
					if dist <= radius_a + radius_b:
						breaker = True
						break
				if breaker:
					break
			if not breaker:
				molecule_ok = True

			else: # binding molecule is too close to the mof
				print "binding molecule doesn't fit. Rotating by 30 degrees"
				mol_coords = np.array([a.pos for a in new_mol.atomlist])
				rotation_angle += np.pi/4.0
				rotated_mol_coords = st.rotate_coords( mol_coords, site.n_hat, np.pi/4.0 )
				for j, atom in enumerate(new_mol.atomlist):
					atom.newpos( rotated_mol_coords[j] )
				
			

			# Check distance from A in mof to B in new_mol
			# If short enough distance found,
				# Rotate new_mol about site.n_hat by 15 degrees
			# Else, molecule_ok=True


		# Add atoms to a mof object with these coordinates
		new_mof = copy.deepcopy(mof)
		for j, atom in enumerate(new_mol.atomlist):
			new_atom = copy.deepcopy( atom )
			new_mof.atomlist.append( new_atom )


		# Print a new poscar for this structure
		new_mof.atomtypes = new_mof.update_atomtypes()
		new_mof.atompop = new_mof.update_atompop()
		subdir = "docking/{0}{1}".format(site.id,binding_mol.face_names[i])
		if not os.path.exists(subdir):
			os.system("mkdir {0}".format(subdir))

		newposcar = "{0}/POSCAR".format(subdir)
		st.write_poscar( new_mof, newposcar )
		print ""




