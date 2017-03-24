This repository contains a few samples of scripts that I have written for use in my research.
I will briefly describe what they do and how I use them in my research.

#### structure_tools.py

This is a module containing classes for objects describing crystal structures, atoms, and molecules as well as methods commonly used for other tasks.

#### QEtoVASP.py

Make a VASP-style POSCAR file from a Quantum Espresso input file.

#### vasp_d2c.py

Convert a VASP-style structure file from direct/fractional coordinates to cartesian coordinates.

#### neb_axsf

Used to create an animated XCrySDen structure file from a VASP nudged-elastic band calculation. The written file can be viewed in XCrySDen or VMD for easy viewing of the image structures.

#### vasp_xray

Given a VASP-style structure file, simulate an Xray diffraction pattern. Several parameters are modifiable in the script header (either to match experimental conditions or for fitting purposes). I was unable to find an existing tool to do this in a scriptable way from the command line, so I wrote one myself. As such, the default parameters may not be suitable for all cases.

#### vasp_vib 

Use the vibrational modes from a VASP OUTCAR file to calculate the zero-point energy contribution. Also calculate the vibrational contribution to enthalpy and entropy at T=300 Kelvin.

#### triangles.py

Given a VASP-style XDATCAR file, identify three given atoms and calculate the area of a triangle that they form for each frame in the file.

This has been useful for monitoring the size of the opening to a pore over the course of a molecular dynamics simulation.

#### binding.py

This script is used to automate the process of intelligently trying to bind small molecules to the inner surfaces of porous materials (in particular Metal-Organic Framework materials). The general procedure is to treat a small molecule as a rectangular prism, which can be placed on any of its six faces above a set of binding sites. If the small molecule doesn't fit, try rotating the box about the face normal. The user must specify the coordinates of the desired small molecule as well as the desired binding sites in the preamble. If The script will create subdirectories for each of the possible structures, which can then be optimized separately.
