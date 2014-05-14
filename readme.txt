Open COMAR
==========

Description
-----------
This is a modified implementation of the COMAR algorithm as given here:
https://dl.acm.org/citation.cfm?id=1435572

It attempts to recreate 3D protein structure from a given contact map.

Files
-----
  * readme.txt: this file
  * atom.py: Atom module
  * contact_map.py: ContactMap module
  * distance_matrix.py: DistanceMatrix module
  * protein_structure.py: ProteinStructure module and executable script
  * 1HXR.pdf: sample PDB file

Quick start
-----------

Install:

numpy
scipy
matplotlib
networkx

Run:

python protein_structure.py 1HXR.pdb

Code Usage
----------
Requirements:
numpy
scipy
matplotlib
networkx

Once requirements are installed, main entry into the code is handled by the
protein_structure.py module. Usage is fairly straight forward, and is given
when the the script is called with no command line parameters:

Usage: python protein_structure.py <PDB file> <chain id (e.g. "A")> <threshold (e.g. "12")>

A sample PDB file has been included, so that a fully functioning demo can run with
the command:

python protein_structure.py 1HXR.pdb

This will parse the PDB file, generate a contact map of the protein, recover the
structure of the protein, and spit out two PDB files: one prefixed with "true_"
which is just the alpha carbons and their coordinates from the given PDB file 
(minus all the other atoms), and one prefixed with "guessed_" which is the
reconstructed protein coordinates.

A number of commented lines are present in the bottom of the protein_structure.py
script that, when uncommented, will generate a few useful figures.

Status
------

It currently only takes a PDB file. It would need slight modification to take an
actual contact map file.

Bugs exist.

License
-------
GPLv3
