# protein_structure.py
# Justin Clark - May 2014

from math import sqrt, log10
from itertools import combinations
import atom
import matplotlib.pyplot as plt
import numpy as np
from scipy.spatial.distance import pdist, squareform
import cProfile
import os.path
import sys

class ProteinStructure:
    def __init__(self):
        self.atoms = []
        self.sequence = []
        self.debug = False

    def from_pdb_file(self, filename, chain_id = 'A'):
        f = open(filename, 'r')
        self.name = os.path.splitext(os.path.basename(filename))[0]
        self.chain = chain_id
        for line in f:
            if line.startswith('SEQRES') and line[11] == chain_id:
                self.sequence += line[19:].split()
            elif line.startswith('ATOM') and line[21] == chain_id:
                aid = int(line[6:11].strip())
                name = line[12:16].strip()
                residue = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[39:46].strip())
                z = float(line[46:54].strip())
                element = line[76:78].strip()
                self.atoms.append(atom.Atom(aid, name, residue, x, y, z, element))

    def filter_atoms(self, selected_atom_name):
        self.atoms = [atom for atom in self.atoms if atom.name == selected_atom_name]
        for i, atom in enumerate(self.atoms):
            self.atoms[i].aid = i

    def to_pdb_file(self, filename, residue_sequence = None):
        if residue_sequence == None: residue_sequence = self.sequence
        f = open(filename, 'w')
        f.write('CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n')
        f.write('ORIGX1      1.000000  0.000000  0.000000        0.00000\n')
        f.write('ORIGX2      0.000000  1.000000  0.000000        0.00000\n')
        f.write('ORIGX3      0.000000  0.000000  1.000000        0.00000\n')
        f.write('SCALE1      1.000000  0.000000  0.000000        0.00000\n')
        f.write('SCALE2      0.000000  1.000000  0.000000        0.00000\n')
        f.write('SCALE3      0.000000  0.000000  1.000000        0.00000\n')
        for i, atom in enumerate(self.atoms):
            line = [' '] * 80
            line[0:5] = 'ATOM'
            line[7:11] = str(atom.aid).rjust(4)
            line[13:16] = atom.name.ljust(3)
            line[17:20] = residue_sequence[i]
            line[21:22] = 'A'
            line[22:26] = str(atom.aid).rjust(4)
            line[31:38] = '{0:.3f}'.format(atom.x).rjust(7)
            line[39:46] = '{0:.3f}'.format(atom.y).rjust(7)
            line[47:54] = '{0:.3f}'.format(atom.z).rjust(7)
            line[77:78] = atom.element
            last_line = ''.join(line) + '\n'
            f.write(last_line)
        f.write('TER    ' + str(len(self.atoms)).rjust(4) + (' ' * 6) + last_line[17:27] + '\n')
        f.close()

    def euclidean_distance(self, coords_1, coords_2):
        return sqrt(sum([(i - j)**2 for i, j in zip(coords_1, coords_2)]))

    def distance_between_atoms(self, atom1, atom2):
        return self.euclidean_distance(atom1.get_coords(), atom2.get_coords())

    def as_distance_matrix(self, selected_atom_name):
        import distance_matrix as dm
        coords = [a.get_coords() for a in self.atoms if a.name == selected_atom_name]
        matrix = squareform(pdist(coords))
        return dm.DistanceMatrix(self, matrix, selected_atom_name)

    def as_distance_matrix_old(self, selected_atom_name):
        import distance_matrix as dm
        atoms = [a for a in self.atoms if a.name == selected_atom_name]
        num_atoms = len(atoms)
        matrix = np.zeros((num_atoms, num_atoms))
        atoms_with_indexes = zip(range(num_atoms), atoms)
        atom_combos = combinations(atoms_with_indexes, 2)
        for atom1, atom2 in atom_combos:
            distance = self.distance_between_atoms(atom1[1], atom2[1])
            matrix[atom1[0], atom2[0]] = distance
            matrix[atom2[0], atom1[0]] = distance
        return dm.DistanceMatrix(self, matrix, selected_atom_name)

    def get_error_map(self, true_contact_map):
        distance_matrix = self.as_distance_matrix(self.atoms[0].name)
        my_contact_map = distance_matrix.as_contact_map(true_contact_map.threshold)
        return true_contact_map - my_contact_map

    def correct(self, true_contact_map):
        distance_matrix = self.as_distance_matrix(self.atoms[0].name)
        error_map = self.get_error_map(true_contact_map)
        for i in range(len(self.atoms)):
            F = [0.0,0.0,0.0]
            i_coords = self.atoms[i].get_coords()
            if self.debug: print 'Old coords:', i_coords
            r = self.get_mobility_radius(i, distance_matrix, true_contact_map, error_map)
            if self.debug: print 'Mobility radius:', r
            error_indices = np.nonzero(error_map[i])[0]
            if len(error_indices) == 0: continue
            if self.debug: print 'Num not well placed:', len(error_indices)
            for j in error_indices:
                j_coords = self.atoms[j].get_coords()
                dist_ij = distance_matrix.matrix[i,j]
                Fj = (np.array(i_coords) - np.array(j_coords)) / dist_ij
                if true_contact_map.matrix[i,j]:
                    F -= Fj
                else:
                    F += Fj
            mag_F = self.euclidean_distance(F, i_coords)
            i_coords += F * (r / mag_F)
            if self.debug: print 'New coords:', i_coords
            self.atoms[i].set_coords(i_coords)

    def perturb(self, true_contact_map, epsilon):
        epsilon_divisor = 10.0
        error_map = self.get_error_map(true_contact_map)
        distance_matrix = self.as_distance_matrix(self.atoms[0].name)
        under_but_close = ((distance_matrix.matrix <= true_contact_map.threshold)
                & (distance_matrix.matrix > true_contact_map.threshold - epsilon)
                & true_contact_map.matrix)
        over_but_close = ((distance_matrix.matrix > true_contact_map.threshold)
                & (distance_matrix.matrix < true_contact_map.threshold + epsilon)
                & ~true_contact_map.matrix)
        under_but_close_i = np.transpose(np.nonzero(under_but_close & (error_map == 0)))
        over_but_close_i = np.transpose(np.nonzero(over_but_close & (error_map == 0)))
        if len(under_but_close_i) > 0:
            under_i = under_but_close_i[under_but_close_i[:,0] < under_but_close_i[:,1]]
            for coords in under_i:
                self.change_interatom_distance(self.atoms[coords[0]], self.atoms[coords[1]], -1 * epsilon / epsilon_divisor)
        if len(over_but_close_i) > 0:
            over_i = over_but_close_i[over_but_close_i[:,0] < over_but_close_i[:,1]]
            for coords in over_i:
                self.change_interatom_distance(self.atoms[coords[0]], self.atoms[coords[1]], epsilon / epsilon_divisor)

    def change_interatom_distance(self, atom1, atom2, by_distance):
        # Negative by_distance is closer
        # Each atom moves half by_distance
        atom1_coords = atom1.get_coords()
        atom2_coords = atom2.get_coords()
        atom1_towards_atom2 = np.array(atom2_coords) - np.array(atom1_coords)
        atom2_towards_atom1 = -1 * atom1_towards_atom2
        magnitude = self.euclidean_distance(atom1_coords, atom2_coords)
        move_atom1_by = atom1_towards_atom2 * (-by_distance / 2.0 / magnitude)
        move_atom2_by = atom2_towards_atom1 * (-by_distance / 2.0 / magnitude)
        atom1.set_coords(np.array(atom1.get_coords()) + np.array(move_atom1_by))
        atom2.set_coords(np.array(atom2.get_coords()) + np.array(move_atom2_by))
        if self.debug: print '    ------------------    '
        if self.debug: print 'Move distance:', by_distance
        if self.debug: print 'Old Atom 1:', atom1_coords, '- Old Atom 2:', atom2_coords
        if self.debug: print 'Difference:', atom1_towards_atom2
        if self.debug: print 'Shift by:', move_atom1_by
        if self.debug: print 'New Atom 1:', atom1.get_coords(), '- New Atom 2:', atom2.get_coords()
        if self.debug: print 'Old Distance between:', magnitude
        if self.debug: print 'New distance between:', self.distance_between_atoms(atom1, atom2)


    def get_mobility_radius(self, atom_i, distance_matrix, contact_map, error_map):
        true_pairs_above_threshold = contact_map.matrix[atom_i] == False
        dist_above_threshold = distance_matrix.matrix[atom_i] > contact_map.threshold
        true_and_above = true_pairs_above_threshold & dist_above_threshold

        true_pairs_below_threshold = contact_map.matrix[atom_i] == True
        dist_below_threshold = distance_matrix.matrix[atom_i] <= contact_map.threshold
        true_and_below = true_pairs_below_threshold & dist_below_threshold

        mins = distance_matrix.matrix[atom_i,true_and_above]
        maxes = distance_matrix.matrix[atom_i,true_and_below]

        if len(mins) == 0 and len(maxes) == 0: return 0.0
        if len(mins) > 0:
            D0 = np.amin(mins)
        else:
            D0 = float('infinity')

        if len(maxes) > 0:
            D1 = np.amax(maxes)
        else:
            D1 = 0.0

        r = min(D0 - contact_map.threshold, contact_map.threshold - D1)
        return r

    def plot_error(self, true_contact_map):
        my_contact_map = self.as_distance_matrix(self.atoms[0].name).as_contact_map(true_contact_map.threshold)
        error_map = self.get_error_map(true_contact_map)
        f, (ax1, ax2, ax3) = plt.subplots(ncols=3)
        ax1.matshow(true_contact_map.matrix, cmap=plt.cm.binary)
        ax1.set_title('True contact map')
        ax2.matshow(my_contact_map.matrix, cmap=plt.cm.binary)
        ax2.set_title('Contact map of guessed coords')
        ax3.matshow(error_map, cmap=plt.cm.binary)
        ax3.set_title('Not well placed: ' + str(self.num_not_well_placed(true_contact_map, error_map)))
        plt.show()

    def num_not_well_placed(self, true_contact_map, error_map = None):
        if error_map == None: error_map = self.get_error_map(true_contact_map)
        return len(np.nonzero(error_map)[0]) / 2

    def set_origin(self, origin):
        translation = -1 * np.array(origin)
        for atom in self.atoms:
            atom.set_coords(np.array(atom.get_coords()) + translation)

    def move_origin_to_first_atom(self):
        self.set_origin(self.atoms[0].get_coords())

def epsilon_delta(epsilon):
    max_delta = 0.010
    min_delta = 0.0001
    max_loops = 200
    delta_scale = np.logspace(log10(max_delta), log10(min_delta), max_loops)
    while sum(delta_scale) <= epsilon:
        max_loops += 1
        delta_scale = np.logspace(log10(max_delta), log10(min_delta), max_loops)
    return delta_scale

def comar(contact_map):
    epsilon = round(contact_map.n / 300.0, 1) + 1
    epsilon_delta_scale = epsilon_delta(epsilon)
    distance_matrix = contact_map.guess_distance_matrix()
    protein_structure = distance_matrix.guess_protein_structure()
    protein_structure.correct(contact_map)
    not_well_placed = protein_structure.num_not_well_placed(contact_map)
    print 'Error after initial guess:', not_well_placed
    loop_count = 0
    history_of_improvement = []
    while not_well_placed > 0 and epsilon >= 0:
        protein_structure.perturb(contact_map, epsilon)
        protein_structure.correct(contact_map)
        prev_not_well_placed = not_well_placed
        not_well_placed = protein_structure.num_not_well_placed(contact_map)
        newly_well_placed = prev_not_well_placed - not_well_placed
        history_of_improvement.append(newly_well_placed)
        print (str(loop_count) + '/' + str(len(epsilon_delta_scale))).rjust(9),\
        '- Error:', str(not_well_placed).rjust(4),\
        '- Epsilon:', str(round(epsilon, 3)).ljust(5),\
        '- Decreased:', newly_well_placed
        epsilon -= epsilon_delta_scale[loop_count]
        loop_count += 1

    return protein_structure

def compare(true_protein_struture, guessed_protein_structure):
    true_protein_struture.filter_atoms('CA')
    true_protein_struture.move_origin_to_first_atom()
    true_protein_struture.to_pdb_file('true_' +
            true_protein_struture.name + '_' +
            true_protein_struture.chain + '.pdb')
    guessed_protein_structure.move_origin_to_first_atom()
    guessed_protein_structure.to_pdb_file('guessed_' + 
            true_protein_struture.name + '_' +
            true_protein_struture.chain + '.pdb',
            true_protein_struture.sequence)

def run(filename, chain = 'A', distance_threshold = 12):
    print "Running for ", filename, '- Chain', chain, '- Distance threshold', distance_threshold
    distance_threshold = float(distance_threshold)
    selected_atom_name = 'CA'
    true_protein_struture = ProteinStructure()
    true_protein_struture.from_pdb_file(filename, chain)
    distance_matrix = true_protein_struture.as_distance_matrix(selected_atom_name)
    true_contact_map = distance_matrix.as_contact_map(distance_threshold)

    guessed_protein_structure = comar(true_contact_map)
    compare(true_protein_struture, guessed_protein_structure)

    # Uncomment to show generated distance matrix
    #distance_matrix.plot().show()

    # Uncomment to show generated contact map
    #true_contact_map.plot().show()

    # Uncomment to see independent components highlighted in gray boxes
    #true_contact_map.plot_with_components().show()

    # Uncomment to show contact maps of given PDB file and generated structure
    # side by side, with any remaining errors shown in difference matrix
    #guessed_protein_structure.plot_error(true_contact_map)

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print 'Usage: python', sys.argv[0], '<PDB file> <chain id (e.g. "A")> <threshold (e.g. "12")>'
        exit()
    run(*sys.argv[1:])
