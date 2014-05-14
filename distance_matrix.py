# distance_matrix.py
# Justin Clark - May 2014

import contact_map as cm
import matplotlib.pyplot as plt
import numpy as np
import copy
import atom

class DistanceMatrix:
    def __init__(self, protein, matrix, selected_atom_name = 'CA'):
        self.protein = protein
        self.matrix = np.array(matrix, dtype=float)
        self.n = len(self.matrix)
        self.selected_atom_name = selected_atom_name

    def as_coordinates(self):
        # Using Wu's algorithm
        dist = self.matrix
        coords = np.zeros((3, self.n))
        coords[0,1] = dist[0,1]

        # I think I have to check if these are in the same plane as the previous
        coords[0,2] = (dist[0,2]**2 - dist[1,2]**2
                ) / (2 * coords[0,1]) + coords[0,1] / 2
        coords[1,2] = sqrt(dist[0,2]**2 - coords[0,2]**2)
        coords[0,3] = (dist[0,3]**2 - dist[1,3]**2
                ) / (2 * coords[0,1]) + coords[0,1]/2

        coords[1,3] = (dist[1,3]**2 - dist[2,3]**2 - (
                coords[0,3] - coords[0,1])**2 + (
                coords[0,3] - coords[0,2])**2) / (
                2 * coords[1,2]) + coords[1,2] / 2

        #print dist[0,3], coords[0,3], coords[1,3]
        coords[2,3] = sqrt(dist[0,3]**2 - coords[0,3]**2 - coords[1,3]**2)
        print coords[:,0:4]
        self.triangulate(coords[:,0:4], 4)

    def as_contact_map(self, threshold):
        matrix = np.copy(self.matrix)
        matrix[matrix <= threshold] = True
        matrix[matrix > threshold] = False
        return cm.ContactMap(self.protein, matrix, threshold)

    def triangulate(self, four_coords, i):
        dist = self.matrix
        u1, v1, w1 = four_coords[:,0]
        u2, v2, w2 = four_coords[:,1]
        u3, v3, w3 = four_coords[:,2]
        u4, v4, w4 = four_coords[:,3]
        x1 = u1**2 + v1**2 + w1**2
        x2 = u2**2 + v2**2 + w2**2
        x3 = u3**2 + v3**2 + w3**2
        x4 = u4**2 + v4**2 + w4**2
        A = 2 * np.matrix([
                [u1 - u2, v1 - v2, w1 - w2],
                [u1 - u3, v1 - v3, w1 - w3],
                [u1 - u4, v1 - v4, w1 - w4]
                ])
        # I'm assuming the four points provided are the first four points
        b = np.matrix([
            (x1 - x2) - (dist[0,i]**2 - dist[1,i]**2),
            (x1 - x3) - (dist[0,i]**2 - dist[2,i]**2),
            (x1 - x4) - (dist[0,i]**2 - dist[3,i]**2)
            ])

        print A.I * b

    def guess_protein_structure(self):
        import protein_structure as ps
        dist = self.matrix
        D = np.zeros((self.n, self.n))
        X = np.zeros((self.n, 3))
        for i in range(self.n):
            for j in range(i, self.n):
                dij = (dist[i,0]**2 - dist[i,j]**2 + dist[j,0]**2) / 2
                D[i,j] = dij
                D[j,i] = dij
        u, s, v = np.linalg.svd(D)
        X = np.asmatrix(u[:,0:3]) * np.asmatrix(np.sqrt(np.diag(s[0:3])))
        protein = ps.ProteinStructure()
        for i in range(X.shape[0]):
            coords = X[i].getA1()
            a = atom.Atom(i, 'CA', '', coords[0], coords[1], coords[2], 'C')
            protein.atoms.append(a)
        return protein

    def plot(self):
        plt.imshow(self.matrix, cmap=plt.cm.PuBu)
        plt.colorbar()
        plt.title('Distance matrix for protein ' +
                self.protein.name +
                ' (' + self.selected_atom_name + ' atoms)')
        return plt

