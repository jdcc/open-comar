# contact_map.py
# Justin Clark - May 2014

import matplotlib.patches as mpatches
import matplotlib.collections as mcoll
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx
import random
import distance_matrix as dm

class ContactMap:
    def __init__(self, protein, matrix, threshold):
        self.protein = protein
        self.matrix = np.array(matrix, dtype=bool)
        self.threshold = threshold
        self.n = len(self.matrix)
        #self.thickness = self.get_thickness()
        self.min_split_size = 13
        self.debug = False

    def __sub__(self, other):
        return np.array(self.matrix, dtype=int) - np.array(other.matrix, dtype=int)

    def component_addresses(self):
        component_addresses = []
        offset_diag = self.offset_diag()
        coordinates = np.array([[x, y] for x, y in zip(offset_diag[0], offset_diag[1])])
        comp_top_left = np.array(coordinates[0] * 2)
        comp_bottom_right = coordinates[self.min_split_size - 1]
        while True:
            if self.debug: print 'This top left:', comp_top_left
            if self.debug: print 'This bottom right:', comp_bottom_right
            upper_right_box = self.matrix[0:comp_bottom_right[0],comp_bottom_right[1]:]
            contacts_in_box = np.transpose(np.nonzero(upper_right_box))
            if contacts_in_box.size == 0:
                if self.debug: print 'Empty upper right'
                if comp_bottom_right[0] - comp_top_left[0] >= self.min_split_size:
                    if self.debug: print 'New component:', comp_top_left[0], comp_bottom_right[0]
                    component_addresses.append([
                        comp_top_left[0],
                        comp_bottom_right[0] + 1
                        ])
                    comp_top_left = np.array([comp_bottom_right[0] + 1] * 2)
                else:
                    if self.debug: print 'Continue adding on'
                    comp_bottom_right += 1
            else:
                last_contact = contacts_in_box[-1]
                last_contact[1] += comp_bottom_right[1]
                if self.debug: print 'Last contact:', last_contact
                # Jump ahead to next possible coord
                coords_outside_upper_right = \
                        coordinates[coordinates[:,1] == last_contact[1] + 1]
                if len(coords_outside_upper_right) == 0:
                    # We're done here
                    break
                else:
                    comp_bottom_right = coords_outside_upper_right[0]
                if self.debug: print 'New bottom right:', comp_bottom_right

        # If the last component we're creating isn't large enough to stand on its
        # own, append to previous component
        if (len(component_addresses) > 0) and (
                    self.n - component_addresses[-1][1] < self.min_split_size):
            component_addresses[-1][1] = self.n - 1
        else:
            component_addresses.append([comp_top_left[0], self.n])
        if self.debug: print component_addresses
        return component_addresses

    def split(self):
        components = []
        component_addresses = self.component_addresses()
        components = map(
                lambda a: self.matrix[a[0]:a[1],a[0]:a[1]],
                component_addresses
                )
        return components

    def get_thickness(self):
        if hasattr(self, 'thickness'):
            return self.thickness

        col_thicknesses = []
        for i in range(self.n):
            col = self.matrix[:(i+1),i]
            last_false_index = self.index_of_last_false(col)
            if last_false_index:
                col_thicknesses.append(i - last_false_index)
            else:
                col_thicknesses.append(1)
        self.thickness = int(round(np.mean(col_thicknesses)))
        return self.thickness

    def offset_diag(self, offset=None):
        offset = offset or self.thickness
        rows = self.n - offset + 1
        return [range(rows), [i + offset - 1 for i in range(rows)]]

    def index_of_last_true(self, array):
        indices = np.nonzero(array)[0]
        if len(indices) == 0:
            return None
        return indices[-1]

    def index_of_last_false(self, array):
        indices = np.nonzero(np.invert(array))[0]
        if len(indices) == 0:
            return None
        return indices[-1]

    def guess_interatom_dist(self, i, j):
        if i == j: return 0.0
        if abs(i - j) == 1: return 3.8
        if abs(i - j) == 2: return 6 + random.uniform(-1.5, 1.5)
        if abs(i - j) == 3: return max(
                0, 7.5 + random.uniform(7.5 - self.threshold, self.threshold - 7.5))
        if abs(i - j) > 3:
            percent = (0.91 - self.threshold / 100) * self.threshold
            return percent + random.uniform(
                    -1 * self.threshold + percent, self.threshold - percent)

    def to_cmap_file(self, filename):
        f = open(filename, 'w')

    def guess_distance_matrix(self):
        distances = np.zeros((self.n, self.n))
        for i in range(self.n):
            for j in range(i, self.n):
                if self.matrix[i, j] == 1:
                    dist = self.guess_interatom_dist(i, j)
                    distances[i,j] = dist
                    distances[j,i] = dist
                else:
                    distances[i,j] = float('infinity')
        return dm.DistanceMatrix(self.protein, self.shortest_paths(distances))

    def shortest_paths(self, guessed_dist):
        graph = nx.Graph()
        graph.add_nodes_from(range(self.n))
        for i in range(self.n):
            for j in range(i, self.n):
                weight = guessed_dist[i, j]
                if weight < float('infinity'): graph.add_edge(i, j, weight=weight)
        shortest_paths = nx.all_pairs_dijkstra_path_length(graph)
        shortest = np.copy(guessed_dist)
        for row, cols in shortest_paths.iteritems():
            for col, weight in cols.iteritems():
                shortest[row][col] = weight
        return shortest

    def plot(self):
        plt.matshow(self.matrix, cmap=plt.cm.binary)
        plt.title('Protein: ' + self.protein.name + ' - Chain: ' + self.protein.chain + ' - Threshold: ' + str(self.threshold))
        return plt

    def plot_with_thickness(self):
        plot = self.plot()
        self.thickness = self.get_thickness()
        offset_diag = self.offset_diag()
        plot.plot([offset_diag[1][0], offset_diag[1][-1]], [offset_diag[0][0], offset_diag[0][-1]])
        return plot

    def plot_with_components(self):
        plot = self.plot()
        self.thickness = self.get_thickness()
        plot.title('Protein: ' + self.protein.name + ' - Chain: ' + self.protein.chain + ' - Threshold: ' + str(self.threshold))
        component_addresses = self.component_addresses()
        components = []
        for address in component_addresses:
            plot.gca().add_patch(plot.Rectangle(
                (address[0], address[0]),
                address[1] - address[0],
                address[1] - address[0],
                facecolor="#aaaaaa",
                alpha=0.4))
        return plot
