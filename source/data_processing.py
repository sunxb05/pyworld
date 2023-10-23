import numpy as np
from tensorflow.keras import utils
from scipy.spatial import distance

def compute_adjacency(coordinates, max_neighbors=2):
    num_samples, flat_dims = coordinates.shape
    num_atoms = flat_dims // 3
    coordinates = coordinates.reshape(num_samples, num_atoms, 3)
    adjacency_matrices_distance = []
    adjacency_matrices_angle = []

    for sample_idx in range(num_samples):
        coords = coordinates[sample_idx]
        dists = distance.cdist(coords, coords)  # 使用 scipy 计算距离
        sorted_neighbors = np.argsort(dists, axis=-1)[:, 1:max_neighbors+1]

        adjacency_distance = np.zeros((num_atoms, num_atoms))
        adjacency_angle = np.zeros((num_atoms, num_atoms))

        for i in range(num_atoms):
            adjacency_distance[i, sorted_neighbors[i]] = 1.0 / (1.0 + dists[i, sorted_neighbors[i]])

            # Compute angle for each pair of neighbors
            for j in range(len(sorted_neighbors[i]) - 1):
                for k in range(j+1, len(sorted_neighbors[i])):
                    v1 = coords[sorted_neighbors[i][j]] - coords[i]
                    v2 = coords[sorted_neighbors[i][k]] - coords[i]
                    
                    epsilon = 1e-7  # 避免除数为0
                    cosine_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + epsilon)
                    cosine_angle = np.clip(cosine_angle, -1, 1)  # 防止精度问题导致的cosine_angle>1或< -1
                    
                    adjacency_angle[i, sorted_neighbors[i][j]] = cosine_angle
                    adjacency_angle[sorted_neighbors[i][j], i] = cosine_angle
                    adjacency_angle[i, sorted_neighbors[i][k]] = cosine_angle
                    adjacency_angle[sorted_neighbors[i][k], i] = cosine_angle

        adjacency_matrices_distance.append(adjacency_distance)
        adjacency_matrices_angle.append(adjacency_angle)

    return np.array(adjacency_matrices_distance), np.array(adjacency_matrices_angle)


class DataProcessing:
    def __init__(self, coordinates_file, atom_types_file, electron_density_file, energy_file, force_file, n):
        self.coordinates = np.load(coordinates_file)
        self.atom_types_data = np.load(atom_types_file)
        self.electron_density = np.load(electron_density_file)
        self.energy = np.load(energy_file)
        self.force = np.load(force_file)
        self.n = n

    def process_data(self):
        # 1. Compute adjacency matrices
        adjacency_distance, adjacency_angle = compute_adjacency(self.coordinates, 2)

        # 2. Convert atom types to one-hot encodings
        ELEMENT_TO_INDEX = {
            "H": 0,
            "O": 1,
            # ... 可以添加其他元素的索引
        }
        num_samples, flat_dims = self.coordinates.shape
        num_atoms = flat_dims // 3
        num_atom_types = len(ELEMENT_TO_INDEX)

        atom_types = np.array([ELEMENT_TO_INDEX[at] for at in self.atom_types_data[0]])
        atom_features_onehot = utils.to_categorical(atom_types, num_classes=num_atom_types)
        atom_features = np.tile(atom_features_onehot[None, :], (num_samples, 1, 1))

        # 3. Split data into training and test sets
        train_percentage = 0.8
        train_samples = int(num_samples * train_percentage)

        x_train = [
            self.coordinates[:train_samples],
            adjacency_distance[:train_samples],
            adjacency_angle[:train_samples],
            atom_features[:train_samples]
        ]
        y_train_density = self.electron_density[:train_samples]
        y_train_density = y_train_density.reshape([-1, num_atoms, self.n, self.n, self.n])
        y_train_energy = self.energy[:train_samples]
        y_train_force = self.force[:train_samples]
    
        x_test = [
            self.coordinates[train_samples:],
            adjacency_distance[train_samples:],
            adjacency_angle[train_samples:],
            atom_features[train_samples:]
        ]
        y_test_density = self.electron_density[train_samples:]
        y_test_density = y_test_density.reshape([-1, num_atoms, self.n, self.n, self.n])
        y_test_energy = self.energy[train_samples:]
        y_test_force = self.force[train_samples:]

        return (x_train, y_train_density, y_train_energy, y_train_force, 
                x_test, y_test_density, y_test_energy, y_test_force, 
                num_atoms, num_atom_types)

    def get_train_test_data(self):
        return self.process_data()