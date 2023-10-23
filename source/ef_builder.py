import numpy as np
import tensorflow as tf
from tensorflow.keras import layers, models, Input

class EfModelBuilder:
    def __init__(self, num_atoms, num_atom_types, grid_size):
        self.num_atoms = num_atoms
        self.num_atom_types = num_atom_types
        self.grid_size = grid_size
        self.atom_types = np.load('data/atom_types.npy')
        
    def _build_subnetwork(self, input_tensor, idx):
        # 为每种原子类型创建一个简化的3D CNN子网络
        x = layers.Conv3D(32, (3, 3, 3), activation='relu')(input_tensor)
        x = layers.MaxPooling3D((2, 2, 2))(x)
        x = layers.Flatten()(x)
        
        # 添加索引到层名称以确保名称的唯一性
        energy = layers.Dense(1, name=f'energy_{idx}')(x)
        force = layers.Dense(3, name=f'force_{idx}')(x)  # 预测x, y, z方向上的力
        return energy, force
        
    def build_energy_force_model(self):
        # 修改输入层以接受整体的电子密度张量
        input_tensor = Input(shape=(self.num_atoms, self.grid_size, self.grid_size, self.grid_size, 1), name='density_input')
        
        energy_outputs = []
        force_outputs = []
        
        for idx, atom_type in enumerate(self.atom_types[0]):  # 注意这里使用self.atom_types[0]来遍历atom_types
            # 使用Lambda层拆分整体的电子密度张量
            atom_density = layers.Lambda(lambda x: x[:, idx])(input_tensor)
            
            energy, force = self._build_subnetwork(atom_density, idx)
            energy_outputs.append(energy)
            force_outputs.append(force)
        
        # 将每个原子的能量输出相加以得到总能量
        total_energy = layers.Add(name='energy')(energy_outputs)
        
        # 将每个原子的力量输出连接在一起
        total_force = layers.Concatenate(name='forces')(force_outputs)
        
        model = models.Model(inputs=input_tensor, outputs=[total_energy, total_force])
        
        return model

