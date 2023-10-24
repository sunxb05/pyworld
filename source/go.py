from data_processing import DataProcessing
from ed_builder import EdModelBuilder
from ef_builder import EfModelBuilder  
from model_trainer import ModelTrainer
import sys
import numpy as np
def main():
    """
    comments1
    """

    # 加载和处理数据
    data_processor = DataProcessing("data/coordinates.npy", "data/atom_types.npy", "data/ed.npy", "data/ener.npy", "data/forces.npy", 25)
    x_train, y_train_density, y_train_energy, y_train_force, x_test, y_test_density, y_test_energy, y_test_force, num_atoms, num_atom_types = data_processor.get_train_test_data()

    # 构建电子密度模型
    ed_builder = EdModelBuilder(num_atoms, num_atom_types, 25, 25, 8)
    model = EdModelBuilder(num_atoms, num_atom_types, 25, 25, 8).build_model()
    print(model.summary())
    ed_model = ed_builder.build_model()

    # 训练电子密度模型
    ed_trainer = ModelTrainer(ed_model)
    ed_trainer.train_ed_model(x_train, y_train_density, x_test, y_test_density)

    # 保存电子密度模型
    ed_trainer.save_trained_model("ed_model.h5")
    
    ef_builder = EfModelBuilder(num_atoms, num_atom_types,25)
    ef_model = ef_builder.build_energy_force_model()

    # 训练能量和力模型
    ef_trainer = ModelTrainer(ef_model)
    ef_trainer.train_ef_model(y_train_density, [y_train_energy, y_train_force], y_test_density, [y_test_energy, y_test_force])

if __name__ == "__main__":
    main()