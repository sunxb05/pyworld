import numpy as np
import tensorflow as tf
from tensorflow.keras.callbacks import LambdaCallback
import sys
from logger import Logger

class MinLRCallback(tf.keras.callbacks.Callback):
    def __init__(self, min_lr):
        super(MinLRCallback, self).__init__()
        self.min_lr = min_lr

    def on_epoch_end(self, epoch, logs=None):
        current_lr = tf.keras.backend.get_value(self.model.optimizer.lr)
        if current_lr < self.min_lr:
            tf.keras.backend.set_value(self.model.optimizer.lr, self.min_lr)
            print(f"Learning rate reached {current_lr}, resetting to {self.min_lr}")

class ModelTrainer:
    def __init__(self, model):
        self.model = model

    def custom_print_ed(self, epoch, logs):
        if epoch == 0:
            print(f"{'Epoch':<5} {'Density_Loss_Val':<15} {'Density_RMSE_Val':<17} {'Density_Loss_Train':<19} {'Density_RMSE_Train':<21} {'LR':<10}")
            print("-" * 90)  # 打印分割线
        
        loss_val_density = logs['val_loss']
        rmse_val_density = np.sqrt(loss_val_density)

        loss_train_density = logs['loss']
        rmse_train_density = np.sqrt(loss_train_density)

        current_lr = float(tf.keras.backend.get_value(self.model.optimizer.learning_rate))
        print(f"{epoch + 1:<5} {loss_val_density:<15.4f} {rmse_val_density:<17.4f} {loss_train_density:<19.4f} {rmse_train_density:<21.4f} {current_lr:<10.6f}")


    def custom_print_ef(self, epoch, logs):
        if epoch == 0:
            print(f"{'Epoch':<5} {'Energy_RMSE_Val':<15} {'Energy_RMSE_Train':<17} {'Force_RMSE_Val':<15} {'Force_RMSE_Train':<17} {'LR':<10}")
            print("-" * 75)  # 打印分割线

        rmse_val_energy = np.sqrt(logs['val_energy_loss'])
        rmse_train_energy = np.sqrt(logs['energy_loss'])

        rmse_val_force = np.sqrt(logs['val_forces_loss'])
        rmse_train_force = np.sqrt(logs['forces_loss'])

        current_lr = float(tf.keras.backend.get_value(self.model.optimizer.learning_rate))
        print(f"{epoch + 1:<5} {rmse_val_energy:<15.4f} {rmse_train_energy:<17.4f} {rmse_val_force:<15.4f} {rmse_train_force:<17.4f} {current_lr:<10.6f}")

    def train_ed_model(self, x_train, y_train_density, x_test, y_test_density, batch_size=20):
        original_stdout = sys.stdout
        sys.stdout = Logger("train_ed.log")

        initial_learning_rate = 0.001
        min_learning_rate = 1e-6
        steps_per_epoch = x_train[0].shape[0] // batch_size
        epochs_to_decay = 10
        lr_schedule = tf.keras.optimizers.schedules.ExponentialDecay(
            initial_learning_rate,
            decay_steps=epochs_to_decay * steps_per_epoch,
            decay_rate=0.999,
            staircase=True
        )

        optimizer = tf.keras.optimizers.Adam(learning_rate=lr_schedule)
        self.model.compile(optimizer=optimizer, 
                   loss={'density_reshaped': 'mse'}, 
                   metrics={'density_reshaped': ['mae']})
    
        custom_callback = LambdaCallback(on_epoch_end=lambda epoch, logs: self.custom_print_ed(epoch, logs))
        min_lr_callback = MinLRCallback(min_learning_rate)  # Create the minimum learning rate callback

        self.model.fit(x_train, y_train_density, epochs=10, batch_size=batch_size, validation_data=(x_test, y_test_density), callbacks=[custom_callback, min_lr_callback], verbose=0)

        sys.stdout = original_stdout



    def train_ef_model(self, x_train, y_train, x_test, y_test, batch_size=20):
        original_stdout = sys.stdout
        sys.stdout = Logger("train_ef.log")

        # You might need to modify the optimizer or its parameters for the EF model
        optimizer = tf.keras.optimizers.Adam(learning_rate=0.001)
        self.model.compile(optimizer=optimizer, loss={'energy': 'mse', 'forces': 'mse'}, metrics={'energy': 'mae', 'forces': 'mae'})
        custom_callback = LambdaCallback(on_epoch_end=lambda epoch, logs: self.custom_print_ef(epoch, logs))
        self.model.fit(x_train, y_train, epochs=5, batch_size=batch_size, validation_data=(x_test, y_test), callbacks=[custom_callback], verbose=0)

        sys.stdout = original_stdout

    def save_trained_model(self, filename="model.h5"):
        self.model.save(filename)