from tensorflow import keras
from tensorflow.keras import layers
from spektral.layers import GATConv
from tensorflow.keras.layers import Lambda
import tensorflow as tf
import numpy as np


"""A typical **MOMAP** control file momap.inp can be as follows:



.. code-block:: bash


    do_evc                 = 1
    do_spec_tvcf_ft        = 0
    do_spec_tvcf_spec      = 0
    do_ic_tvcf_ft          = 1
    do_ic_tvcf_spec        = 1
    do_isc_tvcf_ft         = 0
    do_isc_tvcf_spec       = 0
    do_spec_sums           = 0  

    &evc 
    ...
    /   

    &spec_tvcf
    ...
    /   

    &ic_tvcf 
    ...
    /   

    &isc_tvcf 
    ...
    /   

    &spec_sums
    ...
    /



The typical **evc** control with momap.inp can be as follows:


.. code-block:: bash

    do_evc          = 1                   #1 表示开启dushin计算，0 表示关闭

    &evc
      ffreq(1)      = "s0-freq.log"       #基态结果的日志文件
      ffreq(2)      = "s1-freq.log"       #激发态结果的日志文件
    /

or


.. code-block:: bash

    do_evc          = 1                   #1 表示开启dushin计算，0 表示关闭

    &evc
      ffreq(1)      = "s0-freq.log"       #基态结果的日志文件
      ffreq(2)      = "s1-freq.log"       #激发态结果的日志文件
      fnacme        = "nacme.log"         #nacme计算输出文件
      sort_mode     = 1                   #模式
    /

or


.. code-block:: bash

    do_evc          = 1                   #1 表示开启dushin计算，0 表示关闭

    &evc
      ffreq(1)      = "s0-freq.log"       #基态结果的日志文件
      ffreq(2)      = "s1-freq.log"       #激发态结果的日志文件
      ftdipd        = "numfreq-es.out"    #dip计算输出文件
      sort_mode     = 1                   #模式
    /


or


.. code-block:: bash

    do_evc              = 1               #1 表示开启dushin计算，0 表示关闭

    &evc
      ffreq(1)          = "s0-freq.log"   #基态结果的日志文件
      ffreq(2)          = "s1-freq.log"   #激发态结果的日志文件
      if_add_int_coord  = .t.             #是否增加化学键
      def_del_int_coord = 1 3 0 0         #定义 1 3 原子间成键
                        = 1 4 0 0         #定义 1 4 原子间成键

    /

or


.. code-block:: bash

    do_evc                = 1               #1 表示开启dushin计算，0 表示关闭

    &evc
      ffreq(1)            = "s0-freq.log"   #基态结果的日志文件
      ffreq(2)            = "s1-freq.log"   #激发态结果的日志文件
      if_add_int_coord    = .t.             #是否增加化学键
      def_del_int_coord_c = "a" 20 37 39 0  #定义原子 20 37 39 键角
                            "d" 21 20 37 39 #定义原子 21 20 37 39 二面角
    /



"""

class StructuredDropout(layers.Layer):


    """定义读取的基态或激发态结果的日志文件，可由Gaussian等软件计算得到

    Args:
        “s0-freq.log” (Default)

    Example:
        >>> &evc
        >>>     ffreq(1) = "s0-freq.log"
        >>>     ffreq(2) = "s1-freq.log"        
    """


    def __init__(self, rate):
        super(StructuredDropout, self).__init__()
        self.rate = rate

    def call(self, inputs, training=None):
        if training:
            # Create a binary mask for the entire atom
            mask = tf.random.uniform(shape=(tf.shape(inputs)[0], tf.shape(inputs)[1], 1)) > self.rate
            return inputs * tf.cast(mask, tf.float32)
        return inputs

class EdModelBuilder:
    def __init__(self, num_atoms, num_atom_types, feature_dim, n, n_attn_heads): 
        self.num_atoms = num_atoms
        self.num_atom_types = num_atom_types
        self.feature_dim = feature_dim
        self.n = n
        self.n_attn_heads = n_attn_heads
        self.atom_types_unique = np.unique(np.load('data/atom_types.npy'))
        self.atom_type_to_index = {atom: idx for idx, atom in enumerate(self.atom_types_unique)}

    def build_model(self):  
        X_in = layers.Input(shape=(self.num_atoms * 3,))
        A_distance_in = layers.Input(shape=(self.num_atoms, self.num_atoms))
        A_angle_in = layers.Input(shape=(self.num_atoms, self.num_atoms))
        I_in = layers.Input(shape=(self.num_atoms, self.num_atom_types))
        X_reshaped = layers.Reshape((self.num_atoms, 3))(X_in)

        # Element-specific GNNs
        gnn_outputs = []
        for atom_type in self.atom_types_unique:
            idx = self.atom_type_to_index[atom_type]
            element_mask = layers.Lambda(lambda x: x[:,:,idx:idx+1])(I_in)
            x_element = layers.Multiply()([X_reshaped, element_mask])
            
            x_distance = GATConv(self.feature_dim, activation='relu', attn_heads=self.n_attn_heads)([x_element, A_distance_in])
            x_angle = GATConv(self.feature_dim, activation='relu', attn_heads=self.n_attn_heads)([x_element, A_angle_in])
            
            gnn_outputs.append(layers.Concatenate(axis=-1)([x_distance, x_angle]))

        x = layers.Add()(gnn_outputs)

        x = layers.TimeDistributed(layers.Dense(128, activation='relu'))(x)
        x = StructuredDropout(0.5)(x)
    
        electron_density = layers.TimeDistributed(layers.Dense(self.n*self.n*self.n))(x)
        electron_density = layers.Reshape((self.num_atoms, self.n, self.n, self.n), name='density_reshaped')(electron_density)

        combined_model = keras.Model(inputs=[X_in, A_distance_in, A_angle_in, I_in], outputs=electron_density)
        return combined_model

