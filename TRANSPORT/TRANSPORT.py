"""Documentation about the MOMAP module."""

def do_transport_prepare(int):
    """是否生成预备文件, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_prepare      = 1      
    """

def do_transport_submit_HL_job(int):
    """是否开启计算转移积分, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_submit_HL_job      = 1      
    """

def do_transport_get_transferintegral(int):
    """计算计算转移积分, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_get_transferintegral      = 1      
    """

def do_transport_submit_RE_job(int):
    """计算重整能, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_submit_RE_job      = 1      
    """


def do_transport_get_re_evc(int):
    """使用 evc 程序分析重整能, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_get_re_evc      = 1      
    """

def do_transport_run_MC(int):
    """Monte_Carlo 模拟, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_run_MC      = 1      
    """


def do_transport_get_mob_MC(int):
    """计算迁移率, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_get_mob_MC      = 1      
    """



def do_transport_run_MC_temp(int):
    """不同温度下的Monte_Carlo 模拟, 1表示开启，0表示关闭

    Args:
        0 (Default): 1

    Example:
        >>> &transport
        >>>     do_transport_run_MC_temp      = 0      
    """

def do_transport_get_mob_MC_temp(int):
    """计算不同温度下的迁移率, 1表示开启，0表示关闭

    Args:
        0 (Default): 1

    Example:
        >>> &transport
        >>>     do_transport_get_mob_MC_temp      = 0      
    """

def do_transport_run_ME(int):
    """ME 方法模拟, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_run_ME      = 0   
    """


def do_transport_get_mob_ME(int):
    """计算迁移率, 1表示开启，0表示关闭

    Args:
        0 (Default): 1

    Example:
        >>> &transport
        >>>     do_transport_get_mob_ME      = 0     
    """


def do_transport_run_ME_temp(int):
    """不同温度下的 ME 模拟, 1表示开启，0表示关闭

    Args:
        0 (Default): 1

    Example:
        >>> &transport
        >>>     do_transport_run_ME_temp      = 0     
    """

def do_transport_get_mob_ME_temp(int):
    """计算不同温度下的迁移率, 1表示开启，0表示关闭

    Args:
        0 (Default): 1

    Example:
        >>> &transport
        >>>     do_transport_get_mob_ME_temp      = 0      
    """

def do_transport_gather_momap_data(int):
    """收集计算的相关数据, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     do_transport_gather_momap_data      = 1      
    """

def crystal(str):
    """定义晶体文件

    Args:
        naphthalene.cif (Default)

    Example:
        >>> &transport
        >>>     crystal = naphthalene.cif      
    """

def molecule(str):
    """定义晶体文件, default  2 mol1.mol mol2.mol


    Example:
        >>> &transport
        >>>     mol = 2 mol1.mol mol2.mol     
    """


def lat_cutoff(float):
    """计算相邻转移积分的截断半径(单位:Å)，这意味着如果两个分子的最近原子距离小于 lat_cutoff，则考虑计算它们之间的转移积分。

    Args:
        4 (Default)

    Example:
        >>> &transport
        >>>     lat_cutoff = 4     
    """

def neighbor_scell(float):
    """for neighbor cell search, default 3 3 3


    Example:
        >>> &transport
        >>>     neighbor_scell = 3 3 3
    """

def super_cell(float):
    """for ME calculations, default 5 5 5


    Example:
        >>> &transport
        >>>     super_cell = 5 5 5
    """


def ratetype(str):
    """定义电子空穴迁移速率计算方法，可选经典marcus方法或者量子修正的quantum方法

    Args:
        marcus (Default)：quantum

    Example:
        >>> &transport
        >>>     ratetype = quantum     
    """


def chargetype(str):
    """计算电子或者空穴，可选 e: electron, h: hole or a: all

    Args:
        a (Default)

    Example:
        >>> &transport
        >>>     chargetype = e     
    """



def temp(float):
    """定义模拟温度

    Args:
        300 (Default)

    Example:
        >>> &transport
        >>>     temp = 300     
    """


def start_temp(float):
    """计算不同温度下的电荷迁移率时，定义模拟初始温度

    Args:
        200 (Default)

    Example:
        >>> &transport
        >>>     start_temp = 200     
    """

def end_temp(float):
    """计算不同温度下的电荷迁移率时，定义模拟最终温度

    Args:
        300 (Default)

    Example:
        >>> &transport
        >>>     end_temp = 300     
    """

def delta_temp(float):
    """定义模拟温度间隔，例如，若 Start_Temp，End_Temp 和 delta_Temp 分别为 200，300，50，那么将进行 200K，250K，300 K 下的蒙特卡罗模拟

    Args:
        50 (Default)

    Example:
        >>> &transport
        >>>     delta_temp = 50     
    """

def nsimu(int):
    """定义模拟次数

    Args:
        2000 (Default)

    Example:
        >>> &transport
        >>>     nsimu = 2000
    """

def tsimu(int):
    """定义总模拟时间（in ns）

    Args:
        1000 (Default)

    Example:
        >>> &transport
        >>>     tsimu = 1000
    """

def tsnap(int):
    """定义记录输出文件中的载流子位置的时间间隔

    Args:
        5 (Default)

    Example:
        >>> &transport
        >>>     tsnap = 5
    """


def bond_dis_scale(float):
    """help separating molecules in cif file, users may tune this parameter for abnormal cases.

    Args:
        1.15 (Default)

    Example:
        >>> &transport
        >>>     bond_dis_scale = 1.15     
    """

def HL_unique_ctrl_ratio(float):
    """used for judging unique dipoles, |eigval[i] – eigval[j]| / max(|eigval[i]|, |eigval[j]|) < 0.05

    Args:
        0.05 (Default)

    Example:
        >>> &transport
        >>>     HL_unique_ctrl_ratio = 0.05     
    """


def RE_use_neutral_chk(int):
    """Calculate anion and cation state reorganization energies by using neutral state chk file, can be 0 or 1, default to 0.

    Args:
        0 (Default)

    Example:
        >>> &transport
        >>>     RE_use_neutral_chk = 1    
    """


def RE_calc_lambda_4P(int):
    """If calculate reorganization energies by using the Nelson four- point method, can be 0 or 1, default to 1.

    Args:
        1 (Default)

    Example:
        >>> &transport
        >>>     RE_calc_lambda_4P = 1  
    """

def Thinfilm(int):
    """Format: thinfilm = dir nuc, here dir can be 0(vector_a), 1(vector_b), and 2(vector_c), 
    nuc is the number of repeating unit cell in dir direction, e.g., thinfilm = 0 2


    Example:
        >>> &transport
        >>>     Thinfilm =  0 2   
    """

def V_dynamic_disorder(int):
    """Default to 0, used to control if we need to calculate dynamic disorder of transfer integrals, 
    the per-molecular files VH*- dyn.dat or VL*-dyn.dat are to be provided under data directory.

    Args:
        0 (Default)

    Example:
        >>> &transport
        >>>     V_dynamic_disorder = 0
    """

def HOMOLUMO_dynamic_disorder(int):
    """Default to 0, used to control if we need to calculate dynamic disorder of HOMO/LUMO, files HOMO-dyn.dat or
       LUMO-dyn.dat are to be provided under data directory.

    Args:
        0 (Default)

    Example:
        >>> &transport
        >>>     HOMOLUMO_dynamic_disorder = 0   
    """