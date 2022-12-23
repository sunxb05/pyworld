"""Documentation about the MOMAP module."""

def compute_engine(int):
    """定义使用的计算引擎，1 = Gaussian, 2 = ORCA, 3 = QCHEM， 4 = BDF

    Args:
        1 (Default): 0

    Example:
        >>> &transport
        >>>     compute_engine      = 1      
    """

def module_momap(str):
    """MOMAP程序

    Args:
        momap (Default)

    Example:
        >>> &transport
        >>>     module_momap =  momap/2022B-openmpi      
    """

def module_qc(str):
    """Use module to setup QC running environment

    Args:
        empty (Default)

    Example:
        >>> &transport
        >>>     module_qc = gaussian/g09.e01      
    """

def qc_exe(str):
    """计算引擎可执行程序，g09/g16 or full path/orca or qchem or BDF

    Args:
        g09 (Default)

    Example:
        >>> &transport
        >>>     qc_exe = g09      
    """

def qc_method(str):
    """计算引擎所用方法

    Args:
        b3lyp (Default)

    Example:
        >>> &transport
        >>>     qc_method = b3lyp      
    """


def qc_basis(str):
    """计算引擎所用基组

    Args:
        b3lyp cc-pvdz (Default)

    Example:
        >>> &transport
        >>>     qc_basis = b3lyp cc-pvdz      
    """

def qc_basis_re(str):
    """计算重组能所用基组,无设置的话与qc_basis相同

    Args:
        b3lyp cc-pvdz (Default)

    Example:
        >>> &transport
        >>>     qc_basis_re = b3lyp cc-pvdz      
    """

def qc_memory(int):
    """所用内存（in MB）

    Args:
        4096 (Default)

    Example:
        >>> &transport
        >>>     qc_memory = 4096      
    """

def qc_nodes(int):
    """申请使用节点数

    Args:
        1 (Default)

    Example:
        >>> &transport
        >>>     qc_nodes = 1      
    """

def qc_ppn(int):
    """每节点并行运行核数

    Args:
        8 (Default)

    Example:
        >>> &transport
        >>>     qc_ppn = 20      
    """

def queue_name(str):
    """计算任务提交队列名

    Args:
        workq (Default)

    Example:
        >>> &transport
        >>>     queue_name = workq      
    """

def sched_type(str):
    """作业管理系统类别，pbs, slurm, lsf, or local

    Args:
        local (Default)

    Example:
        >>> &transport
        >>>     sched_type = local      
    """

def job_slot_limit(int):
    """number of jobs at the same time, default to 0, means no limit

    Args:
        0 (Default)

    Example:
        >>> &transport
        >>>     job_slot_limit = 0   
    """

def walltime_hour(int):
    """setting walltime in job submission scripts

    Args:
        1000 (Default)

    Example:
        >>> &transport
        >>>     walltime_hour = 1000   
    """ 

def scratch(str):
    """/tmp, only valid for ORCA calculation.

    Args:
        /tmp (Default)

    Example:
        >>> &transport
        >>>     scratch = /tmp      
    """
