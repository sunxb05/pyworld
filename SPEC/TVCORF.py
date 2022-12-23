"""A typical **spec_tvcf** control with momap.inp can be as follows:

.. code-block:: bash


    do_spec_tvcf_ft   = 1
    do_spec_tvcf_spec = 1    

    &spec_tvcf      
        DUSHIN        = .t.
        Temp          = 300 K
        tmax          = 5000 fs
        dt            = 0.001 fs
        Ead           = 0.075091878 au
        EDMA          = 0.92694 debye
        EDME          = 0.64751 debye 
        FreqScale     = 1.0
        DSFile        = "evc.cart.dat"
        Emax          = 0.3 au
        dE            = 0.00001 au
        logFile       = "spec.tvcf.log" 
        FtFile        = "spec.tvcf.ft.dat" 
        FoFile        = "spec.tvcf.fo.dat" 
        FoSFile       = "spec.tvcf.spec.dat"  

    /


A typical **ic_tvcf** control with momap.inp can be as follows:


.. code-block:: bash

    do_ic_tvcf_ft     =1  
    do_ic_tvcf_spec   =1  

    &ic_tvcf    
        DUSHIN        = .t.
        Temp          = 300 K
        tmax          = 5000 fs
        dt            = 0.001 fs
        Ead           = 0.075091878 au
        DSFile        = "evc.cart.dat"
        CoulFile      = "evc.cart.nac"
        Emax          = 0.3 au
        logFile       = "ic.tvcf.log"
        FtFile        = "ic.tvcf.ft.dat"
        FoSFile       = "ic.tvcf.fo.dat"  

    /



A typical **isc_tvcf** control with momap.inp can be as follows:


.. code-block:: bash

    do_isc_tvcf_ft    = 1
    do_isc_tvcf_spec  = 1    

    &isc_tvcf   
        DUSHIN        = .t.
        Temp          = 298 K
        tmax          = 5000 fs
        dt            = 0.001 fs
        Ead           = 0.0941829 au
        Hso           = 176.08 cm-1 
        DSFile        = "evc.cart.dat" 
        Emax          = 0.3 au
        logFile       = "isc.tvcf.log" 
        FtFile        = "isc.tvcf.ft.dat" 
        FoFile        = "isc.tvcf.fo.dat"
    /


"""



def do_spec_tvcf_ft(int):
    """是否开启计算辐射关联函数, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_spec_tvcf_ft      = 1      
    """

def do_spec_tvcf_spec(int):
    """是否开启计算辐射光谱, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_spec_tvcf_spec      = 1      
    """

def do_ic_tvcf_ft(int):
    """是否开启计算内转换关联函数, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_ic_tvcf_ft      = 1      
    """

def do_ic_tvcf_spec(int):
    """是否开启计算内转换光谱, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_ic_tvcf_spec      = 1      
    """


def do_isc_tvcf_ft(int):
    """是否开启计算系间窜越关联函数, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_isc_tvcf_ft      = 1      
    """

def do_isc_tvcf_spec(int):
    """是否开启计算系间窜越光谱, 1表示开启，0表示关闭

    Args:
        1 (Default): 0

    Example:
        >>>     do_isc_tvcf_spec      = 1      
    """


def DUSHIN(logic):
    """是否开启 Duschinsky t表示开启，f表示关闭

    Args:
        .t. (Default): .f. 

    Example:
        >>> &spec_tvcf
        >>>     DUSHIN      = .t.        
    """

def HERZ(logic):
    """是否开启 Herzberg-Teller 效应的计算

    Args:
        .f. (Default): .t. 

    Example:
        >>> &spec_tvcf
        >>>     HERZ      = .f.        
    """

def Ead(float):
    """ 绝热激发能，可由Gaussian等软件计算得到

    Args:
        0.097374 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     Ead      = 0.097374 au   
    """

def EDMA(float):
    """ 吸收跃迁偶极矩，可由Gaussian等软件计算得到

    Args:
        0.7887 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     EDMA      = 0.7887 debye   
    """

def EDME(float):
    """ 发射跃迁偶极矩，可由Gaussian等软件计算得到

    Args:
        1.0368  (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     EDME      = 1.0368  debye   
    """

def Hso(float):
    """ 旋轨耦合常数，可由Dalton等软件计算得到，用于系间窜越计算

    Args:
        10 (Default)
 

    Example:
        >>> &isc_tvcf
        >>>     Hso      = 116.877376 cm-1   
    """


def Temp(float):
    """ 定义温度

    Args:
        10 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     Temp      = 300 K    
    """

def tmax(float):
    """ 定义积分时间

    Args:
        5000 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     tmax      = 3000 fs    
    """

def dt(float):
    """ 定义积分步长

    Args:
        0.001 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     dt      = 0.001 fs    
    """


def Emax(float):
    """ 定义光谱频率范围上限

    Args:
        0.3 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     Emax      = 0.3 au  
    """

def dE(float):
    """ 定义输出能量间隔

    Args:
        0.00001 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     dE      = 0.00001 au   
    """

def isgauss(logic):
    """是否开启 spectrum broadening， t表示开启，f表示关闭

    Args:
        .f. (Default): .t. 

    Example:
        >>> &spec_tvcf
        >>>     isgauss      = .f.        
    """


def FWHM(float):
    """ Broadening width

    Args:
        0.001 (Default)
 

    Example:
        >>> &spec_tvcf
        >>>     FWHM      = 0.001 au   
    """


def DSFile(str):
    """定义读取的 evc 文件名

    Args:
        "evc.cart.dat" (Default)

    Example:
        >>> &spec_tvcf
        >>>     DSFile = "evc.cart.dat"      
    """

def CoulFile(str):
    """定义读取的nacme文件名，

    Args:
        "evc.cart.nac" (Default)

    Example:
        >>> &ic_tvcf
        >>>     CoulFile = "evc.cart.nac"      
    """

def DDplFile(str):
    """定义读取的 dip 文件名,用于Herzberg-Teller效应计算

    Args:
        "evc.cart.dip" (Default)

    Example:
        >>> &spec_tvcf
        >>>     evc.cart.dip = "evc.cart.dip"      
    """

def logFile(str):
    """定义输出 log 文件名

    Args:
        "spec.tvcf.log" (Default)

    Example:
        >>> &spec_tvcf
        >>>     logFile = "spec.tvcf.log"      
    """


def FtFile(str):
    """定义输出的关联函数文件名

    Args:
        "spec.tvcf.ft.dat" (Default)

    Example:
        >>> &spec_tvcf
        >>>     FtFile = "spec.tvcf.ft.dat"      
    """


def FoFile(str):
    """定义输出的谱函数输出文件名

    Args:
        "spec.tvcf.fo.dat" (Default)

    Example:
        >>> &spec_tvcf
        >>>     FoFile = "spec.tvcf.fo.dat"      
    """


def FoSFile(str):
    """定义输出的归一化的光谱输出文件名

    Args:
        "spec.tvcf.spec.dat" (Default)

    Example:
        >>> &spec_tvcf
        >>>     FoSFile = "spec.tvcf.spec.dat"      
    """

def GFile(str):
    """Output file for broadening info

    Args:
        "spec.tvcf.gauss.dat" (Default)

    Example:
        >>> &spec_tvcf
        >>>     GFile = "spec.tvcf.gauss.dat"      
    """

