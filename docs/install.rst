软件安装
*******

Download and Installation 
=========================


First download the MOMAP package from MOMAP website_, by simply clicking the related link or DOWNLOAD button to download the file.


Install MOMAP for Linux
-----------------------

The MOMAP package for Linux is downloaded as a single zipped installable run file, e.g.,

    ``momap-2022B-linux-el7-openmpi.run.gz``

Unzip the file and add executable attribute to the file, and install MOMAP by running:

	``gunzip momap-2022B-linux-el7-openmpi.run.gz``	

	``chmod a+x momap-2022B-linux-el-openmpi.run``	

	``./momap-2022B-linux-el7-openmpi.run``

The default installation folder is ``$HOME/MOMAP-2022B``, however, one can change the target installation folder in the installation process.

Before using MOMAP, one should first set up the running environment by adding the following line to ``~/.bashrc`` if Bash is used

	``. <installed_momap_folder>/env.sh``

Log out and log back in again for it to take effect, then we can proceed to run MOMAP.




Install MOMAP for macOS
-----------------------
For the macOS case, the installation process is similar to that in the Linux case, however, we need to install the necessary libraries used in the MOMAP package, and here the Homebrew package manager for macOS is used. The installation procedure is as follows:

1. Install Homebrew if not installed already

	``/usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"``

2. Install GCC, Open-mpi, FFTW3 and Lapack

	``brew install gcc@12``	

	``brew install open-mpi``	

	``brew install fftw``	

	``brew install lapack``

Other settings are similar to those in the Linux situation. Once the above installations are finished, we can use the MOMAP for macOS as in the Linux’s case.




Install MOMAP for Windows
-------------------------
For the Windows case, MSYS2 can be used as a terminal to install MOMAP. The installation procedure is as follows:

1. To install MOMAP for Windows, we need first to install MSYS2. Download MSYS2 at the link: ``https://www.msys2.org``, and install MSYS2 to location for example C:\\msys64. Once the installation finishes, launch MSYS2.

2. First, update the package database and core system package. Then use pacman to install Python and vim etc.:

	``pacman -Syu``

	``pacman -S python``

	``pacman -S vim``

	``pacman -S openssh``

3. Download Windows version of MPI program MSMPI at the link: ``https://docs.microsoft.com/en-us/message-passing-interface/microsoft-mpi``, and install to for instance C:\\msys64\\msmpi.

4. Once the above MOMAP running environment settings are done, we can install the MOMAP for Windows.
For example, put the downloaded file to D:\\, and extract the archived files to the directory D:\\MOMAP-2022B. Check the file env-mingw.sh to see if the environment variable MOMAP_ROOT points to the installed directory.

Edit ``~/.bashrc`` and add the following line:

	``. /d/MOMAP-2022B/env-mingw.sh``



Install MOMAP using Environment Modules
---------------------------------------


If the Environment Modules is installed, one can use the command module to manage the MOMAP running environment. 

1. Suppose the Environment Modules is installed to ``/opt/Modules``, and the module files are put under ``/opt/Modules/modulefiles``, one can create a directory and copy MOMAP to that newly created directory:

	``mkdir -p /opt/Modules/modulefiles/momap``

2. Make it available to Environment Modules by using the following command:

	``module use ~/.modulefiles``

3. Now one can use the module files as usual. For example, we can load the MOMAP environment by simply running the command:

	``module load momap/2022B-openmpi``


Licensing and Test
==================

Before we can use MOMAP package, we need a license. There exist four license options, that is, Trial PC version, All-in-one PC version, Intranet cluster version and Super HPC cluster version.

	+ Trial PC version: Choose this option if you’d like to try the software with a desktop/laptop PC.	


	+ All-in-one PC version: Choose this option if you plan to use the software in a workstation.	


	+ Intranet cluster version: Choose this option if you plan to use the software in a small group-wise intranet computing cluster.	


	+ Super HPC cluster version: Choose this option if you plan to use the software in the public domain super HPC cluster.


The program to collect the license data is get_LicenseNumber.exe, the program is called automatically at the end of MOMAP installation.
For the MOMAP for Windows case, we need to run get_LicenseNumber.exe, that is,

	``$MOMAP_ROOT/bin/get_LicenseNumber.exe``

The generated license data file is located at directory ``$MOMAP_ROOT/license/``, named as LicenseNumber.txt. 
One should send this file to HZW Co. Ltd.
Later on, a MOMAP license file, hzwtech.lic, will be sent to you by a sales representative from HZW Co. Ltd., you can simply copy the license file to ``$MOMAP_ROOT/license`` directory.


To verify that the MOMAP package has been properly installed, and the license is correctly configured and installed, users can run a short test to verify the installation.



Troubleshooting
===============

1. In some supercomputing centers, the SSH port may not be the default 22, in that case, we need to setup the SSH environment variable, for example:
``export MOMAP_SSH_PORT=5577``

2. If MOMAP is to be run under the Ubuntu Linux system, before we start to install MOMAP, we need first to promote the user rights and make the user to be an administrator.



.. _website: http://www.momap.net.cn/index.php/download