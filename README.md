
		 Mixed Element Applications with Nodal Discontinuous Galerkin (MEANDG) 

* Copyright (c) 2018, Subodh Joshi, Tapan Mankodi, Shivasubramanian Gopalakrishnan. All rights reserved.
* MEANDG is licensed under the BSD-3 Clause license. See LICENSE for details. 

* Eigen library (see ./include/Eigen/) is licensed under MPL2 (see ./include/Eigen/COPYING.README). 
* Use preprocessor flag -DEIGEN_MPL2_ONLY to ensure that only MPL2 (and possibly BSD) licensed code from Eigen is compiled. 
* Refer (http://eigen.tuxfamily.org/index.php?title=Main_Page///License) for details. 
* No files of the Eigen library are modified for use with MEANDG. */



File structure:

Directories:
1. src: contains all the source (.cpp) files as well as SConscript. Add the new cpp files to the Sconscript. 
2. obj: contains all the compiled files (.o). Can be added to .gitignore as it will be freshly compiled on each system.
3. include: contains the header (.h) files
4. app: contains application directories. Each directory should in turn contain:
            Time folder (such as 0): Contains variable (p, U, etc.) data in OpenFOAM format
            constant folder        : should contain polymesh dir with geometry data inside (in OpenFOAM format)
            system:                  Standard OpenFOAM files. 
5. tests: contains output results of all the tests (tests are described in tests.cpp inside src dir). 


Files:
1. SConstruct:    Scons script. 
2. out.exe:       executable file generated after compiling the code


To run the code:
scons
./out.exe

Note:
1. Each dir contains a READE file for specifics of that dir. Add more information as required.
2. Template functions need to be completely described in the header files. 
3. Add tests to the tests.cpp. Describe clearly what the tests are doing, what arguments are required and what the user is expected to do.
4.  Add comments whenever appropriate.
5. Install scons on your system to compile the code. On Ubuntu 16.04, this translates to:
   $sudo apt-get install scons */


Resources:
1. Git tutorial and help: http://www.vogella.com/tutorials/Git/article.html 
2. Book: Nodal discontinuous Galerkin methods: Algorithms, Analysis and Applications, J. Hesthaven and T. Warburton.
3. dealII library: http://www.dealii.org

