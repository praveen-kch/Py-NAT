# Py-NAT

Details of Application:
-----------------------------------
Name 		: PyNAT.exe
Size 		: 63.5 MB
Language 	: Python
Libraries used  : NumPy, Matlibplot and SciPy
Packaged using 	: Pyinstaller
Paradigm 	: Object Oriented Model
Platform 	: Developed for Windows 10, Needs re compilation of source code to work cross platfom, Stand alone application

PyNAT is a Python based Naval Architecture tool written to compute

 1. Hydro-statics at even keel and and at given trim angles
 2. Cross curves of stability (KN Curves)
 3. GZ curves at a given Load Cases (i.e., Mass and VCG)

The following hull forms can be analysed, 

 1. Only Mono Hulls with port / starboard symmetry
 2. Hull form must be definable using station curves along the length
 2. Hull forms with chines and knuckles
 3. Hull forms with discontinuous curves (near bow and stern)

No installation is required to use this CUI application. 
The application comes as an executable program, that only reads, computes and writes data to its parent directory.

A preliminary validation showed reasonable accuracy. A more in depth validation is to be performed.

The Script is written in object oriented paradigm.
One can experiment with the script and improve the speed and accuracy. 
Typical parameters that one can experiment with,
 1. The fit type / interpolation techniques
 2. Order of curve fitting / interpolation
 3. Spline smoothing
 3. Integration techniques
 4. Extending to perform advanced calculations
 5. Consider advanced hull forms like catamarans

Link: https://github.com/praveen-kch/Py-NAT

The GitHuB Repository contains the following :
 1. The Python Source Scripts
 2. The Application "PyNAT.exe"
 3. User Guide
 4. Demo video links
 5. Test Cases for KVLCC2 and CAD/Excel templates to help prepare the input files


-Praveen Kumar Ch
praveench1888@gmail.com
