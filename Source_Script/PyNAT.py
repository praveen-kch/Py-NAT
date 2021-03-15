'''
****************************************************************
    File        : main.py

    Purpose     : The Base Scipt from where the execution begins
                  Class Objects are Created to perform analysis

    App         : PyNAT - Python Based Naval Architecture Tools 
    Version     : 2021a
    Paradigm    : Object Oriented Model
    Author      : Praveen Kumar Ch
*****************************************************************
'''

from HullGeometry import HullGeometry
from HydroStatics import HydroStatics
from FileSystem import FileSystem
from LAS_KN_1A import LAS_KN_1A
from LAS_KN_0B import LAS_KN_0B
from LAS_GZ import LAS_GZ
from ENV_PROP import ENV_PROP
import matplotlib.pyplot as plt 
import os
import msvcrt as m
import time

#-----------------------------------------------------------------------------
# Creation of 'FileSystem' Object
#       Object Contains The File Names and paths to be utilized by program
#       See "FileSystem.py"
#-----------------------------------------------------------------------------
t0=time.time()
print("\n\n*********************************************************************")
print("                        Py-NAT 2021a")
print("*********************************************************************")
fs=FileSystem()
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Creation of 'ENV_PROP' Object
#       Contains Environment properties like water density
#-----------------------------------------------------------------------------
print("\n\nThe Environemnt Properties are being read...\t...\t...\t...\t...")
ep=ENV_PROP(fs)
print("\nEnvironment Properties File Read")
#-----------------------------------------------------------------------------

# Creation of  'HullGeometry' Object
#    Object Contains Data defining the ship hull geometry 
#    object Contains Script to read and prepare Hull Geometry
#           Contians Script to fit Splines
#    See the Script file "HullGeometry.py"
#-----------------------------------------------------------------------------
print("\n\nThe Hull Geometry is being read and prepared...\t...\t...\t...\t...")
givenShip=HullGeometry(fs)
tend=time.time()
print("\nThe Hull Geometry Data Preperation Complete.")
print("\nCheck the Prepared Hull Geometry Data in Directory:\t...\\",fs.sec_dir)
print("\nTotal Time Elapsed for Hull Geometry Preperations :\t",tend-t0)
print("\nClose Body Plan Plot to Continue.")
plt.show()
#----------------------------------------------------------------------------------
iLoops=0
maxLoops=50
while(iLoops<maxLoops):
    # Loop Main Analysis
    print("\n\n-----------------------------------------------------------------------")
    print("\nAnalysis Loop %d of %d"%(iLoops,maxLoops))
    print("\nChoose the Analysis to be Performed:\n")
    print("\n\tPress  '1'  for Hydrostatics")
    print("\n\tPress  '2'  for Large Angle Stability (KN Curves +GZ : Using Fast Technique) ")
    print("\n\tPress  '3'  for Large Angle Stability (KN Curves +GZ : Using Direct Technique) ")
    print("\n\tPress  '4'  for Performing All Calculations")
    print("\n\tPress  '0' or 'E' or 'e' or 'esc' to Exit")
    char=m.getch()
    print("\n\tPressed Key =",char)
    #-----------------------------------------------------------------------------
    # Creation of  'HydroStatics' Object
    #    Contains Data related to Hydrostatic Calculations
    #    Contains Code to perform hydrostatic calculations
    #    See the Script file "HydroStatics.py"
    #-----------------------------------------------------------------------------
    if char==b'1' or char ==b'4':
        t0=time.time()
        print("\nComputing Hydrostatics... \t...\t...\t...\t...")
        hydrostatics = HydroStatics(fs,givenShip,ep)
        print("\nHydrostatic Analysis Complete")
        print("\nCheck Results in Directory:\t...\\",fs.hs_dir)
        tend=time.time()
        print("\nTotal Time Elapsed for Hydrostatic Calculations:\t",tend-t0)
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Creation of  'LAS_KN' Object
    #    Contains Data related to KN curves
    #    Contains Code to perform KN calculations
    #    See the Script file "LAS_KN.py"
    #-----------------------------------------------------------------------------
    if char==b'2'  or char==b'3' or char==b'4':
        t0=time.time()
        print("\nComputing KN Curves Data...\t...\t...\t...\t...")
        if char==b'2':
            knData = LAS_KN_1A(fs,givenShip,ep)
        else:
            knData = LAS_KN_0B(fs,givenShip,ep)
        print("\nCross Curves of Stability Calculation Complete")
        print("\nCheck Results in Directory:\t...\\",fs.las_dir)
        tend=time.time()
        print("\nTotal Time Elapsed for KN Curve Calculations:\t",tend-t0)
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Creation of  'LAS_GZ' Object
    #    Contains Data related to KN curves
    #    Contains Code to perform KN calculations
    #    See the Script file "LAS_GZ.py"
    #-----------------------------------------------------------------------------
    if char==b'2' or char==b'3' or char==b'4':
        t0=time.time()
        print("\nComputing GZ Curves Data...\t...\t...\t...\t...")
        gzData = LAS_GZ(fs,givenShip,knData,ep)
        print("\nGZ Curve Calculation Complete")
        print("\nCheck Results in Directory:\t...\\",fs.las_dir)
        tend=time.time()
        print("\nTotal Time Elapsed for GZ Curve Calculations=\t",tend-t0)
    #-----------------------------------------------------------------------------

    #-----------------------------------------------------------------------------
    # Prompt User to Close the terminal
    #-----------------------------------------------------------------------------
    if char==b'0' or char==27 or char==b'E' or char==b'e' or char==b'\x1b':
        exit()

    #-----------------------------------------------------------------------------
    # Counting Loops
    #-----------------------------------------------------------------------------
    iLoops+=1
