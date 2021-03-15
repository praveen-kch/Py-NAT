'''
************************************************************************************
    File            : LAS_GZ.py
    
    File Purpose    : Define the Class "LAS_GZ"
    
    Class Purpose   : 
                    -> Read the Load Case Data from the File set_LC.dat
                    -> Interpolate KN values at various heel angles for the given displacement
                    -> Compute GZ values
                    -> write GZ data to file
                    -> Plot GZ data in curves
    
    App             : PyNAT - Python Based Naval Architecture Tools
    Version         : 2021.0
    Paradigm        : Object Based Framework
    Author          : Praveen Kumar Ch
***************************************************************************************
'''
import numpy as np
from scipy import interpolate 
import matplotlib.pyplot as plt
from HullGeometry import HullGeometry
import math as math
import os
# Py-NAT 2021a : PRAVEEN KUMAR CH
class LAS_GZ:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def __init__(self,fs,GivenShip,knData,ep):
        self.readFromFile(fs)
        self.compute(GivenShip,knData,ep)
        self.writeToFile(fs,knData)

    # --------------------------------------------------------------------------------------
    # Function to Read Load Case Data from File
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def readFromFile(self,fs):
        
        print("\t\tReading Load Condition from file...")
        # Reading the Solver properties
        fp=open(fs.ilc_path,"rt")
        lines=fp.readlines()
        fp.close()

        self.nCases=len(lines)
        self.Mass=list()
        self.KG=list()

        for i in range(self.nCases):
            nums=lines[i].split("\t")
            self.Mass.append(float(nums[0].rstrip("\n")))
            self.KG.append(float(nums[1].rstrip("\n")))
    # --------------------------------------------------------------------------------------
    

    # --------------------------------------------------------------------------------------
    # Function to Compute GZ Data
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    
    def compute(self,GivenShip,knData,ep):

        print("\n\t\tComputing GZ Data")
        
        # Interpolate KN at each Heel
        self.KN=np.ones((self.nCases,knData.nHeel),dtype=np.float64)
        self.GZ=np.ones((self.nCases,knData.nHeel),dtype=np.float64)
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        for iha in range(knData.nHeel):
            x=np.copy(knData.MASS[iha,:])
            y=np.copy(knData.KN[iha,:])
            ind=np.argsort(x,axis=0)
            x=np.take_along_axis(x,ind,axis=0)
            y=np.take_along_axis(y,ind,axis=0)
            kn_fMass=interpolate.splrep(x,y,k=3)

            for ics in range(self.nCases):
                self.KN[ics][iha]=interpolate.splev(self.Mass[ics],kn_fMass,der=0)
                phi=abs(knData.rangeHeels[iha])
                self.GZ[ics][iha]=self.KN[ics][iha]-self.KG[ics]*math.sin(phi)
                #print("\n\t\t\tKN=",self.KN[ics][iha])
                print("\n\t\t\tKG= %10.4f    Disp=%10.4f    Heel=%6.4f    GZ=%10.4f"%(self.KG[ics],self.Mass[ics],knData.rangeHeels[iha],self.GZ[ics][iha]))
    # --------------------------------------------------------------------------------------
    

    # --------------------------------------------------------------------------------------
    # Function to Write GZ Data to File
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def writeToFile(self,fs,knData):

        print("\n\t\tPyNAT: Writing GZ Data to File:\t...\\",fs.ogz_path)

        if not os.path.exists(fs.las_dir):
            os.makedirs(fs.las_dir)
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        if(knData.techq=="DIRECT"):
            filepath=fs.ogz_path
            plotpath=os.path.join(fs.las_dir,"LAS_GZ.png")
        else:
            filepath=fs.ogz_fast_path
            plotpath=os.path.join(fs.las_dir,"LAS_GZ_FAST.png")

        fp=open(filepath,'wt')
        fp.write("\nPyNAT: Righting Lever (GZ) Curve Data")

        for ics in range(self.nCases):
            fp.write("\nMass Displacement=%f"%(self.Mass[ics]))
            fp.write("\nVCG (KG) = %f"%(self.KG[ics]))
            fp.write('\nHeel Angle\tGZ        \t')

            for i in range(knData.nHeel):
                fp.write("\n%10.3E\t%10.3E\t"%(knData.rangeHeels[i]*180/math.pi,self.GZ[ics][i]))

            plt.figure("GZ Curve")
            hst=format(ics,"2d")
            plt.plot(knData.rangeHeels*180/math.pi,self.GZ[ics,:],label="Case No"+hst)
            plt.title("Righting Lever (GZ) Curve")
            plt.xlabel("Heel Angle")
            plt.ylabel("GZ")
            lg=plt.legend(bbox_to_anchor=(1.05, 1))
            plt.savefig(plotpath,            
                dpi=300, 
                format='png', 
                bbox_extra_artists=(lg,), 
                bbox_inches='tight')

        fp.close()
        plt.close('all')
        # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------