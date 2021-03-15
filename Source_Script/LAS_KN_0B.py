'''
************************************************************************************
    File            : LAS_KN.py
    
    File Purpose    : Define the Class "LAS_KN"
                        with Direct Integration technique to compute heeled section prop.
    
    Class Purpose   : 
                    -> Compute and store KN values over a gien range of heel angles for each given draft
                    -> Function to Compute Intersection of Water Line with Section Curves
                        at various floating conditions
                    -> Computations are done for Floating Conditions Defined by 
                        -> Draft 
                        -> Heel
                    -> The range of Drafts and Heel are read from the File "Set_LAS_KN.dat"
                    -> computations are performed for each draft and Heel combination
                    -> Functions to write the Large angle Stability Results to File "Res_LAS_KN.dat"
                    -> Plots and Saves KN curves and  as image files of format ".png"
    
    App             : PyNAT - Python Based Naval Architecture Tools
    Version         : 2021.0
    Paradigm        : Object Based Framework
    Author          : Praveen Kumar Ch
***************************************************************************************
'''
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import quad
import matplotlib.pyplot as plt
from HullGeometry import HullGeometry
import math as math
import os
# Py-NAT 2021a :  PRAVEEN KUMAR CH
class LAS_KN_0B:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # Py-NAT 2021a :  PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def __init__(self,fs,GivenShip,ep):
        self.techq="DIRECT"
        self.readFromFile(fs)
        self.compute(GivenShip,ep)
        self.writeToFile(fs)
    # --------------------------------------------------------------------------------------
    

    # --------------------------------------------------------------------------------------
    # The Method to Read the Data from File
    # Py-NAT 2021a :  PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def readFromFile(self,fs):

        print("\n\t\tPyNAT: Reading File:\t...\\",fs.ikn_path)

        # Reading the Solver properties
        fp=open(fs.ikn_path,"rt")
        lines=fp.readlines()
        fp.close()
       
        nums=lines[0].split(",")
        self.nHeel=len(nums)
        self.rangeHeels=np.array(range(self.nHeel),dtype=np.float64)
        for i in range(self.nHeel):
            self.rangeHeels[i] =float(nums[i])*math.pi/180

    # --------------------------------------------------------------------------------------
    

    # --------------------------------------------------------------------------------------
    # The Method to Compute KN Data
    # Py-NAT 2021a :  PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    
    def compute(self,GivenShip,ep):

        print("\n\t\tPyNAT: Computing KN Curves Data")

        # Finding the CL-WL intercept range for each heel angle           
        self.nDraft=11
        self.rangeDrafts=np.zeros((self.nHeel,self.nDraft))

        for iha in range(self.nHeel):

            if(abs(self.rangeHeels[iha])%(math.pi/2)<=0.0001 and abs(self.rangeHeels[iha])>=0.0001):
                self.rangeHeels[iha]-=0.017453292

            # Finding the Limits of Draft at CL
            dclmin=0
            dclmax=0
            for ist in range(GivenShip.nst):
                for ipt in range(GivenShip.npts[ist]):
                    dcl=abs(GivenShip.y[ist][ipt]*math.tan(self.rangeHeels[iha]))
                    db=-dcl+GivenShip.z[ist][ipt]
                    da=dcl+GivenShip.z[ist][ipt]
                    if(db<dclmin):
                        dclmin=db
                    if(da>dclmax):
                        dclmax=da

            delDraft=(dclmax-dclmin)/(self.nDraft-1)
            dclmin=dclmin+delDraft/2
            dclmax=dclmax-delDraft/2

            # creating a range of drafts
            self.rangeDrafts[iha,:]=np.linspace(dclmin,dclmax,self.nDraft,endpoint=True)

        #print(self.rangeDrafts)

        self.VOL=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.MASS=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.LCB=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.TCB=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.VCB=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.WPA=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.WSA=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.KN=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.m=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.c=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.zN=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)
        self.yN=np.zeros([self.nHeel,self.nDraft],dtype=np.float64)

        # Compute the Water Line Intersections for each heel angle
        for iha in range(self.nHeel):

            # For each Draft 
            for idr in range(self.nDraft):

                zd=self.rangeDrafts[iha][idr]

                phi=self.rangeHeels[iha]

                BR=np.zeros(GivenShip.nst,dtype=np.float64)
                Ar=np.zeros(GivenShip.nst,dtype=np.float64)
                LM=np.zeros(GivenShip.nst,dtype=np.float64)
                TM=np.zeros(GivenShip.nst,dtype=np.float64)
                VM=np.zeros(GivenShip.nst,dtype=np.float64)
                CL=np.zeros(GivenShip.nst,dtype=np.float64)

                #print("------\t------\t------\t------\t------\t------\t------\t------\t------\t------\t")
                #print("\t\tZd=%12.6f\tphi=%12.6f"%(zd,phi*180/math.pi))
                # Py-NAT 2 021a :  PRAVEEN KUMAR CH
                #print("------\t------\t------\t------\t------\t------\t------\t------\t------\t------\t")

                for ist in range(GivenShip.nst):

                    if(phi<=0):
                        side=1
                    else:
                        side=-1
                    secData=self.compImmersedSectionProperties(GivenShip.secCrv[ist],zd,phi,side)               

                    BR[ist]+=secData[0]
                    Ar[ist]+=secData[1]
                    LM[ist]+=secData[2]
                    TM[ist]+=secData[3]
                    VM[ist]+=secData[4]
                    CL[ist]+=secData[5]
                 
                    #print("iSta=%5d, BRW=%10.5f, Ar=%10.5f, LM=%10.5f, TM=%10.5f, VM=%10.5f, CL=%10.5f"%(ist,BR[ist],Ar[ist],LM[ist],TM[ist],VM[ist],CL[ist]))


                # Integrating the Areas and Moments over the length
                Ar_fX=interp1d(GivenShip.x,Ar,kind='cubic')
                sq=quad(Ar_fX,GivenShip.x[0],GivenShip.x[-1])
                self.VOL[iha,idr]=sq[0]
                self.MASS[iha,idr]=self.VOL[iha,idr]*ep.waterDensity

                Sa_fX=interp1d(GivenShip.x,CL,kind='cubic')
                sq=quad(Sa_fX,GivenShip.x[0],GivenShip.x[-1])
                self.WSA[iha,idr]=sq[0]

                if self.VOL[iha,idr]!=0:
                    lm_fX=interp1d(GivenShip.x,LM,kind='cubic')
                    sq=quad(lm_fX,GivenShip.x[0],GivenShip.x[-1])
                    self.LCB[iha,idr]=sq[0]/self.VOL[iha,idr]

                    tm_fX=interp1d(GivenShip.x,TM,kind='cubic')
                    sq=quad(tm_fX,GivenShip.x[0],GivenShip.x[-1])
                    self.TCB[iha,idr]=sq[0]/self.VOL[iha,idr]

                    vm_fX=interp1d(GivenShip.x,VM,kind='cubic')
                    sq=quad(vm_fX,GivenShip.x[0],GivenShip.x[-1])
                    self.VCB[iha,idr]=sq[0]/self.VOL[iha,idr]

                else:
                    self.LCB[iha,idr]=0
                    self.TCB[iha,idr]=0
                    self.VCB[iha,idr]=0

                wpa_fX=interp1d(GivenShip.x,BR,kind='cubic')
                sq=quad(wpa_fX,GivenShip.x[0],GivenShip.x[-1])
                self.WPA[iha,idr]=sq[0]
                
                # Py-NAT 2021 a :  PRAVEEN KUMAR CH  
                # Finding KN
                self.m[iha,idr]=math.tan(math.pi/2-phi)
                self.c[iha,idr]=self.VCB[iha,idr]-self.m[iha,idr]*self.TCB[iha,idr]
                self.zN[iha,idr]=self.c[iha,idr]/(1+self.m[iha,idr]**2)
                self.yN[iha,idr]=-self.m[iha,idr]*self.zN[iha,idr]
                self.KN[iha,idr]=math.sqrt(self.zN[iha,idr]**2+self.yN[iha,idr]**2)

                print("\n\t\t\tVol. Disp=%12.6f\tHeel=%12.6f\tKN=%12.6f"%(self.VOL[iha,idr],phi*180/math.pi,self.KN[iha,idr]))
    # --------------------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------------------
    # Function to compute heeled and immersed half section properties
    # Py-NAT  2021a :  PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def compImmersedSectionProperties(self,sec,zd,ang,side):

        # Initialise return Variables with Default Values
        hBw=0
        Ar=0
        LM=0
        VM=0
        TM=0
        CL=0
        YC=0
        ZC=0

        dz=0.25
        z=np.arange(sec.Z[0],sec.Z[-1],dz)
        da=np.arange(z.size,dtype=np.float64)
        dtm=np.arange(z.size,dtype=np.float64)
        dvm=np.arange(z.size,dtype=np.float64)

        # print("\nThe Heeled Section Calculations for %f"%(sec.X))
        # print("Zd=",zd)
        # print("phi=",ang*180/math.pi)
        # print("side=",side)

        if ang!=0:
            for i in range(z.size):
                yw=(z[i]-zd)*math.tan(ang-math.pi/2)
                ycp=sec.getHalfBreadth(z[i])
                ycs=-ycp
                if (side>=0):
                    # Measure strip ends Left Side of WL
                    ya=max(ycp,yw)
                    yb=max(ycs,yw)
                else:
                    # Measure strip ends Right Side of WL
                    ya=min(ycp,yw)
                    yb=min(ycs,yw)
                da[i]=max(ya-yb,0)
                dtm[i]=da[i]*(ya+yb)/2
                dvm[i]=da[i]*z[i]
                # Py-NAT 2021a :  PRAVEEN KUMAR CH
                #print("yw=",yw,"yc=",yc,"ya=",ya,"yb=",yb,"da=",da[i])
            
            a_fz=interp1d(z,da,kind='cubic')
            sq=quad(a_fz,z[0],z[-1])
            Ar=sq[0]

            vm_fz=interp1d(z,dvm,kind='cubic')
            sq=quad(vm_fz,z[0],z[-1])
            VM=sq[0]

            tm_fz=interp1d(z,dtm,kind='cubic')
            sq=quad(tm_fz,z[0],z[-1])
            TM=sq[0]

            LM=Ar*sec.X
            if Ar!=0:
                YC=TM/Ar
                ZC=VM/Ar
            secData=list([hBw,Ar,LM,TM,VM,CL,YC,ZC])
        else :
            secData=sec.sectionHydroStatics(zd)
            secData[0]=2*secData[0]
            secData[1]=2*secData[1]
            secData[2]=2*secData[2]
            secData[3]=0
            secData[4]=2*secData[4]
            secData[5]=2*secData[5]
            secData[6]=0
            if secData[0]!=0:
                secData[7]=secData[4]/secData[1]
            else:
                secData[7]=0

        # print("SecData=",secData)
        return secData
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # The Method to Read the Data from File
    # Py-NAT 2021a :  PRAVEEN  KUMAR CH
    # --------------------------------------------------------------------------------------
    def writeToFile(self,fs):

        print("\n\t\tPyNAT: Writing KN Curves Data to File:\t ...\\",fs.okn_path)

        if not os.path.exists(fs.las_dir):
            os.makedirs(fs.las_dir)

        fp=open(fs.okn_path,"wt")
        fp.write('\nPy-NAT: Large angle Stability - KN Curves   :\t')
        fp.write('\nDraft     \tHeel      \tKN        \t')
        fp.write(  'zN        \tyN        \t')
        fp.write(  'm         \tc         \tWPA       \t')
        fp.write(  'LCB       \tTCB       \tVCB       \t')
        fp.write(  'Volume    \tMass      \tWet Area  \t')


        for i in range(self.nHeel):
            for k in range(self.nDraft):
                # Py-NAT 2021a :  PRAVEEN KUMAR  CH
                fp.write("\n%10.3E\t%10.3E\t%10.3E\t"%(self.rangeDrafts[i][k],self.rangeHeels[i],self.KN[i][k]))
                fp.write(  "%10.3E\t%10.3E\t"%(self.zN[i][k],self.yN[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.m[i][k],self.c[i][k],self.WPA[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.LCB[i][k],self.TCB[i][k],self.VCB[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.VOL[i][k],self.MASS[i][k],self.WSA[i][k]))             
        fp.close()
        
        plt.figure("KN CURVES")
        # PyNAT 2021 a : PraveenKumarCH
        for i in range(self.nHeel):
            hstr = format(self.rangeHeels[i]*180/math.pi,"0.2f")
            plt.plot(self.MASS[i,:],self.KN[i,:],label='KN@Heel='+hstr)
        plt.title("Cross Curves of Stability")
        plt.xlabel("Displacement")
        plt.ylabel("KN")
        lg=plt.legend(bbox_to_anchor=(1.05, 1))
        plt.savefig(os.path.join(fs.las_dir,"LAS_KN.png"),            
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')
        plt.close('all')
    # ----------------------------------------------------------------------------------------