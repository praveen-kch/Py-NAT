'''
************************************************************************************
    File            : LAS_KN_1A.py
    
    File Purpose    : Define the Class "LAS_KN" 
                        with hybrid technique to compute heeled section parameters
    
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

class LAS_KN_1A:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def __init__(self,fs,GivenShip,ep):
        self.techq="FAST"
        self.readFromFile(fs)
        self.compute(GivenShip,ep)
        self.writeToFile(fs)
    # --------------------------------------------------------------------------------------
    

    # --------------------------------------------------------------------------------------
    # The Method to Read the Data from File
    # PyNAT 2021 a : PraveenKumarCH
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
    # PyNAT 2021 a : PraveenKumarCH
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
            # print(iha,self.rangeDrafts[iha,:])
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


        # Compute for each heel angle
        for iha in range(self.nHeel):

            # For each Draft in the draft range
            for idr in range(self.nDraft):

                zd=self.rangeDrafts[iha][idr]

                phi=self.rangeHeels[iha]

                BR=np.zeros(GivenShip.nst,dtype=np.float64)
                Ar=np.zeros(GivenShip.nst,dtype=np.float64)
                LM=np.zeros(GivenShip.nst,dtype=np.float64)
                TM=np.zeros(GivenShip.nst,dtype=np.float64)
                VM=np.zeros(GivenShip.nst,dtype=np.float64)
                CL=np.zeros(GivenShip.nst,dtype=np.float64)

                for ist in range(GivenShip.nst):

      
                    pwli=self.computeIntersections(GivenShip.secCrv[ist],zd,phi)

                    swli=self.computeIntersections(GivenShip.secCrv[ist],zd,-phi)

                    # print("Station no =%5d Heel =%6.3f draft =%6.3f"%(ist,phi,zd),":\n")           
                    # print("Port Side Intersection Points:",pwli)
                    # print("Starboard Side Intersection Points:",swli)
                    # print("Praveen Kumar CH PyNAT : HAULT")

                    # When  WL intersects closed section at only two locations
                    if (len(pwli)==1 and len(swli)==1):

                        secData=GivenShip.secCrv[ist].sectionHydroStatics(pwli[0])
                        #print("The Section Prop at Draft Zp",pwli[0],"\n",secData)

                        ht=zd-pwli[0]
                        bt=abs(ht*math.tan(math.pi/2-(phi)))
                        Atr=0.5*ht*bt
                        hbw=math.sqrt(ht**2+bt**2)
                        tl=Atr*bt/3
                        vl=Atr*(zd-ht*2/3)

                        secData[0]=hbw
                        secData[1]+=Atr
                        secData[2]+=Atr*GivenShip.x[ist]
                        secData[3]+=tl
                        secData[4]+=vl
                        #print("ht=%10.5f\tbt=%10.5f\tAtr=%10.5f\tHb=%10.5f\ttl=%10.5f\tvl=%10.5f\t"%(ht,bt,Atr,hbw,tl,vl))

                        secData2=GivenShip.secCrv[ist].sectionHydroStatics(swli[0])
                        #print("The Section Prop at Draft Zp",swli[0],"\n",secData2)

                        ht=zd-swli[0]
                        bt=abs(ht*math.tan(math.pi/2-(phi)))
                        Atr=0.5*ht*bt
                        hbw=math.sqrt(ht**2+bt**2)
                        tl=Atr*bt/3
                        vl=Atr*(zd-ht*2/3)

                        secData2[0]=hbw
                        secData2[1]+=Atr
                        secData2[2]+=Atr*GivenShip.x[ist]
                        secData2[3]+=tl
                        secData2[4]+=vl

                        #print("ht=%10.5f\tbt=%10.5f\tAtr=%10.5f\tHb=%10.5f\ttl=%10.5f\tvl=%10.5f\t"%(ht,bt,Atr,hbw,tl,vl))

                        secData[0]+=secData2[0]             # Breadth on WL
                        secData[1]+=secData2[1]             # Area
                        secData[2]+=secData2[2]             # LM
                        secData[3]-=secData2[3]             # TM
                        secData[4]+=secData2[4]             # VM
                        secData[5]+=secData2[5]             # CL                           

                        if(abs(phi)>=math.pi/2):
                            # In case if the ship is at a capsizing heel angle...
                            fullSecData=GivenShip.secCrv[ist].sectionHydroStatics(GivenShip.secCrv[ist].Z[-1])
                            secData[1]=2*fullSecData[1]-secData[1]
                            secData[2]=2*fullSecData[2]-secData[2]
                            secData[3]=-secData[3]
                            secData[4]=2*fullSecData[4]-secData[4]
                            secData[5]=2*fullSecData[5]-secData[5]
                            
                        if secData[1]!=0:
                            secData[6]=secData[3]/secData[1]    # YC
                            secData[7]=secData[4]/secData[1]    # ZC
                        else:
                            secData[6]=0    # YC
                            secData[7]=0    # ZC 

                        #print("Section Data Computed using Plane Geometry:\n",secData)

                    else:

                        if(phi<=0):
                            side=1
                        else:
                            side=-1
                        secData=self.compImmersedSectionProperties(GivenShip.secCrv[ist],zd,phi,side) 

                        # This portion of code for computation of Girth length of immersed section 
                        # This is improper and needs to be improved...
                        # Since not a part of KN calculations so ignored for now
                        secData[5]=0
                        if(len(pwli)!=0):
                            secData[5]+=GivenShip.secCrv[ist].getGirthLength(pwli[-1])
                        if(len(swli)!=0):
                            secData[5]+=GivenShip.secCrv[ist].getGirthLength(swli[-1])

                        # This portion of code for computation of section Breadth on water plane 
                        # This is improper and needs to be improved...
                        # Since not a part of KN calculations so ignored for now
                        secData[0]=0
                        pwli.extend(swli)
                        for i in range(0,len(pwli)-1,2):
                            dy=(pwli[i]-pwli[i+1])/abs(math.cos(phi-math.pi/2))
                            secData[0]+=dy
                        
                        #print("Section Data Computed using Numerical Integration:\n",secData)
                    

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

                # Finding KN
                self.m[iha,idr]=math.tan(math.pi/2-phi)
                self.c[iha,idr]=self.VCB[iha,idr]-self.m[iha,idr]*self.TCB[iha,idr]
                self.zN[iha,idr]=self.c[iha,idr]/(1+self.m[iha,idr]**2)
                self.yN[iha,idr]=-self.m[iha,idr]*self.zN[iha,idr]
                self.KN[iha,idr]=math.sqrt(self.zN[iha,idr]**2+self.yN[iha,idr]**2)

                print("\n\t\t\tVol. Disp=%12.6f\tHeel=%12.6f\tKN=%12.6f"%(self.VOL[iha,idr],phi*180/math.pi,self.KN[iha,idr]))
                print("\n\t\t\t\tyn=%10.5f \t zn=%10.5f"%(self.yN[iha,idr],self.zN[iha,idr]))
    # --------------------------------------------------------------------------------------
    
    # --------------------------------------------------------------------------------------
    # Function to compute heeled and immersed half section properties
    # PyNAT 2021 a : PraveenKumarCH
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
        # PyNAT 2021 a : PraveenKumarCH
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
    # The Method to predict water line intersection of a heeled section
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def computeIntersections(self,sec,zd,ang):

        # Tolerance to check intersection
        tol=0.001

        # List of Intersection Z coordinates
        inPts=list()

        # If angle is zero i.e., zero heel angle
        if (math.tan(ang)==0):
            inPts.append(zd)
            return inPts

        # pre calculation of slope of water line
        tanang=math.tan(ang-math.pi/2)

        # Check for Water line cutting the bottom base line of section
        z1=sec.Z[0]
        yC=sec.getHalfBreadth(z1)
        yW=(z1-zd)*tanang
        if(yW<(yC-tol) and yW>=0):
            inPts.append(z1) 

        # Check for water line cutting the section curve 
        for i in range(sec.N-1):
            
            # take each pair of consecutive points 
            z1=sec.Z[i]
            z2=sec.Z[i+1]

            dy1=sec.getHalfBreadth(z1)-(z1-zd)*tanang
            dy2=sec.getHalfBreadth(z2)-(z2-zd)*tanang

            #print("z1=%10.5f,yc1=%10.5f,yw1=%10.5f,dy1=%10.5f,z2=%10.5f,yc2=%10.5f,yw2=%10.5f,dy2=%10.5f"%(z1,yc1,yw1,dy1,z2,yc2,yw2,dy2))

            if(dy1*dy2<0):

                for _ in range(25):

                    if abs(dy1)<=tol:
                        inPts.append(z1)
                        #print("Appended z1",z1)
                        break

                    if abs(dy2)<=tol:
                        inPts.append(z2)
                        #print("Appended z2",z2)
                        break

                    zm=(z1+z2)/2
                    dym=sec.getHalfBreadth(zm)-(zm-zd)*tanang

                    if(dy1*dym<0):
                        z2=zm
                        dy2=dym
                    if(dy2*dym<0):
                        z1=zm
                        dy1=dym
                    
                    #print("z1=%10.5f,dy1=%10.5f,z2=%10.5f,dy2=%10.5f"%(z1,dy1,z2,dy2))

        # Check for Water line cutting the top deck line
        z1=sec.Z[-1]
        yC=sec.getHalfBreadth(z1)
        yW=(z1-zd)*tanang
        if(yW<(yC-tol) and yW>=0):
            inPts.append(z1)

        # Return List of Points
        return inPts

    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # The Method to Read the Data from File
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def writeToFile(self,fs):

        print("\n\t\tPyNAT: Writing KN Curves Data to File:\t ...\\",fs.okn_fast_path)

        if not os.path.exists(fs.las_dir):
            os.makedirs(fs.las_dir)

        fp=open(fs.okn_fast_path,"wt")
        fp.write('\nPy-NAT: Large angle Stability - KN Curves   :\t')
        fp.write('\nDraft     \tHeel      \tKN        \t')
        fp.write(  'zN        \tyN        \t')
        fp.write(  'm         \tc         \tWPA       \t')
        fp.write(  'LCB       \tTCB       \tVCB       \t')
        fp.write(  'Volume    \tMass      \tWet Area  \t')


        for i in range(self.nHeel):
            for k in range(self.nDraft):
                    # PyNAT 2021 a : PraveenKumarCH
                fp.write("\n%10.3E\t%10.3E\t%10.3E\t"%(self.rangeDrafts[i][k],self.rangeHeels[i],self.KN[i][k]))
                fp.write(  "%10.3E\t%10.3E\t"%(self.zN[i][k],self.yN[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.m[i][k],self.c[i][k],self.WPA[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.LCB[i][k],self.TCB[i][k],self.VCB[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.VOL[i][k],self.MASS[i][k],self.WSA[i][k]))             
        fp.close()
        

        plt.figure("KN CURVES")
        for i in range(self.nHeel):
            hstr = format(self.rangeHeels[i]*180/math.pi,"0.2f")
            plt.plot(self.MASS[i,:],self.KN[i,:],label='KN@Heel='+hstr)
        plt.title("Cross Curves of Stability")
        plt.xlabel("Displacement")
        plt.ylabel("KN")
        lg=plt.legend(bbox_to_anchor=(1.05, 1))
        plt.savefig(os.path.join(fs.las_dir,"LAS_KN_FAST.png"),            
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')
        plt.close('all')
    # ----------------------------------------------------------------------------------------
