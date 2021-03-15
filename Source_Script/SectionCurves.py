
'''
************************************************************************************
    File            : SectionCurves.py
    
    File Purpose    : Define the Class "SectionCurves"
    
    Class Purpose   : 
                    -> Create and Store Curve Fitting Data for a Section Curve 
                    -> Assumes symmetry: Stores and Computes for half sections only.
                    -> Section Data Stored
                        -> X location of the stations 
                        -> Curve point coordinates (Y,Z)
                        -> Continuity Type at each point C
                        -> No. of Curvelets and Ends/Knuckle points
                        -> Spline Curve Fit Data of each Curvelet
                    -> Functions to compute and Return the Section Properties at a given draft
                    -> Follwing Section Properties computed at a given draft
                        -> Half Breadth
                        -> Area
                        -> Moments : Long., Trans. , Vertical.
                        -> Girth Lengths
    
    App             : PyNAT - Python Based Naval Architecture Tools
    Version         : 2021.0
    Paradigm        : Object Based Framework 
    Author          : Praveen Kumar Ch
*************************************************************************************
'''
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate 
from scipy.integrate import quad
from scipy.misc import derivative

class SectionCurves:

    # --------------------------------------------------------------------------------------
    # The Constructor for a Section
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def __init__(self,X,Y,Z,C):

        # Set Data from Arguments
        self.N=Y.size
        self.X=X
        self.Y=Y
        self.Z=Z
        self.C=C
        self.chined=0
        self.fittype=1
        self.sf=0.25
        self.maxOrd=2

        # Compute no of Curvelets and their end points
        # temp=np.where(self.C==0)
        # self.iCE=temp[0]
        # self.nCP=self.iCE.size-1

        # 
        self.nCrv=0
        self.sCrv=list()
        self.eCrv=list()
        self.sCrv.append(0)
        flag=1
        for i in range(1,self.N):
            if(flag==1 and self.C[i]==0):
                self.nCrv+=1
                self.eCrv.append(i)
                self.sCrv.append(i)
                flag=1
            elif(flag ==1 and self.C[i]==-1):
                self.nCrv+=1
                self.eCrv.append(i)
                flag=0
            elif(flag==0):
                self.sCrv.append(i)
                flag=1

        # Set Type to Chined Hull Section is all points are knuckles
        if self.nCrv==(self.N-1):
            self.chined=1
        
        # Fit Splines Function for each curvelet
        self.y_fz=list()
        self.zy_fz=list()
        self.yy_fz=list()
        self.L_fz=list()
        for i in range(self.nCrv):
            s=self.sCrv[i]
            e=self.eCrv[i]

            if(self.fittype==0):
                # Interpolation Function Cubic Spline Fit
                e=e+1
                self.y_fz.append(interpolate.interp1d(self.Z[s:e],self.Y[s:e],kind='cubic'))
                self.zy_fz.append(interpolate.interp1d(self.Z[s:e],self.Z[s:e]*self.Y[s:e],kind='cubic'))
                self.yy_fz.append(interpolate.interp1d(self.Z[s:e],self.Y[s:e]*self.Y[s:e]/2,kind='cubic'))
                df=np.arange(e-s,dtype=np.float64)
                df[0]=(self.y_fz[i](self.Z[s+1])-self.y_fz[i](self.Z[s]))/(self.Z[s+1]-self.Z[s])
                df[1:-1]=(self.y_fz[i](self.Z[s+2:e])-self.y_fz[i](self.Z[s:e-2]))/(self.Z[s+2:e]-self.Z[s:e-2])
                df[-1]=(self.y_fz[i](self.Z[e-1])-self.y_fz[i](self.Z[e-2]))/(self.Z[e-1]-self.Z[e-2])
                self.L_fz.append(interpolate.interp1d(self.Z[s:e],np.sqrt(1+df**2),kind='cubic'))
 
            if(self.fittype==1):
                # Smoothing Function Cubic Spline Fit
                if(e-s>self.maxOrd):
                    order=self.maxOrd
                else:
                    order=e-s
                e=e+1
                self.y_fz.append(interpolate.splrep(self.Z[s:e],self.Y[s:e],k=order,s=self.sf))
                self.zy_fz.append(interpolate.splrep(self.Z[s:e],self.Z[s:e]*self.Y[s:e],k=order,s=self.sf))
                self.yy_fz.append(interpolate.splrep(self.Z[s:e],self.Y[s:e]*self.Y[s:e]/2,k=order,s=self.sf))
                df=interpolate.splev(self.Z[s:e],self.y_fz[i],der=1)
                self.L_fz.append(interpolate.splrep(self.Z[s:e],np.sqrt(1+df**2),k=order,s=self.sf))

    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Compute half breadths at a given draft
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def getHalfBreadth(self,draft):

        hBw=0.0
        for i in range(self.nCrv):
            if(self.Z[self.sCrv[i]]<=draft and self.Z[self.eCrv[i]]>=draft):

                if(self.fittype==0):
                    tmp=self.y_fz[i](draft)
                elif(self.fittype==1):
                    tmp=interpolate.splev(draft,self.y_fz[i],0)

                hBw=tmp
                break

        return hBw
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Compute half breadths at a given draft
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def getGirthLength(self,draft):
        CL=0
        for i in range(self.nCrv):
            if(self.Z[self.sCrv[i]]<=draft):
                zs=self.Z[self.sCrv[i]]
                ze=self.Z[self.eCrv[i]]
                if draft<ze:
                    ze=draft
                if (self.fittype==0): 
                    sq=quad(self.L_fz[i],zs,ze)
                    CL+=sq[0]
                if (self.fittype==1):
                    CL+=interpolate.splint(zs,ze,self.L_fz[i])
            else:
                break
        
        return CL
    # --------------------------------------------------------------------------------------
    # --------------------------------------------------------------------------------------
    # Function to Compute sectional hydrostatics at a given draft
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def sectionHydroStatics(self,draft):

        # Initialise return Variables with Default Values
        hBw=0.0
        Ar=0
        LM=0
        VM=0
        TM=0
        CL=0
        YC=0
        ZC=0

        # Compute Area/Moment/Girth Length for each curvelets till given draft 
        # Sum the Area/Mom./Girth for all relavent Curvelets
        for i in range(self.nCrv):
            if(self.Z[self.sCrv[i]]<=draft):
                zs=self.Z[self.sCrv[i]]
                ze=self.Z[self.eCrv[i]]
                if draft<ze:
                    ze=draft
                
                if self.fittype==0:
                    # #Use the below integration for interp1d fit
                    sq=quad(self.y_fz[i],zs,ze)
                    Ar+=sq[0]
                    sq=quad(self.zy_fz[i],zs,ze)
                    VM+=sq[0]
                    sq=quad(self.yy_fz[i],zs,ze)
                    TM+=sq[0]
                    sq=quad(self.L_fz[i],zs,ze)
                    CL+=sq[0]

                elif self.fittype==1:
                    #Use the below integration for Spline Fit
                    Ar+=interpolate.splint(zs,ze,self.y_fz[i])
                    VM+=interpolate.splint(zs,ze,self.zy_fz[i])
                    TM+=interpolate.splint(zs,ze,self.yy_fz[i])
                    CL+=interpolate.splint(zs,ze,self.L_fz[i])

            else:
                break
        
        # Compute the Half Breadth at the given draft
        for i in range(self.nCrv):
            if(self.Z[self.sCrv[i]]<=draft and self.Z[self.eCrv[i]]>=draft):
                if(self.fittype==0):
                    tmp=self.y_fz[i](draft)
                elif(self.fittype==1):
                    tmp=interpolate.splev(draft,self.y_fz[i],0)
                hBw=tmp
                break

        # Longitudinal Area moment of the Section till draft 
        LM=Ar*self.X

        # Check against occurance of NaN 
        # PyNAT 2021 a : Praveen KumarCH
        if(Ar!=0):
            YC=TM/Ar
            ZC=VM/Ar
        
        # Create and return a list of computed section parameters till given draft
        secData=list([hBw,Ar,LM,TM,VM,CL,YC,ZC])

        return secData
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to interpolate offset at given spacing of drafts
    # PyNAT 2021 a : PraveenKumarCH
    # --------------------------------------------------------------------------------------
    def getPoints(self,sp):
        pts=list()

        zts=np.arange(self.Z[0],self.Z[-1],sp,dtype=np.float64)
        yts=np.arange(zts.size,dtype=np.float64)

        for i in range(zts.size):
            yts[i]=self.getHalfBreadth(zts[i])
        
        pts.append(zts)
        pts.append(yts)

        return pts

    # --------------------------------------------------------------------------------------

