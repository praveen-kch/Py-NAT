'''
************************************************************************************
    File            : HydroStatics.py
    
    File Purpose    : Define the Class "HydroStatics"
    
    Class Purpose   : 
                    -> Contains Function to Compute Hydrostatics
                    -> Computations are done for Floating Conditions Defined by 
                        -> Draft 
                        -> Trim
                    -> The range of Drafts and Trim are read from the File "Set_HS.dat"
                    -> Computations are performed over a range of drafts at all Individual Sections
                    -> Section wise drafts are computed for each Floating Condition
                    -> Computations are performed for each trim and draft combination Entire Hull Volume
                    -> Functions to write the sectional and Hull Hydrostatics to File "Res_HS.dat"
                    -> Plots and Saves Hydrostatics as image files of format ".png"
    
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
# ------------------------------------------------------------------------------------------
# The Class to Compute and Store Hydrostatic Data
# Py-NAT 2021a : PRAVEEN KUMAR CH
# ------------------------------------------------------------------------------------------

class HydroStatics:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def __init__(self,fs,shipGrid,ep):
        self.readFromFile(fs,shipGrid)
        self.computeSectionalHydroStatics(shipGrid)
        self.computeHullHydroStatics(shipGrid,ep)
        self.writeSectionalHydrostaticsToFile(fs,shipGrid)
        self.writeHullHydroStaticsToFile(fs,shipGrid)
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Read The Hydrostatic Solver Settings File 
    # create list of arrays of sectional drafts for each floating condition
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def readFromFile(self,fs,shipGrid):

        print("\n\t\tPyNAT : Reading File:\t...\\",fs.ihs_path)
        # Reading the Solver properties
        fp=open(fs.ihs_path,"rt")
        lines=fp.readlines()
        fp.close()

        nums=lines[0].split(",")
        self.nDraft=len(nums)
        self.rangeDrafts = np.array(range(self.nDraft),dtype=np.float64)
        for i in range(self.nDraft):
            self.rangeDrafts[i] =float(nums[i]) 

        nums=lines[1].split(",")
        self.nTrims=len(nums)
        self.rangeTrims = np.array(range(self.nTrims),dtype=np.float64)
        for i in range(self.nTrims):
            self.rangeTrims[i] =float(nums[i])*math.pi/180 

        # Compute Sectional Drafts for each FLoating condition
        # For each floating condition the section drafts are stored in an array
        # A list of all such arrays is created
        self.FC=list()
        for j in range(self.nTrims):
            temp=list()
            for i in range(self.nDraft):
                drafts=np.array(range(shipGrid.nst),dtype=np.float64)
                for k in range(shipGrid.nst):
                    drafts[k]=self.rangeDrafts[i]-shipGrid.xms[k]*math.tan(self.rangeTrims[j])
                temp.append(drafts)
            self.FC.append(temp)
        self.nFC=self.nTrims*self.nDraft

    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to compute sectional hydrostatics over a range of drafts
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def computeSectionalHydroStatics(self,shipGrid):

        print("\n\t\tPyNAT : Computing Section Wise Hydrostatic Properties")

        # Initializing the Arrays to Hold the Sectional Data
        # Each Row is for a particular station
        # Each Column is for a particular draft
        self.Bw=np.zeros((shipGrid.nst, self.nDraft))
        self.Area=np.zeros((shipGrid.nst, self.nDraft))
        self.VM=np.zeros((shipGrid.nst, self.nDraft))
        self.TM=np.zeros((shipGrid.nst, self.nDraft))
        self.LM=np.zeros((shipGrid.nst, self.nDraft))
        self.CL=np.zeros((shipGrid.nst, self.nDraft))
        self.YC=np.zeros((shipGrid.nst, self.nDraft))
        self.ZC=np.zeros((shipGrid.nst, self.nDraft))

        # Loop to Compute for Every section in the shipGrid
        for i in range(shipGrid.nst):
            print("\t\t\t\tStation:",i)
            # Loop to compute for each draft in the range
            for k in range(self.nDraft):
                secData=shipGrid.secCrv[i].sectionHydroStatics(self.rangeDrafts[k])
                self.Bw[i][k]=2*secData[0]
                self.Area[i][k]=2*secData[1]
                self.LM[i][k]=2*secData[2]
                self.TM[i][k]=secData[3]
                self.VM[i][k]=2*secData[4]
                self.CL[i][k]=2*secData[5]
                self.YC[i][k]=secData[5]
                self.ZC[i][k]=secData[7]

    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Compute Hull Volume Hydrostatics
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def computeHullHydroStatics(self,shipGrid,ep):
        
        print("\n\t\tPyNAT : Computing the Hull Hydrostatic Properties")
        # Initializing the arrays to hold hydrostatics
        self.VOL=np.zeros((self.nTrims,self.nDraft))
        self.WSA=np.zeros((self.nTrims,self.nDraft))
        self.WPA=np.zeros((self.nTrims,self.nDraft))
        self.LCF=np.zeros((self.nTrims,self.nDraft))
        self.MASS=np.zeros((self.nTrims,self.nDraft))
        self.VCB=np.zeros((self.nTrims,self.nDraft))
        self.TCB=np.zeros((self.nTrims,self.nDraft))
        self.LCB=np.zeros((self.nTrims,self.nDraft))
        self.IXX0=np.zeros((self.nTrims,self.nDraft))
        self.IXXG=np.zeros((self.nTrims,self.nDraft))
        self.IYY0=np.zeros((self.nTrims,self.nDraft))
        self.BMT=np.zeros((self.nTrims,self.nDraft))
        self.BML=np.zeros((self.nTrims,self.nDraft))
        self.KMT=np.zeros((self.nTrims,self.nDraft))
        self.KML=np.zeros((self.nTrims,self.nDraft))
        self.TPC=np.zeros((self.nTrims,self.nDraft))
        self.MCT=np.zeros((self.nTrims,self.nDraft))
        self.draftAP=np.zeros((self.nTrims,self.nDraft))
        self.draftFP=np.zeros((self.nTrims,self.nDraft))

        # Loop for each Floating Condition [List of Array of Sectional Drafts] 
        for itr in range(self.nTrims):
            for idr in range(self.nDraft):
                # Py-NAT 2021a : PRAVEEN KUMAR CH
                Area=np.zeros(shipGrid.nst)
                VM=np.zeros(shipGrid.nst)
                TM=np.zeros(shipGrid.nst)
                LM=np.zeros(shipGrid.nst)
                CL=np.zeros(shipGrid.nst)
                YC=np.zeros(shipGrid.nst)
                ZC=np.zeros(shipGrid.nst)

                y=np.arange(shipGrid.nst)
                yx=np.arange(shipGrid.nst)
                yxx=np.arange(shipGrid.nst)
                yyy=np.arange(shipGrid.nst)

                for j in range(shipGrid.nst):
                    secData=shipGrid.secCrv[j].sectionHydroStatics(self.FC[itr][idr][j])
                    y[j]=secData[0]
                    Area[j]=secData[1]
                    LM[j]=secData[2]
                    TM[j]=secData[3]
                    VM[j]=secData[4]
                    CL[j]=secData[5]
                    YC[j]=secData[6]
                    ZC[j]=secData[7]
                    yx[j]=y[j]*shipGrid.x[j]
                    yxx[j]=2*yx[j]*shipGrid.x[j]
                    yyy[j]=8*(y[j]**3)/12

                Ar_fX=interp1d(shipGrid.x,Area,kind='cubic')
                sq=quad(Ar_fX,shipGrid.x[0],shipGrid.x[-1])
                self.VOL[itr][idr]=2*sq[0]
                self.MASS[itr][idr]=self.VOL[itr][idr]*ep.waterDensity

                Sa_fX=interp1d(shipGrid.x,CL,kind='cubic')
                sq=quad(Sa_fX,shipGrid.x[0],shipGrid.x[-1])
                self.WSA[itr][idr]=2*sq[0]

                lm_fX=interp1d(shipGrid.x,LM,kind='cubic')
                sq=quad(lm_fX,shipGrid.x[0],shipGrid.x[-1])
                self.LCB[itr][idr]=2*sq[0]/self.VOL[itr][idr]

                tm_fX=interp1d(shipGrid.x,TM-TM,kind='cubic')
                sq=quad(tm_fX,shipGrid.x[0],shipGrid.x[-1])
                self.TCB[itr][idr]=sq[0]/self.VOL[itr][idr]

                vm_fX=interp1d(shipGrid.x,VM,kind='cubic')
                sq=quad(vm_fX,shipGrid.x[0],shipGrid.x[-1])
                self.VCB[itr][idr]=2*sq[0]/self.VOL[itr][idr]

                wpa_fX=interp1d(shipGrid.x,y,kind='cubic')
                sq=quad(wpa_fX,shipGrid.x[0],shipGrid.x[-1])
                self.WPA[itr][idr]=2*sq[0]

                wpm_fX=interp1d(shipGrid.x,yx,kind='cubic')
                sq=quad(wpm_fX,shipGrid.x[0],shipGrid.x[-1])
                self.LCF[itr][idr]=2*sq[0]/self.WPA[itr][idr]

                IXX0_fx=interp1d(shipGrid.x,yxx,kind='cubic')
                sq=quad(IXX0_fx,shipGrid.x[0],shipGrid.x[-1])
                self.IXX0[itr][idr]=sq[0]

                IYY0_fx=interp1d(shipGrid.x,yyy,kind='cubic')
                sq=quad(IYY0_fx,shipGrid.x[0],shipGrid.x[-1])
                self.IYY0[itr][idr]=sq[0]

                self.IXXG[itr][idr]=self.IXX0[itr][idr]-self.WPA[itr][idr]*(self.LCF[itr][idr]**2)
                
                self.draftAP[itr][idr]=self.rangeDrafts[idr]+(shipGrid.MS-shipGrid.AP)*math.tan(self.rangeTrims[itr])
                self.draftFP[itr][idr]=self.rangeDrafts[idr]-(shipGrid.FP-shipGrid.MS)*math.tan(self.rangeTrims[itr])
                # Py-NAT 2021a : PRAVEEN KUMAR CH
                print("\n\t\t\tTrim=%f\tDraft=%f\tVolume=%f"%(self.rangeTrims[itr],self.rangeDrafts[idr],self.VOL[itr][idr]))

        self.BML=self.IXXG/self.VOL
        self.BMT=self.IYY0/self.VOL
        self.KML=self.BML+self.VCB
        self.KMT=self.BMT+self.VCB
        self.TPC=self.WPA*ep.waterDensity/100000
        self.MCT=self.MASS*self.BML/(shipGrid.L*100)
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Write Section wise Hydrostatics to file
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def writeSectionalHydrostaticsToFile(self,fs,shipGrid):

        print("\n\t\tPyNAT : Writing Section wise Hydrostatic Data to File:\t...\\",fs.ohs_path)

        if not os.path.exists(fs.hs_dir):
            os.makedirs(fs.hs_dir)

        fp=open(fs.ohs_path,"wt")
        fp.write("\nPy-NAT: Hydrostatic Properties")
        for i in range(shipGrid.nst):
            fp.write('\n\n\nSection   :\t'+str(i))
            fp.write('\nDraft     \tB/2          \tArea      \tLong. Mom \tVert. Mom \tTrans. Mom\tCurv Leng ')
            for k in range(self.nDraft):
                fp.write("\n%10.3E\t%10.3E\t%10.3E\t%10.3E\t%10.3E\t%10.3E\t%10.3E"%(self.rangeDrafts[k],self.Bw[i,k],self.Area[i,k],self.LM[i,k],self.VM[i,k],self.TM[i,k],self.CL[i,k]))
        fp.close()
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        for i in range(shipGrid.nst):

            plt.figure("Sec Ar vs Draft")
            plt.plot(self.rangeDrafts,self.Area[i],label=str(i))
            plt.title("Sectional Areas")
            plt.xlabel("Draft")
            plt.ylabel("Area")         

            plt.figure("Sec VM vs Draft")
            plt.plot(self.rangeDrafts,self.VM[i],label=str(i))
            plt.title("Vertical Moment abt BL")
            plt.xlabel("Draft")
            plt.ylabel("Moment")

            plt.figure("Sec TM vs Draft")
            plt.plot(self.rangeDrafts,self.TM[i],label=str(i))
            plt.title("Transverse Moment abt CL (for Half Section only)")
            plt.xlabel("Draft")
            plt.ylabel("Moment")

            plt.figure("Sec Girth vs Draft")
            plt.plot(self.rangeDrafts,self.CL[i],label=str(i))
            plt.title("Girth Length")
            plt.xlabel("Draft")
            plt.ylabel("Girth")

            plt.figure("Sec Zc vs Draft")
            plt.plot(self.rangeDrafts,self.ZC[i],label=str(i))
            plt.title("Z Centroid from BL ")
            plt.xlabel("Draft")
            plt.ylabel("Zc")

            plt.figure("Sec Yc vs Draft")
            plt.plot(self.rangeDrafts,self.YC[i],label=str(i))
            plt.title("Y Centroid from CL (for Half Section only)")
            plt.xlabel("Draft")
            plt.ylabel("Yc")

        plt.figure("Sec Ar vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_Areas.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Sec VM vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_VertMom.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Sec TM vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_TransMom.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Sec Girth vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_Girth.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Sec Zc vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_ZCent.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Sec Yc vs Draft")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,'Sect_YCent.png'),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')
        plt.close('all')
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Write HULL Hydrostatics to file
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def writeHullHydroStaticsToFile(self,fs,shipGrid):
        # Py-NAT 2021a : PRAVEEN KUMAR CH   
        print("\n\t\tPyNAT : Writing Hull Hydrostatic Data to File:\t...\\",fs.ohs_path)
        fp=open(fs.ohs_path,"at")
        fp.write('\n\n\nHull Hydrostatics   :\t')
        fp.write('\nDraft MS  \tTrim      \t')
        fp.write(' Draft AP  \tDraft FP  \t')
        fp.write(  'Volume    \tMass      \tWet Hull Ar.\t')
        fp.write(  'LCB       \tTCB       \tVCB       \t')
        fp.write(  'WPA       \tLCF       \t')
        fp.write(  'IXX0      \tIXXG      \tIYY0        \t')        
        fp.write(  'BM_L      \tBM_T      \tKM_L        \tKM_T        \t')
        fp.write(  'TPC       \tMCT       \t')
        fp.write(  '')

        for i in range(self.nTrims):
            for k in range(self.nDraft):
                fp.write("\n%10.3E\t%10.3E\t"%(self.rangeDrafts[k],self.rangeTrims[i]))
                fp.write(  "%10.3E\t%10.3E\t"%(self.draftAP[i][k],self.draftFP[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.VOL[i][k],self.MASS[i][k],self.WSA[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.LCB[i][k],self.TCB[i][k],self.VCB[i][k]))
                fp.write(  "%10.3E\t%10.3E\t"%(self.WPA[i][k],self.LCF[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t"%(self.IXX0[i][k],self.IXXG[i][k],self.IYY0[i][k]))
                fp.write(  "%10.3E\t%10.3E\t%10.3E\t%10.3E\t"%(self.BML[i][k],self.BMT[i][k],self.KML[i][k],self.KMT[i][k]))
                fp.write(  "%10.3E\t%10.3E\t"%(self.TPC[i][k],self.MCT[i][k]))
                
        fp.close()

        plt.figure("Vol Disp vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.VOL[i],label='Volume@Trim='+hstr)
        plt.title("Vol Disp vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("Volume")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_VolDisp.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Mass Disp vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.MASS[i],label='Mass@Trim='+hstr)
        plt.title("Mass Disp vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("Mass")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_MassDisp.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("Wet Area vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.WSA[i],label='WSA@Trim='+hstr)
        plt.title("Wetted Surface Area vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("Wet Surf Area")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_WetArea.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("WPA vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.WPA[i],label='Aw@Trim='+hstr)
        plt.title("Water Plane Area vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("WPA")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_WPA.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("LCB vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.LCB[i]-shipGrid.MS,label='LCB@Trim='+hstr)
        plt.title("LCB (Frm Mid Ship) vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("LCB")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_LCB.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("LCF vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.LCF[i]-shipGrid.MS,label='LCF@Trim='+hstr)
        plt.title("LCF (Frm Mid Ship) vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("LCF")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_LCF.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("VCB vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.VCB[i],label='VCB@Trim='+hstr)
        plt.title("VCB (Frm Base) vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("VCB")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_VCB.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("TCB vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.TCB[i],label='TCB@Trim='+hstr) 
        plt.title("TCB From CL vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("TCB")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_TCB.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("TPC vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.TPC[i],label='TPC@Trim='+hstr)
        plt.title("Tons Per CM immersion vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("TPC")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_TPC.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("MCT vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.MCT[i],label='MCT@Trim='+hstr)
        plt.title("Mom. to Change Trim vs Draft")
        plt.xlabel("Draft@MS")
        plt.ylabel("MCT")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_MCT.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("KMT vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.KMT[i],label='KMT@Trim='+hstr)
        plt.title("KM Trans.")
        plt.xlabel("Draft@MS")
        plt.ylabel("KM_T")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_KMT.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("KML vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.KML[i],label='KML@Trim='+hstr)
        plt.title("KM Long.")
        plt.xlabel("Draft@MS")
        plt.ylabel("KM_L")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_KML.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("BMT vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.BMT[i],label='BMT@Trim='+hstr)
        plt.title("Trans. Metacentric Radius")
        plt.xlabel("Draft@MS")
        plt.ylabel("BM_T")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_BMT.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')

        plt.figure("BML vs Draft")
        for i in range(self.nTrims):
            hstr = format(self.rangeTrims[i],"6.4f")
            plt.plot(self.rangeDrafts,self.BML[i],label='BML@Trim='+hstr)
        plt.title("Long. Metacentric Radius")
        plt.xlabel("Draft@MS")
        plt.ylabel("BM_L")
        lg=plt.legend(bbox_to_anchor=(1.05, 1)) 
        plt.savefig(os.path.join(fs.hs_dir,"HS_BML.png"),
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')
        
        plt.close('all')
    # --------------------------------------------------------------------------------------
    # Py-NAT 2021a : PRAVEEN KUMAR CH
# ------------------------------------------------------------------------------------------