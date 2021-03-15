
"""
************************************************************************************
    File            : HullGeometry.py
    
    File Purpose    : Define the Class "HullGeometry"
    
    Class Purpose   : 
                    -> Read the Hull Geometry definition data from the File "Set_HG.dat" 
                    -> Create the data variables defining the ship Hull Geometry
                    -> Create and prepare a "SectionCurves" object for each station
                    -> Print the read Data to console
                    -> Write the prepared data to "Res_HG.dat"
    
    App             : PyNAT - Python Based Naval Architecture Tools 
    Version         : 2021a
    Paradigm        : Object Oriented Model
    Author          : Praveen Kumar Ch
*************************************************************************************
"""
import numpy as np
import scipy.interpolate as spi
import matplotlib.pyplot as plt 
from SectionCurves import SectionCurves
import os

class HullGeometry:

    # --------------------------------------------------------------------------------------
    # The Constructor    
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def __init__(self,fs):
        self.tol=0.001
        self.readFromFile(fs)
        self.writeToFile(fs)
        #self.printToConsol(fs)
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Read Hull Geometry Data from file
    # Creates and Prepares SectionCurve Objects
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def readFromFile(self,fs):

        print("\n\t\t PyNAT: Reading File : ...\\",fs.ihg_path)

        fp=open(fs.ihg_path,"rt")
        lines=fp.readlines()
        fp.close()
        
        # Reading the Main Dimensions
        nums=lines[1].split("\t")
        self.L=float(nums[0])

        nums=lines[2].split("\t")
        self.B=float(nums[0])

        nums=lines[3].split("\t")
        self.D=float(nums[0])

        nums=lines[4].split("\t")
        self.AP=float(nums[0])
        self.FP=float(nums[1])
        self.MS=(self.AP+self.FP)/2

        x=list()
        y=list()
        z=list()
        c=list()

        nums=lines[5].rstrip("\n").split("\t")
        tx=float(nums[0])
        ty=float(nums[1])
        tz=float(nums[2])
        tc=1
        if(len(nums)>=4):
            if nums[3].count("k"):
                tc=0
            if nums[3].count("b"):
                tc=-1

        x.append(tx)
        y.append(list([ty]))
        z.append(list([tz]))
        c.append(list([tc]))

        # Reading the points and segregating into stations
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        for i in range(6,len(lines),1):

            nums=lines[i].rstrip("\n").split("\t")
            tx=float(nums[0])
            ty=float(nums[1])
            tz=float(nums[2])
            tc=1
            if(len(nums)>=4):
                if nums[3].count("k"):
                    tc=0
                if nums[3].count("b"):
                    tc=-1

            found=-1
            loc=0
            # check which station to add this point
            for j in range(len(x)):
                delx=tx-x[j]
                if delx>0:
                    loc=j+1
                if abs(delx)<=self.tol:
                    found=j
                    loc=j
                    break

            if found<=-1:
                # If the point dosent belong to exisitng station then add new station
                x.insert(loc,tx)
                y.insert(loc,list([ty]))
                z.insert(loc,list([tz]))
                c.insert(loc,list([tc]))
            else:
                # If the point belongs to existing station check were to place it
                for j in range(len(z[loc]),0,-1):
                    delz=tz-z[loc][j-1]
                    if(delz>self.tol):
                        y[loc].insert(j,ty)
                        z[loc].insert(j,tz)
                        c[loc].insert(j,tc)
                        break
                    elif(abs(delz)<=self.tol):
                        dely=ty-y[loc][j-1]
                        if(dely>self.tol):
                            y[loc].insert(j,ty)
                            z[loc].insert(j,tz)
                            c[loc].insert(j,tc)
                            break
                        elif(abs(dely)<=self.tol):
                            break

        # Dealing with Horizontal Lines due to points with same z values
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        for i in range(len(x)):
            # In each station
            for j in range(len(z[i])-1):
                # At each point
                if abs(z[i][j]-z[i][j+1])<=self.tol:
                    # If two adjacent points have same z
                    z[i][j+1]=z[i][j]+self.tol

        self.nst=len(x)
        self.x=np.array(x)
        self.y=list()
        self.z=list()
        self.c=list()
        self.xms=np.array(range(self.nst))
        self.npts=np.array(range(self.nst))

        for i in range(self.nst):
            self.y.append(np.array(y[i]))
            self.z.append(np.array(z[i]))
            self.c.append(np.array(c[i]))
            self.npts[i]=self.y[i].size
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        print("\n\t\t PyNAT : Preparing Section Curve Spline Fit Data")

        self.secCrv=list()
        for i in range(self.nst):    
            ind=np.argsort(self.z[i],axis=0)
            self.z[i]=np.take_along_axis(self.z[i],ind,axis=0)
            self.y[i]=np.take_along_axis(self.y[i],ind,axis=0)
            self.c[i]=np.take_along_axis(self.c[i],ind,axis=0)
            self.c[i][0]=-1
            self.c[i][-1]=-1
            self.secCrv.append(SectionCurves(self.x[i],self.y[i],self.z[i],self.c[i]))
        
        # Compute X locations form Mid Ship
        self.xms=self.x-self.MS

    # --------------------------------------------------------------------------------------- 

    # --------------------------------------------------------------------------------------
    # Function to write prepared hull geometry data to files
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def writeToFile(self, fs):

        print("\n\t\tPy-NAT: Writing Prepared Hull Geometry Data to File : ...\\",fs.ohg_path)

        if not os.path.exists(fs.sec_dir):
            os.makedirs(fs.sec_dir)

        fp=open(fs.ohg_path,'wt')
        fp.write("\nPy-NAT: Hull Geometry Data")
        # Printing the read data to Console
        fp.write("\nLength           =\t%12.4E"%(self.L))
        fp.write("\nBreadth          =\t%12.4E"%(self.B))
        fp.write("\nDepth            =\t%12.4E"%(self.D))
        fp.write("\nX location of AP =\t%12.4E"%(self.AP))
        fp.write("\nX location of FP =\t%12.4E"%(self.FP))
        fp.write("\nX location of MS =\t%12.4E"%(self.MS))
        fp.write("\nNo of Stations   =\t%12.4E"%(self.nst))
        fp.write("\n\n Printing Station Wise Data... \n\n")
        for i in range(self.nst):
            fp.write("\n-----------------------------------------------------------")
            fp.write("\nStation no       =\t%12.4E"%(i))
            fp.write("\nX                =\t%12.4E"%(self.x[i]))
            fp.write("\nX frm MS         =\t%12.4E"%(self.xms[i]))            
            fp.write("\nNo of pts        =\t%12.4E"%(self.npts[i]))
            fp.write("\n  i \t    Z        \t    Y       \t           C")
            for j in range(self.npts[i]):
                fp.write("\n%4d\t%12.4E\t%12.4E\t%12d"%(j,self.z[i][j],self.y[i][j],self.c[i][j]))

            fp.write("\n\n Total No of Curvelets  =%d"%(self.secCrv[i].nCrv))
            for k in range(self.secCrv[i].nCrv):

                fp.write("\n\n The Curve No: %d \t Start index= %10.5f \t End Index= %10.5f \n"%(k,self.secCrv[i].sCrv[k],self.secCrv[i].eCrv[k]))
                fp.write("\nThe Spline curve Fit Y=F(z)")
                tck=self.secCrv[i].y_fz[k]
                fp.write("\nKnots, t=")
                for l in range(len(tck[0])):
                    fp.write("%10.5f,"%(tck[0][l]))
                fp.write("\nB-Spline Coeff, c=")
                for l in range(len(tck[1])):
                    fp.write("%10.5f,"%(tck[1][l]))
                fp.write("\nDegree, k=")
                fp.write("%10.5f,"%(tck[2]))

                fp.write("\nThe Spline curve Fit ZY=F(z)")
                tck=self.secCrv[i].zy_fz[k]
                fp.write("\nKnots, t=")
                for l in range(len(tck[0])):
                    fp.write("%10.5f,"%(tck[0][l]))
                fp.write("\nB-Spline Coeff, c=")
                for l in range(len(tck[1])):
                    fp.write("%10.5f,"%(tck[1][l]))
                fp.write("\nDegree, k=")
                fp.write("%10.5f,"%(tck[2]))

                fp.write("\nThe Spline curve Fit YY=F(z)")
                tck=self.secCrv[i].yy_fz[k]
                fp.write("\nKnots, t=")
                for l in range(len(tck[0])):
                    fp.write("%10.5f,"%(tck[0][l]))
                fp.write("\nB-Spline Coeff, c=")
                for l in range(len(tck[1])):
                    fp.write("%10.5f,"%(tck[1][l]))
                fp.write("\nDegree, k=")
                fp.write("%10.5f,"%(tck[2]))

                fp.write("\nThe Spline curve Fit L=F(z)")
                tck=self.secCrv[i].L_fz[k]
                fp.write("\nKnots, t=")
                for l in range(len(tck[0])):
                    fp.write("%10.5f,"%(tck[0][l]))
                fp.write("\nB-Spline Coeff, c=")
                for l in range(len(tck[1])):
                    fp.write("%10.5f,"%(tck[1][l]))
                fp.write("\nDegree, k=")
                fp.write("%10.5f,"%(tck[2]))
        fp.close()

        #-------------------------------------------------------------------------------------------------
        # Activate this Section to Plot Each individual Section
        # Py-NAT 2021a : PRAVEEN KUMAR CH
        # for i in range(self.nst):
        #     plt.figure(i)
        #     zs,ys=self.secCrv[i].getPoints(0.1)
        #     plt.plot(ys,zs)
        #     plt.plot(self.y[i],self.z[i],'o')
        #     xstr = format(self.x[i],"6.3f")
        #     istr = format(i,"3d")
        #     plt.title("Section @ X="+xstr)
        #     plt.xlabel("Y")
        #     plt.ylabel("Z")
        #     plt.grid(True)
        #     plt.axis('equal')
        #     #lg=plt.legend(bbox_to_anchor=(1.05, 1))
        #     plt.savefig(os.path.join(fs.sec_dir,"Sec_Curv_"+istr+".png"),            
        #         dpi=300, 
        #         format='png', 
        #         bbox_extra_artists=(), 
        #         bbox_inches='tight')
        #     plt.close()
        #-------------------------------------------------------------------------------------------------

        # Plot Section Curves from the Read Points
        for i in range(self.nst):
            plt.figure("Body Plan")
            if(i<self.nst/2):
                zs,ys=self.secCrv[i].getPoints(0.1)
                hstr = format(i,"2d")
                plt.plot(ys,zs,label='Sec:'+hstr)
                plt.plot(self.y[i],self.z[i],'o')
            else:
                zs,ys=self.secCrv[i].getPoints(0.1)
                hstr = format(i,"2d")
                plt.plot(-ys,zs,label='Sec:'+hstr)
                plt.plot(-self.y[i],self.z[i],'o')
        plt.title("Body Plan")
        plt.xlabel("Y")
        plt.ylabel("Z")
        plt.grid(True)
        plt.axis('equal')
        lg=plt.legend(bbox_to_anchor=(1.05, 1))
        plt.savefig(os.path.join(fs.sec_dir,"Body_Plan.png"),            
            dpi=300, 
            format='png', 
            bbox_extra_artists=(lg,), 
            bbox_inches='tight')
    
    # --------------------------------------------------------------------------------------

    # --------------------------------------------------------------------------------------
    # Function to Print the data read from files
    # Py-NAT 2021a : PRAVEEN KUMAR CH
    # --------------------------------------------------------------------------------------
    def printToConsol(self, fs):

        # Printing the read data to Console
        print("Length           =\t%12.4E"%(self.L))
        print("Breadth          =\t%12.4E"%(self.B))
        print("Depth            =\t%12.4E"%(self.D))
        print("X location of AP =\t%12.4E"%(self.AP))
        print("X location of FP =\t%12.4E"%(self.FP))
        print("X location of MS =\t%12.4E"%(self.MS))
        print("No of Stations   =\t%12.4E"%(self.nst))
        for i in range(self.nst):
            print("-------------------------------------------------")
            print("Station no       =\t%12.4E"%(i))
            print("X                =\t%12.4E"%(self.x[i]))
            print("No of pts        =\t%12.4E"%(self.npts[i]))
            print("    Z        \t    Y       \t           C")
            for j in range(self.npts[i]):
                print("%12.4E\t%12.4E\t%12d"%(self.z[i][j],self.y[i][j],self.c[i][j]))
            for k in range(self.secCrv[i].nCrv):
                print("The Spline curve Fit Y=F(z)")
                print(self.secCrv[i].y_fz[k])
                print("The Spline curve Fit ZY=F(z)")
                print(self.secCrv[i].zy_fz[k])
                print("The Spline curve Fit YY=F(z)")
                print(self.secCrv[i].yy_fz[k])
                print("The Spline curve Fit L=F(z)")
                print(self.secCrv[i].L_fz[k])


        # PLot Section Curves from the Read Points
        plt.figure("Section Curves")
        for i in range(self.nst):
            if(i<self.nst/2):
                zs,ys=self.secCrv[i].getPoints(0.5)
                plt.plot(ys,zs,'g-')
                plt.plot(self.y[i],self.z[i],'bo')
            else:
                zs,ys=self.secCrv[i].getPoints(0.5)
                plt.plot(-ys,zs,'g-')
                plt.plot(-self.y[i],self.z[i],'bo')
        plt.title("The Section Curves")
        plt.xlabel("Y")
        plt.ylabel("Z")
        plt.grid(True)
        plt.savefig(os.path.join(fs.out_dir,"SectionCurves.png"))
    # --------------------------------------------------------------------------------------