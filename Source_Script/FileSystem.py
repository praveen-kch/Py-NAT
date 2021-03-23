"""
************************************************************************************
    File            : FileSystems.py
    
    File Purpose    : Define the Class "FileSystems"
    
    Class Purpose   : 
                    -> Stores all the File system related data
                    -> Directory names
                    -> File names
    
    App             : PyNAT - Python Based Naval Architecture Tools 
    Version         : 2021.0
    Paradigm        : Object Based Framework
    Author          : Praveen Kumar Ch
*************************************************************************************
"""
import shutil
import os

class FileSystem:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # --------------------------------------------------------------------------------------
    def __init__(self):
        
        # The Directories for input and output files
        self.in_dir = "in"
        self.out_dir = "out"
        self.sec_dnm = "sec"
        self.sec_dir =os.path.join(self.out_dir,self.sec_dnm)
        self.hs_dnm="hs"
        self.hs_dir =os.path.join(self.out_dir,self.hs_dnm)
        self.las_dnm="las"
        self.las_dir =os.path.join(self.out_dir,self.las_dnm)


        # Check if Output Directory exists and delete if exist
        # if(os.path.isdir(self.out_dir)):
        #         shutil.rmtree(self.out_dir)

        # Create the out directory
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        if not os.path.exists(self.sec_dir):
            os.makedirs(self.sec_dir)
        if not os.path.exists(self.hs_dir):
            os.makedirs(self.hs_dir)
        if not os.path.exists(self.las_dir):
            os.makedirs(self.las_dir)

        self.iep_file = "Set_EP.dat"
        self.iep_path = os.path.join(self.in_dir,self.iep_file)

        self.ihg_file = "Set_HG.dat"
        self.ihg_path = os.path.join(self.in_dir,self.ihg_file)

        self.ihs_file = "Set_HS.dat"
        self.ihs_path = os.path.join(self.in_dir,self.ihs_file)

        self.ikn_file = "Set_KN.dat"
        self.ikn_path = os.path.join(self.in_dir,self.ikn_file)

        self.ilc_file = "Set_LC.dat"
        self.ilc_path = os.path.join(self.in_dir,self.ilc_file)

        self.ohg_file = "Res_HG.dat"
        self.ohg_path = os.path.join(self.sec_dir,self.ohg_file)

        self.ohs_file = "Res_HS.dat"
        self.ohs_path = os.path.join(self.hs_dir,self.ohs_file)

        self.okn_file = "Res_KN.dat"
        self.okn_path = os.path.join(self.las_dir,self.okn_file)

        self.ogz_file = "Res_GZ.dat"
        self.ogz_path = os.path.join(self.las_dir,self.ogz_file)

        self.okn_fast_file = "Res_KN_FAST.dat"
        self.okn_fast_path = os.path.join(self.las_dir,self.okn_fast_file)

        self.ogz_fast_file = "Res_GZ_FAST.dat"
        self.ogz_fast_path = os.path.join(self.las_dir,self.ogz_fast_file)
