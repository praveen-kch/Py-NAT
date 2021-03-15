"""
************************************************************************************
    File            : ENV_PROP.py
    
    File Purpose    : Define the Class "ENV_PROP"
    
    Class Purpose   : 
                    -> Reads and stores Environment properties from file "Set_EP.dat"
                    -> Environment Properties are defined in this class
    
    App             : PyNAT - Python Based Naval Architecture Tools 
    Version         : 2021.0
    Paradigm        : Object Based Framework
    Author          : Praveen Kumar Ch
*************************************************************************************
"""

class ENV_PROP:

    # --------------------------------------------------------------------------------------
    # The Constructor
    # --------------------------------------------------------------------------------------
    def __init__(self,fs):
        
        print("\n\t\tPyNAT: Reading File:\t...\\",fs.iep_path)

        # Reading the Solver properties
        fp=open(fs.iep_path,"rt")
        lines=fp.readlines()
        fp.close()

        nums=lines[0].split('\t')
        self.waterDensity = float(nums[0])