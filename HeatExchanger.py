# -*- coding: utf-8 -*-
"""
Created on Tue May 17 09:28:12 2016

@author: Group 1
Project link here: http://apmonitor.com/che263/uploads/Main/heat_exchanger_project.pdf
"""
class HeatExchanger():
    # Variables that are given in getInput() and then used in another function are declared here
    units=""    # aes or si
    U=0         # Heat Transfer Coefficient
    Tci=0       # T cold fluid, in
    Thi=0       # T hot fluid, in
    mc=0        # flow rate cold fluid
    mh=0        # flow rate hot fluid
    cold=""     # options given below in input
    hot=""      # options given below in input
    given=""    # Tco or Tho

# Marcus    
    def runAll():
    # Source functions, runs all the other functions
    
    
    ### Pseduocode shell
# Gabriel
    def getInput():
        # All variables are set here.
        # Note that cp is the average of inlet and outlet temperatures
    """user can specify:
            aes or si units
            U, Tci, Thi, mc, mh, Tco or Tho
            Cold: water or 1,1,1,2-Tetrafluoroethane or ethanol or 2,2,4-trimethylpentane
            Hot: same options as cold
    """
    
    def isValid(value):
        # checks if a number is negative. If so, returns error message and exits program.    
        
    def convertUnits(value):
        # convert AES to SI units
# Dan the man     
    
    #creat a class for table values
    
    import numpy as np
import pandas as pd

#Importing the File

thermo_file = pd.read_csv('Thermophysical_Properties.csv')



#Assigning Variables to the Columns

compound_name = thermo_file['Compound Name']
mol_wt = thermo_file['Molecular Weight']
melt_pt = thermo_file['Melting Point (K)']
boil_pt = thermo_file['Normal Boiling Point (K)']
cp_A = thermo_file['A']
cp_B = thermo_file['B']
cp_C = thermo_file['C']
cp_D = thermo_file['D']
cp_E = thermo_file['E']


#Creating the Class

class Fluid:
    name = []
    def c_name(self, name):
        self.name.append(name)
        
    Mol_Wt = []
    def c_mw(self, Mol_Wt):
        self.name.append(Mol_Wt)
    
    MP = []
    def c_mp(self, MP):
        self.name.append(MP)
        
    BP = []
    def c_bp(self, BP):
        self.name.append(BP)
        
    A = []
    def c_a(self, A):
        self.name.append(A)
        
    B = []
    def c_b(self, B):
        self.name.append(B)
        
    C = []
    def c_c(self, C):
        self.name.append(C)
        
    D = []
    def c_d(self, D):
        self.name.append(D)
        
    E = []
    def c_e(self, E):
        self.name.append(E)
        
    
C = np.empty(4)
C0 = Fluid()
C1 = Fluid()
C2 = Fluid()
C3 = Fluid()


C0.name = compound_name[0]
C0.Mol_Wt = mol_wt[0]
C0.MP = melt_pt[0]
C0.BP = boil_pt[0]
C0.A = cp_A[0]
C0.B = cp_B[0]
C0.C = cp_C[0]
C0.D = cp_D[0]
C0.E = cp_E[0]
    
C1.name = compound_name[1]
C1.Mol_Wt = mol_wt[1]
C1.MP = melt_pt[1]
C1.BP = boil_pt[1]
C1.A = cp_A[1]
C1.B = cp_B[1]
C1.C = cp_C[1]
C1.D = cp_D[1]
C1.E = cp_E[1]
    
C2.name = compound_name[2]
C2.Mol_Wt = mol_wt[2]
C2.MP = melt_pt[2]
C2.BP = boil_pt[2]
C2.A = cp_A[2]
C2.B = cp_B[2]
C2.C = cp_C[2]
C2.D = cp_D[2]
C2.E = cp_E[2]
    
C3.name = compound_name[3]
C3.Mol_Wt = mol_wt[3]
C3.MP = melt_pt[3]
C3.BP = boil_pt[3]
C3.A = cp_A[3]
C3.B = cp_B[3]
C3.C = cp_C[3]
C3.D = cp_D[3]
C3.E = cp_E[3]
    
    
    
    
    
    
    def solveT():
        # calculates Tco or Tho (whichever one was NOT given in getInput() )
        
        #import numpy
        import numpy as np
        
        #To do!! define mh, mc, Cph, Cpc
        
        #try importing data
            #from xlwings import Workbook, Range
            #from pandas import DataFrame
        
        
        
        
        
        #Solve for Tco or Tho
        if getInput() == "Tho":
            Tco = mh*Cph*(Thi-Tho)/(mc*Cpc) + Tci
            return Tco
        else:
            Tho = -mc*Cpc(Tco-Tci)/(mh*Cph) + Thi
            return Tho
        
    #Solve for Change in Temp 1 and 2
    dT1 = Thi - Tco
    dT2 = Tho - Tci
        
        
    #Solve for Heat Transfer Rate 'q'
    q = mh*Cph*(Tco-Tci)
        
    #Solve for the Surface Area
    def solveArea():
        
        #Correction Factor 'F' as a function of 'R' and 'P'
        R = (Thi-Tho)/(Tco-Tci)
        P = (Tco-Tci)/(Thi-Tci)
        
        F = (np.sqrt(R**2 +1)/(R-1)) * \
            np.log((1-P)/(1-P*R))/ \
            np.log((2-P*(R+1-np.sqrt(R**2 +1)))/(2-P*(R+1+np.sqrt(R**2 +1))))
        
        #Log Mean Temperature Difference
        T_logmean = (dT1 - dT2)/np.log(dT2/dT1)
        
        #Caluclate and Return Area
        area = q/(F*U*T_logmean) #m^2
        return (area)
            
    #Solve for the Cost
    def solveCost():
        # cost = $1000 * area (m^2)
        cost = 1000*area  # $ in USD
        return cost
        
        
# Marcus        
    def output(area, T, cost):
        # Displays solutions all in SI units
        
