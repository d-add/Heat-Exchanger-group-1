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
    
import numpy as np
import pandas as pd
import scipy.optimize as opt
import sympy as sp

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
        
    
C = [0,1,2,3]
i = 0
for i in C:
    C[i] = Fluid()
    C[i].name = compound_name[i]
    C[i].Mol_Wt = mol_wt[i]
    C[i].MP = melt_pt[i]
    C[i].BP = boil_pt[i]
    C[i].A = cp_A[i]
    C[i].B = cp_B[i]
    C[i].C = cp_C[i]
    C[i].D = cp_D[i]
    C[i].E = cp_E[i]
    
    i += 1
    
    
#Calculate Heat Capacity determined by Temperature out (hot or cold) that was selceted at getInput()
if getInput() == "Tho":
    Th_avg = (Thi+Tho)/2
    def Cph(j,Th_avg): #j = compound number o hot fluid, and T is change in temperature in Kelvin
        heat_capacity_h = C[j].A + C[j].B*Th_avg + C[j].C*Th_avg**2 + C[j].D*Th_avg**3 + C[j].E*Th_avg**4
        return heat_capacity_h
else:
    Tc_avg = (Tci + Tco)/2
    def Cpc(k,Tc_avg): #k = compound number of cold fluid, and T is change in temperature in Kelvin
        heat_capacity_c = C[k].A + C[k].B*Tc_avg + C[k].C*Tc_avg**2 + C[k].D*Tc_avg**3 + C[k].E*Tc_avg**4
        return heat_capacity_c
    


if getInput() == "Tho":
    #Solve for Heat Transfer Rate 'q'
    q_h = mh*Cph(j,Th_avg)*(Thi-Tho)  #q_c = mc*Cpc(k,Tc_avg)*(Tco-Tci)    note that q_h = q_c
   

    #solve polynomial equation with one variable (variable = Tc_avg) by setting equation equal to zero
    #0 = 2*mc*(C[k].A + C[k].B*Tc_avg + C[k].C*Tc_avg**2 + C[k].D*Tc_avg**3 + C[k].E*Tc_avg**4)*(Tc_avg - Tci) - q_h
    Tc_avg = sp.Symbol('Tc_avg')
    y = sp.solve(2*mc*(C[k].A + C[k].B*Tc_avg + C[k].C*Tc_avg**2 + C[k].D*Tc_avg**3 + C[k].E*Tc_avg**4)*(Tc_avg - Tci) - q_h, Tc_avg)
    Tc_avg = y[0]

else:
    #Solve for Heat Transfer Rate 'q'
    q_c = mc*Cpc(k,Tc_avg)*(Tco-Tci) #q_h = mh*Cph(j,Th_avg)*(Thi-Tho    note that q_h = q_c

    #solve polynomial equation with one variable (variable = Th_avg) by setting equation equal to zero
    #0 = 2*mc*(C[j].A + C[j].B*Th_avg + C[j].C*Th_avg**2 + C[j].D*Th_avg**3 + C[j].E*Th_avg**4)*(Th_avg - Tci) - q_h
    Th_avg = sp.Symbol('Th_avg')
    z = sp.solve(2*mh*(C[j].A + C[j].B*Th_avg + C[j].C*Th_avg**2 + C[j].D*Th_avg**3 + C[j].E*Th_avg**4)*(Thi-Th_avg) - q_c, Th_avg)
    Th_avg = z[0]





#Solve for Cph or Cpc (heat capacity). Whichever hasn't already been solved.
if getInput() == "Tho":
    def Cpc(k,Tc_avg): #k = compound number of cold fluid, and T is change in temperature in Kelvin
        heat_capacity_c = C[k].A + C[k].B*Tc_avg + C[k].C*Tc_avg**2 + C[k].D*Tc_avg**3 + C[k].E*Tc_avg**4
        return heat_capacity_c
    
else:
    def Cph(j,Th_avg): #j = compound number of hot fluid, and T is change in temperature in Kelvin
        heat_capacity_h = C[j].A + C[j].B*Th_avg + C[j].C*Th_avg**2 + C[j].D*Th_avg**3 + C[j].E*Th_avg**4
        return heat_capacity_h
    
    
        
# calculate Tco or Tho (whichever one was NOT given in getInput() )
        
#Solve for Tco or Tho
if getInput() == "Tho":
    Tco = mh*Cph(j,Th_avg)*(Thi-Tho)/(mc*Cpc(k,Tc_avg)) + Tci
    
else:
    Tho = -mc*Cpc(k,Tc_avg)*(Tco-Tci)/(mh*Cph(j,Th_avg)) + Thi
    


    


    

        
#Solve for the Surface Area
def solveArea():
    #Correction Factor 'F' as a function of 'R' and 'P'
    R = (Thi-Tho)/(Tco-Tci)
    P = (Tco-Tci)/(Thi-Tci)
    
    F = (((R**2 +1.)**(0.5))/(R-1)) * ((np.log((1.-P)/(1.-P*R)))/ np.log((2.-P*(R+1.-((R**2 +1.)**(0.5))))/(2.-P*(R+1+((R**2 +1.)**(0.5))))))
    print(F)

    #Log Mean Temperature Difference
    dT1 = Thi - Tco
    dT2 = Tho - Tci
    T_logmean = (dT1 - dT2)/np.log(dT2/dT1)
        
    #Caluclate and Return Area
    area = q_h/(F*U*T_logmean) #m^2
    return (area)
            
#Solve for the Cost
def solveCost(area):
    # cost = $1000 * area (m^2)
    cost = 1000*area  # $ in USD
    return cost


    
        
# Marcus        
    def output(area, T, cost):
        # Displays solutions all in SI units
        
