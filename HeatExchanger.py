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
    
    def solveT():
        # calculates Tco or Tho (whichever one was NOT given in getInput() )
        
        #import numpy
        import numpy as np
        
        #To do!! define mh, mc, Cph, Cpc
        
        
        
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
        
