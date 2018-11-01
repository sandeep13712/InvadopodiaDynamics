from PlayerPython import * 
import CompuCellSetup

from PySteppables import *
import CompuCell
import sys
import random
from math import exp 
from random import uniform

class InitializeMMPAndECM(SteppableBasePy):
    def start(self):   
       
        random.seed();
        #initialize MMP
        mmpField = self.getConcentrationField("MMP")
        for x,y,z in self.everyPixel():
            mmpField[x,y,z] = 0.0
        
        #initialize ECM
        matrixField = self.getConcentrationField("Matrix")
        for x,y,z in self.everyPixel():
            r = uniform(0,1)
            if(r <= 0.0):
                matrixField[x,y,z]=1
                
class InvadopodiaSteppable(SteppableBasePy):    
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        # Setup invadopodia
        widthOfInvadopodia = 5
        totalInvadopodia = 2
        
#Case1: For 2 invadopodia 
    # alpha= 5 pxls, {(8,8), (14,8)} 
    # alpha= 7 pxls, {(7,8), (15,8)} 
    # alpha= 10 pxls, {(6,8), (17,8)}
    # alpha= 15 pxls, {(3,8), (19,8)}
    
        xMinListDiagonal = [8, 14]
        yMinListDiagonal = [8, 8]

        for i in range(0,totalInvadopodia): 
        
            xMin = xMinListDiagonal[i]
            yMin = yMinListDiagonal[i]
            zMin = 0
            
            invadopodia = self.newCell(1)
            self.cellField[xMin:xMin+widthOfInvadopodia-1,yMin:yMin+widthOfInvadopodia-1,zMin:zMin+1] = invadopodia
            
            invadopodia.dict["zMax"]=1
            invadopodia.dict["xMin"]=xMin
            invadopodia.dict["yMin"]=yMin
        
    def step(self,mcs):
        
        random.seed();
        PPEN_THREASHOLD = 1.0
        
        MT1MMP_THREASHOLD = 2.5
        MT1MMP_BASAL = 0.005
        MT1MMP_MAXVAL = 0.05
        
        totalCollagenAtStartOfMCS = 0
        matrixField = self.getConcentrationField("Matrix")
        for x,y,z in self.everyPixel():
            totalCollagenAtStartOfMCS = totalCollagenAtStartOfMCS + round(matrixField[x,y,z],6);
            
        widthOfInvadopodia = 5
        
        MT1MMPActivity = 0
        totalECMDegradedBYMt1MMP = 0
        for invadopodia in self.cellListByType(1):

            invadopodiaFront = invadopodia.dict["zMax"]
            xMin = invadopodia.dict["xMin"]
            yMin = invadopodia.dict["yMin"]
            
            if(invadopodiaFront >= 20):
                self.stopSimulation()
                
            newInvadopodiaFront = invadopodiaFront+1;
            
            # Calculate total ECM below invadopodia
            totalECMBelowInvadopodia = 0
            for x in range(xMin, xMin+widthOfInvadopodia):
                for y in range(yMin, yMin+widthOfInvadopodia):
                    totalECMBelowInvadopodia = totalECMBelowInvadopodia + round(matrixField[x,y,newInvadopodiaFront],6)
            
            # if no fiber below invadopodia
            if(totalECMBelowInvadopodia == 0):
                # Find mechanical support index
                totalLateralECMDensity = 0.0
                pixelList=self.getCellBoundaryPixelList(invadopodia)
                totalLateralPoint = 0
                for boundaryPixelTrackerData in pixelList:
                    x = boundaryPixelTrackerData.pixel.x
                    y = boundaryPixelTrackerData.pixel.y
                    z = boundaryPixelTrackerData.pixel.z
                    if(self.cellField[x-1,y,z] and not(self.cellField[x+1,y,z])):
                        totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x+1,y,z],6)
                        totalLateralPoint = totalLateralPoint+1
                        
                    if(self.cellField[x+1,y,z] and not(self.cellField[x-1,y,z])):
                        totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x-1,y,z],6)
                        totalLateralPoint = totalLateralPoint+1
                        
                    if(self.cellField[x,y-1,z] and not(self.cellField[x,y+1,z])):
                        totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x,y+1,z],6)
                        totalLateralPoint = totalLateralPoint+1
                        
                    if(self.cellField[x,y+1,z] and not(self.cellField[x,y-1,z])):
                        totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x,y-1,z],6)
                        totalLateralPoint = totalLateralPoint+1
                        
                fractionCovered = totalLateralECMDensity/totalLateralPoint
                prob_Pene = 1.58*(1-exp(-fractionCovered))
                fileName='PPen_New.csv'
                try:
                    fileHandle,fullFileName = self.openFileInSimulationOutputDirectory(fileName,"a")
                except IOError:
                    # print "Could not open file ", fileName," for writing. "
                    return
                print mcs,",",fractionCovered,",",prob_Pene
                print>>fileHandle, mcs,",",fractionCovered,",",prob_Pene
                fileHandle.close()
                
                r = uniform(0,1);
                if(r < prob_Pene):
                    invadopodia.dict["zMax"] = invadopodia.dict["zMax"]+1
                    self.cellField[xMin:xMin+widthOfInvadopodia-1,yMin:yMin+widthOfInvadopodia-1,invadopodiaFront:invadopodiaFront+1] = invadopodia
            
            # if fibers are below invadopodia
            else:    
                # Logic to implement MT1MMP mediated ECM degradation
                MT1MMPActivity = MT1MMP_BASAL + (MT1MMP_MAXVAL-MT1MMP_BASAL)*(1+MT1MMP_THREASHOLD/25.0)*(totalECMBelowInvadopodia/(totalECMBelowInvadopodia+MT1MMP_THREASHOLD))                         
                for x in range(xMin, xMin+widthOfInvadopodia):
                    for y in range(yMin, yMin+widthOfInvadopodia):
                        ecmAtPixel = round(matrixField[x,y,newInvadopodiaFront],6)
                        if(MT1MMPActivity >= ecmAtPixel):
                            matrixField[x,y,newInvadopodiaFront] = 0.0;
                            totalECMDegradedBYMt1MMP = totalECMDegradedBYMt1MMP + ecmAtPixel
                        else:
                            matrixField[x,y,newInvadopodiaFront] = ecmAtPixel - MT1MMPActivity
                            totalECMDegradedBYMt1MMP = totalECMDegradedBYMt1MMP + MT1MMPActivity
            updatedInvadopodiaFront = invadopodia.dict["zMax"]
            
            #Clear the MMP field in invadopodia region
            mmpField = self.getConcentrationField("MMP")
            for x in range(xMin, xMin+widthOfInvadopodia):
                for y in range(yMin, yMin+widthOfInvadopodia):
                    for z in range(1,updatedInvadopodiaFront):
                        mmpField[x,y,z] = 0

class SolubleMMPSecretionSteppable(SteppableBasePy):
    def step(self,mcs):
        SOLUBLE_MMP_SECRETION_RATE = 0.01
        matrixField = self.getConcentrationField("Matrix")
        mmpField = self.getConcentrationField("MMP")
        
        totalMMPSecretedString = ""
        for invadopodiaCell in self.cellListByType(1):
            # Calculate total ECM in the surrounding
            totalLateralECMDensity = 0.0
            pixelList=self.getCellBoundaryPixelList(invadopodiaCell)
            totalLateralPoint = 0
            for boundaryPixelTrackerData in pixelList:
                x = boundaryPixelTrackerData.pixel.x
                y = boundaryPixelTrackerData.pixel.y
                z = boundaryPixelTrackerData.pixel.z
                if(self.cellField[x-1,y,z] and not(self.cellField[x+1,y,z])):
                    totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x+1,y,z],6)
                    totalLateralPoint = totalLateralPoint+1
                    
                if(self.cellField[x+1,y,z] and not(self.cellField[x-1,y,z])):
                    totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x-1,y,z],6)
                    totalLateralPoint = totalLateralPoint+1
                    
                if(self.cellField[x,y-1,z] and not(self.cellField[x,y+1,z])):
                    totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x,y+1,z],6)
                    totalLateralPoint = totalLateralPoint+1
                    
                if(self.cellField[x,y+1,z] and not(self.cellField[x,y-1,z])):
                    totalLateralECMDensity = totalLateralECMDensity +round(matrixField[x,y-1,z],6)
                    totalLateralPoint = totalLateralPoint+1
                    
            mmpToSecrete = round(totalLateralECMDensity*SOLUBLE_MMP_SECRETION_RATE)
            probabilityOfPlacingMMP = mmpToSecrete/totalLateralPoint
            
            
            # Secrete MMP
            totalMMPSecreted = 0
            pixelList=self.getCellBoundaryPixelList(invadopodiaCell)
            for boundaryPixelTrackerData in pixelList:
                ## print "I am to working"
                x = boundaryPixelTrackerData.pixel.x
                y = boundaryPixelTrackerData.pixel.y
                z = boundaryPixelTrackerData.pixel.z
                
                if(self.cellField[x-1,y,z] and not(self.cellField[x+1,y,z])):
                    r = uniform(0,1);
                    if(r < probabilityOfPlacingMMP):
                    # # print "R:", round(matrixField[x+1,y,z],6)," ", 
                        mmpField[x+1,y,z] = round(mmpField[x+1,y,z],6)+ 1.0
                        totalMMPSecreted = totalMMPSecreted + 1.0
        
                if(self.cellField[x+1,y,z] and not(self.cellField[x-1,y,z])):
                    # print "L:", round(matrixField[x-1,y,z],6)
                    r = uniform(0,1);
                    if(r < probabilityOfPlacingMMP):
                        mmpField[x-1,y,z] = round(mmpField[x-1,y,z],6)+1.0
                        totalMMPSecreted = totalMMPSecreted + 1.0
                    
                if(self.cellField[x,y-1,z] and not(self.cellField[x,y+1,z])):
                    # # print "F:",round(matrixField[x,y+1,z],6)
                    r = uniform(0,1);
                    if(r < probabilityOfPlacingMMP):
                        mmpField[x,y+1,z] = round(mmpField[x,y+1,z],6)+1.0
                        totalMMPSecreted = totalMMPSecreted + 1.0
                    
                if(self.cellField[x,y+1,z] and not(self.cellField[x,y-1,z])):
                    # print "B:",round(matrixField[x,y-1,z],6)
                    r = uniform(0,1);
                    if(r < probabilityOfPlacingMMP):
                        mmpField[x,y-1,z] = round(mmpField[x,y-1,z],6)+ 1.0
                        totalMMPSecreted = totalMMPSecreted + 1.0
                    
                if(self.cellField[x,y,z-1] and not(self.cellField[x,y,z+1])):
                    # print "D:",round(matrixField[x,y,z+1],6)
                    r = uniform(0,1);
                    if(r < probabilityOfPlacingMMP):
                        mmpField[x,y,z+1] = round(mmpField[x,y,z+1],6)+1.0
                        totalMMPSecreted = totalMMPSecreted + 1.0
                        
            totalMMPSecretedString = totalMMPSecretedString + ","+ str(totalMMPSecreted)
            
        fileName='SolubleMMPSecretion.csv'
        try:
                fileHandle,fullFileName = self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
                # print "Could not open file ", fileName," for writing. "
                return
        print mcs,",",totalMMPSecreted
        ## print "MaxZ = ", maxz
        print>>fileHandle, mcs,totalMMPSecretedString
        fileHandle.close()
        
class ECMDegradation(SteppableBasePy):
    def start(self):
        pass
    def step(self,mcs):
        
        
        matrixField = self.getConcentrationField("Matrix")
        mmpField = self.getConcentrationField("MMP")
        
        totalECMDegradadedBySolubleMMP = 0
        for x,y,z in self.everyPixel():
            if(mmpField[x,y,z]>= matrixField[x,y,z]):
                mmpField[x,y,z] = round(mmpField[x,y,z],6)- round(matrixField[x,y,z],6)
                totalECMDegradadedBySolubleMMP = totalECMDegradadedBySolubleMMP+round(matrixField[x,y,z],6)
                matrixField[x,y,z] = 0.0
            else:
                matrixField[x,y,z] = round(matrixField[x,y,z],6) - round(mmpField[x,y,z],6)
                mmpField[x,y,z] = 0.0
                
        fileName='ECMDegradation.csv'
        try:
                fileHandle,fullFileName = self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
                # print "Could not open file ", fileName," for writing. "
                return
        
        print mcs,",",totalECMDegradadedBySolubleMMP
        print>>fileHandle, mcs,",",totalECMDegradadedBySolubleMMP

        fileHandle.close()
                   
class LogDataSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        totalInvadopodia = 2
        
        # Quantify total collagen concentration at the end of the MCS
        matrixField = self.getConcentrationField("Matrix")
        totalCollagenConcentration = 0
        for x,y,z in self.everyPixel():
            totalCollagenConcentration = totalCollagenConcentration + round(matrixField[x,y,z],6);
        
        totalMMPAtMCS = 0
        
        updatedInvadopodiaFront = ""
        
        for i in (1,totalInvadopodia+1):
            updatedInvadopodiaFront = updatedInvadopodiaFront +",0"
            
        # Write data to file
        fileName='TotalCollagenAndMMP.csv'
        try:
                fileHandle,fullFileName = self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
                # print "Could not open file ", fileName," for writing. "
                return
        
        print 0,",",totalCollagenConcentration,",",totalMMPAtMCS,updatedInvadopodiaFront
        print>>fileHandle, 0,",",totalCollagenConcentration,",",totalMMPAtMCS,updatedInvadopodiaFront

        fileHandle.close()
        
    def step(self,mcs):
        totalInvadopodia = 2
        
        # Quantify total numnber of soluble MMPs
        mmpField = self.getConcentrationField("MMP")
        totalPixelsHavingMMP = 0
        totalMMPAtMCS = 0.0
        for x,y,z in self.everyPixel():
            totalMMPAtMCS = totalMMPAtMCS+round(mmpField[x,y,z],6) 
            if(mmpField[x,y,z]>0):
                totalPixelsHavingMMP = totalPixelsHavingMMP+1
        
        # Quantify total collagen concentration at the end of the MCS
        matrixField = self.getConcentrationField("Matrix")
        totalCollagenConcentration = 0
        for x,y,z in self.everyPixel():
            totalCollagenConcentration = totalCollagenConcentration + round(matrixField[x,y,z],6);
        
        # Get current invadopodia position
        updatedInvadopodiaFront = ""
        for invadopodia in self.cellListByType(1):
            updatedInvadopodiaFront = updatedInvadopodiaFront + "," + str(invadopodia.dict["zMax"])
     
        # Write data to file
        fileName='TotalCollagenAndMMP.csv'
        try:
                fileHandle,fullFileName = self.openFileInSimulationOutputDirectory(fileName,"a")
        except IOError:
                # print "Could not open file ", fileName," for writing. "
                return
        
        print mcs,",",totalCollagenConcentration,",",totalMMPAtMCS,updatedInvadopodiaFront,",",totalPixelsHavingMMP
        print>>fileHandle, mcs,",",totalCollagenConcentration,",",totalMMPAtMCS,updatedInvadopodiaFront,",",totalPixelsHavingMMP

        fileHandle.close()
        
class ExtraFieldVisualizationSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        self.scalarCellIDField=self.createScalarFieldCellLevelPy("CellID")
    def step(self,mcs):
        self.scalarCellIDField.clear()
        for cell in self.cellList:
            self.scalarCellIDField[cell]=cell.type
