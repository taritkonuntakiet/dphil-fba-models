from core import *
import numpy as np
import pandas as pd

#Initialise model
plantMainCopy = WholePlantModel()
plantMainCopy.initialisation()

#Simulate plant growth
plantMainCopy.outputCycle = np.arange(10,100,5) #Use "outputCycle" to specify days at which to record flux predictions for model output (in DAS)
simulationTime, leafBiomassOutput, rootBiomassOutput, stemBiomassOutput, fluxes = plantMainCopy.simulateGrowth() 

#Read model outputs to csv files
outputDict = dict()
outputDict["DAS"] = simulationTime
outputDict["Leaf_Dry_Weight"] = leafBiomassOutput
outputDict["Stem_Dry_Weight"] = stemBiomassOutput
outputDict["Root_Dry_Weight"] = rootBiomassOutput
df = pd.DataFrame.from_dict(outputDict)
df.to_csv("Output/GrowthPredictions.csv") #Growth predictions

columns = [reaction.id for reaction in plantMainCopy.stoichiometricModel.reactions]
index = [d.split("solution")[1] for d in fluxes.keys()]
fluxArray = np.array([np.arange(len(fluxes))]*len(columns), dtype=np.float32).T
for das in range(len(index)):
    for reac in range(len(columns)):
        fluxArray[das][reac] = fluxes["solution"+str(index[das])][columns[reac]]
df2 = pd.DataFrame(fluxArray, index = index, columns=columns)
df2.index.names = ['DAS']
df2.to_csv("Output/Fluxes.csv") #Metabolic flux predictions