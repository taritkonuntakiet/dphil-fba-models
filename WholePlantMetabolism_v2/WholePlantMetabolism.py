from core import *
import numpy as np
import pandas as pd

#Initialise model
plantMainCopy = WholePlantModel()
plantMainCopy.initialisation()

#Simulate plant growth
simulationTime, leafBiomassOutput, stemBiomassOutput, rootBiomassOutput, fruitBiomassOutput = plantMainCopy.simulateGrowth() 

#Read model outputs to csv files
outputDict = dict()
outputDict["DAS"] = simulationTime
outputDict["Leaf_Dry_Weight"] = leafBiomassOutput
outputDict["Stem_Dry_Weight"] = stemBiomassOutput
outputDict["Root_Dry_Weight"] = rootBiomassOutput
outputDict["Fruit_Dry_Weight"] = fruitBiomassOutput
df = pd.DataFrame.from_dict(outputDict)
df.to_csv("Output/GrowthPredictions.csv") #Growth predictions


