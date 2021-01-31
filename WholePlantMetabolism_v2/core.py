class WholePlantModel:

    from cobra.core import DictList

    """Create whole plant model with leaves, stems, fruits, and root"""

    def __init__(self,id = "tomatoPlant"):

        from cobra.core import DictList
        from cobra.core import Model

        self.id = id
        self.growthCycle = 1
        self.leaves = DictList()
        self.stems = DictList()
        self.fruits = DictList()
        self.root = DictList()
        self.parameters = {}

    def addRoot(self, newOrgan):
        self.root += newOrgan

    def addLeaf(self, newOrgan):
        self.leaves += newOrgan

    def addFruit(self, newOrgan):
        self.fruits += newOrgan

    def addStem(self, newOrgan):
        self.stems += newOrgan

    def loadParams(self, datafile = "./Parameters/parameters.csv"):
        """
        Function to load parameters from CSV file

        Parameters
        ----------
        - datafile : comma-seperated CSV file
            A CSV file with parameters

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import pandas as pd
        import numpy as np
        import math

        df = pd.read_csv(datafile, delimiter=",")

        self.parameters["leafMaint"] = df["leafMaint"][0] #leaf maintenance requirements mmol/gDW/day
        self.parameters["rootMaint"] = df["rootMaint"][0] #root maintenance requirements mmol/gDW/day
        self.parameters["stemMaint"] = df["stemMaint"][0] #stem maintenance requirements mmol/gDW/day 
        self.parameters["fruitMaint"] = df["fruitMaint"][0] #fruit maintenance requirements mmol/gDW/day
        self.parameters["NO3UptakeVm"] = df["NO3UptakeVm"][0]
        self.parameters["NO3UptakeKm"] = df["NO3UptakeKm"][0]
        self.parameters["NO3UptakeVmLow"] = df["NO3UptakeVmLow"][0]
        self.parameters["NO3UptakeKmLow"] = df["NO3UptakeKmLow"][0]
        self.parameters["nitrateConc"] = df["nitrateConc"][0] #Soil/growth medium nitrate concentration mM
        self.parameters["Pmax0"] = df["Pmax0"][0] #umol/m2/s
        self.parameters["RbInit"] = df["RbInit"][0] #initial root mass
        self.parameters["LbInit"] = df["LbInit"][0] #initial leaf mass
        self.parameters["SbInit"] = df["SbInit"][0] #initial stem mass
        self.parameters["temperatureThrsh"] = df["temperatureThrsh"][0]
        self.parameters["temperatureAverageDaily"] = df["temperatureAverageDaily"][0]
        self.parameters["temperatureBaseline"] = df["temperatureBaseline"][0]
        self.parameters["photoperiod"] = df["photoperiod"][0] #photoperiod in hours
        self.parameters["leafAreatoBiomassVeg"] = df["LeafAreatoBiomassVeg"][0] #For conversion between leaf area and mass during vegetative growth
        self.parameters["leafAreatoBiomassRep"] = df["LeafAreatoBiomassRep"][0] #For conversion between leaf area and mass during reproductive growth
        self.gL = df["gL"][0] #Scaling exponent to leaf relative to root
        self.gR = df["gR"][0] #Scaling exponent to root relative to root
        self.gS = df["gS"][0] #Scaling exponent to stem relative to root
        self.parameters["stemGrowthDuration"] = df["stemGrowthDuration"][0]
        self.parameters["leafVarA"] = df["leafVarA"][0]
        self.parameters["leafVarB"] = df["leafVarB"][0]
        self.parameters["leafVarC"] = df["leafVarC"][0]
        self.parameters["stemVarA"] = df["stemVarA"][0]
        self.parameters["stemVarB"] = df["stemVarB"][0]
        self.parameters["stemVarC"] = df["stemVarC"][0]
        self.parameters["fruitVarA"] = df["fruitVarA"][0]
        self.parameters["fruitVarB"] = df["fruitVarB"][0]
        self.parameters["fruitVarC"] = df["fruitVarC"][0]
        self.parameters["m"] = df["m"][0]
        self.parameters["k"] = df["k"][0]
        self.parameters["CotArea"] = df["CotArea"][0] #Fully expanded cotyledon area
        self.parameters["CotAge"] = df["CotAge"][0] #Lifespan of cotyledons in DAS
        self.parameters["InitDAS"] = df["InitDAS"][0] #Start of daily simulations in DAS  
        self.parameters["floweringTime"] = df["floweringTime"][0] #Time from reproductive phytomer appearance to anthesis
        self.parameters["fruitMassAtAnthesis"] = df["fruitMassAtAnthesis"][0] #g DW per fruit at anthesis
        self.parameters["fruitRipenedAge"] = df["fruitRipenedAge"][0] #Age at which fruit reaches ripened stage
        self.parameters["fruitNumber"] = dict() #Dictionary of number of fruits per truss (phytomer rank as key)
        self.parameters["fruitNumber"][10] = df["fruitNumberPhyt10"][0]
        self.parameters["fruitNumber"][13] = df["fruitNumberPhyt13"][0]
        
        #Calculations
        self.parameters["PhyllochronLength"] = float(self.parameters["temperatureThrsh"]/(self.parameters["temperatureAverageDaily"]-self.parameters["temperatureBaseline"]))
        self.cycleThreshold = self.parameters["PhyllochronLength"]
        self.parameters["stemGrowthDuration"] = self.parameters["stemGrowthDuration"]*self.parameters["PhyllochronLength"]
        self.parameters["leafLifeSpan"] = df["leafLifeSpan"][0]*self.parameters["PhyllochronLength"]
        self.parameters["leafGrowthDuration"] = df["leafGrowthDuration"][0]*self.parameters["PhyllochronLength"]

        self.parameters["Pmax"] = (self.parameters["Pmax0"]*60*60*self.parameters["photoperiod"]*self.parameters["leafAreatoBiomassVeg"])/1000 #mmol/gDW/day

        self.parameters["CotWeightAdj"] = self.parameters["CotArea"]/self.parameters["leafAreatoBiomassVeg"]
        
    def initialisation(self):
        """
        Function to initialise model

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """
        
        self.loadParams()
        self.generateStoichiometricModels()
        self.calculateBiomassWeight()
        self.generateModelAtEmergence()
        self.estimateVcFromPPFD()
        self.readStrucutrefromMask()
        
    def simulateGrowth(self, simulationTime = range(1,89)):
        """
        Function to simulate growth

        Parameters
        ----------
        - simulationTime : list specifying the length of growth simulations

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """
        
        import cobra
        from cobra.flux_analysis import pfba
        from cobra import Reaction, Metabolite
        import numpy as np

        leafBiomassOutput = [self.parameters["LbInit"]]
        stemBiomassOutput = [self.parameters["SbInit"]]
        rootBiomassOutput = [self.parameters["RbInit"]]
        fruitBiomassOutput =[0] 
        fluxOutputs = dict()
        biomassOutputs = dict()

        for t in simulationTime:

            #Update age of each organ in days
            for tissue in [self.leaves, self.root, self.stems, self.fruits]:
                for organ in tissue:
                    organ.age = organ.age + 1

            #Update organ developmental stage
            for leaf in self.leaves:
                if leaf.age >= self.parameters["leafLifeSpan"]*self.parameters["leafExpScale"][leaf.position-1]:
                    leaf.stage = 1 #Leaf senescence
            for fruit in self.fruits:
                if fruit.stage == 0 and fruit.age > self.parameters["floweringTime"]:
                    fruit.stage = 1 #Post-anthesis to mature green developmental stage
                    #Initiate fruit growth
                    self.reproductiveSwitch = 1
                    self.estimateVcFromPPFD() #Recalculate carbon assimilation constraint as change in specific leaf area over development
                    fruit.mass = self.parameters["fruitMassAtAnthesis"]*self.parameters["fruitNumber"][fruit.position] #Initiate fruit mass at anthesis
                elif fruit.stage == 1 and fruit.age == self.parameters["fruitRipenedAge"]:
                    fruit.stage = 2 #Red ripened - end fruit growth
                    for rxn in self.stoichiometricModel.reactions:
                        if rxn.id.endswith("FR"+str(fruit.position)):
                            rxn.upper_bound = rxn.lower_bound = 0

            #Check if a complete growth cycle has elapsed
            if t >= self.cycleThreshold:
                self.growthCycle = self.growthCycle + 1
                self.cycleThreshold = self.cycleThreshold + self.parameters["PhyllochronLength"]

                #Check if new organs have appeared and add to stoichiometricModel
                self.createNewOrgansAlongAxis()

            #Update fruit biomass composition
            for fruit in self.fruits:
                if fruit.stage == 1:
                    self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).remove_from_model()
                    self.stoichiometricModel.add_reaction(Reaction("Biomass_tx_FR" + str(fruit.position), name = "Biomass_tx", upper_bound = 1000000, lower_bound = 0))
                    for met in self.fruitBiomassEquations:
                        if met != "DPA":
                            self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).add_metabolites({self.stoichiometricModel.metabolites.get_by_id(
                                met + "_FR" + str(fruit.position)): -self.fruitBiomassEquations[met][(fruit.age-self.parameters["floweringTime"])-1]})

            #Update constraints
            self.updateConstraints(t)

            #Clear biomass partitioning constraints from last cycle
            self.stoichiometricModel.reactions.Total_biomass_tx.remove_from_model()
            self.stoichiometricModel.add_reaction(Reaction("Total_biomass_tx", name = "Total_biomass_tx", lower_bound = 0, upper_bound = 1000000))
            self.stoichiometricModel.reactions.Total_biomass_tx.objective_coefficient = 1

            #Constraint fruit growth rate
            for fruit in self.fruits:
                if fruit.stage == 1:
                    fruitGrowthRate = self.sigmoidalDerivative(age=fruit.age, position=fruit.position, organ = "fruit")/self.fruitBiomassWeight[fruit.age-self.parameters["floweringTime"]]
                    self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR"+str(fruit.position)).upper_bound = fruitGrowthRate
                    self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR"+str(fruit.position)).lower_bound = fruitGrowthRate

            self.totSinkStrengthLeaf = sum([self.sigmoidalDerivative(age=leaf.age, position=leaf.position, organ = "leaf") for leaf in self.leaves])
            self.totSinkStrengthStem = sum([self.sigmoidalDerivative(age=stem.age, position=stem.position, organ = "stem") for stem in self.stems])

            #Biomass partitioning to leaves, roots, and stems
            for leaf in self.leaves:
                self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("Biomass_LD" + str(leaf.position)):
                                                                                      -(self.gL*self.leafBiomass*0.75)*(round(self.sigmoidalDerivative(age=leaf.age, position=leaf.position, organ = "leaf"),5)/self.totSinkStrengthLeaf),
                                                                                     self.stoichiometricModel.metabolites.get_by_id("Biomass_LN" + str(leaf.position)):
                                                                                      -(self.gL*self.leafBiomass*0.25)*(round(self.sigmoidalDerivative(age=leaf.age, position=leaf.position, organ = "leaf"),5)/self.totSinkStrengthLeaf)})

            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.Biomass_R1: -(self.gR*self.rootBiomass), 
                                                                                  self.stoichiometricModel.metabolites.Biomass_ST1: -(self.gS*self.stemBiomass)})


            #Optimise model with new constraints
            #self.stoichiometricModel.solver = "cplex"
            self.stoichiometricModel.slim_optimize()
            self.updateBiomassFromCobraModel()

            #Update organ biomass
            leafBiomassOutput.append(self.leafBiomass)
            stemBiomassOutput.append(self.stemBiomass)
            rootBiomassOutput.append(self.rootBiomass)
            fruitBiomassOutput.append(self.fruitBiomass)
            
            print ("DAS=%d   leaf biomass=%.3f   root biomass=%.3f   stem biomass=%.3f   fruit biomass=%.3f" %((t+self.parameters["InitDAS"]), self.leafBiomass, self.rootBiomass, self.stemBiomass, self.fruitBiomass))

        simulationTime = [self.parameters["InitDAS"]]+[i+self.parameters["InitDAS"] for i in simulationTime]

        return simulationTime, leafBiomassOutput, stemBiomassOutput, rootBiomassOutput, fruitBiomassOutput
        
    def updateBiomassFromSolutionObject(self, sol):

        """Update mass of each organ after a period growth and update the total mass of each organ type

        Parameters
        ----------
        - Optimisation solution
        - Previous organ biomass

        Returns
        -------
        - Mass of organ after a period of growth

        """

        biomassCounter = 0
        for leaf in self.leaves:
            leaf.mass = leaf.mass + self.leafBiomassWeight*(abs(sol["Biomass_tx_LD" + str(leaf.position)]) + abs(sol["Biomass_tx_LN" + str(leaf.position)]))
            biomassCounter = biomassCounter + leaf.mass
        self.leafBiomass = biomassCounter

        for stem in self.stems:
            stem.mass = stem.mass + self.stemBiomassWeight*abs(sol["Biomass_tx_ST1"])*(self.sigmoidalDerivative(age=stem.age, position = stem.position, organ = "stem")/self.totSinkStrengthStem)
        self.stemBiomass = self.stemBiomass + self.stemBiomassWeight*abs(sol["Biomass_tx_ST1"])

        biomassCounter = 0
        for fruit in self.fruits:
            if fruit.stage == 1:
                fruit.mass = fruit.mass + self.fruitBiomassWeight[fruit.age-self.parameters["floweringTime"]]*(abs(sol["Biomass_tx_FR" + str(fruit.position)]))
            biomassCounter = biomassCounter + fruit.mass
        self.fruitBiomass = biomassCounter

        self.root.root1.mass = self.root.root1.mass + (self.rootBiomassWeight*abs(sol["Biomass_tx_R1"]))
        self.rootBiomass = self.root.root1.mass

    def updateBiomassFromCobraModel(self):

        """Update mass of each organ after a period growth and update the total mass of each organ type

        Parameters
        ----------
        - Cobra model
        - Previous organ biomass

        Returns
        -------
        - Mass of organ after a period of growth

        """

        biomassCounter = 0
        for leaf in self.leaves:
            leaf.mass = leaf.mass + self.leafBiomassWeight*(abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_LD" + str(leaf.position)).flux) + abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_LN" + str(leaf.position)).flux))
            biomassCounter = biomassCounter + leaf.mass
        self.leafBiomass = biomassCounter

        self.stemBiomass = self.stemBiomass + self.stemBiomassWeight*abs(self.stoichiometricModel.reactions.Biomass_tx_ST1.flux)
        for stem in self.stems:
            stem.mass = stem.mass + (self.stemBiomassWeight*abs(self.stoichiometricModel.reactions.Biomass_tx_ST1.flux)*(self.sigmoidalDerivative(age=stem.age, position = stem.position, organ = "stem")/self.totSinkStrengthStem))

        self.root.root1.mass = self.root.root1.mass + (self.rootBiomassWeight*abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_R1").flux))
        self.rootBiomass = self.root.root1.mass

        biomassCounter = 0
        for fruit in self.fruits:
            if fruit.stage == 1:
                fruit.mass = fruit.mass + (self.fruitBiomassWeight[fruit.age-self.parameters["floweringTime"]]*abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).flux))
            biomassCounter = biomassCounter + fruit.mass
        self.fruitBiomass = biomassCounter

    def generateModelAtEmergence(self):

        """
        Generate whole plant model for tomato at emergence (1 leaf, 1 stem, 1 root)

        Parameters
        ----------
        - Organ-specific stoichiometric models (leaf, root, stem)
        - Initial organ biomasses (leaf, root, stem)

        Returns
        -------
        - Whole-plant model at emergence

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import cobra
        from cobra.core import Model
        from cobra import Reaction, Metabolite

        plantModel = Model()
        t = 1
        plantModel.add_reaction(Reaction("Total_biomass_tx", name = "Total_biomass_tx", upper_bound = 1000000, lower_bound= 0))

        for model in [self.leafModel, self.rootModel, self.stemModel]:
            if model == self.leafModel:
                organTag = "LD"
                model = self.leafModel.copy()
            if model == self.rootModel:
                organTag = "R"
                model = self.rootModel.copy()
            if model == self.stemModel:
                organTag = "ST"
                model = self.stemModel.copy()
            for rxn in model.reactions:
                rxn.id = rxn.id + str(t)
            for met in model.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(t)
                    if met.compartment != "f":
                        met.compartment = met.compartment + str(t)
            BiomassMet = Metabolite("Biomass_" + organTag + str(t), name = "Biomass_" + organTag + str(t), compartment="b_" + organTag + str(t))
            model.reactions.get_by_id("Biomass_tx_" + organTag + str(t)).add_metabolites({BiomassMet: 1})
            if organTag == "LD":
                BiomassMetLN = Metabolite("Biomass_" + "LN" + str(t), name = "Biomass_" + "LN" + str(t), compartment="b_" + "LN" + str(t))
                model.reactions.get_by_id("Biomass_tx_" + "LN" + str(t)).add_metabolites({BiomassMetLN: 1})
            plantModel = plantModel + model
            plantModel.reactions.Total_biomass_tx.add_metabolites({BiomassMet: -1})
            if organTag == "LD":
                plantModel.reactions.Total_biomass_tx.add_metabolites({BiomassMetLN: -1})

        plantModel.reactions.Total_biomass_tx.objective_coefficient = 1

        self.stoichiometricModel = plantModel
        self.addLeaf([Organ("leaf" + str(t), position=t, mass = self.parameters["LbInit"], age = self.parameters["InitDAS"])])
        self.addRoot([Organ("root" + str(t), position=t, mass = self.parameters["RbInit"], age = self.parameters["InitDAS"])])
        self.addStem([Organ("stem" + str(t), position=t, mass = self.parameters["SbInit"], age = self.parameters["InitDAS"])])

        #Initiate total biomass for organ type attribute
        self.leafBiomass = self.leaves.leaf1.mass
        self.rootBiomass = self.root.root1.mass
        self.stemBiomass = self.stems.stem1.mass
        self.reproductiveSwitch = 0
        self.fruitBiomass = 0

    def canopyEffect(self):
        import math

        a = self.parameters["Pmax"]*self.parameters["k"]*(math.e**(-self.parameters["k"]*self.F))
        b = 1-self.parameters["m"]
        return (a/b)/self.parameters["Pmax"]

    def sigmoidalDerivative(self, age, position, organ):

        import math
        import numpy as np
    
        if organ != "fruit":
            a = self.parameters[organ+"VarA"]
            b = self.parameters[organ+"VarB"]
            c = self.parameters[organ+"VarC"]
            x = age/(self.parameters[organ + "GrowthDuration"]*self.parameters[organ+"ExpScale"][position-1])
            top = (a*b)*(np.exp(-b * (x - c)))
            bottom = (np.exp(-b * (x - c))+1)**2
            ans = top/bottom
        else:
            a = self.parameters[organ+"VarA"]
            b = self.parameters[organ+"VarB"]
            c = self.parameters[organ+"VarC"]
            x = age
            top = (a*b)*(np.exp(-b * (x - c)))
            bottom = (np.exp(-b * (x - c))+1)**2
            ans = (top/bottom)*self.parameters["fruitNumber"][position]        
        return ans
        
    def estimateVcFromPPFD(self, Vc_ID = "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD", CO2in_ID = "CO2_tx_LD"):
        """
        Function to estimate rubisco carboxylase flux based on light intensity based
        on tomato leaf PPFD-assimilation rate reported in Nunes-Nesi et al 2005

        Parameters
        ----------
        - Photon uptake rate
        - Vc_ID : String
            reaction ID of day time rubisco carboxylase
        - CO2in_ID : String
            reaction ID of daytime CO2 uptake
        - modelfile : cobra.core.Model
            SBML file of leaf model

        Returns
        -------
        - Rubisco carboxylase flux based on PPFD-A curve from Nunes-Nesi et al 2005

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import cobra
        from cobra import Reaction

        #Copy model and set objective as sucrose export reactions (3:1 day/night ratio)
        model = self.leafModel.copy()
        Suc = Reaction("SUCROSE_EXPORT", name="SUCROSE_EXPORT", upper_bound=1000000, lower_bound=0)
        Suc.add_metabolites({model.metabolites.SUCROSE_c_LD: -3, model.metabolites.SUCROSE_c_LN: -1})
        model.add_reaction(Suc)
        model.reactions.SUCROSE_EXPORT.objective_coefficient = 1

        for rxn in model.reactions:
            if "_ph" in rxn.id:
                rxn.lower_bound = rxn.upper_bound = 0

        #Calculate net CO2 assimilation from PPFD according to empirical data (Nunes-Nesi, 2005)
        netCO2uptake = (40.9962482*self.parameters["Pmax0"])/(self.parameters["Pmax0"] + 644.48886704)
        model.reactions.get_by_id(Vc_ID).lower_bound = netCO2uptake
        model.reactions.get_by_id(Vc_ID).upper_bound = netCO2uptake

        #perform pFBA
        opt = model.slim_optimize()
        model.reactions.SUCROSE_EXPORT.upper_bound = model.reactions.SUCROSE_EXPORT.lower_bound = opt
        model.slim_optimize()

        #set loop counter
        i=0

        #Use a while loop to increase Vc flux until net CO2 rate is similar to given value (or loop counter hits 10)
        while((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x)/netCO2uptake > 0.01 and i<10):
            i=i+1
            prev = model.reactions.get_by_id(Vc_ID).x
            # Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
            now = prev + ((netCO2uptake - model.reactions.get_by_id(CO2in_ID).x))
            model.reactions.get_by_id(Vc_ID).lower_bound = now
            model.reactions.get_by_id(Vc_ID).upper_bound = now
            model.reactions.SUCROSE_EXPORT.upper_bound = 1000000
            model.reactions.SUCROSE_EXPORT.lower_bound = 0
            opt = model.slim_optimize()
            model.reactions.SUCROSE_EXPORT.upper_bound = model.reactions.SUCROSE_EXPORT.lower_bound = opt
            model.slim_optimize()
        if i==10:
            print("Warning : Loop counter hits maximum in EstimateVcFromPPFD")

        #Calculate Vc in mmol/gDW/day
        if self.reproductiveSwitch == 0:
            self.parameters["VcMax"] = (prev*self.parameters["leafAreatoBiomassVeg"]*60*60*self.parameters["photoperiod"])/1000
        else:
            self.parameters["VcMax"] = (prev*self.parameters["leafAreatoBiomassRep"]*60*60*self.parameters["photoperiod"])/1000
            self.parameters["Pmax"] = (self.parameters["Pmax0"]*60*60*self.parameters["photoperiod"]*self.parameters["leafAreatoBiomassRep"])/1000

    def createNewOrgansAlongAxis(self):

        """
        Function to add organs to PlantModel based on modeling time

        Parameters
        ----------
        - A list of time points of emergence for organ of interest
        - The whole plant model
        - The cobra.core.Model object of organ of interest
        - Thermal time point

        Returns
        -------
        - Updated whole plant model

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import cobra
        from cobra import Reaction, Metabolite

        t = self.growthCycle

        if self.parameters["leafPosition"][t-1]:
            modelNew = self.leafModel.copy()
            organTag = "LD"
            for met in modelNew.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(t)
                    if met.compartment != "f":
                        met.compartment = met.compartment + str(t)
            BiomassMetLD = Metabolite("Biomass_" + organTag + str(t), name = "Biomass_" + organTag + str(t), compartment="b_" + organTag + str(t))
            BiomassMetLN = Metabolite("Biomass_" + "LN" + str(t), name = "Biomass_" + "LN" + str(t), compartment="b_" + "LN" + str(t))
            for rxn in modelNew.reactions:
                rxn.id = rxn.id + str(t)
                if rxn.id == "Biomass_tx_" + organTag + str(t):
                    rxn.add_metabolites({BiomassMetLD: 1})
                if rxn.id == "Biomass_tx_" + "LN" + str(t):
                    rxn.add_metabolites({BiomassMetLN: 1})
                self.stoichiometricModel.add_reaction(rxn)
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({BiomassMetLD: -1, BiomassMetLN: -1})
            
            ###Sucrose valve###
            #rxn = Reaction("HYPO_SucroseOnlyImport_tx_LD"+str(t))
            #rxn.add_metabolites({self.stoichiometricModel.metabolites.get_by_id("SUCROSE_c_LD"+str(t)):-3,
             #                    self.stoichiometricModel.metabolites.get_by_id("PROTON_e_LD"+str(t)):-3,
              #                   self.stoichiometricModel.metabolites.get_by_id("SUCROSE_c_LN"+str(t)):-1,
               #                  self.stoichiometricModel.metabolites.get_by_id("PROTON_e_LN"+str(t)):-1,
                #                 self.stoichiometricModel.metabolites.SUCROSE_ph:4,
                 #                self.stoichiometricModel.metabolites.get_by_id("PROTON_c_LD"+str(t)):3,
                  #               self.stoichiometricModel.metabolites.get_by_id("PROTON_c_LN"+str(t)):1})
            #rxn.lower_bound = 0
            #rxn.upper_bound = 1000000
            #self.stoichiometricModel.add_reaction(rxn)
           #################### 

            self.addLeaf([Organ("leaf" + str(t), position=t)])

        if self.parameters["stemPosition"][t-1]:
            self.addStem([Organ("stem" + str(t), position=t)])

        if self.parameters["fruitPosition"][t-1]:

            modelNew = self.fruitModel.copy()
            organTag = "FR"
            for met in modelNew.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + str(t)
                    if met.compartment != "f":
                        met.compartment = met.compartment + str(t)
            BiomassMet = Metabolite("Biomass_" + organTag + str(t), name = "Biomass_" + organTag + str(t), compartment="b_" + organTag + str(t))
            for rxn in modelNew.reactions:
                rxn.id = rxn.id + str(t)
                if rxn.id == "Biomass_tx_" + organTag + str(t):
                    rxn.add_metabolites({BiomassMet: 1})
                self.stoichiometricModel.add_reaction(rxn)
            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({BiomassMet: -1})

            self.addFruit([Organ("fruit" + str(t), position=t, mass = 0)])

    def calculateBiomassWeight(self, rxnID = "Biomass_tx_"):

        """Calculate the molecular weight of biomass to convert flux from mols to g biomass

        Parameters
        ----------
        - Organ specific stoichiometric model
        - Reaction ID of biomass reaction

        Returns
        -------
        - Molar mass of total biomass components for the organ of interest (g/mmol)

        """

        import cobra
        from cobra import Reaction, Metabolite
        import pandas as pd

        self.leafBiomassWeight = sum([abs(met.formula_weight*self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites[met]) for met in self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites])/1000

        self.stemBiomassWeight = sum([abs(met.formula_weight*self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites[met]) for met in self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites])/1000

        self.rootBiomassWeight = sum([abs(met.formula_weight*self.rootModel.reactions.get_by_id(rxnID + "R").metabolites[met]) for met in self.rootModel.reactions.get_by_id(rxnID + "R").metabolites])/1000

        #Fruit biomass molar mass
        model = cobra.io.read_sbml_model("PlantCoreMetabolism.xml")
        df = pd.read_csv("Parameters/fruitBiomassComposition.csv")
        model.reactions.Biomass_tx.remove_from_model()

        for met in ["ASN","FRU", "FRUCTOSE_6P", "GLC", "GLY", "ILE", "LYS", "MALTOSE", "MET", "PHE", "Pi", "PRO", "PYRUVATE", "THR", "TRP", "TYR", "VAL"]:
            met_b = model.metabolites.get_by_id(met + "_c").copy()
            met_b.id = "s" + met + "_b"
            met_b.compartment = "b"
            met_b.name = "s" + met + "_b"
            rxn = Reaction("s" + met + "_biomass", upper_bound = 0, lower_bound = -1000, name = met + "_biomass")
            model.add_reaction(rxn)
            model.reactions.get_by_id("s" + met + "_biomass").add_metabolites({met_b: -1, model.metabolites.get_by_id(met + "_c"): 1})

        fruitBiomassWeight = dict()
        for row in range(len(df)):
            model.add_reaction(Reaction("Biomass_tx", name = "Biomass_tx", upper_bound=1000, lower_bound=-1000))
            for met in df.keys():
                if met != "DPA":
                    model.reactions.Biomass_tx.add_metabolites({model.metabolites.get_by_id(met): -(df[met][row])})
            fruitBiomassWeight[row+1] = sum([met.formula_weight*(-model.reactions.Biomass_tx.metabolites[met]) for met in model.reactions.Biomass_tx.metabolites])/1000
            model.reactions.Biomass_tx.remove_from_model()

        self.fruitBiomassWeight = fruitBiomassWeight
        self.fruitBiomassEquations = df

    def readStrucutrefromMask(self, datafile = "./Parameters/tomMask.csv"):

        """
        Function to import and read mask file for plant structure, specifying organ
        position along the developmental axis.

        Parameters
        ----------
        - datafile : GREENLAB mask file
            A GREENLAB mask file containing information on organ emergence in thermal time

        Returns
        -------
        - A list of time points of leaf emergence
        - A list of time points of stem emergence
        - A list of time points of fruit emergence
        - A list of expansion time scaling of leaves at each phytomer unit
        - A list of expansion time scaling of stems at each phytomer unit
        - A list of expansion time scaling of fruits at each phytomer unit

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import pandas as pd
        import numpy as np

        #Import positions
        df = pd.read_csv(datafile, delimiter=",")

        self.parameters["leafPosition"] = [int(df["leaf"][n]) for n in range(0, len(df))]
        self.parameters["stemPosition"] = [int(df["stem"][n]) for n in range(0, len(df))]
        self.parameters["fruitPosition"] = [int(df["fruit"][n]) for n in range(0, len(df))]

        self.parameters["leafExpScale"] = [float(df["leafExpScale"][n]) for n in range(0, len(df))]
        self.parameters["stemExpScale"] = [float(df["stemExpScale"][n]) for n in range(0, len(df))]
        self.parameters["fruitExpScale"] = [float(df["fruitExpScale"][n]) for n in range(0, len(df))]

    def updateConstraints(self, t):

        """
        Function to calculate maximum nitrate uptake rate by root, update maintenance respiration constraints,
        and maximum light uptake and carbon assimilation rate by leaf
        """

        #Root nitrate uptake
        if self.parameters["nitrateConc"] >= 1:
            num = float(self.parameters["NO3UptakeVm"]*self.parameters["nitrateConc"])
            denom = float(self.parameters["NO3UptakeKm"]+self.parameters["nitrateConc"])
        else:
            num = float(self.parameters["NO3UptakeVmLow"]*self.parameters["nitrateConc"])
            denom = float(self.parameters["NO3UptakeKmLow"]+self.parameters["nitrateConc"])

        MichMen = float(num/denom)
        MichMen = float(MichMen*self.root.root1.mass)
        self.stoichiometricModel.reactions.Nitrate_ec_R1.upper_bound = MichMen

        #Update maintenace respiration requirements
        #Leaves
        for leaf in self.leaves:
            if leaf.stage == 0:
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*self.parameters["photoperiod"]
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*(24-self.parameters["photoperiod"])
            else: #Leaf has senseced
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).lower_bound = 0
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).lower_bound = 0

        #Stem
        self.stoichiometricModel.reactions.ATPase_tx_ST1.upper_bound = self.stoichiometricModel.reactions.ATPase_tx_ST1.lower_bound = self.parameters["stemMaint"]*self.stemBiomass

        #Fruits
        for fruit in self.fruits:
            if fruit.stage == 1:
                self.stoichiometricModel.reactions.get_by_id("ATPase_tx_FR" + str(fruit.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_FR" + str(fruit.position)).lower_bound = self.parameters["fruitMaint"]*fruit.mass

        #Root
        self.stoichiometricModel.reactions.ATPase_tx_R1.upper_bound = self.stoichiometricModel.reactions.ATPase_tx_R1.lower_bound = self.parameters["rootMaint"]*self.root.root1.mass
        
        #Leaf photon uptake rate and RuBisCO carboxylation
        if t<8:
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomassVeg"]*15
        elif t>=8 and t<45:
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomassVeg"]*3.75
        elif t>=45 and t<70:
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomassRep"]*2.1
        else:
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomassRep"]*2.1
        canopyReductionFactor = self.canopyEffect()
        cotLeaf = min([l.position for l in self.leaves if l.stage == 0])
        for leaf in self.leaves:
            if leaf.stage == 0: #Leaf is alive
                if leaf.position == cotLeaf and t<(self.parameters["CotAge"]-self.parameters["InitDAS"]):
                    self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD" + str(leaf.position)).upper_bound = self.parameters["VcMax"]*(leaf.mass+self.parameters["CotWeightAdj"])*canopyReductionFactor
                    self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD" + str(leaf.position)).upper_bound = self.parameters["Pmax"]*(leaf.mass+self.parameters["CotWeightAdj"])*canopyReductionFactor
                else:
                    self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD" + str(leaf.position)).upper_bound = self.parameters["VcMax"]*leaf.mass*canopyReductionFactor
                    self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD" + str(leaf.position)).upper_bound = self.parameters["Pmax"]*leaf.mass*canopyReductionFactor
            else: #Leaf has senesced
                self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD" + str(leaf.position)).upper_bound = 0
                self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD" + str(leaf.position)).upper_bound = 0

    def generateStoichiometricModels(self, phloemComposition = "Shameer2018", sbmlFile = "PlantCoreMetabolism.xml", biomassEquations = "./Parameters/tomatoBiomassEquations.csv", biomassEquationsFruit = "Parameters/fruitBiomassComposition.csv"):

        """
        Generate organ-specific stoichiometric models (diel leaf, root, stem, and fruit) from PlantCoreMetabolism_v1_2.xml

        Parameters
        ----------
        - Core stoichiometric model (PlantCoreMetabolism_v1_2_3.xml)
        - Organ-specific biomass compositions
        - Phloem composition option
            - Shameer2018: tomato composition reported in Shameer et al., 2018
            - free: unfixed composition

        Returns
        -------
        - Organ-specific models for fruit, diel leaf, stem, and root metabolism

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import cobra
        from cobra.core import Model
        from cobra import Reaction, Metabolite
        from cobra.flux_analysis import pfba
        import csv
        import pandas as pd

        model = cobra.io.read_sbml_model(sbmlFile)

        #Turn off exchange reactions
        reactionsOff = ["GLC_tx", "Sucrose_tx", "NH4_tx", "NADPH_Dehydrogenase_p", "Plastoquinol_Oxidase_p", "ATP_pc",
                        "Photon_tx", "Pi_tx", "SO4_tx", "Nitrate_tx", "Ca_tx", "Mg_tx", "K_tx"]
        for rxn in reactionsOff:
                model.reactions.get_by_id(rxn).upper_bound = model.reactions.get_by_id(rxn).lower_bound = 0
        model.reactions.Biomass_tx.remove_from_model()
        model.reactions.AraCore_Biomass_tx.remove_from_model()

        #Increase flux bounds for all models
        for rxn in model.reactions:
            if rxn.upper_bound == 1000:
                rxn.upper_bound = 1000000
            if rxn.lower_bound == -1000:
                rxn.lower_bound = -1000000

        #Add reactions to complete branch-chain amino acid catabolism pathways
        model.add_reaction(Reaction("ACETOACETATE_COA_LIGASE_RXN_m", name = "ACETOACETATE_COA_LIGASE_RXN", upper_bound = 1000000, lower_bound = 0))
        met = model.metabolites.ACETOACETYL_COA_x.copy()
        met.compartment = "m"
        met.id = "ACETOACETYL_COA_m"
        model.reactions.ACETOACETATE_COA_LIGASE_RXN_m.add_metabolites({model.metabolites.get_by_id("3_KETOBUTYRATE_m"): -1,
                                                                            model.metabolites.CO_A_m: -1,
                                                                            model.metabolites.ATP_m: -0.9,
                                                                            model.metabolites.aATP_m: -0.1,
                                                                            model.metabolites.AMP_m: 1,
                                                                            model.metabolites.PPI_m: 0.55,
                                                                            model.metabolites.bPPI_m: 0.45,
                                                                            met: 1})
        met = model.metabolites.ACETOACETYL_COA_x.copy()
        met.compartment = "c"
        met.id = "ACETOACETYL_COA_c"
        model.add_reaction(Reaction("ACETOACETATE_COA_COA_mc", name = "ACETOACETATE_COA_COA_mc", upper_bound = 1000000, lower_bound = -1000000))
        model.reactions.ACETOACETATE_COA_COA_mc.add_metabolites({model.metabolites.ACETOACETYL_COA_m: -1,
                                                                      model.metabolites.CO_A_c: -1,
                                                                      model.metabolites.CO_A_m: 1,
                                                                      met: 1})
        model.add_reaction(Reaction("ACETYL_COA_ACETYLTRANSFER_RXN_c", name = "ACETYL_COA_ACETYLTRANSFER_RXN", upper_bound = 1000000, lower_bound = -1000000))
        model.reactions.ACETYL_COA_ACETYLTRANSFER_RXN_c.add_metabolites({model.metabolites.ACETYL_COA_c: -2,
                                                                            model.metabolites.CO_A_c: 1,
                                                                            met: 1})

        model.reactions.PROPIONYL_COA_xc.add_metabolites({model.metabolites.CO_A_c: -1,
                                                               model.metabolites.CO_A_x: 1})
        model.reactions.PROPIONYL_COA_mc.add_metabolites({model.metabolites.CO_A_c: -1,
                                                               model.metabolites.CO_A_m: 1})

        #Fix ratio of ATPase and NAD(P)Hox energy consumption reactions to a fixed ratio of 3:1.
        ATP_NADPH_f = Metabolite('ATP_NADPH_f', name='ATP NADPH pseudo', compartment='f')
        model.reactions.ATPase_tx.add_metabolites({ATP_NADPH_f: 1.0})
        model.reactions.NADPHoxc_tx.add_metabolites({ATP_NADPH_f: -3.0})
        model.reactions.NADPHoxm_tx.add_metabolites({ATP_NADPH_f: -3.0})
        model.reactions.NADPHoxp_tx.add_metabolites({ATP_NADPH_f: -3.0})

        #Adding biomass reactions for new biomass components
        for met in ["ASN","FRU", "FRUCTOSE_6P", "GLC", "GLY", "ILE", "LYS", "MALTOSE", "MET", "PHE", "Pi", "PRO", "PYRUVATE", "THR", "TRP", "TYR", "VAL"]:
            met_b = model.metabolites.get_by_id(met + "_c").copy()
            met_b.id = "s" + met + "_b"
            met_b.compartment = "b"
            met_b.name = "s" + met + "_b"
            rxn = Reaction("s" + met + "_biomass", upper_bound = 0, lower_bound = -1000000, name = met + "_biomass")
            model.add_reaction(rxn)
            model.reactions.get_by_id("s" + met + "_biomass").add_metabolites({met_b: -1, model.metabolites.get_by_id(met + "_c"): 1})

        #Establishing the xylem (reactions defined as uptake)
        xylem_metabolites = ["KI", "CAII", "MGII", "NITRATE", "Pi", "SULFATE"]
        for met in xylem_metabolites:
            rxn = Reaction(met + "_xy")
            rxn.name = met + "_xy"
            rxn.subsystem = "Xylem"
            rxn.lower_bound = 0
            rxn.upper_bound = 1000000
            new_met = model.metabolites.get_by_id(met + "_c").copy()
            new_met.name = met + "_xy"
            new_met.id = met + "_xy"
            new_met.compartment = "xy"
            if met == "Pi":
                rxn.add_metabolites({model.metabolites.get_by_id(met + "_c"): 0.7, model.metabolites.get_by_id("a" + met + "_c"): 0.3,
                                     new_met: -1, model.metabolites.PROTON_e: -3, model.metabolites.PROTON_c: 4.7})
            else:
                rxn.add_metabolites({model.metabolites.get_by_id(met + "_c"): 1, new_met: -1})
            if met == "NITRATE":
                rxn.add_metabolites({model.metabolites.PROTON_e: -2, model.metabolites.PROTON_c: 2})
            elif met == "SULFATE":
                rxn.add_metabolites({model.metabolites.PROTON_e: -3, model.metabolites.PROTON_c: 3})
            elif met == "KI":
                rxn.add_metabolites({model.metabolites.PROTON_e: -1, model.metabolites.PROTON_c: 1})

            model.add_reaction(rxn)

        #Adding molecular weights to undefined metabolites
        model.metabolites.FRUCTAN_v.formula = model.metabolites.FRU_v.formula
        model.metabolites.Heteroglycans_c.formula = model.metabolites.GLC_c.formula

        #Adjust phloem to include cytosolic sucrose
        sucrose_stoic = abs(model.reactions.Phloem_output_tx.metabolites[model.metabolites.sSUCROSE_b])
        model.reactions.Phloem_output_tx.add_metabolites({model.metabolites.sSUCROSE_b: sucrose_stoic,
                                                          model.metabolites.SUCROSE_c: -sucrose_stoic})

        #Make separate reaction for each phloem component and disable defined phloem reaction
        if phloemComposition == "free":
            for met in model.reactions.Phloem_output_tx.metabolites:
                if not((met.id.startswith("PROTON"))):
                    met_name = met.id[:-2]
                    rxn = Reaction(met_name + "_ph", name = met_name + "_ph", subsystem = "Phloem", lower_bound = 0, upper_bound = 1000000)
                    new_met = model.metabolites.get_by_id(met.id).copy()
                    new_met.name = met_name + "_ph"
                    new_met.id = met_name + "_ph"
                    new_met.compartment = "ph"
                    model.add_reaction(rxn)
                    model.reactions.get_by_id(met_name + "_ph").add_metabolites({model.metabolites.PROTON_c: 1,
                                                                                 model.metabolites.PROTON_e: -1,
                                                                                 met: -1,
                                                                                 new_met: 1})
            #Remove phloem reaction
            model.reactions.Phloem_output_tx.remove_from_model()

        #Establishing fixed phloem composition
        if phloemComposition == "Shameer2018":

            #Removing amino acids without degradation pathways from phloem composition
            for met in ["MET", "HIS", "CYS", "PHE", "TRP", "TYR"]: #["HIS", "CYS", "PHE", "TRP", "TYR"]
                metStoic = abs(model.reactions.Phloem_output_tx.metabolites[model.metabolites.get_by_id(met+"_c")])
                model.reactions.Phloem_output_tx.add_metabolites({model.metabolites.get_by_id(met+"_c"): metStoic, model.metabolites.PROTON_c: -metStoic, model.metabolites.PROTON_e: metStoic})

            for met in model.reactions.Phloem_output_tx.metabolites:
                if "PROTON" not in met.id:
                    newMet = met.copy()
                    newMet.compartment = "ph"
                    newMet.name = met.id[:-2] + "_ph"
                    newMet.id = met.id[:-2] + "_ph"
                    model.reactions.Phloem_output_tx.add_metabolites({newMet: -model.reactions.Phloem_output_tx.metabolites[met]})

        self.PlantCoreMetabolism_v1_2 = model

        #Create diel source leaf model
        diel_model = Model('Diel_model')
        model = self.PlantCoreMetabolism_v1_2.copy()
        for met in model.metabolites:
            if not(met.id.endswith("_xy") or (met.id.endswith("_ph"))):
                met.id = met.id + "_LD"
                if met.compartment != "f":
                    met.compartment = met.compartment + "_LD"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_LD"
            diel_model.add_reaction(rxn_copy)
        model = self.PlantCoreMetabolism_v1_2.copy()
        for met in model.metabolites:
            if not(met.id.endswith("_xy") or (met.id.endswith("_ph"))):
                met.id = met.id + "_LN"
                if met.compartment != "f":
                    met.compartment = met.compartment + "_LN"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_LN"
            diel_model.add_reaction(rxn_copy)

        #Adding linker reactions between day and night metabolites that are allowed to accumulate.
        linker_reactions = ["STARCH", "NITRATE", "GLC", "FRU"]
        for linker in linker_reactions:
            if linker == "STARCH":
                met_day = diel_model.metabolites.get_by_id(linker + "_p_LD")
                met_night = diel_model.metabolites.get_by_id(linker + "_p_LN")
                rxn = Reaction(linker + "_p_linker_LD")
            else:
                met_day = diel_model.metabolites.get_by_id(linker + "_v_LD")
                met_night = diel_model.metabolites.get_by_id(linker + "_v_LN")
                rxn = Reaction(linker + "_v_linker_LD")
            rxn.name = linker + " day to night linker"
            rxn.subsystem = "Day-Night linker reaction"
            rxn.lower_bound = -1000000
            rxn.upper_bound = 1000000
            rxn.add_metabolites({met_day: -1.0, met_night: 1.0})
            diel_model.add_reaction(rxn)
        rxn = Reaction("MAL_v_linker_LD", name = "MAL day to night linker", lower_bound = -1000000, upper_bound = 1000000)
        rxn.add_metabolites({diel_model.metabolites.MAL_v_LD: -0.7, diel_model.metabolites.aMAL_v_LD: -0.3, diel_model.metabolites.MAL_v_LN: 0.7, diel_model.metabolites.aMAL_v_LN: 0.3})
        diel_model.add_reaction(rxn)
        rxn = Reaction("CIT_v_linker_LD", name = "CIT day to night linker", lower_bound = -1000000, upper_bound = 1000000)
        rxn.add_metabolites({diel_model.metabolites.CIT_v_LD: -0.5, diel_model.metabolites.aCIT_v_LD: -0.5, diel_model.metabolites.CIT_v_LN: 0.5, diel_model.metabolites.aCIT_v_LN: 0.5})
        diel_model.add_reaction(rxn)
        Amino_acids = ["L_ALPHA_ALANINE", "L_ASPARTATE", "ARG", "ASN", "CYS", "GLN", "GLT", "GLY", "ILE",
                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "bHIS"]
        for aa in Amino_acids:
            rxn = Reaction(aa + "_v_linker_LD", name = aa + " day to night linker", subsystem = "Day-Night linker reaction", lower_bound = 0, upper_bound = 1000000)
            rxn.add_metabolites({diel_model.metabolites.get_by_id(aa + "_v_LD"): -1.0, diel_model.metabolites.get_by_id(aa + "_v_LN"): 1.0})
            diel_model.add_reaction(rxn)

        #Fix the nitrate uptake ratio to 3:2 (day:night).
        Nitrate_f = Metabolite('Nitrate_f', name='Nitrate pseudo-metabolite', compartment='f')
        diel_model.reactions.NITRATE_xy_LD.add_metabolites({Nitrate_f: 2.0})
        diel_model.reactions.NITRATE_xy_LN.add_metabolites({Nitrate_f: -3.0})

        #Fix the Vc/Vo ratio to 3:1 using pseudo-metabolites for both day and night.
        Rubisco_f = Metabolite('Rubisco_f', name='Rubisco Vc/Vo day pseudo-metabolite', compartment='f')
        diel_model.reactions.RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD.add_metabolites({Rubisco_f: 1.0})
        diel_model.reactions.RXN_961_p_LD.add_metabolites({Rubisco_f: -3.0})

        #Add biomass reaction
        biomass_file = csv.reader(open(biomassEquations,"r"))
        diel_model.add_reaction(Reaction("Biomass_tx_LD", name = "Biomass_tx_LD", lower_bound = 0, upper_bound = 1000000))
        for row in biomass_file:
            met_id = row[0]+"_LD"
            stoich = float(row[1])
            diel_model.reactions.get_by_id("Biomass_tx_LD").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})
        biomass_file = csv.reader(open(biomassEquations,"r"))
        diel_model.add_reaction(Reaction("Biomass_tx_LN", name = "Biomass_tx_LN", lower_bound = 0, upper_bound = 1000000))
        for row in biomass_file:
            met_id = row[0]+"_LN"
            stoich = float(row[1])
            diel_model.reactions.get_by_id("Biomass_tx_LN").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})

        #Constrain starch phosphorylase reaction in leaf
        diel_model.reactions.G6P_Pi_pc_LD.upper_bound = diel_model.reactions.G6P_Pi_pc_LD.lower_bound = 0
        diel_model.reactions.G6P_Pi_pc_LN.upper_bound = diel_model.reactions.G6P_Pi_pc_LN.lower_bound = 0

        #Allow photon uptake during the day
        diel_model.reactions.Photon_tx_LD.upper_bound = 1000000

        #Create reverse phloem reactions for uptake from the phloem
        if phloemComposition == "free":
            revRxnDict = dict()
            for rxn in diel_model.reactions:
                if "_ph" in rxn.id:
                    revRxn = Reaction(rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], name = rxn.id.split("_ph")[0] + "_rev" + rxn.id.split(rxn.id.split("_ph")[0])[1], upper_bound = 1000000, lower_bound = 0, subsystem = "Phloem")
                    for met in rxn.metabolites:
                        if "PROTON" not in met.id:
                            revRxn.add_metabolites({met: -rxn.metabolites[met]})
                        else:
                            revRxn.add_metabolites({met: rxn.metabolites[met]})
                    revRxnDict[revRxn] = revRxn
            for rxn in revRxnDict:
                diel_model.add_reaction(revRxnDict[rxn])

            #Phloem day-night 3:1 constraint
            #Phloem_output_day = Metabolite("Phloem_LD_f", name = "Phloem output LD", compartment = "f")
            #Phloem_output_night = Metabolite("Phloem_LN_f", name = "Phloem output LN", compartment = "f")
            #Phloem_output_rxn = Reaction("Phloem_LD_LN_constraint_LD", name = "Phloem_LD_LN_constraint", upper_bound = 1000000, lower_bound = -1000000)
            #diel_model.add_reaction(Phloem_output_rxn)
            #diel_model.reactions.get_by_id("Phloem_LD_LN_constraint_LD").add_metabolites({Phloem_output_day: -3, Phloem_output_night: -1})
            
            #for rxn in diel_model.reactions.query("ph_LD"):
             #   if "rev" in rxn.id:
             #       rxn.add_metabolites({diel_model.metabolites.Phloem_LD_f:-1})
             #   if "rev" not in rxn.id:
             #       rxn.add_metabolites({diel_model.metabolites.Phloem_LD_f:1})
            #for rxn in diel_model.reactions.query("ph_LN"):
             #   if "rev" in rxn.id:
             #       rxn.add_metabolites({diel_model.metabolites.Phloem_LN_f:-1})
             #   if "rev" not in rxn.id:
             #       rxn.add_metabolites({diel_model.metabolites.Phloem_LN_f:1})
            
            phloemRxnList = list()
            for rxn in diel_model.reactions:
                if "_rev_ph_LD" in rxn.id:
                    phloemRxnList.append(rxn.id.split("_rev")[0])    
            for metName in phloemRxnList:
                pseudoMet = Metabolite(metName + "_PseudoPhloemLeaf_f", name = metName + "_PseudoPhloemLeaf_f", compartment="f")
                diel_model.reactions.get_by_id(metName + "_ph_LD").add_metabolites({pseudoMet: 1})
                diel_model.reactions.get_by_id(metName + "_ph_LN").add_metabolites({pseudoMet: -3})
                pseudoMet = Metabolite(metName + "_rev_PseudoPhloemLeaf_f", name = metName + "_rev_PseudoPhloemLeaf_f", compartment="f")
                diel_model.reactions.get_by_id(metName + "_rev_ph_LD").add_metabolites({pseudoMet: 1})
                diel_model.reactions.get_by_id(metName + "_rev_ph_LN").add_metabolites({pseudoMet: -3})

        elif phloemComposition == "Shameer2018":
            for rxn in [diel_model.reactions.Phloem_output_tx_LD, diel_model.reactions.Phloem_output_tx_LN]:
                revRxn = Reaction(rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], name = rxn.id.split("_tx")[0] + "_rev" + rxn.id.split(rxn.id.split("_tx")[0])[1], upper_bound = 1000000, lower_bound = 0, subsystem = "Phloem")
                for met in rxn.metabolites:
                    if "PROTON" not in met.id:
                        revRxn.add_metabolites({met: -rxn.metabolites[met]})
                    else:
                        revRxn.add_metabolites({met: rxn.metabolites[met]})
                diel_model.add_reaction(revRxn)

            #Phloem day-night 3:1 constraint
            Phloem_output_day = Metabolite("Phloem_LD_f", name = "Phloem output LD", compartment = "f")
            Phloem_output_night = Metabolite("Phloem_LN_f", name = "Phloem output LN", compartment = "f")
            diel_model.add_reaction(Reaction("Phloem_LD_LN_constraint_LD", name = "Phloem_LD_LN_constraint", upper_bound = 1000000, lower_bound = -1000000))
            diel_model.reactions.Phloem_LD_LN_constraint_LD.add_metabolites({Phloem_output_day: -3, Phloem_output_night: -1})
            diel_model.reactions.Phloem_output_tx_LD.add_metabolites({Phloem_output_day: 1})
            diel_model.reactions.Phloem_output_rev_tx_LD.add_metabolites({Phloem_output_day: -1})
            diel_model.reactions.Phloem_output_tx_LN.add_metabolites({Phloem_output_night: 1})
            diel_model.reactions.Phloem_output_rev_tx_LN.add_metabolites({Phloem_output_night: -1})

        self.leafModel = diel_model

        #Create sink models for root (R), fruit (F), and stem (S)
        sink_tissues = ["R", "FR", "ST"]
        for tissue in sink_tissues:
            model = self.PlantCoreMetabolism_v1_2.copy()
            for rxn in model.reactions:
                if rxn.id.endswith("_ph"): #Reverse phloem reaction direction for uptake
                    rxn.upper_bound = 0
                    rxn.lower_bound = -1000000
                    rxn.add_metabolites({model.metabolites.PROTON_e: 2, model.metabolites.PROTON_c:-2})
                if "Phloem_output_tx" in rxn.id:
                    rxn.upper_bound = 0
                    rxn.lower_bound = -1000000
                    rxn.add_metabolites({model.metabolites.PROTON_e: 1.9207920792, model.metabolites.PROTON_c:-1.9207920792})
                rxn.id = rxn.id + "_" + tissue
            for met in model.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + "_" + tissue
                    if met.compartment != "f":
                        met.compartment = met.compartment + "_" + tissue
            model.add_reaction(Reaction("Biomass_tx_" + tissue, name = "Biomass_tx_" + tissue, lower_bound = 0, upper_bound = 1000000))
            if tissue == "FR":
                biomass_file = pd.read_csv(biomassEquationsFruit)
                for met in biomass_file.keys():
                    if met != "DPA":
                        model.reactions.get_by_id("Biomass_tx_" + tissue).add_metabolites({model.metabolites.get_by_id(met+"_"+tissue): -biomass_file[met][0]})
                #Make biomass metabolite reversible for metabolites that can be degraded by the fruit during ripening
                for met in ['Starch_b', 'sASN_b', 'sFRU_b', 'sFUM_b', 'sGABA_b', 'sGLC_b', 'sGLY_b', 'sMALTOSE_b', 'sMAL_b', 'sPRO_b', 'sPYRUVATE_b', 'sSER_b', 'sSUCROSE_b', 'sTHR_b', 'sVAL_b']:
                    rxn = model.reactions.get_by_id(met + "iomass_" + tissue)
                    rxn.upper_bound = 1000000
                self.fruitModel = model
            elif tissue == "ST":
                biomass_file = csv.reader(open(biomassEquations,"r"))
                for row in biomass_file:
                    model.reactions.get_by_id("Biomass_tx_" + tissue).add_metabolites({model.metabolites.get_by_id(row[0]+"_"+tissue): -float(row[2])})
                self.stemModel = model
            elif tissue == "R":
                biomass_file = csv.reader(open(biomassEquations,"r"))
                for row in biomass_file:
                    model.reactions.get_by_id("Biomass_tx_" + tissue).add_metabolites({model.metabolites.get_by_id(row[0]+"_"+tissue): -float(row[3])})
                #Allow for nutrient uptake
                for rxn in ["Pi", "Nitrate", "SO4", "Ca", "Mg", "K"]:
                    model.reactions.get_by_id(rxn + "_tx_R").upper_bound = 1000000
                    model.reactions.get_by_id(rxn + "_tx_R").lower_bound = -1000000
                #Create new compartment for soil uptake (protons)
                proton_met = model.metabolites.PROTON_e_R.copy()
                proton_met.id = "PROTON_s_R"
                proton_met.compartment = "s_R"
                for rxn in model.reactions:
                    if "_ec_R" in rxn.id and model.metabolites.PROTON_e_R in rxn.metabolites:
                        stoich = rxn.metabolites[model.metabolites.PROTON_e_R]
                        rxn.add_metabolites({model.metabolites.PROTON_e_R: abs(stoich), proton_met: stoich})
                #Create new proton pump for soil
                PROTON_pump = model.reactions.PROTON_ATPase_c_R.copy()
                PROTON_pump.id = "PROTON_ATPase_s_R"
                PROTON_pump.name = "PROTON_ATPase_s"
                stoich = model.reactions.PROTON_ATPase_c_R.metabolites[model.metabolites.PROTON_e_R]
                PROTON_pump.add_metabolites({model.metabolites.PROTON_e_R: -abs(stoich), model.metabolites.PROTON_s_R: stoich})
                model.add_reaction(PROTON_pump)
                #Reverse xylem reaction direction for output into the xylem
                for rxn in model.reactions:
                    if rxn.id.endswith("_xy_R"):
                        rxn.upper_bound = 0
                        rxn.lower_bound = -1000000
                        #Adjust cost of xylem loading
                        for met in rxn.reactants:
                            if "PROTON" in met.id:
                                rxn.add_metabolites({met: abs(rxn.metabolites[met])})
                        for met in rxn.products:
                            if "PROTON" in met.id:
                                rxn.add_metabolites({met: -rxn.metabolites[met]})
                self.rootModel = model

class Organ:

    #"""Create new organ objects with age, mass, and position attributes"""

    def __init__(self,id, position, mass = 0, age = 1, stage = 0):
        self.id = id
        self.age = age
        self.mass = mass
        self.position = position
        self.stage = stage

    def __repr__(self):
        return '%s \nPosition: %d \nAge: %d \nMass: %.4f \nStage: %d' %(self.id, self.position, self.age, self.mass, self.stage)
