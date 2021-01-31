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

    def addStem(self, newOrgan):
        self.stems += newOrgan

    def initialisation(self):

        self.loadParams()
        self.generateStoichiometricModels()
        self.calculateBiomassWeight()
        self.generateModelAtEmergence()
        self.estimateVcFromPPFD()

    def loadParams(self, datafile = "parameters.csv"):
        """
        Function to load parameters from CSV file

        Parameters
        ----------
        - datafile : comma-seperated CSV file
            A CSV file with parameters

        Returns
        -------
        - ATP hydrolysis flux for NGAM
        - NADPH oxidase flux for NGAM
        - Vmax for Nitrate uptake
        - Km for Nitrate uptake
        - [NO3-] for Nitrate uptake
        - Maximum photon uptake rate
        - Initial root biomass
        - Initial leaf biomass
        - Photoperiod length
        - Conversion factor between leaf area and biomass (m2 to g)

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import pandas as pd
        import numpy as np
        import math

        df = pd.read_csv(datafile, delimiter=",")

        self.parameters["leafMaint"] = df["leafMaint"][0]
        self.parameters["rootMaint"] = df["rootMaint"][0]
        self.parameters["stemMaint"] = df["stemMaint"][0]
        self.parameters["NO3UptakeVm"] = df["NO3UptakeVm"][0]
        self.parameters["NO3UptakeKm"] = df["NO3UptakeKm"][0]
        self.parameters["NO3UptakeVmLow"] = df["NO3UptakeVmLow"][0]
        self.parameters["NO3UptakeKmLow"] = df["NO3UptakeKmLow"][0]
        self.parameters["soilNitrate"] = df["InitNO3Con"][0]
        self.parameters["RbInit"] = df["Rb"][0]
        self.parameters["LbInit"] = df["Lb"][0]
        self.parameters["SbInit"] = df["Sb"][0]
        self.parameters["photoperiod"] = df["photoperiod"][0]
        self.parameters["leafAreatoBiomass"] = df["LeafAreatoBiomass"][0]
        self.parameters["m"] = df["m"][0]
        self.parameters["k"] = df["k"][0]
        self.gL = df["gL"][0]
        self.gR = df["gR"][0]
        self.gS = df["gS"][0]
        self.parameters["CotArea"] = df["CotArea"][0]
        self.parameters["CotAge"] = df["CotAge"][0]
        self.parameters["InitDAS"] = df["InitDAS"][0]

        self.parameters["Pmax"] = df["Pmax"][0]*self.parameters["photoperiod"]*self.parameters["leafAreatoBiomass"] #mmol/gDW/day
        self.soilNitrate = df["InitNO3Con"][0] #mM

        self.parameters["CotWeightAdj"] = self.parameters["CotArea"]/self.parameters["leafAreatoBiomass"]
        self.outputCycle = []
        
    def simulateGrowth(self, simulationTime = list(range(1,40))):

        #Import packages
        import cobra
        from cobra.flux_analysis import pfba
        from cobra import Reaction, Metabolite
        import numpy as np
        import pandas as pd
        #from IPython.display import display, clear_output
        import copy

        leafBiomassOutput = [self.parameters["LbInit"]] #list to store leaf growth predictions
        stemBiomassOutput = [self.parameters["SbInit"]] #list to store stem growht predictions
        rootBiomassOutput = [self.parameters["RbInit"]] #list to store root growth predictions
        fluxOutputs = dict() #dictionary to store flux predictions (specify days of interest in outputcycle)
        self.outputCycle = [i - self.parameters["InitDAS"] for i in self.outputCycle]

        for t in simulationTime:

            #Constrain maximum nitrate uptake rate by root
            self.maxNitrateUptake() #mmol/root/day

            #Constraint maximum light uptake and carbon assimilation by leaves
            self.maxLightUptakeCarbonAssimilation(t) #mmol/leaves/day

            #Update maintenace respiration requirements
            self.maintRespirationConstraint() #mmol/organ/day

            self.stoichiometricModel.reactions.Total_biomass_tx.remove_from_model()
            self.stoichiometricModel.add_reaction(Reaction("Total_biomass_tx", name = "Total_biomass_tx", lower_bound = 0, upper_bound = 1000000))
            self.stoichiometricModel.reactions.Total_biomass_tx.objective_coefficient = 1

            self.stoichiometricModel.reactions.Total_biomass_tx.add_metabolites({self.stoichiometricModel.metabolites.Biomass_LD1: -(self.gL)*0.75*self.leafBiomass,
                                                                                  self.stoichiometricModel.metabolites.Biomass_LN1: - (self.gL)*0.25*self.leafBiomass,
                                                                                  self.stoichiometricModel.metabolites.Biomass_ST1: -self.gS*self.stemBiomass,
                                                                                  self.stoichiometricModel.metabolites.Biomass_R1: -self.gR*self.rootBiomass})

            #Optimise model with new constraints
            if t in self.outputCycle:
                sol = cobra.flux_analysis.parsimonious.pfba(self.stoichiometricModel)
                fluxOutputs["solution"+str(t+self.parameters["InitDAS"])] = sol
                self.updateBiomassFromSolutionObject(sol)
            if t not in self.outputCycle:
                self.stoichiometricModel.slim_optimize()
                self.updateBiomassFromCobraModel()

            #Update organ biomass for output lists
            leafBiomassOutput.append(self.leafBiomass)
            stemBiomassOutput.append(self.stemBiomass)
            rootBiomassOutput.append(self.rootBiomass)

            #Display plant growth
            print ("DAS=%d   leaf biomass=%.3f   root biomass=%.3f   stem biomass=%.3f" %((t+self.parameters["InitDAS"]), self.leafBiomass, self.rootBiomass, self.stemBiomass))
            #if t != simulationTime[-1]:
             #   clear_output(wait=True)

        simulationTime = [0] + simulationTime
        for i in range(len(simulationTime)):
            simulationTime[i] = simulationTime[i] + self.parameters["InitDAS"]

        return simulationTime, leafBiomassOutput, rootBiomassOutput, stemBiomassOutput, fluxOutputs

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

        biomassCounter = 0
        for stem in self.stems:
            stem.mass = stem.mass + self.stemBiomassWeight*abs(sol["Biomass_tx_ST" + str(stem.position)])
            biomassCounter = biomassCounter + stem.mass
        self.stemBiomass = biomassCounter

        biomassCounter = 0
        for fruit in self.fruits:
            fruit.mass = fruit.mass + self.fruitBiomassWeight["stage" + str(fruit.stage)]*(abs(sol["Biomass_tx_FR" + str(fruit.position)]))
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

        biomassCounter = 0
        for stem in self.stems:
            stem.mass = stem.mass + self.stemBiomassWeight*abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_ST" + str(stem.position)).flux)
            biomassCounter = biomassCounter + stem.mass
        self.stemBiomass = biomassCounter

        biomassCounter = 0
        for fruit in self.fruits:
            fruit.mass = fruit.mass + (self.fruitBiomassWeight["stage" + str(fruit.stage)]*abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_FR" + str(fruit.position)).flux))
            biomassCounter = biomassCounter + fruit.mass
        self.fruitBiomass = biomassCounter

        self.root.root1.mass = self.root.root1.mass + (self.rootBiomassWeight*abs(self.stoichiometricModel.reactions.get_by_id("Biomass_tx_R1").flux))
        self.rootBiomass = self.root.root1.mass

    def generateModelAtEmergence(self):

        """
        Generate whole plant model for tomato at emergence (1 leaf, 1 stem, 1 root)

        Parameters
        ----------
        - Organ-specific stoichiometric models (leaf, root, stem)
        - Initial organ biomasses (leaf, root, stem)

        Returns
        -------
        - Whole plant model at emergence

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        import cobra
        from cobra.core import Model
        from cobra import Reaction, Metabolite

        plantModel = Model()
        t = 1
        BiomassRxn = Reaction("Total_biomass_tx", name = "Total_biomass_tx", upper_bound = 1000000, lower_bound= 0)
        plantModel.add_reaction(BiomassRxn)

        leafModel = self.leafModel
        rootModel = self.rootModel
        stemModel = self.stemModel

        for model in [leafModel, rootModel, stemModel]:
            if model == leafModel:
                organTag = "LD"
                model = leafModel.copy()
            if model == rootModel:
                organTag = "R"
                model = rootModel.copy()
            if model == stemModel:
                organTag = "ST"
                model = stemModel.copy()
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
        self.addLeaf([Organ("leaf" + str(1), position=1, mass = self.parameters["LbInit"])])
        self.addRoot([Organ("root" + str(1), position=1, mass = self.parameters["RbInit"])])
        self.addStem([Organ("stem" + str(1), position=1, mass = self.parameters["SbInit"])])

        #Initiate total biomass for organ type attribute
        self.leafBiomass = self.leaves.leaf1.mass
        self.rootBiomass = self.root.root1.mass
        self.stemBiomass = self.stems.stem1.mass
        self.fruitBiomass = 0

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

        BiomassWeight = 0
        for met in self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites:
            BiomassWeight = BiomassWeight + abs(met.formula_weight*self.leafModel.reactions.get_by_id(rxnID + "LD").metabolites[met])
        self.leafBiomassWeight = BiomassWeight/1000

        BiomassWeight = 0
        for met in self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites:
            BiomassWeight = BiomassWeight + abs(met.formula_weight*self.stemModel.reactions.get_by_id(rxnID + "ST").metabolites[met])
        self.stemBiomassWeight = BiomassWeight/1000

        BiomassWeight = 0
        for met in self.rootModel.reactions.get_by_id(rxnID + "R").metabolites:
            BiomassWeight = BiomassWeight + abs(met.formula_weight*self.rootModel.reactions.get_by_id(rxnID + "R").metabolites[met])
        self.rootBiomassWeight = BiomassWeight/1000

    def maxNitrateUptake(self):

        """
        Function to calculate maximum nitrate uptake rate.

        Parameters
        ----------
        - Nitrate concentration
        - Nitrate uptake Vmax
        - Nitrate uptake Km
        - Root biomass

        Returns
        -------
        - Maximum nitrate uptake by root (organ biomass adjusted)

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        if self.soilNitrate > 1:
            num = float(self.parameters["NO3UptakeVm"]*self.soilNitrate)
            denom = float(self.parameters["NO3UptakeKm"]+self.soilNitrate)
        else:
            num = float(self.parameters["NO3UptakeVmLow"]*self.soilNitrate)
            denom = float(self.parameters["NO3UptakeKmLow"]+self.soilNitrate)

        MichMen = float(num/denom)
        MichMen = float(MichMen*self.root.root1.mass)
        self.stoichiometricModel.reactions.Nitrate_ec_R1.upper_bound = MichMen

    def generateStoichiometricModels(self, phloemComposition = "Shameer2018", sbmlFile = "PlantCoreMetabolism.xml", biomassEquations = "tomatoBiomassEquations.csv"):

        """
        Generate organ-specific stoichiometric models (diel leaf, root, stem, and fruit) from PlantCoreMetabolism_v1_2.xml

        Parameters
        ----------
        - Core stoichiometric model (PlantCoreMetabolism_v1_2.xml)
        - Organ-specific biomass compositions
        - Phloem composition option
            - Shameer2018: tomato composition reported in Shameer et al., 2018
            - free: unfixed composition
            - sucroseConstrained: only sucrosed fixed

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

        #Add amino acid degradation pathways for BCAA
        model.add_reaction(Reaction("ACETOACETATE_COA_LIGASE_RXN_m", name = "ACETOACETATE_COA_LIGASE_RXN", upper_bound = 1000, lower_bound = 0))
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
        model.add_reaction(Reaction("ACETOACETATE_COA_COA_mc", name = "ACETOACETATE_COA_COA_mc", upper_bound = 1000, lower_bound = -1000))
        model.reactions.ACETOACETATE_COA_COA_mc.add_metabolites({model.metabolites.ACETOACETYL_COA_m: -1,
                                                                      model.metabolites.CO_A_c: -1,
                                                                      model.metabolites.CO_A_m: 1,
                                                                      met: 1})
        model.add_reaction(Reaction("ACETYL_COA_ACETYLTRANSFER_RXN_c", name = "ACETYL_COA_ACETYLTRANSFER_RXN", upper_bound = 1000, lower_bound = -1000))
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
            rxn = Reaction("s" + met + "_biomass", upper_bound = 0, lower_bound = -1000, name = met + "_biomass")
            model.add_reaction(rxn)
            model.reactions.get_by_id("s" + met + "_biomass").add_metabolites({met_b: -1, model.metabolites.get_by_id(met + "_c"): 1})

        #Establishing the xylem (reactions defined as uptake)
        xylem_metabolites = ["KI", "CAII", "MGII", "NITRATE", "Pi", "SULFATE"]
        for met in xylem_metabolites:
            rxn = Reaction(met + "_xy")
            rxn.name = met + "_xy"
            rxn.subsystem = "Xylem"
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
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
        if phloemComposition == "free" or phloemComposition == "sucroseConstrained":
            for met in model.reactions.Phloem_output_tx.metabolites:
                if not((met.id.startswith("PROTON"))):
                    met_name = met.id[:-2]
                    rxn = Reaction(met_name + "_ph", name = met_name + "_ph", subsystem = "Phloem", lower_bound = 0, upper_bound = 1000)
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
            for met in ["MET"]: # ["HIS", "CYS", "PHE", "TRP", "TYR"]
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

        source_model = model.copy()

        #Create diel source leaf model
        diel_model = Model('Diel_model')
        model = source_model.copy()
        for met in model.metabolites:
            if not(met.id.endswith("_xy") or (met.id.endswith("_ph"))):
                met.id = met.id + "_LD"
                if met.compartment != "f":
                    met.compartment = met.compartment + "_LD"
        for rxn in model.reactions:
            rxn_copy = rxn.copy()
            rxn_copy.id = rxn_copy.id + "_LD"
            diel_model.add_reaction(rxn_copy)
        model = source_model.copy()
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
            rxn.lower_bound = -1000
            rxn.upper_bound = 1000
            rxn.add_metabolites({met_day: -1.0, met_night: 1.0})
            diel_model.add_reaction(rxn)
        MAL_v_day_night_linker = Reaction("MAL_v_linker_LD")
        MAL_v_day_night_linker.name = "MAL day to night linker"
        MAL_v_day_night_linker.subsystem = 'Day-Night linker reaction'
        MAL_v_day_night_linker.lower_bound = -1000
        MAL_v_day_night_linker.upper_bound = 1000
        MAL_v_day_night_linker.add_metabolites({diel_model.metabolites.MAL_v_LD: -0.7, diel_model.metabolites.aMAL_v_LD: -0.3,
                                                    diel_model.metabolites.MAL_v_LN: 0.7, diel_model.metabolites.aMAL_v_LN: 0.3})
        diel_model.add_reaction(MAL_v_day_night_linker)
        CIT_v_day_night_linker = Reaction("CIT_v_linker_LD")
        CIT_v_day_night_linker.name = "CIT day to night linker"
        CIT_v_day_night_linker.subsystem = 'Day-Night linker reaction'
        CIT_v_day_night_linker.lower_bound = -1000
        CIT_v_day_night_linker.upper_bound = 1000
        CIT_v_day_night_linker.add_metabolites({diel_model.metabolites.CIT_v_LD: -0.5, diel_model.metabolites.aCIT_v_LD: -0.5,
                                                    diel_model.metabolites.CIT_v_LN: 0.5, diel_model.metabolites.aCIT_v_LN: 0.5})
        diel_model.add_reaction(CIT_v_day_night_linker)
        Amino_acids = ["L_ALPHA_ALANINE", "L_ASPARTATE", "ARG", "ASN", "CYS", "GLN", "GLT", "GLY", "ILE",
                       "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "bHIS"]
        for aa in Amino_acids:
            met_day = diel_model.metabolites.get_by_id(aa + "_v_LD")
            met_night = diel_model.metabolites.get_by_id(aa + "_v_LN")
            rxn = Reaction(aa + "_v_linker_LD")
            rxn.name = aa + " day to night linker"
            rxn.subsystem = "Day-Night linker reaction"
            rxn.lower_bound = 0
            rxn.upper_bound = 1000
            rxn.add_metabolites({met_day: -1.0, met_night: 1.0})
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
        biomass_rxn = Reaction("Biomass_tx_LD", name = "Biomass_tx_LD", lower_bound = 0, upper_bound = 1000)
        diel_model.add_reaction(biomass_rxn)
        for row in biomass_file:
            met_id = row[0]+"_LD"
            stoich = float(row[1])
            diel_model.reactions.get_by_id("Biomass_tx_LD").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})
        biomass_file = csv.reader(open(biomassEquations,"r"))
        biomass_rxn = Reaction("Biomass_tx_LN", name = "Biomass_tx_LN", lower_bound = 0, upper_bound = 1000)
        diel_model.add_reaction(biomass_rxn)
        for row in biomass_file:
            met_id = row[0]+"_LN"
            stoich = float(row[1])
            diel_model.reactions.get_by_id("Biomass_tx_LN").add_metabolites({diel_model.metabolites.get_by_id(met_id): -stoich})

        #Light uptake constraints
        diel_model.reactions.Photon_tx_LD.upper_bound = 1000

        #Constrain starch phosphorylase reaction in leaf
        diel_model.reactions.G6P_Pi_pc_LD.upper_bound = diel_model.reactions.G6P_Pi_pc_LD.lower_bound = 0
        diel_model.reactions.G6P_Pi_pc_LN.upper_bound = diel_model.reactions.G6P_Pi_pc_LN.lower_bound = 0

        #Fix the flux ratio for phloem export day:night to 3:1
        if phloemComposition == "Shameer2018":
            Phloem_output_day = Metabolite("Phloem_LD_f", name = "Phloem output LD", compartment = "f")
            Phloem_output_night = Metabolite("Phloem_LN_f", name = "Phloem output LN", compartment = "f")
            diel_model.add_reaction(Reaction("Phloem_LD_LN_constraint_LD", name = "Phloem_LD_LN_constraint", upper_bound = 1000, lower_bound = 0))
            diel_model.reactions.Phloem_LD_LN_constraint_LD.add_metabolites({Phloem_output_day: -3, Phloem_output_night: -1})
            diel_model.reactions.Phloem_output_tx_LD.add_metabolites({Phloem_output_day: 1})
            diel_model.reactions.Phloem_output_tx_LN.add_metabolites({Phloem_output_night: 1})
        if phloemComposition == "free" or phloemComposition == "sucroseConstrained":
            phloemMetList = [rxn.id.replace("_ph_LD","") for rxn in diel_model.reactions if "ph_LD" in rxn.id]
            for met in phloemMetList:
                pseudoMet = Metabolite(met + "_LD_LN_constraint_LD", name=met + "_LD_LN_constraint", compartment="f")
                diel_model.reactions.get_by_id(met + "_ph_LD").add_metabolites({pseudoMet: 1})
                diel_model.reactions.get_by_id(met + "_ph_LN").add_metabolites({pseudoMet: -3})

        #Additional phloem constraints
        if phloemComposition == "sucroseConstrained":
            for organ in ["_LD", "_LN"]:
                SucroseMet = Metabolite("SUCROSE" + organ + "_f", name = "SUCROSE" + organ + "_f", compartment="f")
                AminoAcidMet = Metabolite("AMINO_ACID" + organ + "_f", name = "AMINO_ACID" + organ + "_f", compartment="f")
                for rxn in diel_model.reactions:
                    if "_ph" + organ in rxn.id:
                        if rxn.id == "SUCROSE_ph" + organ:
                            rxn.add_metabolites({SucroseMet: 1})
                        else:
                            rxn.add_metabolites({AminoAcidMet: 1})
                diel_model.add_reaction(Reaction("AMINO_ACID_VALVE" + organ, name = "AMINO_ACID_VALVE" + organ, upper_bound = 1000, lower_bound = 0))
                diel_model.reactions.get_by_id("AMINO_ACID_VALVE" + organ).add_metabolites({AminoAcidMet: 1})
                diel_model.add_reaction(Reaction("Phloem_SUC_AA_constraint" + organ, name = "Phloem_SUC_AA_constraint" + organ, upper_bound = 1000, lower_bound = 0))
                diel_model.reactions.get_by_id("Phloem_SUC_AA_constraint" + organ).add_metabolites({AminoAcidMet: -0.284, SucroseMet: -0.716})

        ######################################
        #Generate models for all sink organs #
        ######################################

        #Create sink models for root (R) and stem (S)
        sink_tissues = ["R", "ST"]
        for tissue in sink_tissues:
            model = source_model.copy()
            for rxn in model.reactions:
                if rxn.id.endswith("_ph"):
                    rxn.upper_bound = 0
                    rxn.lower_bound = -1000
                    rxn.add_metabolites({model.metabolites.PROTON_e: 2, model.metabolites.PROTON_c:-2})
                if "Phloem_output_tx" in rxn.id:
                    rxn.upper_bound = 0
                    rxn.lower_bound = -1000
                    rxn.add_metabolites({model.metabolites.PROTON_e: 1.9207920792, model.metabolites.PROTON_c:-1.9207920792})
                rxn.id = rxn.id + "_" + tissue
            for met in model.metabolites:
                if not(met.id.endswith("_xy")) and not(met.id.endswith("_ph")):
                    met.id = met.id + "_" + tissue
                    if met.compartment != "f":
                        met.compartment = met.compartment + "_" + tissue
            biomass_rxn = Reaction("Biomass_tx_" + tissue, name = "Biomass_tx_" + tissue, lower_bound = 0, upper_bound = 1000)
            model.add_reaction(biomass_rxn)
            if tissue == "ST":
                biomass_file = csv.reader(open(biomassEquations,"r"))
                for row in biomass_file:
                    model.reactions.get_by_id("Biomass_tx_" + tissue).add_metabolites({model.metabolites.get_by_id(row[0]+"_"+tissue): -float(row[2])})
                model_S = model.copy()
            elif tissue == "R":
                biomass_file = csv.reader(open(biomassEquations,"r"))
                for row in biomass_file:
                    model.reactions.get_by_id("Biomass_tx_" + tissue).add_metabolites({model.metabolites.get_by_id(row[0]+"_"+tissue): -float(row[3])})
                #Allow for nutrient uptake
                exchange_rxn = ["Pi", "Nitrate", "SO4", "Ca", "Mg", "K"]
                for rxn in exchange_rxn:
                    model.reactions.get_by_id(rxn + "_tx_R").upper_bound = 1000
                    model.reactions.get_by_id(rxn + "_tx_R").lower_bound = -1000
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
                        rxn.lower_bound = -1000
                        #Adjust cost of xylem loading
                        for met in rxn.reactants:
                            if "PROTON" in met.id:
                                rxn.add_metabolites({met: abs(rxn.metabolites[met])})
                        for met in rxn.products:
                            if "PROTON" in met.id:
                                rxn.add_metabolites({met: -rxn.metabolites[met]})
                model_R = model.copy()

        #Increase flux bounds for all models
        for organModel in [diel_model, model_R, model_S]:
            for rxn in organModel.reactions:
                if rxn.upper_bound == 1000:
                    rxn.upper_bound = 1000000
                if rxn.lower_bound == -1000:
                    rxn.lower_bound = -1000000

        self.leafModel = diel_model
        self.rootModel = model_R
        self.stemModel = model_S

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
                rxn.lower_bound = 0

        #Calculate net CO2 assimilation from PPFD according to empirical data (Nunes-Nesi, 2005)
        Pmax1 = (self.parameters["Pmax"]/self.parameters["photoperiod"]/self.parameters["leafAreatoBiomass"]/60/60)*1000
        netCO2uptake = (40.9962482*Pmax1)/(Pmax1 + 644.48886704)
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

        self.parameters["VcMax"] = (prev*self.parameters["leafAreatoBiomass"]*60*60*self.parameters["photoperiod"])/1000

    def maintRespirationConstraint(self):

        """
        Update maintenance respiration constraint in each organ

        Parameters
        ----------
        - Maintenance respiration requirements

        Returns
        -------
        - Updated maintenance respiration constraints

        Author: Tarit Konuntakiet
        Email: tarit.konuntakiet@seh.ox.ac.uk
        """

        for leaf in self.leaves:
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LD" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*self.parameters["photoperiod"]
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_LN" + str(leaf.position)).lower_bound = (self.parameters["leafMaint"]*leaf.mass/24)*(24-self.parameters["photoperiod"])
        for stem in self.stems:
            self.stoichiometricModel.reactions.get_by_id("ATPase_tx_ST" + str(stem.position)).upper_bound = self.stoichiometricModel.reactions.get_by_id("ATPase_tx_ST" + str(stem.position)).lower_bound = self.parameters["stemMaint"]*stem.mass
        self.stoichiometricModel.reactions.ATPase_tx_R1.upper_bound = self.stoichiometricModel.reactions.ATPase_tx_R1.lower_bound = self.parameters["rootMaint"]*self.root.root1.mass

    def canopyEffect(self):
        import math

        a = self.parameters["Pmax"]*self.parameters["k"]*(math.e**(-self.parameters["k"]*self.F))
        b = 1-self.parameters["m"]
        return (a/b)/self.parameters["Pmax"]

    def maxLightUptakeCarbonAssimilation(self, t):

        if t<8: #14
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomass"]*15
        elif t>=8 and t<45: #14-45
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomass"]*3.75
        else:
            self.F = self.leafBiomass*self.parameters["leafAreatoBiomass"]*2.1
        canopyReductionFactor = self.canopyEffect()
        #Constrain maximum light uptake rate and RuBisCO carboxylation rate by leaves (mmol/leaves/day)
        if t < self.parameters["CotAge"]-self.parameters["InitDAS"]: #Cotyledons are active
            self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD1").upper_bound = self.parameters["VcMax"]*(self.leafBiomass+self.parameters["CotWeightAdj"])*canopyReductionFactor
            self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD1").upper_bound = self.parameters["Pmax"]*(self.leafBiomass+self.parameters["CotWeightAdj"])*canopyReductionFactor
        else: #Cotyledons have dropped off
            self.stoichiometricModel.reactions.get_by_id("RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_LD1").upper_bound = self.parameters["VcMax"]*(self.leafBiomass)*canopyReductionFactor
            self.stoichiometricModel.reactions.get_by_id("Photon_tx_LD1").upper_bound = self.parameters["Pmax"]*(self.leafBiomass)*canopyReductionFactor

class Organ:

    """Create new organ objects with age, mass, and position attributes"""

    def __init__(self,id, position, mass = 0, age = 1, stage = 0):
        self.id = id
        self.age = age
        self.mass = mass
        self.position = position
        self.stage = stage

    def __repr__(self):
        return '%s \nPosition: %d \nAge: %d \nMass: %.4f \nStage: %d' %(self.id, self.position, self.age, self.mass, self.stage)