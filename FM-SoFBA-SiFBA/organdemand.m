function [rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand(timestep,...
    rsratio,leaf_c,rosette_area, sinkFBA_model)

global p

%Leaf growth
%___________

maxgrowth = p(63); %0.4080
rlc_pt = 0;

Biomass_tx = 'AraCore2014_Biomass_tx';
carbons_conversion = 155.6411422;

%Run leaf sink pFBA
%____________________
leaf_growth = maxgrowth*leaf_c*timestep; %gC/plant/hour
leaf_growth = convert_to_micromol(leaf_growth, carbons_conversion,rosette_area); %Convert to micromol/m2/s
sinkFBA_model = changeRxnBounds(sinkFBA_model, Biomass_tx, leaf_growth, 'u');
sinkFBA_model = changeRxnBounds(sinkFBA_model, Biomass_tx, leaf_growth, 'l');

sol = pFBA(sinkFBA_model, 'max');
Phloem_uptake_index = find(sinkFBA_model.rxns=="sSUCROSE_phloem"); %Phloem_output_tx sSUCROSE_phloem
Phloem_uptake_leaf_flux = sol.x(Phloem_uptake_index); %micromol/m2/s
Phloem_uptake_leaf_flux = convert_to_gramsC(Phloem_uptake_leaf_flux,12, rosette_area);
leaf_growth = Phloem_uptake_leaf_flux; %gC/plant/hour
CO2_tx_index = find(sinkFBA_model.rxns=="CO2_tx");
CO2_respiration_leaf = abs(sol.x(CO2_tx_index));
CO2_respiration_leaf = convert_to_gramsC(CO2_respiration_leaf, 1, rosette_area);

%Root growth
%___________

root_growth = leaf_growth*rsratio; %Root growth (g C/plant/timestep)
rrc_pt = 0;

%Run root sink pFBA
%____________________
root_growth = convert_to_micromol(root_growth, carbons_conversion,rosette_area); %Convert to micromol/m2/s
sinkFBA_model = changeRxnBounds(sinkFBA_model, Biomass_tx, root_growth, 'u');
sinkFBA_model = changeRxnBounds(sinkFBA_model, Biomass_tx, root_growth, 'l');

sol = pFBA(sinkFBA_model, 'max');
Phloem_uptake_root_flux = sol.x(Phloem_uptake_index);
Phloem_uptake_root_flux = convert_to_gramsC(Phloem_uptake_root_flux,12, rosette_area);
root_growth = Phloem_uptake_root_flux; %gC/plant/hour
CO2_tx_index = find(sinkFBA_model.rxns=="CO2_tx");
CO2_respiration_root = abs(sol.x(CO2_tx_index));
CO2_respiration_root = convert_to_gramsC(CO2_respiration_root, 1, rosette_area);

%Total growth demand
%___________________

totdem = leaf_growth + root_growth + rlc_pt + rrc_pt;
rlc_pt = CO2_respiration_leaf;
rrc_pt = CO2_respiration_root;

end
