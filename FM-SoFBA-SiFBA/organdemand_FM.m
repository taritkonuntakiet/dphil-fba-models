function [rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand_FM(timestep,rsratio,leaf_c)

global p

%Leaf growth
%___________

maxgrowth = p(63); 
leaf_growth = maxgrowth*leaf_c*timestep; %Maximal leaf growth (gC/plant/timestep)

rc_coef = p(64); %Growth respiration to total growth allocation

numrlc = leaf_growth*rc_coef;
denomrlc = 1-rc_coef;

rlc_pt = numrlc/denomrlc; % Growth respiration In g C/plant/timestep (required)

%Root growth
%___________

root_growth = leaf_growth*rsratio; %Root growth (g C/plant/time step)

numrrc = root_growth*rc_coef;
denomrrc = 1-rc_coef;

rrc_pt = numrrc/denomrrc; % Growth respiration In g C/plant/timestep (required)

%Total growth demand
%___________________

totdem = leaf_growth + root_growth + rlc_pt + rrc_pt;
