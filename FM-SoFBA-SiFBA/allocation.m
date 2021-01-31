function [rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1] = allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp)

if      totdem < suc_c_disp

    suc_growth = totdem;
    rrc_pt1 = rrc_pt;
    root_gro1 = root_growth;
    rlc_pt1 = rlc_pt;
    leaf_gro1 = leaf_growth;
    suc_sta = 0; %turn off starch overflow mechanism

else
    suc_growth = suc_c_disp; %sucrose for growth (g C/plant)
    rrc_pt1 = rrc_pt*(suc_growth/totdem); % Root growth respiration In g C/plant/timestep (actually used)
    root_gro1 = root_growth*(suc_growth/totdem); %root growth actually achieved
    rlc_pt1 = rlc_pt*(suc_growth/totdem); % Leaf growth respiration In g C/plant/timestep (actually used)
    leaf_gro1 = leaf_growth*(suc_growth/totdem); %leaf growth actually achieved (g C/plant/timestep)
    suc_sta = 0;
end