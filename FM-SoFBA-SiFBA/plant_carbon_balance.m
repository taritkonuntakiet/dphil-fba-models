function [rlc_pt1,rrc_pt1,leaf_res,root_res,leaf_carbon,root_carbon,sucrose_carbon,starch_carbon,rgtotal,...
    rmtotal,totalCarbon,Assim] = plant_carbon_balance(Tleaf,...
    CO2,PAR,sunrise,sunset,is_light, rsratio,rosette_area,leaf_c,root_c,suc_c,sta_c,rgtot,rmtot,timestep,...
    sta_c_endday, Sucrose_biomass_day_flux, Starch_accumulation_flux,...
    Sucrose_biomass_night_flux,net_rate_dielFBA, sinkFBA_model,t)

global p

%Initial values
%______________

vlmax20 = p(38); %Photosynthetic rubisco capacity per unit leaf area at 20degC (micromol CO2 m-2 s-1)
vlmax25 = vlmax20/0.64;

%Initial setting
%_______________

daylength = sunset - sunrise;
      
        
%Calculating photosynthesis
%__________________________
        
if      is_light == 1
        
        [net_rate] = photosynthesis(CO2,Tleaf,PAR,vlmax25,daylength);
else
        net_rate = 0;
end    
        
        
%Calculating maintenance respiration
%___________________________________
        
[leaf_res,root_res] = mainres(Tleaf,leaf_c,root_c,suc_c,rosette_area,timestep);
        
        
%Calculating Assimilatory flux
%_____________________________     
      
[suc_sta_base,sta_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,Assim]...
    = assimilation(daylength,is_light,net_rate,timestep,leaf_res,...
    root_res,suc_c,rosette_area,sta_c_endday, Sucrose_biomass_day_flux,...
    Starch_accumulation_flux, Sucrose_biomass_night_flux, net_rate_dielFBA);

%Calculating organ demand
%________________________    
	
[rrc_pt,totdem,rlc_pt,root_growth,leaf_growth] = organdemand(timestep,...
    rsratio,leaf_c, rosette_area, sinkFBA_model);

%Calculating allocation 
%______________________

[rlc_pt1,suc_sta,root_gro1,rrc_pt1,leaf_gro1]...
= allocation(rrc_pt,totdem,rlc_pt,root_growth,leaf_growth,suc_c_disp);

%Calculating translocation
%_________________________

[leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c);

%Calculating amount in each carbon pool
%______________________________________

leaf_carbon = leaf_c + leaf_gro1 - leaf_trans; % Leaf carbon: g C/plant
starch_carbon = sta_c + suc_sta_base + suc_sta - sta_use; %Starch content per plant: gC/plant
root_carbon = root_c + root_gro1 - root_trans; %Root carbon: gC/plant
sucrose_carbon = suc_c - rlc_pt1 - rrc_pt1 - leaf_res...
                 - root_res - root_gro1 + root_trans...
                 + al_suc - leaf_gro1 + leaf_trans...
                 - suc_sta + sta_use; %Sucrose carbon: gC/plant
rgtotal = rgtot + rlc_pt1 + rrc_pt1; %Total growth respiration: gC/plant
rmtotal = rmtot + leaf_res + root_res; %Total maintenance respiration: gC/plant 
        
totalCarbon = leaf_carbon + root_carbon + rmtotal + rgtotal + starch_carbon + sucrose_carbon; %total C (g/plant)

end


