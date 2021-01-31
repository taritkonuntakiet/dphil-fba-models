function [leaf_res,root_res] = mainres(Tleaf,leaf_c,root_c,suc_c_perplant,rosette_area,timestep)

global p

%Calculation for maintenance respiration
%_______________________________________
%_______________________________________

act_ener_res20 = p(56); %activation energy for rl 20 (kJ/mol)

suc_conc_s = suc_c_perplant/rosette_area;  

rl20_leafres = (p(57)*suc_conc_s + p(58))*24; %Leaf respiration at 20 (g CO2 C/m2/day)

numres = act_ener_res20*(Tleaf - 20);
denomres = 293*p(39)*(Tleaf + 273);

if  suc_conc_s <= 0
    
    rl_leafres = 0;
else
    rl_leafres = rl20_leafres*exp(numres/denomres); %leaf respiration at leaf temperature (g CO2 C/m2/day)
end

leaf_res = rl_leafres*rosette_area*timestep; %leaf respiration per plant per time step (gC per plant per time step)
root_res = leaf_res*root_c/leaf_c; %gC/plant/time step
