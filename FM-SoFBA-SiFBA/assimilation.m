function [suc_sta_base,sta_use,suc_equi,al_suc,suc_c_disp,suc_c_interm,al_pt_plant_assim]...
    = assimilation(daylength,is_light,net_rate,timestep,leaf_res,...
    root_res,suc_c_perplant,rosette_area,sta_c_endday, Sucrose_biomass_day_flux,...
    Starch_accumulation_flux, Sucrose_biomass_night_flux, net_rate_dielFBA)

global p

convert_to_gC = timestep*p(59)*10^(-6)*12; %conversion factor for micromol/m2 leaf/sec to gC/m2 leaf/timestep

%Baseline conversion coefficient (default p(60)=0.125)
sta_base = p(60); %Baseline starch conversion coefficient

%Starch turnover (default p(61)=0.84)
if daylength == 18
    sta_convert_night = 0.6; %Based on observation in Sulpice et al (2013)
else
    sta_convert_night = p(61);
end

if      is_light == 1 %Daytime

        sta_use = 0; %conversion of starch to sugar
        al_pt_plant_assim = net_rate*rosette_area*convert_to_gC; %Assimilatory flux per plant (gC/plant/timestep)

        %Calculate daytime sucrose export rate
        Vc_ratio = net_rate/net_rate_dielFBA;
        sucrose_export_day = 3600*rosette_area*Sucrose_biomass_day_flux*Vc_ratio; %umol/plant/hr
        al_suc = (sucrose_export_day*12*12)/(10^6); %gC/plant/hr

        %Calculate daytime starch accumulation rate
        starch_accumulation_rate = 3600*rosette_area*Starch_accumulation_flux*Vc_ratio; %umol/plant/hr
        suc_sta_base = (starch_accumulation_rate*12*6)/(10^6); %gC/plant/hr

else    %Night-time
        al_pt_plant_assim = 0;
        suc_sta_base = 0;
        al_suc = 0;

        %Calculate night-time sucrose export rate
        Sucrose_biomass_night_flux1 = convert_to_gramsC(Sucrose_biomass_night_flux, 12, rosette_area);
        Starch_accumulation_flux1 = convert_to_gramsC(Starch_accumulation_flux, 6, rosette_area);
        sta_use = (Sucrose_biomass_night_flux1*sta_c_endday)/(Starch_accumulation_flux1*(24-daylength));

end

suc_equi = p(62); %equilibrium sucrose plus hexose concentration in leaves (g C/m2 leaf)
suc_c_interm = suc_c_perplant + sta_use + al_suc - root_res - leaf_res;

current_value = suc_c_interm - (suc_equi*rosette_area);

if      current_value <= 0
        suc_c_disp = 0;
else
        suc_c_disp = current_value; %Amount of sugar available for growth
end

end
