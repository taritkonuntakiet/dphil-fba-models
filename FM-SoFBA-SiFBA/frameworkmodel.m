function [ShootFreshWeight] = frameworkmodel(data_set)
%data_set = 1 for growth conditions simulated in DPhil thesis for
%Konuntakiet, T. (2021), University of Oxford

% genotype = 1 for Ler, 2 for Col
% water content can be specified at line 362
% flowering time can be tuned either directly at line 63 (the photothermal model in line 62 can be optionally switched off) or
%adjusting the MPTU threshold in the phenology function file (Line 23 for Ler and Line 34 for Col in phenology.m)
%dataset = 1 for Chew et al. (2014), 2 for Sulpice et al. (2014), 3 for
%Rubio-Asensio et al. (2015)

%Calling for parameters
%______________________

global p
global v
load('parameter.mat')
load('vpaper.mat')
p=parameter;
v=vpaper;

probability = 1;
genotype = 2;

%Calling for meteorological data
%_______________________________

if data_set == 1
    temperature = 21.3;
    sunrise = 6;
    sunset = 18;
    co2 = 37.5;
    light = 110;
    load('weather.mat')
    hour=weather(:,1);
    T=temperature*weather(:,2);
    sunrise=sunrise*weather(:,3);
    sunset=sunset*weather(:,4);
    CO2=co2*weather(:,5); %CO2 partial pressure (Pa)
    PAR=light*weather(:,6); %total absorbed PAR per unit leaf area (micromol m-2 s-1)
    photoperiod = sunset(1) - sunrise(1);
end

%Parameters for vegetative growth component
%__________________________________________

Tb = p(1); %Base temperature for thermal-time calculation
TT0 = p(19); %degree days required from seed germination until plant emergence (stage 1.0 in Boyes et al. 2001)
PL=p(20); %sink strength of leaf
aL=p(21); %parameter for leaf sink variation (beta law)
bL=p(22); %parameter for leaf sink variation (beta law)
PR=p(23); %sink strength of root
aR=p(24); %parameter for root sink variation (beta law)
bR=p(25); %parameter for root sink variation (beta law)
TS = p(26); %degree days of leaf lifespan (from the point of leaf initiation)
TLcot = p(27); %cotyledon expansion period
TLtrue = p(28); %expansion period of true leaves
Root_LC = p(29); %fold of flowering time for the root expansion period (set to 1.3)
leaf_factor = p(30); %carbon-dry mass conversion in leaves
root_factor = p(31); %carbon-dry mass conversion in the root
seed_input = p(32); %seed biomass
cot_area = p(33); %cotyledon area at emergence
sucrose_initial = p(34); %initial sucrose content per unit area (from Rasse and Tocquin)
starch_initial = p(35); %initial starch content per initial sucrose content (from Rasse and Tocquin)
min_alpha = p(65); %minimum zenithal angle
inc_alpha = p(66); %maximum increase in zenithal angle from the minimum value
juv_TT = p(67); %thermal time since plant emergence for the juvenile-to-adult transition
SLA_cot = p(71); %SLA for cotyledons (in m2/g dry mass) Fig. 3 Christophe et al (2008)
SLA_exp = p(72); %parameter for SLA
juv_phyll = p(82); %phyllochron at the juvenile stage
ad_phyll = p(83); %phyllochron at the adult stage

ShootFreshWeight = [];

%Determining flowering (bolting) time
%____________________________________

%[input_size] = phenology(hour,T,sunrise,sunset,genotype) %Calling for the photothermal model
input_size = 40*24; %hours to flowering


%To calculate thermal time and determine the seedling emergence point
%____________________________________________________________________

for n=1:input_size
    Thrm(n)=max(0,(T(n)-Tb))/24; %in degree days
    CumThrm(n) = sum(Thrm(1:n)); %cumulative thermal time
    
    if CumThrm(n) < TT0
        eme = n+1;
    end
end


%Expansion period for the root system
%____________________________________

if genotype == 1 %Ler
    TR = Root_LC*CumThrm(end);
elseif genotype == 2 %Col
    TR = 963.495;
end

MaxfRpoint = TR/(1+(bR-1)/(aR-1))-0.5; %determining the maximum point by solving the function derivative
MaxfR = ((MaxfRpoint+0.5)/TR)^(aR-1)*(1-(MaxfRpoint+0.5)/TR)^(bR-1);

%Maximum values of leaf sink variation
%_____________________________________

Maxfcotpoint = TLcot/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
Maxfcot = ((Maxfcotpoint+0.5)/TLcot)^(aL-1)*(1-(Maxfcotpoint+0.5)/TLcot)^(bL-1);

MaxfLpoint = TLtrue/(1+(bL-1)/(aL-1))-0.5; %determining the maximum point by solving the function derivative
MaxfL = ((MaxfLpoint+0.5)/TLtrue)^(aL-1)*(1-(MaxfLpoint+0.5)/TLtrue)^(bL-1);


%For the first growth cycle (from germination to seedling emergence at TT0):
%__________________________________________________________________________

Seed_mass = seed_input; %Seed dry mass (g), emptied after the first cycle to reach stage 1.0
Q(eme) = Seed_mass;
fR = ((CumThrm(eme)+0.5)/TR)^(aR-1)*(1-(CumThrm(eme)+0.5)/TR)^(bR-1)/MaxfR;
D(eme) = PR*fR + 2*PL;
Root_mass(eme) = PR*fR*Q(eme)/D(eme); %root dry mass (g)

Si(eme,1:2) = cot_area; %area for EACH cotyledon (m2), cycle 1 for leaf ranked 1 and 2
Leaf_mass(eme,1:2) = Si(eme,1)/SLA_cot; %dry mass for EACH cotyledon (g)
Hypocotyl_mass = Seed_mass - Root_mass(eme) - sum(Leaf_mass(eme,:));
Total_shoot(eme) = Seed_mass - Root_mass(eme);
Destructive_area(eme) = sum(Si(eme,:));
Leaf_no(eme) = 2;
Appear(1:2) = eme; %The timepoint when the cotyledons start expanding using photosynthate

%Initialisation
%______________

S(eme) = sum(Si(eme,1:2));
S_intercept(eme) = S(eme)*cos(10/180*pi);
Func_area(eme)=S(eme);
steps_perhour = 1;
steps_perday = 24*steps_perhour;
timestep = 1/steps_perday; % Time step length used to solve the model

Leaf_carbon(eme) = leaf_factor*sum(Leaf_mass(eme,:));
Root_carbon(eme) = root_factor*Root_mass(eme);
Sucrose_carbon(eme) = sucrose_initial*S(eme); %Initial amount of sucrose carbon (g/plant)
Starch_carbon(eme) = starch_initial*Sucrose_carbon(eme); %Initial amount of starch carbon (g/plant)
rsratio(eme) = Root_carbon(eme)/Leaf_carbon(eme);

%To determine the first end-of-day starch level in the case where emergence
%is modelled to occur at night:
if      hour(eme) < sunrise(eme)
    sta_c_endday = Starch_carbon(eme)/(p(61)*((sunrise(eme)-hour(eme))/(sunrise(eme)+24-sunset(eme-1)))+1-p(61));
elseif  hour(eme) > sunset(eme)
    sta_c_endday = Starch_carbon(eme)/(p(61)*((24-hour(eme)+sunrise(eme+1))/(24-sunset(eme)+sunrise(eme+1)))+1-p(61));
else
    sta_c_endday = 0;
end

%To initialise the Rasse and Tocquin carbon dynamic model:
[rlc_pt1(eme),rrc_pt1(eme),leaf_res(eme),root_res(eme),rgtot(eme),rmtot(eme),totalC(eme),assim(eme)] = ...
    ini_carbon_balance(T(eme),CO2(eme),PAR(eme),...
    sunrise(eme),sunset(eme),rsratio(eme),...
    S_intercept(eme),Leaf_carbon(eme),...
    Root_carbon(eme),Sucrose_carbon(eme),...
    Starch_carbon(eme),timestep);

GC = 1; %number of growth cycle elapsed
GCstart(GC+1) = eme; %The timepoint when the second growth cycle starts

%Load SoFBA, a core C3 diel FBA model with 1.5:1 phloem export ratio over the diel
%cycle
%______________________________________________________________________________________
dielFBA_model = readCbModel('SoFBA.xml');

%Adjust light and net rate to account for the photoperiod
%________________________________________________________
night_length = 24-photoperiod;
light = (light*photoperiod)/12;
vlmax25 = p(37);
[net_rate_adjusted] = photosynthesis(CO2(505),T(505),PAR(505),vlmax25,photoperiod);
net_rate_adjusted = (net_rate_adjusted*photoperiod)/12;

%Initialise SoFBA and run pFBA to acquire flux predictions for RuBisCO
%carboxylase rate, starch accumulation, and phloem output
%_________________________________________________________________________
dielFBA_model = changeRxnBounds(dielFBA_model, "Photon_tx_day", light, 'u');
dielFBA_model = changeRxnBounds(dielFBA_model, "Photon_tx_day", light, 'l');
estimatedVc = estimateVcFromNetRate(dielFBA_model, net_rate_adjusted);
dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", estimatedVc, 'u');
dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", estimatedVc, 'l');
sol = pFBA(dielFBA_model, 'max');

%Record fluxes
%______________
%Net rate
CO2_tx_index = find(dielFBA_model.rxns=="CO2_tx_day");
net_rate_dielFBA = sol.x(CO2_tx_index);
net_rate_dielFBA = (net_rate_dielFBA*photoperiod)/12;
%Phloem output (day)
Sucrose_biomass_day_index = find(dielFBA_model.rxns=="sSUCROSE_biomass_day");
Sucrose_biomass_day_flux = sol.x(Sucrose_biomass_day_index);
Sucrose_biomass_day_flux = (Sucrose_biomass_day_flux*12)/photoperiod;
%Starch accumulation
Starch_accumulation_index = find(dielFBA_model.rxns=="Starch_p_day_night_transfer");
Starch_accumulation_flux = sol.x(Starch_accumulation_index);
Starch_accumulation_flux = (Starch_accumulation_flux*12)/photoperiod;
%Phloem output (night)
Sucrose_biomass_night_index = find(dielFBA_model.rxns=="sSUCROSE_biomass_night");
Sucrose_biomass_night_flux = sol.x(Sucrose_biomass_night_index);
Sucrose_biomass_night_flux = (Sucrose_biomass_night_flux*12)/night_length;

%Load SiFBA: core C3 sink model
%_________________________
sinkFBA_model = readCbModel('SiFBA.xml');

%To calculate growth at each time point t
%_________________________________________

for t=eme+1:input_size
    
    %To determine the rosette area at the previous time point for light interception
    %_________________________________________________________________________________
    
    
    %Determining the status (functioning/senesced) for each leaf rank i at the previous timepoint:
    
    for i = 1:Leaf_no(t-1)
        if CumThrm(t-1) < (CumThrm(Appear(i))+TS)
            fS(i) = 1;
        else
            fS(i) = 0;
        end
    end
    
    %To determine rosette structure (at the previous time point):
    
    [Max_area,i_max] = max(Si(t-1,:));
    i_max = max(i_max,2);
    
    for i=1:Leaf_no(t-1) %for each leaf rank i
        
        if      i<=i_max
            
            alpha(i)=min_alpha; %zenithal angle
        else
            alpha(i)=min_alpha+inc_alpha*(i-i_max)/(Leaf_no(t-1)-i_max);
        end
        
        S(t-1,i)=Si(t-1,i)*cos(alpha(i)/180*pi)*fS(i); %Projected functioning leaf area
        Sf(t-1,i)=Si(t-1,i)*fS(i); %Functioning leaf area
        
    end
    
    
    %Rosette area (at the previous time point):
    
    if  Leaf_no(t-1) <= 15
        
        S_intercept(t-1) = sum(S(t-1,1:Leaf_no(t-1))); %Projected rosette area
        
    else
        maxblade = sort(S(t-1,:),'descend');
        S_intercept(t-1) = sum(maxblade(1:13));
    end
    
    %Total functioning area (at the previous time point):
    
    Func_area(t-1)=sum(Sf(t-1,1:Leaf_no(t-1)));
    Destructive_projected_area(t-1) = sum(S(t-1,1:Leaf_no(t-1)));
    
    
    %To determine whole-plant carbon balance (at the current time point):
    %___________________________________________________________________
    
    
    if CumThrm(t) < (juv_TT+TT0)
        phyllochron = juv_phyll*probability;
    else
        phyllochron = ad_phyll*probability;
    end
    
    %To determine current leaf number:
    if  (CumThrm(t) - CumThrm(GCstart(end))) >= phyllochron % A growth cycle has been reached
        
        GC = GC+1;
        GCstart(GC+1) = t; % the timepoint when a new growth cycle starts
        
        Leaf_no(t) = Leaf_no(t-1) + binornd(1,probability); % A new leaf may/may not be initiated
        
        if Leaf_no(t) ~= Leaf_no(t-1) % A new leaf is initiated
            Appear(Leaf_no(t)) = t;
            Leaf_mass(t-1,Leaf_no(t))=0;
            Si(t-1,Leaf_no(t))=0;
        end
        
    else
        Leaf_no(t) = Leaf_no(t-1); %No new leaf initiated
    end
    
    
    %To calculate sink variation for leaf:
    for i = 1:Leaf_no(t)
        
        if i <= 2
            TL = TLcot;
            Maxf = Maxfcot;
            
            fL(t,i) = ((CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(aL-1)...
                *(1-(CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(bL-1);
            
            if real(fL(t,i))~=fL(t,i);
                fL(t,i) = 0;
            end
            
            fL(t,i) = fL(t,i)/Maxf;
            
        else
            TL =TLtrue;
            Maxf = MaxfL;
            
            fL(t,i) = ((CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(aL-1)...
                *(1-(CumThrm(t)-CumThrm(Appear(i))+0.5)/TL)^(bL-1);
            
            if real(fL(t,i))~=fL(t,i)
                fL(t,i) = 0;
            end
            
            fL(t,i) = fL(t,i)/Maxf;
            
        end
        
        
    end
    
    %To calculate sink variation for root:
    fR(t) = ((CumThrm(t)+0.5)/TR)^(aR-1)*(1-(CumThrm(t)+0.5)/TR)^(bR-1);
    
    fR(t) = fR(t)/MaxfR; %sink strength
    
    
    %Calculating root-to-shoot ratio (in g Carbon):
    leafconversion = leaf_factor;     % Gorsuch et al 2010a and 2010b
    
    rootconversion = root_factor; %g Carbon per g dry mass
    %from Kumar et al (Agroforest Syst 2010) around 30-35%
    %from UK BioChar Research Centre (UoE) around 35 % in
    %wheat, cited as Prendergast-Miller M and Sohi SP 2010.
    %Investigating biochar impacts on plant roots and root carbon.
    %Poster presented in the Organic Matter Stabilization
    %and Ecosystem Functions session at Soil Organic
    %Matter Conference, Cote d'Azur, France (Sept 2010)
    
    num   = fR(t)*PR*rootconversion;
    denom = sum(fL(t,:))*PL*leafconversion;
    rsratio(t) = num/denom;
    
    % To determine light status:
    
    if  hour(t) > sunrise(t) && hour(t) <= sunset(t)
        is_light(t) = 1;
    else
        is_light(t) = 0;
    end
    
    %To determine starch status:
    
    if      is_light(t) == 1
        sta_c_endday = 0;
        
    elseif  is_light(t) == 0 && is_light(t-1) == 1
        sta_c_endday = Starch_carbon(t-1);
    else
        sta_c_endday = sta_c_endday;
    end
    
    %Calling for the Rasse and Tocquin carbon dynamic model:
    [rlc_pt1(t),rrc_pt1(t),leaf_res(t),root_res(t),...
        Leaf_carbon(t),Root_carbon(t),Sucrose_carbon(t),Starch_carbon(t),rgtot(t),rmtot(t),totalC(t),assim(t),...
        ] = plant_carbon_balance(T(t),CO2(t),PAR(t),...
        sunrise(t),sunset(t),is_light(t),rsratio(t),S_intercept(t-1),Leaf_carbon(t-1),Root_carbon(t-1),...
        Sucrose_carbon(t-1),Starch_carbon(t-1),rgtot(t-1),rmtot(t-1),timestep,sta_c_endday,...
        Sucrose_biomass_day_flux, Starch_accumulation_flux, Sucrose_biomass_night_flux,...
        net_rate_dielFBA, sinkFBA_model,t);
    
    %Root mass:
    Root_mass(t) = Root_carbon(t)/rootconversion; %root dry mass (g)
    %Actual root C growth
    CarbonR = Root_carbon(t) - Root_carbon(t-1);
    
    %Total dry biomass accumulation for leaves (g):
    CarbonL = Leaf_carbon(t) - Leaf_carbon(t-1);
    BiomassL = CarbonL/leafconversion;
    
    %Leaf demand:
    for     i = 1:Leaf_no(t)
        Ld(t,i) = PL*fL(t,i);
    end
    
    Totalleafdemand = sum(Ld(t,:));
    
    
    %To calculate SLA:
    SLA(t) = SLA_cot*exp(-SLA_exp*(CumThrm(t)-TT0)); %SLA (in m2/g dry mass) Fig. 3 Christophe et al (2008)
    

    %Leaf mass and area:
    for     i = 1:Leaf_no(t)
        Leaf_mass(t,i) = Leaf_mass(t-1,i) + Ld(t,i)/Totalleafdemand*BiomassL;
        Si(t,i) = max(Leaf_mass(t,i)*SLA(t),Si(t-1,i));
    end
    
    Destructive_area(t) = sum(Si(t,1:Leaf_no(t)));
    Total_leaf_mass(t) = sum(Leaf_mass(t,:));
    Total_shoot(t) = Total_leaf_mass(t) + Hypocotyl_mass;
    Total_plant(t) = Total_shoot(t) + Root_mass(t);
    if t == input_size && data_set == 3
        disp(Total_plant(t))
    end
    
    Shoot_FW_gFW(t) = Total_shoot(t)/0.08;
    ShootFreshWeight = add_to_list(ShootFreshWeight, Shoot_FW_gFW(t));

end

%-------------------------------End-of-model-------------------------------

w = 0.92; %water content
d = 1-w; %dry matter
FW = Total_shoot/d; %fresh weight in g
Shoot_fresh_weight_at_flowering = FW(end)

