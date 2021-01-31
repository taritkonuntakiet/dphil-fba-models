%Function to convert gC/plant/hour to umol/m2/s

function [new_value] = convert_to_micromol(value, carbons, surface_area)

%Convert per plant to per m2

new_value = value/surface_area;

%Convert per hour to per second

new_value = new_value/3600;

%Convert gC to mol

total_mass = carbons*12;
new_value = new_value/total_mass;

%Convert mol to umol

new_value = new_value*(10^(6));

end
