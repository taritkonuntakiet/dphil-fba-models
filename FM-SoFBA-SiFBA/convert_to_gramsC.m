%Function to convert umol/m2/s to gC/plant/hour

function [new_value] = convert_to_gramsC(value, carbons, surface_area)

%Convert per m2 to per plant

new_value = value*surface_area;

%Convert per second to per hour

new_value = new_value*3600;

%Convert umol to ug

new_value = new_value*carbons*12;

%Convert ug to g

new_value = new_value/(10^(6));

end