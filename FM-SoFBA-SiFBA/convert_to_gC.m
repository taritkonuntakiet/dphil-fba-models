%Function to convert umol/.../... to gC/.../...

function [new_value] = convert_to_gC(value, carbons)

%Convert umol to ug

new_value = value*carbons*12;

%Convert ug to g

new_value = new_value/(10^(6));

end