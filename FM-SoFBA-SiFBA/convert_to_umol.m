%Function to convert umol/.../... to gC/.../...

function [new_value] = convert_to_umol(value, carbons)

%Convert gC to mol

new_value = value/carbons/12;

%Convert mol to umol

new_value = new_value*(10^(6));

end