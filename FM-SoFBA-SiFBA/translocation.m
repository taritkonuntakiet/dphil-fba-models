function [leaf_trans,root_trans] = translocation(rosette_area,suc_c_interm,suc_equi,leaf_c,root_c)


suc_equi_plant = suc_equi*rosette_area; %Equilibrium sucrose plus hexose mass for whole plant (g C/plant)
root_and_leaf_c = root_c + leaf_c; %Total root and leaf C (g C/plant)

if      suc_c_interm <= suc_equi_plant %translocation needed
        
        scalingl = leaf_c/root_and_leaf_c;
        leaf_trans = (suc_equi_plant - suc_c_interm)*scalingl; %Translocation from leaves (gC/plant)
        scalingr = root_c/root_and_leaf_c;
        root_trans = (suc_equi_plant - suc_c_interm)*scalingr; %Translocation from root (g C/plant)
else
        leaf_trans = 0;
        root_trans = 0;
end        
    
