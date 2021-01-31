%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function estimates Rubisco carboxylase flux %
% at which the net CO2 uptake rate is equal to the %
% user defined value                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [estimatedVc] = estimateVcFromNetRate(dielFBA_model, net_rate)

CO2_index = find(dielFBA_model.rxns=="CO2_tx_day");
Vc_index = find(dielFBA_model.rxns=="RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day");

%Initally constraint Vc flux to net CO2 uptake rate
dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", net_rate, 'u');
dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", net_rate, 'l');

%Perform pFBA
sol = pFBA(dielFBA_model, 'max');

%Set loop counter
i = 0;

%Use a while loop to increase Vc until net CO2 rate is similar to given value (or loop counter hits 10)
while (net_rate - sol.x(CO2_index))/net_rate > 0.0001 && i<10
    i = i+1;
    prev = sol.x(Vc_index);
    %Increment in Vc flux is set by given netCo2 uptake - model predicted CO2 uptake rate in previous pFBA run
    now = prev + (net_rate - sol.x(CO2_index));
    dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", now, 'u');
    dielFBA_model = changeRxnBounds(dielFBA_model, "RIBULOSE_BISPHOSPHATE_CARBOXYLASE_RXN_p_day", now, 'l');
    sol = pFBA(dielFBA_model, 'max');
end

estimatedVc = prev;

end


