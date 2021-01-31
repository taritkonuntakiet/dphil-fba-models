function [FBAsoln] = pFBA(model,max_or_min)
%This function performs simple pFBA.
%!Warning only flux values are updated!

FBAsoln = optimizeCbModel(model,max_or_min);
if FBAsoln.f == 0
    disp(FBAsoln.f)
    error('INFEASIBLE 1!')
end

indx_obj = find(model.c == 1);
[model.lb(indx_obj), model.ub(indx_obj)] = deal(FBAsoln.f);

[irr_model,matchRev,rev2irrev,irrev2rev] = convertToIrreversible(model);

irr_model.c = irr_model.c * 0;
irr_model.c = irr_model.c + 1;

FBAsoln2 = optimizeCbModel(irr_model,'min');

if FBAsoln2.f == 0
    disp(FBAsoln2.f)
    error('INFEASIBLE 2!')
end

if length(rev2irrev) ~= length(model.rxns)
    error('PROBLEM!!')
end

for i = 1:length(rev2irrev)
    temp = rev2irrev(i);
    if length(temp{1}) == 2
        
        FBAsoln.x(i) = FBAsoln2.x(temp{1}(1)) - FBAsoln2.x(temp{1}(2));
    else
        FBAsoln.x(i) = FBAsoln2.x(temp{1}(1));
    end
end