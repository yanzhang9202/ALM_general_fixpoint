function [ u ] = calc_proj(v, lb, ub)
global verbose
u = v;
if ~isempty(lb)
    u(v < lb) = lb(v < lb);
end
if ~isempty(ub)
    u(v > ub) = ub(v > ub);
end
if any(u ~= v) && verbose
   fprintf('Warning: Projection on multiplier is effective!\n') 
end
end