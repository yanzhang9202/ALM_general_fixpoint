function [flag] = check_stopcriter(g, ind_lb, ind_ub, delta1, delta2)
% Check optimality condition;
if any(g(ind_lb) <= -delta1)
    flag = 0;
    return
end

if any(g(ind_ub) >= delta1)
    flag = 0;
    return
end

ind_act = [ind_lb(g(ind_lb > delta1)); ind_ub(g(ind_ub) < -delta1)];
g(ind_act) = [];
for i = 1 : length(g)
   if abs(g(i)) > delta2
       flag = 0;
       return
   end 
end

flag = 1;

end
