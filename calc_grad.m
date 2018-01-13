function [g] = calc_grad(x, data, lambda)
global prob_type ALMparam
switch prob_type
    case 'waterfilling'
        g = -1./(data.a + x) + lambda + ALMparam.rho*(sum(x)-1);
        
    case 'mpc'
        
        
    otherwise
        error('Error:Undefined gradient expression for the problem type!')
end
end