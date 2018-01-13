function [ sol ] = calc_closeform( obj )
global prob_type
switch prob_type
    case 'waterfilling'
        calc_cf_waterfilling;
        
    case 'mpc'
        calc_cf_mpc;
        
    otherwise
        error('Error: Undefined closed form solution for this problem type!')
end
end