function [ data ] = calc_scale(data, FP)
global prob_type P1
switch prob_type
    case 'waterfilling'
        
    case 'mpc'
        H = double(data.H);
        f = double(data.f);
        data.H = fi(H*P1, FP.T, FP.F);
        data.f = fi(f*P1, FP.T, FP.F);
    otherwise
        error('Error: In scaling data!')
end