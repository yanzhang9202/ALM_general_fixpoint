switch prob_type
    case 'waterfilling'
        
    case 'mpc'
        n = data.N;
        H = double(data.H);
        A = double(data.Aeq);
        rho = double(ALMparam.rho);
        Lp = GPMparam.L;
        T = FPparam.T;
        F = FPparam.F;
        
        delta1 = U2/Lp/sqrt(n);
        GPMparam.delta1 = fi(delta1, T, F);
        sigma = min(svd(H + rho*A'*A));
        delta2 = sqrt(sigma*B_in/2/n/Lp) - delta1;
%         delta2 = delta2/2;
        if delta2 <= 0
            GPMparam.delta2 = 0;
            fprintf('Warning: No stopping criteria applied!\n')
        else
            GPMparam.delta2 = fi(delta2, T, F);
        end
        
        if delta2 < delta1
            fprintf('Warning: delta2 < delta1, stopping criteria is never satisfied!')
        end
        
        clear n H A rho Lp T F delta1 sigma delta2
end