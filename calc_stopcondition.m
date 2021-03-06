switch prob_type
    case 'waterfilling'
        n = data.N;
        a = double(data.a);
        b = double(data.b);
        Lp = GPMparam.L;
        T = FPparam.T;
        F = FPparam.F;
        
        delta1 = U2/Lp/sqrt(n);
        GPMparam.delta1 = fi(delta1, T, F);
        sigma = 1/(max(a)+b)^2;
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
        
        delta3 = sqrt(sigma*B_in/2/Lp) - delta1*sqrt(n);
        GPMparam.delta3 = fi(delta3, T, F);
        
        clear n a b Lp T F delta1 delta2 sigma B_in U2 delta3
        
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
        
        delta3 = 0; % Unused in the MPC problem
        GPMparam.delta3 = fi(delta3, T, F);
        
        clear n H A rho Lp T F delta1 sigma delta2 B_in U2 delta3
end