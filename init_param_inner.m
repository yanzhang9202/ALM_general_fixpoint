switch ALMparam.inner_solver
    case 'GPM'
        GPMparam.stopcriter = 3; % 1 - stop by criteria, 2 - stop by iter number, 3 - stop by embedded criteria
        GPMparam.iter_max = 1e3;
        GPMparam.epsilon = 1e-6;
        GPMparam.sol_type = 1; % 1 - direct iterate, 2 - ergodic average
        switch prob_type
            case 'waterfilling'
                Hess = ALMparam.rho*ones(data.N) + diag(1./(double(data.a).^2));   % Hessian of the objective function
                if cond(Hess) > 1e3
                    fprintf('Warning: The inner problem is ill conditioned!\n')
                end
                GPMparam.L = max(eig(Hess));
                GPMparam.alpha = 1/GPMparam.L;
                
                rho = ALMparam.rho;
                Aeq = double(data.A);
                beq = double(data.b);
                L = ALMparam.L;
                
                fl = PreFPparam.fl;
                wl = PreFPparam.wl;
                init_fixpoint_param;                
                M1 = Aeq/L;
                M2 = beq/L;
                M3 = GPMparam.alpha;
                M4 = GPMparam.alpha*rho;
                
                data.M1 = fi(M1, T, F);
                data.M2 = fi(M2, T, F);
                data.M3 = fi(M3, T, F);
                data.M4 = fi(M4, T, F);
                
                clear Hess M1 M2 M3 M4 Aeq beq rho fl wl L 
                
            case 'mpc'
                H = double(data.H);
                Aeq = double(data.Aeq);
                beq = double(data.beq);
                rho = ALMparam.rho;
                L = ALMparam.L;
                GPMparam.L = norm(H + rho*Aeq'*Aeq);
                GPMparam.alpha = 1/GPMparam.L;
                
                M1 = Aeq/L;
                M2 = beq/L;
                M3 = (H + rho*Aeq'*Aeq)/GPMparam.L;
                M4 = Aeq'/GPMparam.L;
                M5 = rho*Aeq'*beq/GPMparam.L;
                fl = PreFPparam.fl;
                wl = PreFPparam.wl;
                init_fixpoint_param;
                data.M1 = fi(M1, T, F);
                data.M2 = fi(M2, T, F);
                data.M3 = fi(M3, T, F);
                data.M4 = fi(M4, T, F);
                data.M5 = fi(M5, T, F);
                
                clear H Aeq beq rho fl wl M1 M2 M3 M4 M5 L
                
            otherwise
                error('Error:Undefined GPM step size for the problem type!')
        end
end