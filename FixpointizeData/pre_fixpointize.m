function [ data, FPparam ] = pre_fixpointize(data, fl)
global prob_type P1
% Transform the generated problem data to fix point data
switch prob_type
    case 'waterfilling'
        % Decide the word length
        a = data.a;
        pw = data.pw;
        lb = data.lb;
        ub = data.ub;
        nm = 2^(-fl); % The minimum number can be represented;
        nM = max([max(a), pw, max(abs(lb)), max(abs(ub))]);
        wl = ceil(log(nM/nm)/log(2)) + 2;   % +1 to has the sign bit
        % Transform data
        init_fixpoint_param;
        data.a = fi(a, T, F);
        data.pw = fi(pw, T, F);
        data.lb = fi(lb, T, F);
        data.ub = fi(ub, T, F);
        data.b = fi(data.b, T, F);
        data.A = fi(data.A, T, F);
        % Output Fixpoint parameter
        FPparam.fl = fl;    FPparam.wl = wl;
        
    case 'mpc'
        H = data.H;
        f = data.f;
        Aeq = data.Aeq;
        beq = data.beq;
        ub = data.ub;
        lb = data.lb;
        % Decide wordlength and fraction length to represent the problem under fixed point arithmetic
        max_entry = zeros(6,1);
        min_entry = zeros(6,1);
        max_entry(1) = max(max(abs(H(H ~= 0))));
        min_entry(1) = min(min(abs(H(H ~= 0))));
        % max_entry(2) = max(abs(f(f ~= 0)));
        % min_entry(2) = min(abs(f(f ~= 0)));
        max_entry(3) = max(max(abs(Aeq(Aeq ~= 0))));
        min_entry(3) = min(min(abs(Aeq(Aeq ~= 0))));
        max_entry(4) = max(abs(beq(beq ~= 0)));
        min_entry(4) = min(abs(beq(beq ~= 0)));
        max_entry(5) = max(abs(ub(ub ~= 0)));
        min_entry(5) = min(abs(ub(ub ~= 0)));
        max_entry(6) = max(abs(lb(lb ~= 0)));
        min_entry(6) = min(abs(lb(lb ~= 0)));
        max_entry(2) = [];
        min_entry(2) = [];
        fl_min = ceil(log(1/min(min_entry))/log(2));    % Minimum fraction bit needed to represent the problem most exactly
        fl = fl_min + ceil(-log(P1)/log(2));
        wl = fl + 1 + ceil(log(max(max_entry))/log(2)); 
        % Transform data
        init_fixpoint_param;
        data.H = fi(H, T, F);
        data.f = fi(f, T, F);
        data.beq = fi(beq, T, F);
        data.Aeq = fi(Aeq, T, F);
        data.lb = fi(lb, T, F);
        data.ub = fi(ub, T, F);        
        % Output Fixpoint parameter
        FPparam.fl = fl;    FPparam.wl = wl;  
        FPparam.T = T; FPparam.F = F;
end
end