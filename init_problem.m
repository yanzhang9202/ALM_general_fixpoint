global P1 P2
%% Add all functions in this directory and subdir to searching path
addpath(genpath('/Users/Yan_Zhang/Documents/MATLAB/Optimization/ALM_general_fixpoint'), '-end'); 
%% Code starts here
rng(1); % Set random generator seed
switch prob_type
    case 'waterfilling'
        % Example 5.2 in "Convex Optimization" by S.Boyd
        data.N = 10; % # of channels, or dimension of decision variables
        data.m = 1;  % Number of equality constraint
        data.a = sort(rand(data.N,1)+1);  % parameter for each channel
        data.pw = 1; % Total power, RHS of the equality constraint        
        data.lb = zeros(data.N,1);  % All decision variables >= 0
        data.ub = data.pw*ones(data.N,1);
        data.b = data.pw;
        data.A = ones(1,data.N);
        fl = 5;
        [data, PreFPparam] = pre_fixpointize(data,fl);
        % Solve the problem using the closed form solution in the textbook
        data.sol = calc_closeform(data);
        
    case 'mpc'
        % Scaling
%         P1 = 1;    % Scale on objective
        P1 = 2^(-4);
        P2 = 2^4;   % Scale on Aeq and beq

        main_pf_mpc;
        % Solve the problem using Matlab QP function
        data.sol = calc_closeform(data);
        data = calc_scale(data, PreFPparam);
        
    otherwise
        error('Error: Undefined problem type and data!')
end