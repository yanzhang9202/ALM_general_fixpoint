%% Set algorithm parameters for ALM
ALMparam.rho = 1;
ALMparam.iter_max = 1e3;    % May be changed to the iteration complexity result
ALMparam.inner_solver = 'GPM';
ALMparam.stopcriter = 2;   % 1 - stop by criteria, 2 - stop by iter number

% ALMparam.epsilon = 1e0;   % 1 - \|Axk - b\|, 2 - user specified accuracy
ALMparam.epsilon = 1e-3 * P1;    % target epsilon * scaling factor

ALMparam.exitflag = 0;  % 0 - fail, 1 - success
ALMparam.L = 2/ALMparam.rho;

%% set parameter for inner solver
init_param_inner;

%% Set error weights for fixed point ALM
alpha = 0.5;
beta = 0.5;
gamma = 100;
calc_FPparam_V2;

%% Decide stopping criteria for inner solver
calc_stopcondition;

%% Assign variable space
% var_fp.x = fi(zeros(data.N, ALMparam.iter_max), FPparam.T, FPparam.F);
var.x_avg = zeros(data.N, ALMparam.iter_max);
% var_fp.lambda = fi(zeros(data.m, ALMparam.iter_max), FPparam.T, FPparam.F);
var.lambda_avg = zeros(data.m, ALMparam.iter_max); % Currently the ergodic average is double precision
var.x = zeros(data.N, ALMparam.iter_max);
var.lambda = zeros(data.m, ALMparam.iter_max);

x_iter = fi(zeros(data.N,1), FPparam.T, FPparam.F);
lambda_prev = fi(zeros(data.m,1), FPparam.T, FPparam.F);
lambda_iter = lambda_prev;

err = zeros(1, ALMparam.iter_max); % Error of GPM
GPMflag = zeros(1, ALMparam.iter_max);