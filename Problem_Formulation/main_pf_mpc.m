% Constrained-MPC parameter
N = 5; % Horizon
nx = 2; % State dimension
nu = 1; % Control dimension
Q = [10, 0; 0, 1];    % Energy cost
R = 1;
Qn = Q;
Ad = [1, 0.01; 0, 1];    % State dynamics
Bd = [-0.0004; -0.0701];
lim_x = [-0.2, 0.01; -0.1, 0.1];    % Box Constraints
lim_u = [-0.0524, 0.0524];
x0 = [-0.15; 0.05];  % Starting state
% x0 = [0.009; 0.09];

% Build QP problem out of MPC
dim_z = (N+1)*nx + N*nu;    % Dimension of the stacked variable
H = zeros(dim_z);
H(1:(N+1)*nx, 1:(N+1)*nx) = kron(eye(N+1), Q);
H(((N+1)*nx+1):end,((N+1)*nx+1):end) = kron(eye(N), R);
f = zeros(dim_z, 1);
Aeq = [zeros(nx, dim_z); kron(-eye(N), Ad), zeros(N*nx, nx), kron(-eye(N), Bd)];
Aeq = Aeq + [eye((N+1)*nx), zeros((N+1)*nx, N*nu)];
beq = [x0; zeros(N*nx, 1)];
ub = [kron(ones(N+1, 1), lim_x(:,2)); kron(ones(N,1), lim_u(2))];
lb = [kron(ones(N+1, 1), lim_x(:,1)); kron(ones(N,1), lim_u(1))];

data.N = dim_z;
data.m = size(Aeq, 1);
data.H = H;
data.f = f;
data.Aeq = Aeq * P2;
data.beq = beq * P2;
data.lb = lb;
data.ub = ub;

clear N nx nu Q R Qn Ad Bd lim_x lim_u x0 dim_z H f Aeq beq lb ub

[data, PreFPparam] = pre_fixpointize(data, []);  % Decide fl inside the function
