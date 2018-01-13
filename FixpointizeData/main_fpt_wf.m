rho = ALMparam.rho;
L = 2/rho;
C1 = L/2*max([16*bd_lambda^2, 4*(bd_lambda + 1)^2]) + ...
    1/2*max([9*bd_lambda^2, (2*bd_lambda + 1)^2]);

% Kout
Ko = ceil(C1/(1-alpha)/epsilon);

% U1 B_out B_in
B_lambda = 2*bd_lambda;
U1 = (sqrt((1+4/L)^2*(B_lambda+gamma)^2+2*(1+1/L)*alpha*epsilon)...
    -(1+4/L)*((B_lambda+gamma)))/(1+1/L);
B_out = U1;
B_in = gamma*U1;

% K_in U2
a = double(data.a);
n = data.N;
%         Lp = 1/min(a)^2+n;
Lp = GPMparam.L;
Bx = sqrt(n);
Ki = ceil(Lp*Bx^2/2/(1-beta)/B_in - 1);
U2 = beta*B_in/2/Bx;  

% bd_1_fl
m = 1;
x_infty = double(data.pw);
C2 = L*sqrt(m)*((1+x_infty)*n+1);
bd_1_fl = ceil(-log(U1/C2)/log(2)-1);

% bd_2_fl
C3 = Lp*sqrt(n)*(4+max(1./a)+bd_lambda+n*x_infty);
bd_2_fl = ceil(-log(U2/C3)/log(2) - 1);

% fl
fl = max([bd_1_fl, bd_2_fl]);

% U3 wl
A = double(data.A);
b = double(data.b);
U3 = max([bd_lambda+A/L*x_infty-b/L, ...
    b+1/Lp*max(1./a)+1/Lp*bd_lambda+rho/Lp, ...
    1/Lp*bd_lambda+rho/Lp*n*b-b-1/Lp*max(1./a)-rho/Lp, ...
    max(a)+b]); % Should use matrix inf norm
wl = ceil(log(U3)/log(2) + fl + 2);

C(4) = 2*bd_lambda^2/rho+1/2*(1/rho+1)*9*bd_lambda^2;
C(5) = 8*bd_lambda^2/rho+1/2*(1/rho+1)*9*bd_lambda^2;
C(6) = 1/2*1/rho*(3*bd_lambda+1)^2+1/2*(1/rho+1)*9*bd_lambda^2; 

% Test fixed point parameters:
%         E = (1+4/L)*B_lambda*B_out+(1+4/L)*B_in+(1/2+1/2/L)*B_out^2;
%         fprintf(['K_out: ', num2str(Ko), ', K_in: ', num2str(Ki),...
%             ', fl: ', num2str(fl), ', wl: ', num2str(wl), ...
%             ', B_in: ', num2str(B_in), ', B_out: ', num2str(B_out), ...
%             ', E: ', num2str(E)]);
%         Ko*Ki*wl
Ki = 400;

%% Apply the fix point setting and iteration numbers
init_fixpoint_param;
GPMparam.alpha = fi(GPMparam.alpha, T, F);
GPMparam.iter_max = Ki;
ALMparam.rho = fi(ALMparam.rho, T, F);
ALMparam.iter_max = Ko;
data.a = fi(data.a, T, F);
data.A = fi(data.A, T, F);
data.b = fi(data.b, T, F);

FPparam.fl = fl;    FPparam.wl = wl;
FPparam.T = T;  FPparam.F = F;
FPparam.E = alpha*epsilon;

clear rho a bd_lambda Lp Bx C1 C2 C3 U1 U2 U3 B_lambda B_out ...
    B_in n A b bd_1_fl bd_2_fl epsilon