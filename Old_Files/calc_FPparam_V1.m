switch prob_type
    case 'waterfilling'
    rho = ALMparam.rho;
    N = data.N;
    a = double(data.a);
    bd_lambda = 1;  % Tested by the code in the folder ~/ALM_general
    Lp = 1/min(a)^2+N;
    Bx = sqrt(N);        
        
    %% Compute the constants in the error expression
    C = zeros(7,1);
    C(1) = 4*rho+3;
    C(2) = 1/2+1/2*rho;
    C(3) = 2*(rho+1)*(1+2*sqrt(Lp)*Bx);
    C(4) = 2*bd_lambda^2/rho+1/2*(1/rho+1)*9*bd_lambda^2;
    C(5) = 8*bd_lambda^2/rho+1/2*(1/rho+1)*9*bd_lambda^2;
    C(6) = 1/2*1/rho*(3*bd_lambda+1)^2+1/2*(1/rho+1)*9*bd_lambda^2;
    C(7) = max([C(4),C(5),C(6)]);
    C(8) = Lp*Bx^2/2;

    %% Compute Ko Ki fl
    e = ALMparam.epsilon;
    Ko = ceil(C(7)/alpha/e);
    fl = ceil(-log((1-alpha)^2*e^2/N/C(3)^2/(1+1/GPMparam.alpha))/log(2));
    Bout = 0;
    delta_gp = 2*(1+1/GPMparam.alpha)*sqrt(N)*2^(-fl-1)*sqrt(N);
    if C(1)*Bout*4*bd_lambda^2+C(2)*Bout^2+C(3)*sqrt(delta_gp) > (1-alpha)*e
        error('Warning: Increase fl so that wanted steady error E can be achieved!')
    end
    Ki = C(8)/((1-alpha)*e - C(1)*Bout*4*bd_lambda^2+C(2)*Bout^2+C(3)*sqrt(delta_gp)) - 1;
    Ki = ceil(Ki);
    %% Compute wl
    gi = 1/min(a)+2*bd_lambda+rho*(N-1);
    rg_x = 1 + GPMparam.alpha*gi;
    rg_lambda = max(2*bd_lambda, bd_lambda+1) + rho*(N-1);
    wl = ceil(log(max([rg_x,rg_lambda])/(2^-fl))/log(2))+1;
    % test the product KoKiwl to optimize the computing time
%     Ko*Ki*wl
    
    %% Apply the fix point setting and iteration numbers
    init_fixpoint_param;
    GPMparam.alpha = fi(GPMparam.alpha, T, F);
    GPMparam.iter_max = Ki;
    ALMparam.rho = fi(ALMparam.rho, T, F);
    ALMparam.iter_max = Ko;
    data.a = fi(data.a, T, F);
    data.A = fi(data.A, T, F);
    data.b = fi(data.b, F, F);
    
    FPparam.fl = fl;    FPparam.wl = wl;
    FPparam.T = T;  FPparam.F = F;
    FPparam.E = (1-alpha)*e;
    
    clear rho N a bd_lambda Lp Bx e Bout delta_gp gi rg_x rg_lambda fl wl T F  
end