%% Use the closed form solution in Example 5.2 in "Convex
%% Optimization" by S.Boyd
a = obj.a; N = obj.N; pw = obj.pw;
% Change fixpoint data back to double floating point data
a = double(a);
pw = double(pw);
% Assume N > 2;
flag = 0; % The optimal multiplier exists outside vector a
for i = 2 : N
   v = a(i) - a; 
   v = v(v>=0);
   if sum(v) > pw
       flag = 1; % The optimal multiplier exists within vector a
       break;
   end
end
% Initialize bisectioning
bl = 0; bh = 0; bm = 0;
if flag
    bl = a(i-1); bh = a(i); bm = (bl+bh)/2;
else
    bl = a(N); bh = bl+pw; bm = bl+pw/2;
end
% Start bisectioning
cnt = 0;
while((bh - bl) > 1e-15)
    v = bm - a; 
    v = v(v>=0);
    if sum(v) > pw
        bh = bm;
        bm = (bl+bh)/2;
    else if sum(v) < pw
            bl = bm;
            bm = (bl+bh)/2;        
        else
            break;
        end
    end
    cnt = cnt + 1;
end
nu = bm;    % This is the optimal multiplier
v = nu - a;
sol = zeros(N,1);
sol(v>=0) = v(v>=0);


