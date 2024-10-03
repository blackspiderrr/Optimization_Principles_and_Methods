f = @(x)8*exp(1-x)+7*log(x);
left = 1;
right = 2;
tol = 0.23;
e=0.05;
F_min = (right-left)*(1+2*e)/tol;
F(1) = 1;
F(2) = 2;
N = 2;
while true
    F(N+1) = F(N) + F(N-1);
    if F(N+1) > F_min
        break;
    end
    N = N + 1;
end
rho = 1 - F(N)/F(N+1);
a = left + rho * (right - left);
b = left + (1 - rho)*(right - left);
for i = 1: N
    %此处与黄金分割法区分，涉及到rho的迭代
    if i == N
        rho = 1/2 - e;
    else
        rho = 1- F(N+1-i)/F(N+2-i);
    end
    if i < N-1
        rho_latter_one = 1- F(N-i)/F(N+1-i);
    else 
        rho_latter_one = 1/2 - e;
    end
    fprintf(['Iteration', num2str(i),' rho:', num2str(rho)]);
    f_a = f(a);
    f_b = f(b);
    fprintf(['\na=',num2str(a),' b=',num2str(b),' f(a)=',num2str(f_a),' f(b)=',num2str(f_b)]);
    if f_a < f_b 
        c = b;
        b = a;
        a = left + rho_latter_one * (c - left);
        right = right - rho * (right - left);
    else 
        c = a;
        a = b;
        b = c + (1 - rho_latter_one) * (right - c);
        left = left + rho * (right - left);
    end
    new_interval = [left, right];
    fprintf(['\nnew_interval: ', num2str(new_interval),'\n']);
end