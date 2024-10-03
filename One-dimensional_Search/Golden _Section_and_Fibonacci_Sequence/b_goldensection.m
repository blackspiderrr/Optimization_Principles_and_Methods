f = @(x)8*exp(1-x)+7*log(x);
left = 1;
right = 2;
tol = 0.23;
rho = (3 - sqrt(5)) / 2;
golden_ratio = 1 - rho;
N = log(tol/ (right - left))/ log(golden_ratio);
a = left + rho * (right - left);
b = left + (1 - rho)*(right - left);
for i = 1: ceil(N)
    fprintf(['Iteration', num2str(i)]);
    f_a = f(a);
    f_b = f(b);
    fprintf(['\na=',num2str(a),' b=',num2str(b),' f(a)=',num2str(f_a),' f(b)=',num2str(f_b)]);
    if f_a < f_b 
        c = b;
        b = a;
        a = left + rho * (c-left);
        right = right - rho * (right - left);
    else 
        c = a;
        a = b;
        b = c + (1 - rho)*(right - c);
        left = left + rho * (right - left);
    end
    new_interval = [left, right];
    fprintf(['\nnew_interval: ', num2str(new_interval),'\n']);
end