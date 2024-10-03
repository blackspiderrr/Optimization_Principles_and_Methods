syms x;
f = 8*exp(1-x) + 7*log(x);
f_prime = diff(f,x);
points = solve(f_prime == 0, x);
points = double(points);
points = points(points >=1 & points <= 2);
if length(points) <= 1
    disp('f在[1,2]上是单峰的');
else
    disp('f在[1,2]上不是单峰的');
end