function f = objfun(x,last_x)
eta1 = 1/((1-last_x(1)*last_x(2))^2);
eta2 = (last_x(2))^2 /2;

f = (x(3)+x(2))/2 + 1/(x(1)) + eta1 * ( eta2 * x(1) + 0.5 * x(2));
end