function y = csimpson(f,a,b,n)
%y = csimpson(f,a,b,n)
%This is an algorithm by Alexander Winkles used to perform composite
%Simpson's rule to compute numerical integral values.
%
%f : the function being integrated
%a : the lower bound of the integral
%b : the upper bound of the integral
%n : the number of subintervals used

h=(b-a)/n;

x = zeros(1,n+1);

x(1) = a;
x(n+1) = b;

p = 0;
q = 0;

for i=2:n
    x(i) = a + h*(i-1);
end;

for i=2:n
    p = p + f(x(i));
end;

for i = 2:n+1
    q = q + f((x(i)+x(i-1))/2);
end;

y = (h/6)*(f(a) + 2*p + 4*q + f(b));