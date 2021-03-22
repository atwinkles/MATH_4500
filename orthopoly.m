function P = orthopoly(w,a,b,n,q)
%P = orthopoly(f,g,w,a,b,n)
%
%This is an algorithm designed by Alexander Winkles that producing
%orthogonal polynomials of degree n+1 based on a weight function w(x) on
%the set [a,b] used for the sake of integration.
%f : a function of f
%g : a function of f
%w : a weight function
%a : lower bound of intergral
%b : the upper bound of the integral
%n : the degree of the polynomial generated
%q : if q == 1 then gives roots for polynomial

syms x

p(1) = 1;
p(2) = x - k(1);

k(1) = (integral(w*x,a,b))/(integral(w,a,b));

