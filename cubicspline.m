function res = cubicspline(t,y,x)

%This is an algorithm for testing a more flexible spline code. Currently
%only for natural cubic splines.

l = length(t);
k = length(y);

if l ~= k
    fprintf('The size of the knots does not match the size of the results')
end;

h = zeros(1, l-1);
for i=1:l-1
    h(1,i) = t(1,i+1) - t(1,i);
end;

A = zeros(

