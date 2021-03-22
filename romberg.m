function y = romberg(f,a,b,n,q1,q2)
%y = romberg(f,a,b,n,r,q)
%
%This is an algorithm written by Alexander Winkles that performs Romberg
%integration.
%
%f  : the function being integrated
%a  : the lower bound of the integral
%b  : the upper bound of the integral
%n  : the 2^n subinterval specification
%q1 : if q == 1, returns an array of all Romberg values computed
%q2 : a two element vector. if q2(1) == 1, then prints the error matrix
%       where q2(2) is the true value

h = b-a;

R = zeros(n+1,n+1);

R(1,1) = (1/2)*h*(feval(f,a) + feval(f,b));

for i = 2 : n+1
    h = h/2;
    sum = 0;
    for u = 1 : (2^(i-2))
        sum = sum + feval(f,(a + (2*u-1)*h));
    end;
    R(i,1) = (1/2)*R(i-1,1) + h*(sum);
    for j = 2 : i
        R(i,j) = R(i,j-1) + (R(i,j-1) - R(i-1,j-1))/(4^(j-1) - 1);
    end;
end;


if q1 == 1
    disp(R)
end;

if q2(1) == 1
    E = zeros(n+1,n+1);
    for i=1:n+1
        for j=1:n+1
            if i >= j
                E(i,j) = q2(2) - R(i,j);
            end;
        end;
    end;
    disp(E)
end;

y = R(n+1,n+1);