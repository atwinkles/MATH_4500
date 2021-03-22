function p = newtondivdiff(xx,yy,q1,q2,q3)
%p(x) = newtondivdiff(xx,yy,h)
%
%This is an algorithm developed by Alexander Winkles that takes data and
%utilizes divided differences to generate the Newton interpolating
%polynomial of the given data.
%
%The output p(x) may be used for evaluation the polynomial in various ways.
%
%xx : the nodes
%yy : the results
%q1 : prints divided difference table if q1 == 1
%q2 : prints the interpolated polynomial if q2 == 1
%q3 : prints figures to compare data to interpolating polynomial if q3 == 1

format rat

%Builds a matrix of the divided difference values, including f(x_i)
%values.

l = length(xx);

k = length(yy);

if l ~= k
    fprintf('The data points do not match in length!')
end;

t = zeros(l,l);

t(:,1) = yy';

for j = 2 : l
    for i = 1 : (l - j + 1)
        t(i,j) = (t(i+1,j-1) - t(i,j-1))./(xx(i+j-1) - xx(i));
    end;
end;

%Uses values from the first row of the divided difference table to
%generate the Newton interpolating polynomial.

q =@(x) 0;
v =@(x) 1;

for j = 1 : l
    if j == 1;
        q =@(x) t(1,j);
    else
        u = t(1,j);
        v =@(x) v(x).*(x - xx(j-1));
        q =@(x) q(x) + u.*v(x);
    end;
end;

%Generates the same polynomial as above for printing purposes. (find out
%how to remove

b = 0;
a = 1;

for j = 1 : l
    if j == 1;
        b = t(1,j);
    else
        syms x;
        u = t(1,j);
        a = a*(x - xx(j-1));
        b = b + u*a;
    end;
end;

%Print results.

if q1 == 1
    fprintf('\nThe divided difference table is as follows:\n\n')
    disp(t)
end;

if q2 == 1
    fprintf('\nThe Newton interpolating polynomial is:\n\n %s\n\nUse the following to further evaluate the found interpolating polynomial.\n',b);
end;

p =@(x) q(x);

%Plots points and the interpolating polynomial.

if q3 == 1
    figure('Name','Data points')
    scatter(xx,yy);

    figure('Name','Data points on polynomial')
    %sorted = sort(xx);
    %ss = sorted(2) : 0.01: sorted(length(sorted)-1);
    rr = min(xx) : 0.1: max(xx);
    %plot(ss,q(ss));
    plot(rr,q(rr));
    hold on
    plot(xx,yy,'+')
    hold off
    
    figure('Name','Interpolating polynomial over a large domain')
    domain = 0 : 1 : 500;
    range = q(domain);
    plot(domain,range)
    hold on
    plot(xx,yy,'+')
    hold off
end;