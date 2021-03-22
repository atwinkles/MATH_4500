function y = naturalspline(t,fv,x,q)

%This is an algorithm by Alexander Winkles that generates a natural cubic
%spline from given data.
%
%t : the knots
%y : the function values
%x : the point being evaluated
%q : if q == 1, plots the cubic splines for the given domain

l = length(t);
k = length(fv);

if l ~= k
    fprintf('The size of the knots does not match the size of the results')
end;

h = zeros(1,l-1);
for i=1:l-1
    h(1,i) = t(1,i+1) - t(1,i);
end;

A = zeros(l-2,l-2);
for i=1:l-3
    A(i+1,i) = h(1,i+1);
    A(i,i+1) = h(1,i+1);
end;

for i=2:l-1
    A(i-1,i-1) = 2.*(h(1,i)+h(1,i-1));
end;

v = zeros(1,l-2);
for i=2:k-1
    v(1,i-1) = 6./h(1,i)*(fv(1,i+1)-fv(1,i)) - 6./h(1,i-1)*(fv(1,i)-fv(i-1));
end;

z = v/A;

disp(l)
disp(length(h))
disp(length(z))

for i=1:l-1
    if (x >= t(1,i)) && (x <= t(1,i+1))
        if (i == 1)
            S = z(1,i+1)./(6*h(1,i)).*(x - t(1,i)).^3 + (fv(1,i+1)./h(1,i)...
                - z(1,i+1)*h(1,i)./6).*(x - t(1,i)) + (fv(1,i)./h(1,i)).*(t(1,i+1) - x);
            break
        elseif (i == l-1)
            S = z(1,i-1)./(6.*h(1,i)).*(t(1,i+1)-x).^3 + (fv(1,i+1)./h(1,i))...
                .*(x - t(1,i)) + (fv(1,i)./h(1,i) - z(1,i-1).*h(1,i)./6).*(t(1,i+1)-x);
            break
        else
            S = z(1,i-1)./(6*h(1,i)).*(t(1,i+1)-x).^3 + z(1,i)./(6.*h(1,i)).*(x - t(1,i)).^3 ...
                + (fv(1,i+1)./h(1,i) - z(1,i).*h(1,i)./6).*(x - t(1,i))...
                + (fv(1,i)./h(1,i) - z(1,i-1).*h(1,i)./6).*(t(i+1)-x);
            break
        end;
    end;
end;

y = S;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% attempting to make a better function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fone =@(u,i) z(1,i+1)./(6*h(1,i)).*(u - t(1,i)).^3 + (fv(1,i+1)./h(1,i)...
                - z(1,i+1)*h(1,i)./6).*(u - t(1,i)) + (fv(1,i)./h(1,i)).*(t(1,i+1) - u);
            
ftwo =@(u,i) z(1,i-1)./(6.*h(1,i)).*(t(1,i+1)-u).^3 + (fv(1,i+1)./h(1,i))...
                .*(u - t(1,i)) + (fv(1,i)./h(1,i) - z(1,i-1).*h(1,i)./6).*(t(1,i+1)-u);

fthree =@(u,i) z(1,i-1)./(6*h(1,i)).*(t(1,i+1)-u).^3 + z(1,i)./(6.*h(1,i)).*(u - t(1,i)).^3 ...
                + (fv(1,i+1)./h(1,i) - z(1,i).*h(1,i)./6).*(u- t(1,i))...
                + (fv(1,i)./h(1,i) - z(1,i-1).*h(1,i)./6).*(t(i+1)-u);

if q == 1
    r = min(t) : 0.001 : max(t);
    s = zeros(1,length(r));
    for i = 1:length(r)
        if (r(1,i) <= t(1,2))
            s(1,i) = fone(r(1,i),1);
        elseif (r(1,i) >= t(1,length(t))-1)
            s(1,i) = ftwo(r(1,i),length(t)-2);
        else
            for j=1:l-1
                if (r(1,i) >= t(1,j)) && (r(1,i) <= t(1,j+1))
                    s(1,i) = fthree(r(1,i),j);
                else
                    continue
                end;
            end;
        end;
    end;
    
    figure('Name', 'Spline interpolation & data points')
    disp(r)
    plot(r,s)
    hold on
    plot(t,fv,'+')
    hold off
end;
            