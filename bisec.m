function bisec(f,a,b,T,S,q)
% bisect(f,a,b,D,E)
%This is a bisection algorithm by Alexander Winkles used to find roots of 
% polynomials over a domain.
%
% f : the function being evaluated
% a : the lower domain value
% b : the upper domain value
% T : the tolerance of the final result
% S : the number of iterations
% q : prints the iterations (1 is yes, 0 is no)


if sign(feval(f,a))==sign(feval(f,b))
    disp('\nError: f(a) and f(b) have the same signs.\n\n');
else
    if q == 1
        fprintf('\nStep \t\t Result\n---\t\t --------\n');
    end;
    i = 0;
    while i <= S
        e = b - a;
        e = e/2;
        c = a + e;
        if q  == 1
            fprintf('%d \t\t %f\n',i,c)
        end;
        if feval(f,c) == 0 || abs(feval(f,c)) < T
            fprintf('\nThe solution is %d. The computation was a success after %d iterations!\n\n',c,i)
            break;
        end;
        i = i+1;
        if sign(feval(f,a))==sign(feval(f,c))
            a = c;
        else
            b = c;
        end;
    end;
    if i == S+1
        fprintf('\nMethod failed after %d iterations.\n\n',S)
    end;
end;
    