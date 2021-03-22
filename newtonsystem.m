function newtonsystem(F,J,x0,y0,T,S,q)
%newtonsystem(F,J,x0,y1,T,S,q)
%
%This is a Newton's method algorithm developed by Alexander Winkles used to
% find roots in two-dimension systems of nonlinear equations. 
%
%F  : the function matrix being studied.
%J  : the Jacobian of the function matrix being studied.
%x0 : the initial vector being used.
%T  : the tolerance vector of the final answer.
%S  : the number of iterations.
%q  : prints the steps (1 for yes, 0 for no).

format long g

if q == 1
    fprintf('\nStep \t\t x \t\t\t y\n---\t\t -------- \t\t --------\n');
end;
if abs(feval(F,x0,y0)) < T
    fprintf('\nYour inital guess (%d %d) was close enough to zero!',x0,y0);
else
    i = 0;
    while i <= S
        H = -inv(J(x0,y0))*F(x0,y0);
        x1 = x0 + H(1);
        y1 = y0 + H(2);
        if q == 1
            fprintf('%d \t\t %f \t\t %f\n',i,x1,y1);
        end;
        if abs(feval(F,x1,y1)) < T
            fprintf('\nThe solution is (%f, %f). The computation was a success after %d iterations!\n\n',x1,y1,i);
            break;
        end;
        i = i + 1;
        x0 = x1;
        y0 = y1;
    end;
    if i == S+1
        fprintf('\nMethod failed after %d iterations.\n\n',S)
    end;
end;