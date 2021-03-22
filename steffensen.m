function steffensen(f,x0,T,S,q)
%steffensen(f,x0,T,S,q)
%
%This is a Steffensen's method algorithm designed by Alexander Winkles to
% find roots of functions.
%
%f  : the function being evaluated.
%x0 : the initial value being used.
%T  : the tolerance of the final answer.
%S  : the number of iterations.
%q  : prints the steps (1 for yes, 0 for no)
%
%*Current issues: does not necessarily converge to the closest zero from the
% starting point.

format long
if q == 1
    fprintf('\nStep \t\t Result\n---\t\t --------\n');
end;
if abs(feval(f,x0)) < T
    fprintf('\n Your inital guess %d was close enough to zero!',x0);
else
    i = 0;
    while i<= S
        x1 = x0 - feval(f,x0)/((feval(f,x0 + feval(f,x0))-feval(f,x0))/feval(f,x0));
        if q == 1
            fprintf('%d \t\t %f\n',i,x1)
        end;
        if abs(feval(f,x1)) < T
            fprintf('\nThe solution is %f. The computation was successful after %d iterations!\n\n',x1,i)
            break;
        end;
        i = i+1;
        x0 = x1;
    end;
    if i == S+1
        fprintf('\nThe computation failed after %d iterations.\n\n',S)
    end;
end;
    