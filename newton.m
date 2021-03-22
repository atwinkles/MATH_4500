function newton(f,df,x0,T,S,q)
%newton(f,df,x0,T,S,q)
%
%This is a Newton's method algorithm developed by Alexander Winkles used to
% find roots in functions.
%
%f  : the function being studied.
%df : the derivative of the function being studied.
%x0 : the initial point being used.
%T  : the tolerance of the final answer.
%S  : the number of iterations.
%q  : prints the steps (1 for yes, 0 for no).

format long
if q == 1
    fprintf('\nStep \t\t Result\n---\t\t --------\n');
end;
if feval(f,x0) == 0 || abs(feval(f,x0)) < T
    fprintf('\n Your inital guess %f was close enough to zero!',x0);
else
    i = 0;
    while i <= S
        x1 = x0 - feval(f,x0)/feval(df,x0);
        if q == 1
            fprintf('%d \t\t %f\n',i,x1);
        end;
        if feval(f,x1) == 0 || abs(feval(f,x1)) < T
            fprintf('\nThe solution is %f. The computation was a success after %d iterations!\n\n',x1,i);
            break;
        end;
        i = i+1;
        x0 = x1;
    end;
    if i == S+1
        fprintf('\nMethod failed after %d iterations.\n\n',S)
    end;
end;
