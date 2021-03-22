function secand(f,x0,x1,T,S,q)
%secand(f,x0,x1,T,S,q)
%
%This is a method developed by Alexander Winkles used to
% find roots in functions that is modelled after a question on my numerical
% anlysis exam.
%
%f  : the function being studied.
%x0 : the first initial point being used.
%x1 : the second initial point being used.
%T  : the tolerance of the final answer.
%S  : the number of iterations.
%q  : prints the steps (1 for yes, 0 for no).

format long
if q == 1
    fprintf('\nStep \t\t Result\n---\t\t --------\n');
end;
    i = 0;
    while i <= S
        x2 = x1 - feval(f,x1)*(x1-x0)/(feval(f,x1)-feval(f,x0));
        if q == 1
            fprintf('%d \t\t %f\n',i,x2);
        end;
        if feval(f,x2) == 0 || abs(feval(f,x2)) < T
            fprintf('\nThe solution is %d. The computation was a success after %d iterations!\n\n',x2,i);
            break;
        end;
        i = i+1;
        x0 = x1;
        x1 = x2;
        
    end;
    if i == S+1
        fprintf('\nMethod failed after %d iterations.\n\n',S)
    end;