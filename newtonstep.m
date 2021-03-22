function newtonstep(f,df,x0,S)
%newtonstep(f,df,x0,T,S)
%
%This is a program derived from Newton's method algorithm developed by Alexander Winkles used to
% find roots in functions for a certain number of steps.
%
%f  : the function being studied.
%df : the derivative of the function being studied.
%x0 : the initial point being used.
%S  : the number of steps required.

format long

i = 0;
fprintf('\nStep \t\t Result\n---\t\t --------\n');
while i <= S
    x1 = x0 - feval(f,x0)/feval(df,x0);
    fprintf('%d \t\t %f\n',i,x1);
    i = i+1;
    x0 = x1;
end;
fprintf('\nAfter %d steps, the solution is %d.\n\n',i-1,x1)

