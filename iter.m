function x = iter(A,x0,b,S,T,q)
%function x = iter(A,x0,b,S,T,q)
%
%This is an algorithm designed by Alexander Winkles that iteratively solves
%systems of equations of the form Ax = b using various standard methods.
%
%A  : the matrix
%x0 : the guess solution to start with
%b  : the vector
%S  : the number of iterations
%T  : the tolerance of the result using an infinity norm
%q  : indicates which method to use
%       q == 1 - Jacobi
%       q == 2 - Gauss-Seidel
%       q == 3 - steepest descent

n = size(A,1);
m = size(A,2);

if n ~= m
    fprintf('\nThis is not an n x n matrix!\n')
elseif det(A) < 1e-14
    fprintf('\nThis matrix is singular!\n\n')
else
    if q == 1
        
        % Jacobi %
        
        x = x0;
        u = zeros(n,1);
        k = 1;
        while k <= S
            for i=1:n
                sum = 0;
                for j=1:n
                    if j ~= i
                        sum = sum + A(i,j)*x(j);
                    end;
                end;
                u(i) = (b(i) - sum)/A(i,i);
            end;
            for i=1:n
                x(i) = u(i);
            end;
            err = norm(A*x-b,Inf);
            if err < T
                fprintf('\n The iteration was successful after %d iterations with ||Ax - b|| = %d!\n', k,err)
                break;
            end;
            k = k+1;
        end;
        if k == S + 1
            err = norm(A*x-b, Inf);
            fprintf('\nMethod failed after %d iterations with ||Ax-b|| = %d.\n\n',S,err)
        end;
    end;

    % Gauss-Seidel %
    
    if q == 2
        x = x0;
        k = 1;
        while k <= S
            for i=1:n
                sum = 0;
                for j=1:n
                    if j ~= i
                        sum = sum + A(i,j)*x(j);
                    end; 
                end;
                x(i) = (b(i)-sum)/A(i,i);
            end;
            err = norm(A*x-b, Inf);
            if err < T
                fprintf('\nThe iteration was successful after %d iterations with ||Ax - b|| = %d!\n', k,err)
                break;
            end;
            k = k+1;
        end;
        if k == S+1
            err = norm(A*x-b, Inf);
            fprintf('\nMethod failed after %d iterations with ||Ax-b|| = %d.\n\n',S,err)
        end;
    end;
    
    % Steepest descent %
    
    if q == 3
        k = 1;
        x = x0;
        while k <= S
            v = b - A*x;
            t = dot(v,v)/dot(v,A*v);
            x = x + t*v;
            err = norm(A*x-b, Inf);
            if err < T
                fprintf('\nThe iteration was successful after %d iterations with ||Ax - b|| = %d!\n',k,err);
                break;
            end;
            k = k+1;
        end;
        if k == S+1
            err = norm(A*x-b, Inf);
            fprintf('\nMethod failed after %d iterations with ||Ax-b|| = %d.\n\n',S,err)
        end;
    end;
    
end;