function x = gaussianelim(q1,A,p,b,q2)
%function x = gaussianelim(A,b,p,q)
%
%This is an algorithm made by Alexander Winkles to manipulate linear 
%systems to solve various problems involving matrices and problems of the
%form Ax = b.
%
%A  : the matrix
%b  : the vector (optional)
%p  : the permutation of the matrix, where (1,2,...,n) is the given matrix
%q1 : controls what method is used to solve problem at hand
%       if b exists:
%           q == 1 - Gaussian elimination
%           q == 2 - LU factorization
%           q == 3 - Scaled row pivoting
%       if b does not exist:
%           q == 1 - returns an inverse matrix
%           q == 2 - computes the determinant (in progress)
%q2 : displays intermediate matrices (optional)

if size(A,1) ~= size(A,2)
    fprintf('\nThis matrix is not square!\n\n')
elseif det(A) < 1e-14
    fprintf('\nThis matrix is singular!\n\n')
else
    n = size(A,1);
    A1 = A;
    
%%% Gaussian elimination %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if exist('b','var')
        if q1 == 1
        
            %Performs the Gaussian elimination to create an upper triangular
            %matrix A
        
            for k=1:n-1
                for i=k+1:n
                    b(i) = b(i) - (A(p(i),k)/A(p(k),k))*b(k);
                    z = A(p(i),k)/A(p(k),k);
                    A(p(i),k) = 0;
                    for j = k+1:n
                        A(p(i),j) = A(p(i),j) - z*A(p(k),j);
                    end;
                end;
            end;
        
            if exist('q2','var')
                fprintf('\nThe row reduced matrix is:\n\n')
                disp(A);
            end;
        
            %Solves the system Ax = b with new matrix A

            for i=n:-1:1
                sum = 0;
                for j=i+1:n
                    sum = sum + A(p(i),j)*x(j);
                end;
                x(i,1) = (b(i) - sum)/A(p(i),i);
            end;
    
        end;
    
%%% LU factorization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if q1 == 2 
            U = zeros(n,n);
            L = zeros(n,n);
            for i=1:n
                L(i,i) = 1;
                for j=i:n
                    sum = 0;
                    for k=1:i-1
                        sum = sum + L(i,k)*U(k,j);
                    end;
                
                    U(i,j) = A(i,j) - sum;
                end;
                for j=i+1:n
                    sum = 0;
                    for k=1:i-1
                        sum = sum + L(j,k)*U(k,i);
                    end;
                
                    L(j,i) = (A(j,i) - sum)/U(i,i);
                end;
            end;
        
        %Solves the system using LU results
        
            x = zeros(n,1);
            for i=1:n
                sum = 0;
                for j=1:i-1
                    sum = sum + L(i,j)*z(j);
                end;
                z(i) = b(p(i)) - sum;
            end;
       
            for i = n:-1:1
                sum = 0;
                for j=i+1:n
                    sum = sum + U(i,j)*x(j);
                end;
                x(i) = (z(i) - sum)/U(i,i);
            end;
        
            if exist('q2','var')
                fprintf('\nThe L matrix is:\n')
                disp(L)
                fprintf('\nThe U matrix is:\n')
                disp(U)
            end;
        end;
    
%%% Scaled row pivoting %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        if q1 == 3
            for i=1:n
                p(i) = i;
                s = max(abs(A'));
            end;
            for k=1:n-1
                r = abs(A(p(k),k)/s(p(k)));
                kp = k;
                for i = (k+1):n
                    t = abs(A(p(i),k)/s(p(i)));
                    if t > r, r = t; kp = i; end;
                end;
                l = p(kp); p(kp) = p(k); p(k) = l;
                for i = (k+1):n
                    A(p(i),k) = A(p(i),k)/A(p(k),k);
                    for j = (k+1):n
                        A(p(i),j) = A(p(i),j)-A(p(i),k)*A(p(k),j);
                    end;
                end;
            end;
        
            y = zeros(n,1);
            y(1) = b(p(1));
            for i = 2:n
                y(i) = b(p(i));
                for j = 1:(i-1)
                    y(i) = y(i)-A(p(i),j)*y(j);
                end;
            end;
        
            x = zeros(n,1);
            x(n) = y(n)/A(p(n),n);
            for i = (n-1):-1:1
                x(i) = y(i);
                for j = (i+1):n
                    x(i) = x(i) - A(p(i),j)*x(j);
                end;
                x(i) = x(i)/A(p(i),i);
            end;
        
            P = zeros(n,n);
            for i=1:length(p)
                for j=1:length(p)
                    if j == p(i)
                        P(i,j) = 1;
                    end;
                end;
            end;

            K = P*A1;
        
            [L,U] = lufactor(K,1);
        
            if exist('q2','var')
                fprintf('\nThe L matrix is:\n')
                disp(L)
                fprintf('\nThe U matrix is:\n')
                disp(U)
            end;
        
        end;
        
%%% Inverse matrix generator %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
    else
        if q1 == 1
            I = eye(n);
            for k=1:n-1
                for i=k+1:n
                    z = A(p(i),k)/A(p(k),k);
                    for j = 1:n
                        A(p(i),j) = A(p(i),j) - z*A(p(k),j);
                        I(p(i),j) = I(p(i),j) - z*I(p(k),j);
                    end;
                end;
            end;
            for k=n:-1:1
                for i=k-1:-1:1
                    z = A(p(i),k)/A(p(k),k);
                    for j = n:-1:1
                        A(p(i),j) = A(p(i),j) - z*A(p(k),j);
                        I(p(i),j) = I(p(i),j) - z*I(p(k),j);
                    end;
                end;
            end;
            for i=1:n
                for j=1:n
                    I(p(i),j) = I(p(i),j)/A(p(i),p(i));
                end;
                A(p(i),p(i)) = 1;
            end;
            x = I;
        end;
        
    end;
end;

