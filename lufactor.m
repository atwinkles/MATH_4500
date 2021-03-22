function [L,U] = lufactor(A,q)
%[L,U] = lufactor(A)
%
%This is an LU factorization algorithm written by Alexander Winkles that
%performs Crout factorization.
%
%A : The matrix to be factorized
%q : Determines what specific LU factorization will be returned:
%       q == 1 : Doolittle factorization
%       q == 2 : Crout factorization
%       q == 3 : UL factorization with L unit lower triangular

if size(A,1) ~= size(A,2)
    fprintf('This matrix is not square!')
else
    n = size(A,1);
    
    L = zeros(n,n);
    U = zeros(n,n);
    if q == 1
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
    end;
        
    if q == 2
        for i=1:n
            for j=1:n
            
                sum = 0;
                m = i - 1;
                for k=1:m
                    sum = sum + (L(i,k)*U(k,j));
                end;
            
                if i == j
                    U(i,j) = 1;
                end;
            
                if i >= j
                    L(i,j) = A(i,j) - sum;
                else
                    U(i,j) = (A(i,j)-sum)/L(i,i);
                end;
            end;
        end;
    end;
    if q == 3 %IN PROGRESS
        for i=n:-1:1
            for j=n:-1:1

                for m=1:n
                    L(m,m) = 1;
                end;
                
                sum = 0;
                for k=1:n-1
                    sum = sum + (U(i,k)*L(k,j));
                end;
                
                if i > j
                    L(i,j) = (A(i,j)-sum)/U(i,i);
                else
                    U(i,j) = (A(i,j)-sum)/L(j,j);
                end;
            end;
        end;
    end;
end;    