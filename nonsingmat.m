function A = nonsingmat(n,q)
%function A = nonsingmat
%
%A quick script that generates nonsingular matricies to solve.
%
%n : the size of the n*n matrix to be generated
%q : the type of matrix desired
%       q == 1 - nonsingular
%       q == 2 - positive definite
%       q == 3 - diagonal dominant
%       q == 4 - diagonally dominant & positive definite

A = randi(10,n,n);
if q == 1
    while det(A) < 1e-14
        A = randi(10,n,n);
    end;
end;
if q == 2
    B = randi(10,n,n);
    A = B'*B + eye(n);
    while det(A) < 1e-14;
        B = randi(10,n,n);
        A = B'*B + eye(n);
    end;
end;
if q == 3
    A = randi(10,n,n);
    A = A + diag(sum(abs(A),2));
end;
if q == 4
    B = randi(10,n,n);
    A = B'*B + eye(n);
    A = A + diag(sum(abs(A),2));
end;