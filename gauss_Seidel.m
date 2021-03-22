function x = gauss_Seidel(A,b,tol,max_iters)
    n = length(A(:,1));
    m = length(A(1,:));
    x = zeros(n,1);
    iterations = 0;
    while iterations < max_iters
        for k = 1:n
            for i = 1:m
                sm = 0;
                for j = 1:n
                    if j == i
                        continue
                    end
                    sm = sm + A(i,j)*x(j);
                end
                x(i) = (b(i) - sm)/A(i,i);
            end
        end
        iterations = iterations + 1;
        if max(abs(A*x-b)) < tol
            break
        end

    end
    
    
    fprintf('Number of Iterations: %.0f \nError : %d',iterations,max(abs(A*x-b)));
end