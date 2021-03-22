function mytester(A,b,p)

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