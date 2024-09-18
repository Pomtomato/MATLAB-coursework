a = [9 4 1; 1 6 0; 1 -2 -6];
b = [-17; 4; 14];
n = 3;
C = [a b];


ans_jacobi = jacobi_method(a,b,n)

ans_seidel = gauss_seidel_method(a,b,n)
ans_gauss_elimination = gauss_elimination(C)
ans_jordan_elimination = gauss_jordan_elimination(C)
ans_jordan_easy = gauss_jordan_easy(C)

function x = jacobi_method(a,b,n)
    for i = 1:n
        x(i) = 0;
        xold(i) = 0;
        xnew(i) = 0;
    end
    
    sum = 0;
    err = 0;
    iteration = 0;
    
    for k = 1:100
            for p =1:n
                x(p) = xnew(p);
            end
            for i = 1:n       
                    for j =1:n                    
                        if j ~= i
                            sum = sum + a(i,j)*x(j);
                        end
                    end
                    xnew(i) = (b(i) - sum)/a(i,i);
                    sum = 0.0;
                    err = err + abs((xnew(i) - xold(i))/xnew(i));
                    xold(i) = xnew(i);              
            end      
        if err <= 0.0000001
            break;
        end
        err = 0.0;
        iteration = iteration +1;   
        values(:,iteration) = x;
    end
    
    values;
    figure(1);
    plot(values(1,:),"Marker","diamond")
    xlabel("No of Iterations")
    ylabel("Functional Values")
    grid on
    hold on
    plot(values(2,:),"Marker","hexagram")
    plot(values(3,:),"Marker","square")
    hold off
    legend("X1","X2","X3")
    
    fprintf("The solution by Jacobi Method after %d itertion is ", iteration);
    x = xnew';

end

function x = gauss_seidel_method(a,b,n)
    for i = 1:n
        x(i) = 0;
        xold(i) = 0;
        
    end
    
    sum = 0;
    err = 0;
    iteration = 0;
    
    for k = 1:100
            for i = 1:n
                    for j =1:n
                        if j ~= i
                            sum = sum + a(i,j)*x(j);
                        end
                    end
                    x(i) = (b(i) - sum)/a(i,i);
                    sum = 0.0;
                    err = err + abs((x(i) - xold(i))/x(i));
                    xold(i) = x(i);
                         
            end
        if err <= 0.0000001
            break;
        end
        err = 0.0;
        iteration = iteration +1; 
        val(:,1) = [0 0 0]';
        val(:,k+1) = x;
        
    end
    
    val;
    figure(2);
    plot(val(1,:),"Marker","diamond")
    xlabel("No of Iterations")
    ylabel("Functional val")
    grid on
    hold on
    plot(val(2,:),"Marker","hexagram")
    plot(val(3,:),"Marker","square")
    hold off
    legend("X1","X2","X3")
    
    fprintf("The solution by Gauss-Seidel Method after %d itertion is ", iteration);
    x = x';
end

function x = gauss_elimination(C)
    n = size(C,1);
    x = zeros(n,1);
    
    for i = 1:n-1
        for j = i+1:n
            factor = C(j,i)/C(i,i);
            C(j,:) = C(j,:) - factor*C(i,:);
        end
    end
    
    x(n) = C(n,n+1)/C(n,n);
    
    for i = n-1:-1:1
        summ = 0;
        for j = i+1 : n
            summ = summ + C(i,j)*x(j,:);
            x(i,:) = (C(i,n+1) - summ) / C(i,i);
        end
    end
    
    disp('The required solution is:')
    x;
end

function x = gauss_jordan_elimination(P)



    [ row col ] = size( P); 
      for i = 1:row-1 % Finding zeros of lower triangular matrix.
            if P(i,i) == 0 % checking wheather diagonal elements are all zeros or not
                disp(' Gauss elimination method can not applicale. Rearrange the equations!!!') 
               return
            end
         a=P(i,i);  
         P(i,:)= P(i,:)/a; % normalizing the diagonal entries
        for j=i+1:row     
          P(j,:)= P(j,:)- P(j,i)* P(i,:);
        end
      end
    
    a=P(row,row);  
    P(row,:)= P(row,:)/a;
    for i=row:-1:2   % Finding zeros of the upper triangular matrix.
        for j=i-1:-1:1    
          P(j,:)= P(j,:)- P(j,i)* P(i,:);
        end
    end 
    disp('The required solution is:')
    x = P(:,col);
end

function x = gauss_jordan_easy(C)
    R = rref(C);
    [m,n] = size(C);
    x = zeros(m,1);
    
    for i = 1:m
        x(i) = R(i,n);
    end
    disp('The required solution is:')
    x;
end