function Q = sparse_chol(A)   
    % Initialising the problem
    m = size(A,1);
    Q = zeros(m);  
    
    % Hard coding the results for row 1
    Q(1,1) = sqrt(A(1,1));
    Q(1,2) = A(1,2) / Q(1,1);
    Q(1,m) = A(1,m) / Q(1,1);
    
    % Exploiting that A is a tridiagonal sparse matrix - rows 2 to m-2
    for i = 2 : m-2
        Q(i,i) = sqrt( A(i,i) - Q(i-1,i)^2 );
        Q(i,i+1) = A(i,i+1) / Q(i,i);
        Q(i,m) = - ( Q(i-1,i) * Q(i-1,m) ) / Q(i,i);
    end
    
    % Hard coding the results for row m-1
    i = m-1;
    Q(i,i) = sqrt( A(i,i) - Q(i-1,i)^2 );
    Q(i,m) = ( A(i,m) - Q(i-1,i) * Q(i-1,m) ) / Q(i,i);
    
    % Calculating the (m,m) point
    s = 0;
    for i = 1 : m-1
        s = s + Q(i,m)^2;
    end
    Q(m,m) = sqrt( A(m,m) - s );
end
