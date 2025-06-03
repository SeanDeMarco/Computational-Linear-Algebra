function A = random_A(m)   
    % Creating the case that m =1
    if m == 1
        % Positive number between 0.1 and 1
        A(1,1) = 0.1 + 0.9*rand();  
        return;
    end    
    % Generating random off-diagonal elements for tridiagonal part 
    % (between -1 and 1)
    b = 2*rand(1, m-1) - 1;  
    c = 2*rand() - 1;        
    
    % Computing the diagonal entries such that there 
    % is strict diagonal dominance
    diag_entries = zeros(1, m);
    for i = 1:m
        if i == 1
            sum_off_diag = abs(b(1)) + abs(c);
        elseif i == m
            sum_off_diag = abs(b(end)) + abs(c);
        else
            sum_off_diag = abs(b(i-1)) + abs(b(i));
        end
        % Ensuring diagonal dominance with a margin of at least 0.1 -
        % ensures diagonal entry is not exactly equal to sum of off
        % diagonal entries
        diag_entries(i) = sum_off_diag + 0.1 + 0.9*rand(); 
    end
    
    % Building the matrix
    A = diag(diag_entries);
    for i = 1:m-1
        A(i, i+1) = b(i);
        A(i+1, i) = b(i);
    end
    A(1, m) = A(1, m) + c;
    A(m, 1) = A(m, 1) + c;   
    
    % Checking if all the eigenvalues are positive - Hint condition
    eigenvalues = eig(A);
    % Adding a tolerance for floating-point errors
    tolerance = 1e-10; 
    if any(eigenvalues <= tolerance)
        warning('Matrix might not be positive definite!');
    end

end