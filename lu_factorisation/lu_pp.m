function [LU, p] = lu_pp(A)
    % Implementing the LU factorisation with partial pivoting based on
    % algorithm 4.10 (pg45 in notes)
    % Initialising size and permutation vector
    [n, ~] = size(A);
    p = 1:n;     
    for k = 1:n-1
        % Finding the pivot row index in column k
        [~, r] = max(abs(A(k:n, k)));
        r = r + k - 1;        
        % Checking if the matrix is invertible
        if A(r, k) == 0
            error("A is not invertible");
        end
        % Swapping rows k and r in A
        if r ~= k
            temp = A(k, :);
            A(k, :) = A(r, :);
            A(r, :) = temp;            
            % Swapping elements in permutation vector p
            temp = p(k);
            p(k) = p(r);
            p(r) = temp;
        end        
        % Gaussian elimination
        for i = k+1:n
            A(i, k) = A(i, k) / A(k, k);  
            A(i, k+1:n) = A(i, k+1:n) - A(i, k) * A(k, k+1:n); 
        end      
    end
    % Return combined LU matrix and permutation vector
    LU = A;
end