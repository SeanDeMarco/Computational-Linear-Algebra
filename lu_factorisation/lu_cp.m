function [LU, p, q] = lu_cp(A)
    % Initialising the size and permutation vectors
    [n, ~] = size(A);
    p = 1:n;  
    q = 1:n;  
    
    for k = 1:n-1
        % Finding the pivot (largest absolute value in submatrix A(k:n, k:n))
        [row_max, row_idx] = max(abs(A(k:n, k:n)), [], 1);  
        [~, col_idx] = max(row_max);  
        row_idx = row_idx(col_idx) + k - 1;
        col_idx = col_idx + k - 1;
        
        % Checking if the matrix is invertible
        if A(row_idx, col_idx) == 0
            error("A is not invertible");
        end

        % Swapping rows in A and updating permutation vector p
        if row_idx ~= k
            temp = A(k, :);
            A(k, :) = A(row_idx, :);
            A(row_idx, :) = temp;
            
            % Swapping elements in row permutation vector p
            temp = p(k);
            p(k) = p(row_idx);
            p(row_idx) = temp;
        end

        % Swapping columns in A and updating permutation vector q
        if col_idx ~= k
            temp = A(:, k);
            A(:, k) = A(:, col_idx);
            A(:, col_idx) = temp;
            
            % Swapping elements in column permutation vector q
            temp = q(k);
            q(k) = q(col_idx);
            q(col_idx) = temp;
        end
        
        % Gaussian elimination
        for i = k+1:n
            A(i, k) = A(i, k) / A(k, k);  % Store L factor in A
            A(i, k+1:n) = A(i, k+1:n) - A(i, k) * A(k, k+1:n);  % Update U
        end      
    end

    % Return LU decomposition and permutation vectors
    LU = A;
end

