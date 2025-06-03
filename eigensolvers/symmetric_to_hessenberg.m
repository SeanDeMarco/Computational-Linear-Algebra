function [A, U] = symmetric_to_hessenberg(A, compute_eigenvectors)
    % Input:
    %   A: Symmetric matrix (will be overwritten by its Hessenberg form)
    %   compute_eigenvectors: (Optional) If true, computes transformation matrix U
    % Output:
    %   A: Tridiagonal matrix (Hessenberg form for symmetric A)
    %   U: (Optional) Orthogonal transformation matrix s.t. A_original = U*A*U'

    n = size(A, 1);
    
    % Initialize U if eigenvectors are needed
    if nargin < 2
        compute_eigenvectors = false;
    end
    
    if compute_eigenvectors
        U = eye(n);
    else
        U = [];
    end

    for k = 1:n-2
        % Extract the column below the diagonal
        x = A(k+1:n, k);
        
        % Generate Householder reflector [v, beta] = householder(x)
        [v, beta] = householder_vector(x);
        
        % Apply P_k from the left: A := P_k * A
        A(k+1:n, k:n) = A(k+1:n, k:n) - beta * v * (v' * A(k+1:n, k:n));
        
        % Apply P_k from the right: A := A * P_k
        A(1:n, k+1:n) = A(1:n, k+1:n) - beta * (A(1:n, k+1:n) * v) * v';
        
        % Accumulate U if eigenvectors are needed
        if compute_eigenvectors
            U(k+1:n, k+1:n) = U(k+1:n, k+1:n) - beta * v * (v' * U(k+1:n, k+1:n));
        end
    end
end

% Helper function: Compute Householder reflector for a vector x
function [v, beta] = householder_vector(x)
    sigma = norm(x(2:end))^2;
    v = [1; x(2:end)];
    
    if sigma == 0
        beta = 0;
    else
        mu = sqrt(x(1)^2 + sigma);
        if x(1) <= 0
            v(1) = x(1) - mu;
        else
            v(1) = -sigma / (x(1) + mu);
        end
        beta = 2 * v(1)^2 / (sigma + v(1)^2);
        v = v / v(1);
    end
end