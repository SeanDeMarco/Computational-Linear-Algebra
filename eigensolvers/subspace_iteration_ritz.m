function [eigenvalues, eigenvectors] = Q2c(K, M, m, max_iter, tol)
    
    n = size(K, 1);
    % Initialising random trial vectors of size n x m
    X = rand(n, m);
    % Mass normalising the initial vectors
    for i = 1:m
        X(:, i) = X(:, i) / sqrt(X(:, i)' * M * X(:, i));
    end

    % Performing LU decomposition on the K matrix - this is a special LU
    % decomposition taking advantage of the tridiagonal nature of K
    [L, U] = tridiag_lu(K);
    
    % Performing the subspace iteration method
    [eigenvalues, eigenvectors] = SSI(K, M, X, L, U, m, max_iter, tol);
    
end

function [eigenvalues, eigenvectors] = SSI(K, M, X, L, U, m, max_iter, tol)
    % Initialising covnergence test eigenvalues
    eigenvalues_prev = zeros(m, 1);
    for iter = 1:max_iter
        % Solving the triagonal system using the LU decomposed K, forward
        % and bacward subsitution. Solving K * X_new = M * X
        X_new = zeros(size(X));
        for i = 1:m
            X_new(:, i) = tridiag_solve(L, U, M * X(:, i));
        end

        % Performing the RITZ projection step - returns m eigenmodes
        [X_new, eigenvalues_ritz] = RITZ(K, M, X_new, m);

        % Checking convergence
        if iter > 1 && max(abs(eigenvalues_ritz - eigenvalues_prev)) < tol
            break;
        end
        eigenvalues_prev = eigenvalues_ritz;
        X = X_new;
    end

    eigenvalues = eigenvalues_ritz(1:m);
    eigenvectors = X_new(:, 1:m);
end

function [X_orthogonal, eigenvalues] = RITZ(K, M, X, m)
    % Transforming K and M to RITZ subspace
    K_hat = X' * K * X;
    M_hat = X' * M * X;
    % Solving the reduced eigenvalue problem. Here any method to solve the
    % problem can be used but must be suitable for the application itself.
    % SOlving K_hat * P = M_hat * P * D
    [P, D] = eig(K_hat, M_hat);
    eigenvalues = diag(D);    
    [eigenvalues, idx] = sort(eigenvalues);
    P = P(:, idx);  
    % Transforming the RITZ eigenvectors back to the original eigenspace
    X_orthogonal = X * P;    
    % Mass normalising the eigenvectors
    for i = 1:size(X_orthogonal, 2)
        X_orthogonal(:, i) = X_orthogonal(:, i) / sqrt(X_orthogonal(:, i)' * M * X_orthogonal(:, i));
    end
    eigenvalues = eigenvalues(1:m);
    X_orthogonal = X_orthogonal(:, 1:m);
end

function [L, U] = tridiag_lu(K)
    n = size(K, 1);
    a = diag(K, -1);
    d = diag(K);
    c = diag(K, 1);

    L = eye(n);
    U = zeros(n);
    U(1,1) = d(1);
    U(1,2) = c(1);

    for i = 2:n
        L(i, i-1) = a(i-1) / U(i-1, i-1);
        U(i, i) = d(i) - L(i, i-1) * c(i-1);
        if i < n
            U(i, i+1) = c(i);
        end
    end
end

function x = tridiag_solve(L, U, b)
    n = length(b);
    y = zeros(n, 1);
    
    y(1) = b(1);
    for i = 2:n
        y(i) = b(i) - L(i, i-1) * y(i-1);
    end
    
    x = zeros(n, 1);
    x(n) = y(n) / U(n, n);
    for i = n-1:-1:1
        x(i) = (y(i) - U(i, i+1) * x(i+1)) / U(i, i);
    end
end