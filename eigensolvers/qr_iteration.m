function [eigenvalues, eigenvectors] = Q2b(K, M, max_iter, tol)
    
    % Transforming from general to standard eigenvalue problem form
    R = cholesky_tridiagonal(M);
    A = R' \ (K / R); 
    n = size(A, 1);
    % Converting A to a hessenberg matrix - will be tridiagonal since A is
    % symmetric, done to make program more efficient. Functions stores the transformation matrix U
    [H, U] = hessenberg(A);
    % Initialising the orthogonal matrix as the identity matrix 
    Q = eye(n);  
    % QR iteration step
    [H, Q] = QRITER(H, Q, max_iter, tol);  % Q accumulates properly
    % Storing eigenvalues
    eigenvalues = diag(H);
    % Transforming back into the original eigenproblem space
    eigenvectors = R \ (U \ Q); 
    % Normalising eigenvectors w.r.t. M
    for i = 1:n
        eigenvectors(:, i) = eigenvectors(:, i) / sqrt(eigenvectors(:, i)' * M * eigenvectors(:, i));
    end
    
    
end

function [H, Q] = QRITER(H, Q, max_iter, tol)  
    n = size(H, 1);
    eigenvalues_old = zeros(n,1);
    for k = 1:max_iter
        % Computing a shift which is equal to the eigenvalue of the
        % bottom-right 2x2 submatric
        delta = (H(n-1, n-1) - H(n, n)) / 2;
        sign_delta = sign(delta) + (delta == 0);
        mu = H(n, n) - sign_delta * H(n, n-1)^2 / (abs(delta) + sqrt(delta^2 + H(n, n-1)^2));
        % Applying the shift only if necessary - no need for a shift if the
        % matrix is nearly diagonal
        if abs(H(n, n-1)) > tol
            sigma = mu;
        else
            sigma = 0;
        end
        % Applying a shfit to accelerate convergence
        H = H - sigma * eye(n);
        % Computing the QR decomposition using the Gram-Schmidt process
        [Q_k, R_k] = gram_schmidt(H);
        % Updating the hessenberg matrix for the next iteration and undoing
        % the shift
        H = R_k * Q_k + sigma * eye(n);
        % Accumulating Q
        Q = Q * Q_k; 
        eigenvalues_new = diag(H);
        % Checking for convergence 
        if max(abs(eigenvalues_new - eigenvalues_old)) / max(abs(eigenvalues_new)) < tol
            break;
        end
        eigenvalues_old = eigenvalues_new;
    end
end

function [Q, R] = gram_schmidt(A)
    n = size(A, 1);
    Q = zeros(n);
    R = zeros(n);
    for j = 1:n
        v = A(:, j);  
        for i = 1:j-1
            R(i, j) = Q(:, i)' * A(:, j);
            v = v - R(i, j) * Q(:, i);
        end
        R(j, j) = norm(v); 
        Q(:, j) = v / R(j, j);
    end
end