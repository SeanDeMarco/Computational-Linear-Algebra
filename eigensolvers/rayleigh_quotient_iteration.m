function [eigenvalues, eigenvectors] = rayleight_quotient_iteration(K,M,max_iter,tol)
    
    % Transforming from general to standard eigenvealue problem form
    R = cholesky_tridiagonal(M);
    A = R' \ (K / R);
    N = size(A, 1);
    % Initialising the eigenvalue storage vector
    eigenvalues = zeros(N, 1);
    % Initilaising general eigenvalue matrix
    Q = eye(N);
    % Initialising the matrix to be deflated
    current_A = A;

    for k = 1:N       
        % Computing eigenvalue and eigenvector using RQI
        [lambda, x] = RQI(current_A,max_iter,tol);       
        eigenvalues(k) = lambda;        
        % Computing Householder reflector and deflated matrix
        [H, A_deflated] = HHD(current_A, x);        
        m = length(x);
        % Updating Q with the Householder reflector
        H_embedded = eye(N);
        cols = k:(k + m - 1);
        H_embedded(cols, cols) = H;
        Q = Q * H_embedded;        
        % Deflating the matrix
        current_A = A_deflated;
    end
    % Recovering the original problems eigenvectors due to similarity
    % transformation from the cholesky decomposition
    eigenvectors = R \ Q;
end

function [lambda, x] = RQI(A, max_iter, tol)

    n = size(A, 1);   
    % Initialising a random eigenvector
    x = rand(n, 1);
    x = x / norm(x);
    % Initialising the  Rayleigh quotient
    lambda = (x' * A * x) / (x' * x);
    sigma = lambda;  
    for iter = 1:max_iter
        % Applying a shift and usign the current eigenvalue as an estimate
        shifted_A = A - sigma * eye(n);        
        % Solving the linear system (shifted_A) * w = x - regularisation
        % was added as the shift made sigma too close to the eigenvalue,
        % making the matrix singular
        w = (shifted_A+ 1e-8 * eye(n)) \ x; 
        % Normalising the new eigenvector estimate
        x_new = w / norm(w);        
        % Updating eigenvalue estimate using Rayleigh Quotient
        lambda_new = (x_new' * A * x_new) / (x_new' * x_new);
        
        % Checking for convergence
        residual = norm(A * x_new - lambda_new * x_new);
        if residual < tol
            break;
        end

        % Updating eigenvector, eigenvalue, and shift
        x = x_new;
        lambda = lambda_new;
        sigma = lambda; % Update the shift
    end
end

function [H, A_deflated] = HHD(A, x)

    m = length(x);
    v = x(:);
    % If problem has been fully deflated, return empty values for H and A
    if m == 0
        H = [];
        A_deflated = [];
        return;
    end
    
    % Computing Householder reflector
    sigma = sign(v(1)) * norm(v);
    u = v - sigma * eye(m, 1);    
    if norm(u) < eps(class(u)) * 1e2
        % No reflection is needed
        H = eye(m);
    else
        H = eye(m) - (2 / (u' * u)) * (u * u');
    end
    
    % Applying a similarity transformaiton - eigenvalues are preserved
    % since H is orthogonal
    A_transformed = H * A * H';
    
    % Deflating the matrix
    if m > 1
        A_deflated = A_transformed(2:end, 2:end);
    else
        A_deflated = [];
    end
end
