function [sens_rho,sens_E] = SENS(m,k,eigenvectors,eigenvalues,rho0,E0,pert_rho,pert_E)
    N = length(rho0);
    % Initialising the differential sensitivities vectors
    S_rho = zeros(length(eigenvalues), N); 
    S_E = zeros(length(eigenvalues), N);  
    eigenvectors = mass_normalise(eigenvectors,m,N);
    for i = 1:length(eigenvectors)
        % Mass normalising the eigenvectors to make the solution unique
        lambda = eigenvalues(i);
        for j = 1:length(E0)
            % Since it is element wise differentiation, only the affected
            % element has a derivative value. dm and dk are therefore sparse.
            dm = zeros(N+1,N+1);
            dk = zeros(N+1,N+1);
    
            % Derivative of M w.r.t. rho0 for element j
            dm(j:j+1, j:j+1) = dm(j:j+1, j:j+1) + (m{j} / rho0(j));
    
            % Derivative of K w.r.t. E0 for element j
            dk(j:j+1, j:j+1) = dk(j:j+1, j:j+1) + (k{j} / E0(j));
    
            % Removing boundary DOF from derivative matrices
            dm = dm(2:end, 2:end);
            dk = dk(2:end, 2:end);
    
            % Sensitivity w.r.t. rho0 for element j
            S_rho(i, j) = eigenvectors(:,i)' * (-lambda * dm) * eigenvectors(:,i);
    
            % Sensitivity w.r.t. E0 for element j
            S_E(i, j) = eigenvectors(:,i)' * dk * eigenvectors(:,i);
        end    
    end
    % perturbing each element a specific amount
    sens_rho = S_rho*pert_rho';
    sens_E = S_E*pert_E';
end

function eigenvectors = mass_normalise(eigenvectors,m,N)
        
        M = zeros(N+1, N+1);
        % Assemble with overlapping nodes
        for i = 1:N
            nodes = [i, i+1]; % Nodes for element i
            M(nodes, nodes) = M(nodes, nodes) + m{i};
        end
        M = M(2:end, 2:end);
        for i = 1:N
            v = eigenvectors(:, i);
            m_i = v' * M * v;
            vnormd = v / sqrt(m_i);
            eigenvectors(:,i) = vnormd;
        end
end