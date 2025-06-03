% Store stiffness submatrices in a cell array
k = cell(4, 1); % 4 elements
k{1} = [5.4931, -5.4931; -5.4931, 5.4931]; % k1
k{2} = [4.0024, -4.0024; -4.0024, 4.0024]; % k2
k{3} = [2.9024, -2.9024; -2.9024, 2.9024]; % k3
k{4} = [2.0932, -2.0932; -2.0932, 2.0932]; % k4

% Store mass submatrices in a cell array
m = cell(4, 1); % 4 elements
m{1} = [2.2484, 1.0087; 1.0087, 1.8052]; % m1
m{2} = [1.4473, 0.6486; 0.6486, 1.1594]; % m2
m{3} = [0.9273, 0.4151; 0.4151, 0.7410]; % m3
m{4} = [0.5910, 0.2642; 0.2642, 0.4709]; % m4

% Original parameters
rho0 = [3, 3, 3, 3]; % Density parameters for each element
E0 = [4, 4, 4, 4];   % Young's modulus parameters for each element

pert_rho = (1/3).*rho0;
pert_E = (1/4).*E0;

% Initialising the global stiffness and mass matrices
dof = 5;
K = zeros(dof, dof);
M = zeros(dof, dof);

% Assembling global stiffness matrix
K(1:2, 1:2) = K(1:2, 1:2) + k{1};
K(2:3, 2:3) = K(2:3, 2:3) + k{2};
K(3:4, 3:4) = K(3:4, 3:4) + k{3};
K(4:5, 4:5) = K(4:5, 4:5) + k{4};

% Assembling global mass matrix
M(1:2, 1:2) = M(1:2, 1:2) + m{1};
M(2:3, 2:3) = M(2:3, 2:3) + m{2};
M(3:4, 3:4) = M(3:4, 3:4) + m{3};
M(4:5, 4:5) = M(4:5, 4:5) + m{4};

% Removing first row and column to exclude boundary DOF - to match K and M
% in coursework sheet
K = K(2:end, 2:end);
M = M(2:end, 2:end);

% Compute eigenvalues and eigenvectors
[eigenvectors, eigenvalues] = eig(K, M);
eigenvalues = diag(eigenvalues); % Extract diagonal elements


%% Test 1: Uniform perturbations (all elements)
pert_magnitudes = [0.2, 0.4, 0.6, 0.8, 1.0];  
num_perturbations = length(pert_magnitudes);
num_eigenvalues = length(eigenvalues);

sens_rho_uniform = zeros(num_eigenvalues, num_perturbations);
sens_E_uniform = zeros(num_eigenvalues, num_perturbations);

for p = 1:num_perturbations
    pert_rho = pert_magnitudes(p) * ones(size(rho0));  
    pert_E = pert_magnitudes(p) * ones(size(E0));      
    
    [sens_rho_uniform(:,p), sens_E_uniform(:,p)] = ...
        SENS(m, k, eigenvectors, eigenvalues, rho0, E0, pert_rho, pert_E);
end

figure;
plot(pert_magnitudes, sens_rho_uniform', 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Density Perturbation Magnitude (Δρ)');
ylabel('Sensitivity to Pertubation');
legend(arrayfun(@(i) sprintf('λ_%d', i), 1:num_eigenvalues, 'UniformOutput', false));
grid on;

figure;
plot(pert_magnitudes, sens_E_uniform', 's-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('Youngs Modulus Perturbation Magnitude (ΔE)');
ylabel('Sensitivity to Pertubation');
legend(arrayfun(@(i) sprintf('λ_%d', i), 1:num_eigenvalues, 'UniformOutput', false));
grid on;
%% Test 2: Single-element perturbations
num_elements = length(rho0);
sens_rho_single = zeros(num_eigenvalues, num_elements);
sens_E_single = zeros(num_eigenvalues, num_elements);

for j = 1:num_elements
    
    pert_rho = zeros(size(rho0));
    pert_rho(j) = 1.0; 
    
    pert_E = zeros(size(E0));
    pert_E(j) = 1.0;   
    
    [sens_rho_single(:,j), sens_E_single(:,j)] = ...
        SENS(m, k, eigenvectors, eigenvalues, rho0, E0, pert_rho, pert_E);
end

figure;
bar(sens_rho_single', 'FaceColor', 'flat');
xlabel('Element Index of Pertubed Element');
ylabel('Sensitivity to a Unit ρ Pertubation in a Single Element');
legend(arrayfun(@(i) sprintf('λ_%d', i), 1:num_eigenvalues, 'UniformOutput', false));
grid on;
colormap(parula(num_eigenvalues));

figure;
bar(sens_E_single', 'FaceColor', 'flat');
xlabel('Element Index of Perturbed Element');
ylabel('Sensitivity to a E Unit Pertubation in a Single Element');
legend(arrayfun(@(i) sprintf('λ_%d', i), 1:num_eigenvalues, 'UniformOutput', false));
grid on;
colormap(parula(num_eigenvalues));