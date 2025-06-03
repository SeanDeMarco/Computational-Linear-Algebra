%Initialising parameters and matrices
T=50;
alp=0.5;
kap=1;
gma=0.25;
u0=1;
N_total = 10:5:200;
% Arrays to store errors at different N values
max_error_pp_values = zeros(size(N_total));
max_error_cp_values = zeros(size(N_total));
tick=1;
for N=N_total
    h=T/N;
    a1=(1-(kap*h/2));
    a0 = (1+(kap*h/2));
    RHS = zeros(N,1);
    A = zeros(N,N);
    
    % Solving for the RHS matrix
    for i=1:N
        RHS(i,1) = a0*(exp(-gma*i*h)*u0);
    end
    
    % Creating the LHS matrix A
    for i = 1:N
        A(i,i) = a1;
        A(i,N) = (alp/gma)*(1-exp(-gma*i*h));
        beta = (alp/gma)*(1-exp(-gma*i*h));
        if (1 < i) && (1< N)
            for j = 1:i-1
                A(i,j)= (-kap*h)*exp(-gma*(i-j)*h);
            end
        if i==N
            A(N,N) = a1+beta;
        end
        end
    
    end
    
    % LU factorisation
    [A_pp,p_pp] = lu_pp(A);
    [A_cp, p_cp, q_cp] = lu_cp(A);
    
    % Solving using partial pivoting - x_pp is the final answer
    RHS_pp = RHS(p_pp);
    [y_pp] = lt_solve(A_pp, RHS_pp);
    [x_pp] = ut_solve(A_pp, y_pp);
    
    % Solving using complete pivoting - x_cp is the final answer
    RHS_cp = RHS(p_cp);
    [y_cp] = lt_solve(A_cp, RHS_cp);
    [x_temp] = ut_solve(A_cp, y_cp);
    x_cp = zeros(N,1);
    % Reordering x_cp using the column permutation vector q_cp
    x_cp(q_cp) = x_temp;
  
    i_vals = (1:N)'; % Column vector of indices
    exact = u0 * ((alp + (kap - gma - alp) * exp((kap - gma) * (i_vals * h - T))) ./ ...
                  (alp + (kap - gma - alp) * exp(-(kap - gma) * T)));
    
    % Compute relative errors in one vectorized step
    error_pp = abs((exact - x_pp) ./ exact);
    error_cp = abs((exact - x_cp) ./ exact);

    % Store results for plotting
    max_error_pp_values(tick) =  max(error_pp);
    max_error_cp_values(tick) = max(error_cp);
    tick=tick+1;
end

% Plotting error profiles
figure;
plot(N_total, log(max_error_pp_values), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6);
hold on;
plot(N_total, log(max_error_cp_values), 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6);
hold off;

xlabel('Number of Discretisation Points, N');
ylabel('log(Maximum Error)');
legend('Partial Pivoting', 'Complete Pivoting');
grid on;