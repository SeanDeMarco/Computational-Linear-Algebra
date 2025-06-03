% Initialising m and timing list
m_list = 6:2:1000;
t = zeros(size(m_list));
tick=1;

% Creating random function A, performing the cholesky algorithm and timing
% for each size of mxm matrix
for m = m_list
    A = random_A(m);
    Q = sparse_chol(A);
    f = @() sparse_chol(A);
    t(tick) = timeit(f);
    tick=tick+1;    
end
%% 
% Fitting a line to the answer
p = polyfit(m_list(270:end), t(270:end), 1);
y_fit_2 = polyval(p, m_list(270:end));
SS_res = sum((t(270:end) - y_fit_2).^2);       
SS_tot = sum((t(270:end) - mean(t(270:end))).^2);     
R2 = 1 - (SS_res / SS_tot);  

% Plotting
figure;
hold on
scatter((m_list(270:end)), (t(270:end)));
plot(m_list(270:end),y_fit_2)
legend('Run time', 'Trendline')
xlabel('Matrix Size (m)');
ylabel('Execution Time');
grid on;
disp(R2)
