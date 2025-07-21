function run_sweep()
% Solve the optimal current sharing problem for a multi-winding transformer.
%
% The optimal current sharing problem is solved with two methods:
%     - Eigenvalue method.
%         - Neglect the mutual resistances.
%         - Minimize the losses.
%     - Numerical solver.
%         - Include the mutual resistances.
%         - Minimize the losses and reactive power.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add the functions to the path
addpath('fct')
close('all')

% number of sweeps
n_sweep = 25;

% scaling factor for minimizing the total losses
k_P_vec = linspace(1.0, 0.3, n_sweep);

% scaling factor for minimizing the reactive power
k_Q_vec = linspace(0.0, 0.7, n_sweep);

% scaling factor for minimizing the apparent power
k_S_vec = linspace(0.0, 0.0, n_sweep);

% numerical tolerances (see definition below)
[tol_eig, tol_opt] = get_tolerance();

% sweep the solution for different objective functions
for i=1:n_sweep
    % load the problem definition (see definition below)
    [trf_mat, op_mat] = get_problem(k_P_vec(i), k_Q_vec(i), k_S_vec(i));

    % check that the transformer design is valid
    get_design_check(trf_mat, op_mat);

    % compute the full multi-winding solution
    [sol_eig, sol_num] = get_main_solve(['sweep : ' num2str(i)], trf_mat, op_mat, tol_eig, tol_opt);

    % extract the results
    pf_num_vec(i) = sol_num.pf_sum;
    pf_eig_vec(i) = sol_eig.pf_sum;
    eta_num_vec(i) = sol_num.eta;
    eta_eig_vec(i) = sol_eig.eta;
end

% plot the results
figure()
hold('on')
plot(1e2.*pf_num_vec, 1e2.*eta_num_vec, 'r', 'LineWidth', 2.0)
plot(1e2.*pf_eig_vec, 1e2.*eta_eig_vec, 'xb', 'MarkerSize', 10.0)
grid('on')
xlabel('Power Factor (%)')
ylabel('Efficiency (%)')
legend('numerical solver', 'eigenvalue')
title('Optimal Current Sharing')

end
