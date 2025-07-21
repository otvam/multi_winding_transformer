function run_single()
% Optimal current sharing problem for a multi-winding transformer.
%
% In a first step, the multi-winding transformer is reduced:
%     - The different primary windings are series connected.
%     - The different secondary windings are series connected.
%     - The transformer is reduced into a two-winding transformer.
%     - The optimal current distribution is computed (fixed power flow).
%
% In a second step, the complete multi-port transformer is considered:
%     - All the primary windings are driven independently.
%     - All the secondary windings are driven independently.
%     - The optimal current distribution is computed (fixed power flow).
%
% For the computation of the optimal current distribution, two methods are used:
%     - Eigenvalue method (neglect the mutual resistances).
%     - Numerical solver (include the mutual resistances).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Guillod - Dartmouth College.
% 2025 - MIT License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add the functions to the path
addpath('fct')
close('all')

% scaling factor for minimizing the total losses
k_P = 1.0;

% scaling factor for minimizing the reactive power
k_Q = 0.0;

% scaling factor for minimizing the apparent power
k_S = 0.0;

% numerical tolerances (see definition below)
[tol_eig, tol_opt] = get_tolerance();

% load the problem definition (see definition below)
[trf_mat, op_mat] = get_problem(k_P, k_Q, k_S);

% check that the transformer design is valid
get_design_check(trf_mat, op_mat);

% reduce the problem into a two-winding transformer
[trf_red, op_red] = get_matrix_reduce(trf_mat, op_mat);

% compute the reduced two-winding solution
get_main_solve('two-winding', trf_red, op_red, tol_eig, tol_opt);

% compute the full multi-winding solution
get_main_solve('multi-winding', trf_mat, op_mat, tol_eig, tol_opt);

end
