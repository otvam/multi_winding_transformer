function [sol_eig, sol_num, I_num] = get_main_solve(tag, trf, op, tol_eig, tol_opt)
% Get the optimal current sharing with the eigenvalue and numerical methods.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the eigenvalue solution
I_eig = get_solve_eig(trf, tol_eig);
I_eig = get_power_flow(I_eig, trf, op);
sol_eig = get_eval_solution(I_eig, trf);

% get the numerical solution
[I_num, solver] = get_solve_num(I_eig, trf, op, tol_opt);
I_num = get_power_flow(I_num, trf, op);
sol_num = get_eval_solution(I_num, trf);

% show the results
fprintf('%s\n', tag)
fprintf('    solver output\n')
get_disp_solver(solver)
fprintf('    eigenvalue method\n')
get_disp_sol(sol_eig)
fprintf('    numerical solver\n')
get_disp_sol(sol_num)

end

function get_disp_solver(solver)
% Display a solver output.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('        fval = %.5f\n', solver.fval)
fprintf('        exitflag = %d\n', solver.exitflag)
fprintf('        n_trial = %d\n', solver.n_trial)
fprintf('        n_converged = %d\n', solver.n_converged)
fprintf('        n_eval = %d\n', solver.n_eval)

end

function get_disp_sol(sol)
% Display a solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('        I_abs = %s A\n', get_array(abs(sol.I), '%+.5f'))
fprintf('        I_angle = %s deg\n', get_array(rad2deg(angle(sol.I)), '%+.5f'))
fprintf('        V_abs = %s V\n', get_array(abs(sol.V), '%+.5f'))
fprintf('        V_angle = %s deg\n', get_array(rad2deg(angle(sol.V)), '%+.5f'))
fprintf('        P = %s W\n', get_array(real(sol.S), '%+.5f'))
fprintf('        Q = %s Var\n', get_array(imag(sol.S), '%+.5f'))
fprintf('        S_1 = %+.5f W / %+.5f VAr / %+.5f VA\n', real(sol.S_1), imag(sol.S_1), abs(sol.S_1))
fprintf('        S_2 = %+.5f W / %+.5f VAr / %+.5f VA\n', real(sol.S_2), imag(sol.S_2), abs(sol.S_2))
fprintf('        S_sum = %+.5f W / %+.5f VAr / %+.5f VA\n', real(sol.S_sum), imag(sol.S_sum), abs(sol.S_sum))
fprintf('        S_abs = %+.5f VA\n', sol.S_abs)
fprintf('        eta = %+.5f %%\n', 1e2.*sol.eta)
fprintf('        pf_sum = %+.5f %%\n', 1e2.*sol.pf_sum)
fprintf('        pf_abs = %+.5f %%\n', 1e2.*sol.pf_abs)

end

function str = get_array(vec, fmt)
% Format a array into a string.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

str = compose(fmt, vec);
str = strjoin(str, ' / ');
str = sprintf('[%s]', str);

end