function [tol_eig, tol_opt] = get_tolerance()
% Definition of the numerical tolerance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% relative tolerance for the detecting common eigenvalues
tol_eig = 1e-6;

% numerical solver algoritm
solver = 'quadratic';

% options for the fmincon solver
solver_opt.MaxIterations = 50e3;
solver_opt.MaxFunctionEvaluations = 50e3;
solver_opt.ConstraintTolerance = 1e-7;
solver_opt.OptimalityTolerance = 1e-7;
solver_opt.StepTolerance = 1e-7;
solver_opt.ScaleProblem = true;

% options for the global/quadratic solver
tol_opt.use_hessian = true; % use (or not) the Hessian for fmincon
tol_opt.use_gradient = true; % use (or not) the gradient for fmincon
tol_opt.fix_equality = false; % set (or not) the equality constraint for fmincon
tol_opt.fix_inequality = false; % set (or not) the inequality constraint for fmincon
tol_opt.k_bounds = 10.0; % factor for determining the bounds (relative to initial guess)
tol_opt.solver_opt = solver_opt; % options for the fmincon solver
switch solver
    case 'global'
        tol_opt.use_quadratic = false; % use (or not) the quadratic function
        tol_opt.use_bounds = true; % use (or not) the bounds for fmincon
        tol_opt.n_start = 10; % number of starting points for global
    case 'quadratic'
        tol_opt.use_quadratic = true; % use (or not) the quadratic function
        tol_opt.use_bounds = false; % use (or not) the bounds for fmincon
        tol_opt.n_start = 1; % number of starting points for global
    otherwise
        error('invalid solver')
end

end
