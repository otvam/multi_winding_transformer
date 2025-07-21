function [I, solver] = get_solve_num(I_init, trf, op, tol_opt)
% Get the optimal currents with the a numerical solver (include the mutual resistances).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract values
R = trf.R;
X = trf.X;
idx_1 = trf.idx_1;
idx_2 = trf.idx_2;

% extract values
k_P = op.k_P;
k_Q = op.k_Q;
k_S = op.k_S;
P_trg = op.P_trg;
P_flow = op.P_flow;
idx_fix = op.idx_fix;

% extract options
n_start = tol_opt.n_start;
k_bounds = tol_opt.k_bounds;
solver_opt = tol_opt.solver_opt;
fix_equality = tol_opt.fix_equality;
fix_inequality = tol_opt.fix_inequality;
use_hessian = tol_opt.use_hessian;
use_gradient = tol_opt.use_gradient;
use_quadratic = tol_opt.use_quadratic;
use_bounds = tol_opt.use_bounds;

% get the optim functions
if use_quadratic==true
    % get the quadratically constrained quadratic problem
    [P, Q, H] = get_matrix_quad(R, X, idx_1, idx_2, P_flow);

    % scaling for the objective function
    P = k_P.*P;
    Q = k_Q.*Q;

    % make all matrices symmetric
    P = (P+P')./2;
    Q = (Q+Q')./2;
    H = (H+H')./2;

    % get the optim functions
    fct_obj = @(x) get_obj_quad(x, P, Q);
    fct_con = @(x) get_con_quad(x, H, P_trg);
    if use_hessian==true
        fct_hes = @(x,lambda) get_hes_quad(lambda, P, Q, H);
    else
        fct_hes = [];
    end
else
    fct_obj = @(x) get_obj_direct(x, R, X, idx_1, idx_2, k_P, k_Q, k_S);
    fct_con = @(x) get_con_direct(x, R, X, idx_1, idx_2, P_flow, P_trg);
    fct_hes = [];
end

% encode the initial guess
x_init = get_encode(I_init, idx_1, idx_2);

% find the value for the bounds
v_bounds = k_bounds.*max(abs(x_init));

% get the linear constraints
[Aineq, bineq, Aeq, beq, lb, ub] = get_matrix_linear(idx_1, idx_2, idx_fix, v_bounds);

% remove constraints if required
if fix_inequality==false
    Aineq = [];
    bineq = [];
end
if fix_equality==false
    Aeq = [];
    beq = [];
end
if use_bounds==false
    lb = [];
    ub = [];
end

% solver base options
options = optimoptions(...
    'fmincon', ...
    'Algorithm','interior-point', ...
    'ScaleProblem', true, ...
    'SpecifyObjectiveGradient', use_gradient&&use_quadratic, ...
    'SpecifyConstraintGradient', use_gradient&&use_quadratic, ...
    'HessianFcn', fct_hes, ...
    'Display', 'off');

% solver additional options
field = fieldnames(solver_opt);
for i=1:length(field)
    options.(field{i}) = solver_opt.(field{i});
end

% create the optimization problem
problem = createOptimProblem(...
    'fmincon', ...
    x0=x_init, ...
    lb=lb, ...
    ub=ub, ...
    Aeq=Aeq, ...
    beq=beq, ...
    Aineq=Aineq, ...
    bineq=bineq, ...
    objective=fct_obj, ...
    nonlcon=fct_con, ...
    options=options);

% get the global optimizer
gs = MultiStart('Display', 'off');

% solve the problem
[x, fval, exitflag, output] = run(gs, problem, n_start);

% decode the solution
I = get_decode(x, idx_1, idx_2);

% assign the solver results
solver.fval = fval;
solver.exitflag = exitflag;
solver.n_eval = output.funcCount;
solver.n_trial = output.localSolverTotal;
solver.n_converged = output.localSolverSuccess;
solver.n_eval = output.funcCount;

end

function y = get_obj_direct(x, R, X, idx_1, idx_2, k_P, k_Q, k_S)

% get the currents
I = get_decode(x, idx_1, idx_2);

% extract the impedances
Z = R+1i.*X;

% compute the voltage
V = Z*I;

% compute the complex power
S = 0.5.*(V.*conj(I));

% compute the power metrics
P = sum(real(S));
Q = sum(imag(S));
S = sum(abs(S));

% compute objective
y = k_P.*P+k_Q.*Q+k_S.*S;

end

function [y, yeq] = get_con_direct(x, R, X, idx_1, idx_2, P_flow, P_trg)
% Computation of the constraint (power flow).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the power
S = get_power(x, R, X, idx_1, idx_2);

% compute the terminal power
S_1 = sum(S(idx_1));
S_2 = sum(S(idx_2));

% inequality constraint
y = [];

% equality constraint
switch P_flow
    case 'primary' % the input power (primary coil) is fixed
        yeq = P_trg-real(S_1);
    case 'secondary' % the output power (secondary coil) is fixed
        yeq = P_trg+real(S_2);
    otherwise
        error('invalid terminal for the power flow')
end

end

function S = get_power(x, R, X, idx_1, idx_2)
% Computation of the complex power.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the currents
I = get_decode(x, idx_1, idx_2);

% extract the impedances
Z = R+1i.*X;

% compute the voltage
V = Z*I;

% compute the complex power
S = 0.5.*(V.*conj(I));

end

function [y, grady] = get_obj_quad(x, P, Q)
% Computation of the objective (total losses and reactive power).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = 0.5.*(x'*P*x)+0.5.*(x'*Q*x);
grady = P*x+Q*x;

end

function [y, yeq, grady, gradyeq] = get_con_quad(x, H, P_trg)
% Computation of the constraint (power flow).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y = [];
grady = [];

yeq = 0.5.*(x'*H*x)-P_trg;
gradyeq = H*x;

end

function hess = get_hes_quad(lambda, P, Q, H)
% Computation of the Hessian matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hess = P+Q+lambda.eqnonlin*H;

end

function I = get_decode(x, idx_1, idx_2)
% Decode a solution vector into a complex current vector.
% The solution vector has the form: [I_1_re, I_2_re, I_1_im, I_2_im].
% The current vector has the form: [I_1, I_2].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[idx_1_r, idx_2_r, idx_1_i, idx_2_i, n_var] = get_indices(idx_1, idx_2);
assert(length(x)==n_var, 'invalid data')

I = zeros(n_var./2, 1);
I(idx_1) = x(idx_1_r)+1i.*x(idx_1_i);
I(idx_2) = x(idx_2_r)+1i.*x(idx_2_i);

end

function x = get_encode(I, idx_1, idx_2)
% Encode a complex current vector into a solution vector.
% The solution vector has the form: [I_1_re, I_2_re, I_1_im, I_2_im].
% The current vector has the form: [I_1, I_2].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[idx_1_r, idx_2_r, idx_1_i, idx_2_i, n_var] = get_indices(idx_1, idx_2);
assert(length(I)==(n_var./2), 'invalid data')

x = zeros(n_var, 1);
x(idx_1_r) = real(I(idx_1));
x(idx_2_r) = real(I(idx_2));
x(idx_1_i) = imag(I(idx_1));
x(idx_2_i) = imag(I(idx_2));

end

function [idx_1_r, idx_2_r, idx_1_i, idx_2_i, n_var] = get_indices(idx_1, idx_2)
% Get the indices of the different components of the solution vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

count = 1;
n_var = 2.*(length(idx_1)+length(idx_2));

idx_1_r = count:(count+length(idx_1)-1);
count = count+length(idx_1);

idx_2_r = count:(count+length(idx_2)-1);
count = count+length(idx_2);

idx_1_i = count:(count+length(idx_1)-1);
count = count+length(idx_1);

idx_2_i = count:(count+length(idx_2)-1);
count = count+length(idx_2);

assert((count-1)==n_var, 'invalid data')

end

function [P, Q, H] = get_matrix_quad(R, X, idx_1, idx_2, P_flow)
% Get the matrices defining the quadratically constrained quadratic problem:
%     - P: quadratic matrix for the total losses
%     - Q: quadratic matrix for the reactive power
%     - H: quadratic matrix for the power flow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the indices of the different components of the solution vector
[idx_1_r, idx_2_r, idx_1_i, idx_2_i, n_var] = get_indices(idx_1, idx_2);

% assign the losses computation matrix
P = zeros(n_var, n_var);
P(idx_1_r, idx_1_r) = R(idx_1, idx_1);
P(idx_2_r, idx_2_r) = R(idx_2, idx_2);
P(idx_1_i, idx_1_i) = R(idx_1, idx_1);
P(idx_2_i, idx_2_i) = R(idx_2, idx_2);
P(idx_1_r, idx_2_r) = R(idx_1, idx_2);
P(idx_2_r, idx_1_r) = R(idx_2, idx_1);
P(idx_1_i, idx_2_i) = R(idx_1, idx_2);
P(idx_2_i, idx_1_i) = R(idx_2, idx_1);

% assign the rectative computation matrix
Q = zeros(n_var, n_var);
Q(idx_1_r, idx_1_r) = X(idx_1, idx_1);
Q(idx_2_r, idx_2_r) = X(idx_2, idx_2);
Q(idx_1_i, idx_1_i) = X(idx_1, idx_1);
Q(idx_2_i, idx_2_i) = X(idx_2, idx_2);
Q(idx_1_r, idx_2_r) = X(idx_1, idx_2);
Q(idx_2_r, idx_1_r) = X(idx_2, idx_1);
Q(idx_1_i, idx_2_i) = X(idx_1, idx_2);
Q(idx_2_i, idx_1_i) = X(idx_2, idx_1);

% assign the power flow computation matrix
H = zeros(n_var, n_var);
H(idx_1_r, idx_2_i) = -X(idx_1, idx_2);
H(idx_1_i, idx_2_r) = +X(idx_1, idx_2);

% add the loss contributions to the power flow
switch P_flow
    case 'primary' % the input power (primary coil) is fixed
        H(idx_1_r, idx_1_r) = +R(idx_1, idx_1);
        H(idx_1_i, idx_1_i) = +R(idx_1, idx_1);
        H(idx_1_r, idx_2_r) = +R(idx_1, idx_2);
        H(idx_1_i, idx_2_i) = +R(idx_1, idx_2);
    case 'secondary' % the output power (secondary coil) is fixed
        H(idx_2_r, idx_2_r) = -R(idx_2, idx_2);
        H(idx_2_i, idx_2_i) = -R(idx_2, idx_2);
        H(idx_1_r, idx_2_r) = -R(idx_1, idx_2);
        H(idx_1_i, idx_2_i) = -R(idx_1, idx_2);
    otherwise
        error('invalid terminal for the power flow')
end

end

function [Aineq, bineq, Aeq, beq, lb, ub] = get_matrix_linear(idx_1, idx_2, idx_fix, v_bounds)
% Get the matrices defining the linear conditions:
%     - Aineq, bineq: inequality condition for fixing the solution
%     - Aeq, beq: equality condition for fixing the solution
%     - lb, ub: bounds for the variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the indices of the different components of the solution vector
[idx_1_r, idx_2_r, idx_1_i, idx_2_i, n_var] = get_indices(idx_1, idx_2);

% empty inequality condition
Aineq = zeros(1, n_var);
bineq = 0;

% empty equality condition
Aeq = zeros(1, n_var);
beq = 0;

% fix the solution for a terminal
%     - real part should be positive
%     - imaginary part should be zero
if any(idx_fix==idx_1)
    Aineq(idx_1_r(idx_fix==idx_1)) = -1;
    Aeq(idx_1_i(idx_fix==idx_1)) = 1;
end
if any(idx_fix==idx_2)
    Aineq(idx_2_r(idx_fix==idx_2)) = -1;
    Aeq(idx_2_i(idx_fix==idx_2)) = 1;
end

% bounds
lb = -v_bounds.*ones(1, n_var);
ub = +v_bounds.*ones(1, n_var);

end