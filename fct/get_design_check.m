function get_design_check(trf, op)
% Check the validity of a transformer design.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Guillod - Dartmouth College.
% 2025 - MIT License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract values
R = trf.R;
X = trf.X;
idx_1 = trf.idx_1;
idx_2 = trf.idx_2;

% extract values
k_P = op.k_P;
k_Q = op.k_Q;
P_trg = op.P_trg;
P_flow = op.P_flow;
idx_fix = op.idx_fix;

% validate types
validateattributes(R, {'double'},{'2d', 'nonempty', 'real','finite'});
validateattributes(X, {'double'},{'2d', 'nonempty', 'real','finite'});
validateattributes(idx_1, {'double'},{'row', 'nonempty' 'real', 'finite', 'integer'});
validateattributes(idx_2, {'double'},{'row', 'nonempty' 'real', 'finite', 'integer'});
validateattributes(idx_fix, {'double'},{'scalar', 'nonempty' 'real', 'finite', 'integer'});
validateattributes(k_P, {'double'},{'scalar', 'nonnan', 'real','finite', 'nonnegative'});
validateattributes(k_Q, {'double'},{'scalar', 'nonnan', 'real','finite', 'nonnegative'});
validateattributes(P_trg, {'double'},{'scalar', 'nonnan', 'real','finite', 'positive'});
validateattributes(P_flow, {'char'}, {'nonempty'});

% compute the size of the transformer
idx_all = sort([idx_1, idx_2]);
n = length(idx_all);

% check the winding indices
assert(all(idx_all==1:n), 'invalid winding indices')
assert(any(idx_fix==1:n), 'invalid winding indices')

% check power flow
assert(any(strcmp(P_flow, {'primary', 'secondary'})), 'invalid terminal for the power flow')

% check that the matrices are have the right size
assert(all(size(R)==[n, n]), 'resistance matrix has an incorrect size')
assert(all(size(X)==[n, n]), 'admittance matrix has an incorrect size')

% check that the matrices are positive definite
assert(issymmetric(R)&all(eig(R)>0), 'resistance matrix is not symmetric positive')
assert(issymmetric(X)&all(eig(X)>0), 'admittance matrix is not symmetric positive')

end