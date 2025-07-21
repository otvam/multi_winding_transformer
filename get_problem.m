function [trf, op] = get_problem(k_P, k_Q)
% Definition of the transformer and the power flow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Guillod - Dartmouth College.
% 2025 - MIT License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% power flow value
P_trg = 1.0;

% coil where the power flow is imposed (primary or secondary)
P_flow = 'secondary';

% terminal where the phase is imposed
idx_fix = 1;

% indices of the primary windings
idx_1 = 1:10;

% indices of the secondary windings
idx_2 = 11:20;

% construct the resistance and admittance matrices
R = get_matrix_fake(idx_1, idx_2, 0.5, 0.5, 1.*[0.4, 0.3, 0.2]);
X = get_matrix_fake(idx_1, idx_2, 10.0, 10.0, 1.*[0.7, 0.4, 0.2]);

% assign transformer
trf.R = R; % resistance matrix
trf.X = X; % admittance matrix
trf.idx_1 = idx_1; % indices of the primary windings
trf.idx_2 = idx_2; % indices of the secondary windings

% assign operation
op.k_P = k_P; % scaling factor for minimizing the total losses (only used for the numerical solver)
op.k_Q = k_Q; % scaling factor for minimizing the reactive power (only used for the numerical solver)
op.P_trg = P_trg; % power flow value (imposed value)
op.P_flow = P_flow; % coil where the power flow is imposed
op.idx_fix = idx_fix; % terminal where the phase is imposed

end

function mat = get_matrix_fake(idx_1, idx_2, v_1, v_2, k)
% Construct a fake matrix with band diagonal value.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the diagonal elements
v_diag = [v_1.*ones(1, length(idx_1)), v_2.*ones(1, length(idx_2))];

% get the off-diagonal couplings
k_off = cell(1, length(k));
for i=1:length(k)
    k_off{i} = k(i).*ones(1, length(idx_1)+length(idx_2)-i);
end

% construct the matrix
mat = get_matrix_coupling(v_diag, k_off);

end

function mat = get_matrix_coupling(v_diag, k_off)
% Construct a matrix from the diagonal values and coupling coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% create the coupling
k_mat = eye(length(v_diag));
for i=1:length(k_off)
    k_mat = k_mat+diag(k_off{i}, i)+diag(k_off{i}, -i);
end

% assign the matrix
[v_mat_1, v_mat_2] = ndgrid(v_diag, v_diag);
mat = k_mat.*sqrt(v_mat_1.*v_mat_2);

end
