function [trf, op] = get_matrix_reduce(trf, op)
% Reduce a multi-winding transformer into a two-winding transformer.
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
P_trg = op.P_trg;
P_flow = op.P_flow;
idx_fix = op.idx_fix;

% reduce matrices
trf.idx_1 = 1;
trf.idx_2 = 2;
trf.R = [
    sum(sum(R(idx_1, idx_1))), sum(sum(R(idx_1, idx_2))) ; ...
    sum(sum(R(idx_2, idx_1))), sum(sum(R(idx_2, idx_2))) ; ...
    ];
trf.X = [
    sum(sum(X(idx_1, idx_1))), sum(sum(X(idx_1, idx_2))) ; ...
    sum(sum(X(idx_2, idx_1))), sum(sum(X(idx_2, idx_2))) ; ...
    ];

% assign power flow
op.P_trg = P_trg;
op.P_flow = P_flow;
if any(idx_fix==idx_1)
    op.idx_fix = 1;
end
if any(idx_fix==idx_2)
    op.idx_fix = 2;
end

end