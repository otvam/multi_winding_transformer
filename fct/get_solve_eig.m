function I = get_solve_eig(trf, tol_eig)
% Get the optimal currents with the eigenvalue method (neglect the mutual resistances).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Guillod - Dartmouth College.
% 2025 - MIT License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract values
R = trf.R;
X = trf.X;
idx_1 = trf.idx_1;
idx_2 = trf.idx_2;

% extract the matrices
X_m = X(idx_1, idx_2);
R_1 = R(idx_1, idx_1);
R_2 = R(idx_2, idx_2);

% get the lamda matrices
Q_1 = R_1\(X_m*(R_2\X_m'));
Q_2 = R_2\(X_m'*(R_1\X_m));

% compute the eigenvalues / eigenvectors
[lbd, I_1, I_2] = get_common_eig(Q_1, Q_2, tol_eig);

% get the maximum efficiency
eta_max = 1+(2./lbd)-sqrt((1+(2./lbd)).^2-1);

% get the scaling between the currents
scl = ((1-eta_max)./2).*((I_1'*X_m*I_2)./(I_2'*R_2*I_2));

% construct the current vector
I = zeros(length(idx_1)+length(idx_2), 1);
I(idx_1) = I_1;
I(idx_2) = I_2;

% add the correct scaling
I(idx_2) = I(idx_2).*scl;

% add the correct phase shift
I(idx_2) = I(idx_2).*exp(-1i.*pi./2);

end

function [lbd, ev_1, ev_2] = get_common_eig(Q_1, Q_2, tol_eig)

% compute the eigenvalues
[ev_1, lbd_1] = eig(Q_1);
[ev_2, lbd_2] = eig(Q_2);
lbd_1 = diag(lbd_1);
lbd_2 = diag(lbd_2);

% absolute tolerance
tol_abs = tol_eig.*max([lbd_1 ; lbd_2]);

% get the largest common eigenvalue
[lbd_mat_1, lbd_mat_2] = ndgrid(lbd_1, lbd_2);
idx = abs(lbd_mat_1-lbd_mat_2) < tol_abs;
lbd = max((lbd_mat_1(idx)+lbd_mat_2(idx))./2);

% get the eigenvectors
idx_lbd_1 = find(abs(lbd_1-lbd)<tol_abs, true, 'last');
idx_lbd_2 = find(abs(lbd_2-lbd)<tol_abs, true, 'last');
assert(isscalar(idx_lbd_1), 'common eigenvalue not found')
assert(isscalar(idx_lbd_2), 'common eigenvalue not found')

% extract the corresponding eigenvectors
ev_1 = ev_1(:, idx_lbd_1);
ev_2 = ev_2(:, idx_lbd_2);

end

