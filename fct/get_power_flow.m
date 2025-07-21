function I = get_power_flow(I, trf, op)
% Fix the power flow direction and magnitude.
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

% get the power flow and correct the direction
P_cmp = get_terminal(I, R, X, idx_1, idx_2, P_flow);
if sign(P_trg)~=sign(P_cmp)
    I(idx_2) = I(idx_2).*exp(+1i.*pi);
end

% get the power flow and the amplitude
P_cmp = get_terminal(I, R, X, idx_1, idx_2, P_flow);
if P_trg~=P_cmp
    I = sqrt(abs(P_trg./P_cmp)).*I;
end

% align the phase
I = I.*exp(-1i.*angle(I(idx_fix)));

end

function P_cmp = get_terminal(I, R, X, idx_1, idx_2, P_flow)
% Compute the power flow at the primary coil.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract the matrices
X_m = X(idx_1, idx_2);
R_m = R(idx_1, idx_2);
R_1 = R(idx_1, idx_1);
R_2 = R(idx_2, idx_2);

% extract the current
I_1 = I(idx_1);
I_2 = I(idx_2);

% compute the power
P_trf = 0.5.*(imag(I_1)'*X_m*real(I_2))-0.5.*(real(I_1)'*X_m*imag(I_2));
P_1 = 0.5.*(real(I_1)'*R_1*real(I_1))+0.5.*(imag(I_1)'*R_1*imag(I_1));
P_2 = 0.5.*(real(I_2)'*R_2*real(I_2))+0.5.*(imag(I_2)'*R_2*imag(I_2));
P_cross = 0.5.*(real(I_1)'*R_m*real(I_2))+0.5.*(imag(I_1)'*R_m*imag(I_2));

% get the right power flow
switch P_flow
    case 'primary' % the input power (primary coil) is fixed
        P_cmp = P_trf+P_1+P_cross;
    case 'secondary' % the output power (secondary coil) is fixed
        P_cmp = P_trf-P_2-P_cross;
    otherwise
        error('invalid terminal for the power flow')
end

end