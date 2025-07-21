function sol = get_eval_solution(I, trf)
% Compute the problem solution (current, voltage, and power).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Guillod - Dartmouth College.
% 2025 - MIT License.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% extract values
R = trf.R;
X = trf.X;
idx_1 = trf.idx_1;
idx_2 = trf.idx_2;

% extract the impedances
Z = R+1i.*X;

% compute the voltage
V = Z*I;

% compute the complex power
S = 0.5.*(V.*conj(I));

% compute the terminal power
S_1 = sum(S(idx_1));
S_2 = sum(S(idx_2));
S_tot = sum(S);

% compute the efficiency
if sign(real(S_1))==sign(real(S_2))
    eta = NaN;
elseif real(S_1) > 0
    eta = -real(S_2)./real(S_1);
elseif real(S_2) > 0
    eta = -real(S_1)./real(S_2);
else
    eta = NaN;
end

% assign solution
sol.eta = eta;
sol.S_1 = S_1;
sol.S_2 = S_2;
sol.S_tot = S_tot;
sol.I = I;
sol.V = V;
sol.S = S;

end
