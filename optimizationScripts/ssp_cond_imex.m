function con = ssp_cond_imex(s, A, b, r1, At, bt, r2, enforce_positivity)

% using Higueras' ARK formulation
% r1 > 0 (SSP-Explicit)
% r2 > 0 (SSP-Implicit)

b = b(:)'; bt = bt(:)';

zero_vector = zeros(s,1);
K = [A zero_vector; b 0];
Kt = [At zero_vector; bt 0];

G = eye(s+1) + r1*K + r2*Kt;
con1 = G\(r1*K);
con2 = G\(r2*Kt);
con3 = G\ones(s+1);

con = -[con1(:); con2(:); con3(:)];

if enforce_positivity
    % enforce all the coefficients are positive
    con_butcher_coef = -[A(:); b(:); At(:); bt(:)];
    
    con = [con(:); con_butcher_coef(:)];
end

end
