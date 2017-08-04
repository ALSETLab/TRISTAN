function F = Estimate_a_b(z)
global PP QQ  

% Unknown variables,
% one delta variable for each operating point.

a  = z(1);
b  = z(2);

% Power Flow Equations
F = [a+b*PP-QQ];