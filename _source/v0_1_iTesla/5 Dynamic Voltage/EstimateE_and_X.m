function F = EstimateE_and_X(z)

global P Q V 

% Unknown variables,
% one delta variable for each operating point.
n=size(P,1);
delta = z(1:n);
E  = z(n+1);
X  = z(n+2);

% Power Flow Equations
F = [P*X - E*V.*sin(delta)
     Q*X - E*V.*cos(delta) + V.^2];