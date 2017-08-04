function F = Estimate(z)

global PP QQ VV 

% Unknown variables,
% one delta variable for each operating point.
n=size(PP,1);
delta = z(1:n);
E  = z(n+1);
X  = z(n+2);

% Power Flow Equations
F = [PP*X - E*VV.*sin(delta)
     QQ*X - E*VV.*cos(delta) + VV.^2];