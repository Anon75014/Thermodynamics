function z= ReducedMatrixFitOsmoticCoefficient21electrolyte(b, Q0, B0, C0, D0, E0)
%---------- Ideal gas constant
R = 8.314; %J/K/mol
%---- Reference temperature
Theta = 298.15;

%----------Ak coefficients from Clarke 1980
A0 = -0.391940*298.15;

I = 3*b;
p=1;
q=2;

%---------- Molality Matrix for osmotic coefficient as in Clarke
b0_phi = -2*I.^(1/2)./(1+1.2*I.^(1/2));
b1_phi = (2*p*q)/(p+q)*b.*exp(-2*I.^(1/2));
b2_phi = (2*p*q)/(p+q)*b;
b3_phi = (2*(p*q)^(3/2))/(p+q)*b.^2;
b4_phi = (2*(p*q)^(2))/(p+q)*b.^3;
b5_phi = (2*(p*q)^(5/2))/(p+q)*b.^4;

B_phi = -1/Theta.*[b0_phi, b1_phi, b2_phi, b3_phi, b4_phi, b5_phi]; %(-Theta is for the -1/Theta term in T)

V0 = [A0; Q0; B0; C0; D0; E0];

z = 1+B_phi*V0;