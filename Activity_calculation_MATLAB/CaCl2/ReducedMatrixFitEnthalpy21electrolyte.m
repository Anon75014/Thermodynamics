function z= ReducedMatrixFitEnthalpy21electrolyte(b, Q1, B1, C1, D1, E1)
%---------- Ideal gas constant
R = 8.314; %J/K/mol
%---- Reference temperature
Theta = 298.15;

%----------Ak coefficients from Clarke 1980
A1 = 0.198653*Theta;
I = 3*b;
p=1;
q=2;


 %---------- Molality Matrix for Lphi, Cpphi and Jphi as in Clarke
b0_LCp = -4*3/1.2*log(1+1.2*(I).^(1/2));
b1_LCp = (4*p*q)/(3*4)*(1-(1+2*(I).^(1/2)).*exp(-2*(I).^(1/2)));
b2_LCp = p*q*2*b;
b3_LCp = (p*q)^(3/2)*b.^2;
b4_LCp = 2/3*(p*q)^(2)*b.^3;
b5_LCp = 1/2*(p*q)^(5/2)*b.^4;

B_LCp = -[b0_LCp, b1_LCp, b2_LCp, b3_LCp, b4_LCp, b5_LCp]; %(- is for the -1 sign in Tref)

V1 = [A1; Q1; B1; C1; D1; E1];

z = R*B_LCp*V1;