%----------Reading the experimental data
Data_gamma = readmatrix('ActivityCoeffT_CaCl2.txt');
Data_phi = readmatrix('OsmoticCoeffCaCl2.txt');
Data_L = readmatrix('EnthalpyMixingCaCl2.txt');
Data_J = readmatrix('HeatCapacityCaCl2.txt');
%----------matrixes of different parameters
%activity
gamma = Data_gamma(:,2:6);
b_gamma = Data_gamma(:,1);
T_data_gamma = [273.15, 283.15, 298.15, 313.15, 333.15];
%osmotic coefficient
phi = Data_phi(:,2:6);
b_phi = Data_phi(:,1);
T_data_phi = [273.15, 283.15, 298.15, 313.15, 333.15];
%enthalpy
Lphi = Data_L(:,2:8);
b_L = Data_L(:,1);
T_data_L = [273.15, 283.15, 293.15, 298.15, 303.15, 323.15, 343.15];
%heat capacity
Jphi = -(Data_J(:,2)+275.7);
b_J = Data_J(:,1);
T_data_J = [298.15];

%-----------Data at 25°C for fitting
phi_25 = phi(:,3);
L_25 = Lphi(:,4);
J_25 = Jphi(:,1);

%----------Ak coefficients from Clarke 1980 (nu_gk is not exactly u_k)
A0 = -0.391940*298.15;
A1 = 0.198653*298.15;
A2 = 0.772533;
A3 = 2/298.15*1.68848;
A4 = 6/(298.15)^2*1.99113;

%-----Enthalpy
%General model:
     %f(b) = ReducedMatrixFitEnthalpy(b, Q1, B1, C1, D1, E1)
     %Coefficients (with 95% confidence bounds):
       B1 =        15.8 ;% (0.6358, 30.97)
       C1 =      -6.714 ;% (-19.22, 5.794)
       D1 =      -3.107 ;% (-6.162, -0.05163)
       E1 =      0.3062 ;% (0.08118, 0.5313)
       Q1 =       278.1 ;% (225.4, 330.8)

%Goodness of fit:
  %SSE: 1.934e+05
  %R-square: 0.9997
  %Adjusted R-square: 0.9997
  %RMSE: 72.3
  
  
  %-----Heat Capacity
  %General model:
     %f(b) = ReducedMatrixFitHeatCapacity21electrolyte(b, Q2, B2, C2, D2, E2)
%Coefficients (with 95% confidence bounds):
       B2 =     -0.6466 ;% (-0.9517, -0.3415)
       C2 =     -0.3304 ;% (-0.5496, -0.1112)
       D2 =     0.07392 ;% (0.02731, 0.1205)
       E2 =   -0.003602 ;% (-0.006586, -0.0006182)
       Q2 =      0.1441 ;% (-1.062, 1.351)

%Goodness of fit:
  %SSE: 101.1
 % R-square: 0.9997
  %Adjusted R-square: 0.9996
  %RMSE: 1.778
 
  %-----Osmotic coefficient
  %General model:
     %f(b) = ReducedMatrixFitOsmoticCoefficient21electrolyte(b, Q0, B0, C0, D0, E0)
%Coefficients (with 95% confidence bounds):
       B0 =       -92.6 ;% (-94.48, -90.73)
       C0 =      0.9724 ;% (0.01978, 1.925)
       D0 =     -0.5565 ;% (-0.7115, -0.4015)
       E0 =     0.05521 ;% (0.04721, 0.0632)
       Q0 =      -497.4 ;% (-516.9, -477.9)

%Goodness of fit:
  %SSE: 0.0003217
  %R-square: 1
  %Adjusted R-square: 1
  %RMSE: 0.003171
 
  %---------MOLALITY COEFFICIENTS

p=1;
q=2;

%---------- Molality Matrix for osmotic coefficient as in Clarke
I_phi = 3*b_phi;
b0_phi = -2*I_phi.^(1/2)./(1+1.2*I_phi.^(1/2));
b1_phi = (2*p*q)/(p+q)*b_phi.*exp(-2*I_phi.^(1/2));
b2_phi = (2*p*q)/(p+q)*b_phi;
b3_phi = (2*(p*q)^(3/2))/(p+q)*b_phi.^2;
b4_phi = (2*(p*q)^(2))/(p+q)*b_phi.^3;
b5_phi = (2*(p*q)^(5/2))/(p+q)*b_phi.^4;

B_phi = [b0_phi, b1_phi, b2_phi, b3_phi, b4_phi, b5_phi];

%---------- Molality Matrix for activity coefficient as in Clarke
I_gamma = 3*b_gamma;
b0_gamma = -2*(I_gamma.^(1/2)./(1+1.2*I_gamma.^(1/2))+2/1.2*log(1+1.2*I_gamma.^(1/2)));
b1_gamma = (2*p*q)/(p+q)*2/(4*3)*(1-(1+2*I_gamma.^(1/2)-2*I_gamma).*exp(-2*I_gamma.^(1/2)));
b2_gamma = (2*p*q)/(p+q)*2*b_gamma;
b3_gamma = 3/2*(2*(p*q)^(3/2))/(p+q)*b_gamma.^2;
b4_gamma = 4/3*(2*(p*q)^(2))/(p+q)*b_gamma.^3;
b5_gamma = 5/4*(2*(p*q)^(5/2))/(p+q)*b_gamma.^4;

B_gamma = [b0_gamma, b1_gamma, b2_gamma, b3_gamma, b4_gamma, b5_gamma];

%---------- Molality Matrix for Lphi
I_L = 3*b_L;
b0_L = -4*3/1.2*log(1+1.2*(I_L).^(1/2));
b1_L = (4*p*q)/(3*4)*(1-(1+2*(I_L).^(1/2)).*exp(-2*(I_L).^(1/2)));
b2_L = p*q*2*b_L;
b3_L = (p*q)^(3/2)*b_L.^2;
b4_L = 2/3*(p*q)^(2)*b_L.^3;
b5_L = 1/2*(p*q)^(5/2)*b_L.^4;

B_L = [b0_L, b1_L, b2_L, b3_L, b4_L, b5_L];

 %---------- Molality Matrix for Jphi
I_J = 3*b_J;
b0_J = -4*3/1.2*log(1+1.2*(I_J).^(1/2));
b1_J = (4*p*q)/(3*4)*(1-(1+2*(I_J).^(1/2)).*exp(-2*(I_J).^(1/2)));
b2_J = p*q*2*b_J;
b3_J = (p*q)^(3/2)*b_J.^2;
b4_J = 2/3*(p*q)^(2)*b_J.^3;
b5_J = 1/2*(p*q)^(5/2)*b_J.^4;

B_J = [b0_J, b1_J, b2_J, b3_J, b4_J, b5_J];

%---------TEMPERATURE COEFFICIENTS
%---Reference temperature
Theta = 298.15;

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for Lphi
T0_L = zeros(size(T_data_L));
T1_L = -ones(size(T_data_L));
T2_L = -(T_data_L-Theta);
T3_L = -(1/2.*(T_data_L-Theta).^2);
T4_L = -(1/6.*(T_data_L-Theta).^3);

T_L = [T0_L; T1_L; T2_L; T3_L; T4_L];

%-----Temperature matrix for Jphi
T0_J = zeros(size(T_data_J));
T1_J = zeros(size(T_data_J));
T2_J = -ones(size(T_data_J));
T3_J = -(T_data_J-Theta);
T4_J = -(1/2.*(T_data_J-Theta).^2);

T_J = [T0_J; T1_J; T2_J; T3_J; T4_J];

%-----Temperature matrix for gamma %---Problem with - signs
x_gamma = (T_data_gamma-Theta)./Theta;

T0_gamma = -1/Theta.*ones(size(T_data_gamma));
T1_gamma = (1/Theta.*x_gamma./(x_gamma+1));
T2_gamma = (x_gamma.^2.*vpa(S2(x_gamma)));
T3_gamma = (1/2*Theta.*x_gamma.^3.*vpa(S3(x_gamma)));
T4_gamma = (1/6*Theta^2.*x_gamma.^4.*vpa(S4(x_gamma)));

T_gamma = [T0_gamma; T1_gamma; T2_gamma; T3_gamma; T4_gamma];

%-----Temperature matrix for phi %---Problem with - signs
x_phi = (T_data_phi-Theta)./Theta;

T0_phi = -1/Theta.*ones(size(T_data_phi));
T1_phi = (1/Theta.*x_phi./(x_phi+1));
T2_phi = (x_phi.^2.*vpa(S2(x_phi)));
T3_phi = (1/2*Theta.*x_phi.^3.*vpa(S3(x_phi)));
T4_phi = (1/6*Theta^2.*x_phi.^4.*vpa(S4(x_phi)));


T_phi = [T0_phi; T1_phi; T2_phi; T3_phi; T4_phi];

%Virial matrix
V_reduced = [A0, A1, A2, A3, A4; Q0, Q1, Q2, 0, 0; B0, B1, B2, 0, 0; C0, C1, C2, 0, 0; D0, D1, D2, 0, 0; E0, E1, E2, 0, 0];

%calculation of osmotic coefficient and activity
phi_calc_reduced = double(1+B_phi*V_reduced*T_phi);
gamma_calc_reduced = double(exp(B_gamma*V_reduced*T_gamma));

%{
%------ Plotting
%Phi as a function of b and T
figure(1)
[x1, y1, z1] = prepareSurfaceData( T_data_gamma, b_gamma, gamma_calc_reduced);
[x2, y2, z2] = prepareSurfaceData( T_data_gamma, b_gamma, gamma);
h1 = scatter3(x1, y1, z1);
hold on
h2 = scatter3(x2, y2, z2);
legend([h1,h2],'Calculated','data', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'none' );
ylabel( 'b (mol/kg)', 'Interpreter', 'none' );
zlabel( 'phi', 'Interpreter', 'none' );
view( -53.2, 26.2 );
hold off

gamma_CaCl2_calc = horzcat(b_gamma, gamma_calc_reduced);
writematrix(gamma_CaCl2_calc,'gammaCaClCalc.csv')
gamma_CaCl2_exp = horzcat(b_gamma, gamma);
writematrix(gamma_CaCl2_exp,'gammaCaClExp.csv') 

phi_CaCl2_calc = horzcat(b_phi, phi_calc_reduced);
writematrix(phi_CaCl2_calc,'phiCaClCalc.csv')
phi_CaCl2_exp = horzcat(b_phi, phi);
writematrix(phi_CaCl2_exp,'phiCaClExp.csv') 

%residuals
residuals_phi_CaCl2 = (phi_calc_reduced-phi)./phi*100;
residuals_gamma_CaCl2 = (gamma_calc_reduced-gamma)./gamma*100;
residuals_phi= horzcat(b_phi, residuals_phi_CaCl2);
residuals_gamma= horzcat(b_gamma, residuals_gamma_CaCl2);
writematrix(residuals_phi,'phiCaClresiduals.csv') 
writematrix(residuals_gamma,'gammaCaClresiduals.csv')
%}

%-------SECTION 2: FREEZING POINT

%-----Bodnar
%     f(x) = 273.15+a1*x^(1/2)+a2*x+a3*x^(3/2)+a4*x^2
%Coefficients (with 95% confidence bounds):
       t1 =       3.389 ;% (2.451, 4.326)
       t2 =      -15.12 ;% (-17.21, -13.03)
       t3 =       11.33 ;% (9.829, 12.82)
       t4 =      -5.539 ;% (-5.885, -5.193)
syms x
T_f(x) = t1*x^(1/2)+t2*x+t3*x^(3/2)+t4*x^2+ 273.15;

%Freezing Point of KCl evaluated at input molalities
T_fus = zeros(1,size(b_phi,1));
m = size(b_phi,1);
for i=1:m
    T_fus(1,i) = T_f(b_phi(i,1));
end
%estimation of properties at the freezing point: each property needs to be
%estimated at a given b and only the corresponding T_f


 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for phi
x_phi_fus = (T_fus-Theta)./Theta;

T0_phi_fus = -1/Theta.*ones(size(T_fus));
T1_phi_fus = (1/Theta.*x_phi_fus./(x_phi_fus+1));
T2_phi_fus = (x_phi_fus.^2.*vpa(S2(x_phi_fus)));
T3_phi_fus = (1/2*Theta.*x_phi_fus.^3.*vpa(S3(x_phi_fus)));
T4_phi_fus = (1/6*Theta^2.*x_phi_fus.^4.*vpa(S4(x_phi_fus)));

T_phi_fus = [T0_phi_fus; T1_phi_fus; T2_phi_fus; T3_phi_fus; T4_phi_fus];


%----43x43 matrix where all parameters have been evaluated on all
%molalities at all temperatures
phi_Tf = 1+B_phi*V_reduced*T_phi_fus;


%----43x1 matrixes where only the value at (b(Tf),Tf) is selected
phi_fp = zeros(size(b_phi));
for i=1:m
    phi_fp(i,1) = phi_Tf(i,i);
end

figure(7)
plot(b_phi, phi_fp, b_phi, phi(:,5))

phi_CaCl_fp = horzcat(b_phi, phi_fp, phi(:,4));
writematrix(phi_CaCl_fp ,'phiCaClFP.csv')