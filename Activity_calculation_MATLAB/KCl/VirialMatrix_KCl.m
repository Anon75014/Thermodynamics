%----------Reading the experimental data
Data_gamma = readmatrix('ActivityCoefficientsKCl_clean.txt');
Data_phi = readmatrix('OsmoticCoefficientsKCl_clean.txt');
Data_L = readmatrix('MixingEnthalpyKCl.txt');
Data_J = readmatrix('ApparentMolarHeatCapacityKCl_clean.txt');
%----------matrixes of different parameters
%activity
gamma = Data_gamma(:,2:7);
b_gamma = Data_gamma(:,1);
T_data_gamma = [283.15, 288.15, 298.15, 313.15, 323.15, 333.15];
%osmotic coefficient
phi = Data_phi(:,2:8);
b_phi = Data_phi(:,1);
T_data_phi = [273.15, 283.15, 288.15, 298.15, 313.15, 323.15, 333.15];
%enthalpy
Lphi = Data_L(:,2:5)*4.184; %conversion cal to J
b_L = Data_L(:,1);
T_data_L = [298.15, 313.15, 333.15, 353.15];
%heat capacity
Jphi = Data_J(:,2:6);
b_J = Data_J(:,1);
T_data_J = [268.15, 298.15, 318.15, 338.15, 358.15];

%-----------Data at 25°C for fitting
phi_25 = phi(:,4);
L_25 = Lphi(:,1);
J_25 = Jphi(:,2);

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
       B1 =        61.8 ;% (58.07, 65.52)
       C1 =      -7.998 ;% (-12.3, -3.698)
       D1 =     -0.1125 ;% (-1.883, 1.658)
       E1 =     0.04836 ;% (-0.1858, 0.2825)
       Q1 =       80.02 ;% (73.16, 86.88)

%Goodness of fit:
  %SSE: 86.94
  %R-square: 1
  %Adjusted R-square: 1
  %RMSE: 2.331

  %-----Heat capacity
  %General model:
     %f(b) = ReducedMatrixFitHeatCapacity(b, Q2, B2, C2, D2, E2)
%Coefficients (with 95% confidence bounds):
       B2 =     -0.5265 ;% (-4.857, 3.804)
       C2 =     -0.4449 ;% (-5.116, 4.226)
       D2 =      0.2421 ;% (-1.527, 2.011)
       E2 =    -0.03133 ;% (-0.2466, 0.184)
       Q2 =      -2.028 ;% (-10.28, 6.226)

%Goodness of fit:
  %SSE: 8.079
  %R-square: 0.999
  %Adjusted R-square: 0.9981
  %RMSE: 1.421
  
  %General model:
     %f(b) = ReducedMatrixFitOsmoticCoefficient(b, Q0, B0, C0, D0, E0)
%Coefficients (with 95% confidence bounds):
       B0 =      -11.97 ;% (-13.25, -10.7)
       C0 =     -0.9969 ;% (-2.143, 0.1494)
       D0 =      0.2107 ;% (-0.1782, 0.5996)
       E0 =   -0.009914 ;% (-0.05489, 0.03506)
       Q0 =      -73.02 ;% (-76.27, -69.76)

%Goodness of fit:
  %SSE: 2.32e-06
  %R-square: 0.9998
  %Adjusted R-square: 0.9998
  %RMSE: 0.0003176

%---------MOLALITY COEFFICIENTS
 
%---------- Molality Matrix for osmotic coefficient as in Clarke
b0_phi = -b_phi.^(1/2)./(1+1.2*b_phi.^(1/2));
b1_phi = b_phi.*exp(-2*b_phi.^(1/2));
b2_phi = b_phi;
b3_phi = b_phi.^2;
b4_phi = b_phi.^3;
b5_phi = b_phi.^4;

B_phi = [b0_phi, b1_phi, b2_phi, b3_phi, b4_phi, b5_phi];

 %---------- Molality Matrix for activity coefficient as in Clarke
b0_gamma = -(b_gamma.^(1/2)./(1+1.2*b_gamma.^(1/2))+2/1.2*log(1+1.2*b_gamma.^(1/2)));
b1_gamma = 1/2*(1-(1+2*b_gamma.^(1/2)-2*b_gamma).*exp(-2*b_gamma.^(1/2)));
b2_gamma = 2*b_gamma;
b3_gamma = 3/2*b_gamma.^2;
b4_gamma = 4/3*b_gamma.^3;
b5_gamma = 5/4*b_gamma.^4;

B_gamma = [b0_gamma, b1_gamma, b2_gamma, b3_gamma, b4_gamma, b5_gamma];

 %---------- Molality Matrix for Lphi
b0_L = -4/1.2*log(1+1.2*b_L.^(1/2));
b1_L = 1-(1+2*b_L.^(1/2)).*exp(-2*b_L.^(1/2));
b2_L = 2*b_L;
b3_L = b_L.^2;
b4_L = 2/3*b_L.^3;
b5_L = 1/2*b_L.^4;

B_L = [b0_L, b1_L, b2_L, b3_L, b4_L, b5_L];

 %---------- Molality Matrix for Lphi
b0_J = -4/1.2*log(1+1.2*b_J.^(1/2));
b1_J = 1-(1+2*b_J.^(1/2)).*exp(-2*b_J.^(1/2));
b2_J = 2*b_J;
b3_J = b_J.^2;
b4_J = 2/3*b_J.^3;
b5_J = 1/2*b_J.^4;

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

%-----Temperature matrix for gamma %---
x_gamma = (T_data_gamma-Theta)./Theta;

T0_gamma = -1/Theta.*ones(size(T_data_gamma));
T1_gamma = (1/Theta.*x_gamma./(x_gamma+1)); %Wolfram alpha
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


%------ Plotting
%Phi as a function of b and T
figure(1)
[x1, y1, z1] = prepareSurfaceData( T_data_phi, b_phi, phi_calc_reduced);
[x2, y2, z2] = prepareSurfaceData( T_data_phi, b_phi, phi);
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


gamma_KCl_calc = horzcat(b_gamma, gamma_calc_reduced);
writematrix(gamma_KCl_calc,'gammaKClCalc.csv')
gamma_KCl_exp = horzcat(b_gamma, gamma);
writematrix(gamma_KCl_exp,'gammaKClExp.csv') 

phi_KCl_calc = horzcat(b_phi, phi_calc_reduced);
writematrix(phi_KCl_calc,'phiKClCalc.csv')
phi_KCl_exp = horzcat(b_phi, phi);
writematrix(phi_KCl_exp,'phiKClExp.csv') 

%residuals
residuals_phi_KCl = (phi_calc_reduced-phi)./phi*100;
residuals_gamma_KCl = (gamma_calc_reduced-gamma)./gamma*100;
residuals_phi= horzcat(b_phi, residuals_phi_KCl);
residuals_gamma= horzcat(b_gamma, residuals_gamma_KCl);
writematrix(residuals_phi,'phiKClresiduals.csv') 
writematrix(residuals_gamma,'gammaKClresiduals.csv')


%{

%-------SECTION 2: FREEZING POINT

%-----Scatchard
%     f(x) = 273.15+a1*x^(1/2)+a2*x+a3*x^(3/2)+a4*x^2
%Coefficients (with 95% confidence bounds):
       TfusH2O =       273.15 ;% (fixed at bound)
       t1 =    -0.02087 ;% (-0.0277, -0.01403)
       t2 =      -3.463 ;% (-3.498, -3.428)
       t3 =      0.3278 ;% (0.2715, 0.3842)
       t4 =     -0.1072 ;% (-0.1352, -0.0792)
syms x
T_f(x) = t1*x^(1/2)+t2*x+t3*x^(3/2)+t4*x^2+ TfusH2O;

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
Scatch

%----43x1 matrixes where only the value at (b(Tf),Tf) is selected
phi_fp = zeros(size(b_phi));
for i=1:m
    phi_fp(i,1) = phi_Tf(i,i);
end

figure(7)
plot(b_phi, phi_fp, b_phi, phi(:,5))

phi_KCl_fp = horzcat(b_phi, phi_fp, phi(:,5));
writematrix(phi_KCl_fp ,'phiKClFP.csv')

%}