%----------Reading the experimental data
Data = readmatrix('Data_NaCl.txt');

%---------- Ideal gas constant
R = 8.314; %J/K/mol
%----------matrixes of different parameters
phi_Clarke = Data(:,2:7);
gamma_Clarke = Data(:,8:13);
Lphi_Clarke = Data(:,14:18);
Cpphi_Clarke = Data(:,19:23);
Jphi_Clarke = Data(:,24:28);
%------------Temperature matrixes
T_Clarke = [273.15, 283.15, 298.15, 313.15, 333.15];
T_phi_Clarke = [273.15, 278.15, 288.15, 298.15, 313.15, 333.15];
T_gamma_Clarke = [273.15, 283.15, 288.15, 298.15, 313.15, 333.15];
%----------Column vectors for individual temperatures
b_Clarke = Data(:,1);
phi_0_Clarke = phi_Clarke(:,1);
phi_5_Clarke = phi_Clarke(:,2);
phi_15_Clarke = phi_Clarke(:,3);
phi_25_Clarke = phi_Clarke(:,4);
phi_40_Clarke = phi_Clarke(:,5);
phi_60_Clarke = phi_Clarke(:,6);
gamma_0_Clarke = gamma_Clarke(:,1);
gamma_10_Clarke = gamma_Clarke(:,2);
gamma_15_Clarke = gamma_Clarke(:,3);
gamma_25_Clarke = gamma_Clarke(:,4);
gamma_40_Clarke = gamma_Clarke(:,5);
gamma_60_Clarke = gamma_Clarke(:,6);
Lphi_0_Clarke = Lphi_Clarke(:,1);
Lphi_10_Clarke = Lphi_Clarke(:,2);
Lphi_25_Clarke = Lphi_Clarke(:,3);
Lphi_40_Clarke = Lphi_Clarke(:,4);
Lphi_60_Clarke = Lphi_Clarke(:,5);
Cpphi_0_Clarke = Cpphi_Clarke(:,1);
Cpphi_10_Clarke = Cpphi_Clarke(:,2);
Cpphi_25_Clarke = Cpphi_Clarke(:,3);
Cpphi_40_Clarke = Cpphi_Clarke(:,4);
Cpphi_60_Clarke = Cpphi_Clarke(:,5);
Jphi_0_Clarke = Jphi_Clarke(:,1);
Jphi_10_Clarke = Jphi_Clarke(:,2);
Jphi_25_Clarke = Jphi_Clarke(:,3);
Jphi_40_Clarke = Jphi_Clarke(:,4);
Jphi_60_Clarke = Jphi_Clarke(:,5);

%----------- COEFFICIENTS FROM CLARKE
%----------Ak coefficients from Clarke 1980 (nu_gk is not exactly u_k)
A0 = -0.391940*298.15;
A1 = 0.198653*298.15;
A2 = 0.772533;
A3 = 2/298.15*1.68848;
A4 = 6/(298.15)^2*1.99113;

%Other coefficients from Clarke 1985
Q0 = -8.290627e1;
 Q1 = 3.680527e1;
 Q2 =-4.135823e-1;
 Q3 =9.928923e-2;
 Q4 =-6.966493e-3;
 Q5 =3.545976e-4;
 Q6 = -7.375379e-6;
 
 B0 = -2.250462e1;
 B1 = 8.178707e1;
 B2 = -1.722838;
 B3 = 3.655402e-2;
 B4 = -8.698438e-4;
 B5 = 1.077619e-5;
 
 C0 = -7.849751e-1;
 C1 = -2.045013e1;
 C2 = 4.449950e-1;
 C3 = -9.817574e-3;
 C4 = 2.042558e-4;
 C5 = -2.065926e-6;
 
 D0 = 9.433827e-02;
 D1 = 2.153966;
 D2 = -5.631073e-2;
 D3= 7.473613e-4;
 D4= -9.246901e-6;
 
 E0 = -4.764496e-3;
 E1 = -1.245283e-1;
 E2 = 2.492125e-3;
 
 %------Virial Matric from Clarke
 V = [A0, A1, A2, A3, A4, 0, 0; Q0, Q1, Q2, Q3, Q4, Q5, Q6; B0, B1, B2, B3, B4, B5, 0; C0, C1, C2, C3, C4, C5, 0; D0, D1, D2, D3, D4, 0, 0; E0, E1, E2, 0, 0, 0, 0];

%---------MOLALITY COEFFICIENTS
 
%---------- Molality Matrix for osmotic coefficient as in Clarke
b0_phi = -b_Clarke.^(1/2)./(1+1.2*b_Clarke.^(1/2));
b1_phi = b_Clarke.*exp(-2*b_Clarke.^(1/2));
b2_phi = b_Clarke;
b3_phi = b_Clarke.^2;
b4_phi = b_Clarke.^3;
b5_phi = b_Clarke.^4;

B_phi = [b0_phi, b1_phi, b2_phi, b3_phi, b4_phi, b5_phi];

 %---------- Molality Matrix for activity coefficient as in Clarke
b0_gamma = -(b_Clarke.^(1/2)./(1+1.2*b_Clarke.^(1/2))+2/1.2*log(1+1.2*b_Clarke.^(1/2)));
b1_gamma = 1/2*(1-(1+2*b_Clarke.^(1/2)-2*b_Clarke).*exp(-2*b_Clarke.^(1/2)));
b2_gamma = 2*b_Clarke;
b3_gamma = 3/2*b_Clarke.^2;
b4_gamma = 4/3*b_Clarke.^3;
b5_gamma = 5/4*b_Clarke.^4;

B_gamma = [b0_gamma, b1_gamma, b2_gamma, b3_gamma, b4_gamma, b5_gamma];

 %---------- Molality Matrix for Lphi, Cpphi and Jphi as in Clarke
b0_LCp = -4/1.2*log(1+1.2*b_Clarke.^(1/2));
b1_LCp = 1-(1+2*b_Clarke.^(1/2)).*exp(-2*b_Clarke.^(1/2));
b2_LCp = 2*b_Clarke;
b3_LCp = b_Clarke.^2;
b4_LCp = 2/3*b_Clarke.^3;
b5_LCp = 1/2*b_Clarke.^4;

B_LCp = [b0_LCp, b1_LCp, b2_LCp, b3_LCp, b4_LCp, b5_LCp];

%---------TEMPERATURE COEFFICIENTS
%---Reference temperature
Theta = 298.15;

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);
S5(x) = symsum((n/(4+n)*(-x)^(n-1)),n,0,Inf);
S6(x) = symsum((n/(5+n)*(-x)^(n-1)),n,0,Inf);
S7(x) = symsum((n/(6+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for Lphi
T0_L = zeros(size(T_Clarke));
T1_L = -ones(size(T_Clarke));
T2_L = -(T_Clarke-Theta);
T3_L = -(1/2.*(T_Clarke-Theta).^2);
T4_L = -(1/6.*(T_Clarke-Theta).^3);
T5_L = -(1/24.*(T_Clarke-Theta).^4);
T6_L = -(1/120.*(T_Clarke-Theta).^5);

T_L = [T0_L; T1_L; T2_L; T3_L; T4_L; T5_L; T6_L];

%-----Temperature matrix for Jphi
T0_J = zeros(size(T_Clarke));
T1_J = zeros(size(T_Clarke));
T2_J = -ones(size(T_Clarke));
T3_J = -(T_Clarke-Theta);
T4_J = -(1/2.*(T_Clarke-Theta).^2);
T5_J = -(1/6.*(T_Clarke-Theta).^3);
T6_J = -(1/24.*(T_Clarke-Theta).^4);

T_J = [T0_J; T1_J; T2_J; T3_J; T4_J; T5_J; T6_J];

%-----Temperature matrix for gamma 
x_gamma = (T_gamma_Clarke-Theta)./Theta;

T0_gamma = -1/Theta.*ones(size(T_gamma_Clarke));
T1_gamma = (1/Theta.*x_gamma./(x_gamma+1));
T2_gamma = (x_gamma.^2.*vpa(S2(x_gamma)));
T3_gamma = (1/2*Theta.*x_gamma.^3.*vpa(S3(x_gamma)));
T4_gamma = (1/6*Theta^2.*x_gamma.^4.*vpa(S4(x_gamma)));
T5_gamma = (1/24*Theta^3.*x_gamma.^5.*vpa(S5(x_gamma)));
T6_gamma = (1/120*Theta^4.*x_gamma.^6.*vpa(S6(x_gamma)));

T_gamma = [T0_gamma; T1_gamma; T2_gamma; T3_gamma; T4_gamma; T5_gamma; T6_gamma];

%-----Temperature matrix for phi
x_phi = (T_phi_Clarke-Theta)./Theta;

T0_phi = -1/Theta.*ones(size(T_phi_Clarke));
T1_phi = (1/Theta.*x_phi./(x_phi+1));
T2_phi = (x_phi.^2.*vpa(S2(x_phi)));
T3_phi = (1/2*Theta.*x_phi.^3.*vpa(S3(x_phi)));
T4_phi = (1/6*Theta^2.*x_phi.^4.*vpa(S4(x_phi)));
T5_phi = (1/24*Theta^3.*x_phi.^5.*vpa(S5(x_phi)));
T6_phi = (1/120*Theta^4.*x_phi.^6.*vpa(S6(x_phi)));

T_phi = [T0_phi; T1_phi; T2_phi; T3_phi; T4_phi; T5_phi; T6_phi];


%-----CORRECTNESS CHECKS

%check of  phi at 25C
figure(1)
T_ref = [-1/Theta; 0; 0; 0; 0; 0; 0];
phi_ref = 1+B_phi*V*T_ref;
plot(b_Clarke,phi_ref);
hold on
scatter(b_Clarke,phi_25_Clarke);
hold off
%WORKS!

%check of gamma at 25C
figure(2)
gamma_ref = B_gamma*V*T_ref;
gamma_not_log = exp(gamma_ref);
plot(b_Clarke,gamma_not_log);
hold on
scatter(b_Clarke,gamma_25_Clarke);
hold off
%WORKS!

%check of Lphi at 25C
%must change position of Tref
figure(3)
T_ref = [0; -1; 0; 0; 0; 0; 0];
Lphi_ref = R*B_LCp*V*T_ref;
plot(b_Clarke,Lphi_ref);
hold on
scatter(b_Clarke,Lphi_25_Clarke);
hold off
%CHECK!

%check of Jphi at 25C
%must change position of Tref again
figure(4)
T_ref = [0; 0; 1; 0; 0; 0; 0];
Jphi_ref = R*B_LCp*V*T_ref;
plot(b_Clarke,Jphi_ref);
hold on
scatter(b_Clarke,Jphi_25_Clarke);
hold off
%CHECK!

%-------Check of Lphi's surface
L_phi_calc = R*B_LCp*V*T_L;
%CHECK

%-------Check of Jphi's surface
J_phi_calc = R*B_LCp*V*T_J;
%CHECK

%-------Check of phi's surface
phi_calc = 1+B_phi*V*T_phi;

%-------Check of gamma's surface
gamma_calc = B_gamma*V*T_gamma;
gamma_not_log = exp(gamma_calc);

%return; % Exit script or return to calling routine.


%------ Plotting
%Phi as a function of b and T
figure(5)
[x1, y1, z1] = prepareSurfaceData( T_gamma_Clarke, b_Clarke, gamma_not_log);
[x2, y2, z2] = prepareSurfaceData( T_gamma_Clarke, b_Clarke, gamma_Clarke);
h1 = scatter3(x1, y1, z1);
hold on
h2 = scatter3(x2, y2, z2);
legend([h1,h2],'Calculated','Clarke data', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'none' );
ylabel( 'b (mol/kg)', 'Interpreter', 'none' );
zlabel( 'phi', 'Interpreter', 'none' );
view( -53.2, 26.2 );
hold off

%-----Reduced matrix fit
%------- Lphi at 25??C
%General model:
     %f(b) = ReducedMatrixFitEnthalpy(b, Q1, B1, C1, D1, E1)
%Coefficients (with 95% confidence bounds):
       B1 =       81.84 ;% (81.76, 81.92)
       C1 =       -20.5 ;% (-20.57, -20.43)
       D1 =        2.17 ;% (2.148, 2.193)
       E1 =     -0.1262 ;% (-0.1284, -0.1239)
       Q1 =        36.7 ;% (36.53, 36.86)

%Goodness of fit:
  %SSE: 0.3733
  %R-square: 1
  %Adjusted R-square: 1
  %RMSE: 0.09911
 
  
  %------- Jphi at 25??C
%General model:
     %f(b) = ReducedMatrixFitHeatCapacity(b, Q2, B2, C2, D2, E2)
%Coefficients (with 95% confidence bounds):
       B2 =      -1.723 ;% (-1.723, -1.722)
       C2 =      0.4449 ;% (0.4443, 0.4455)
       D2 =    -0.05627 ;% (-0.05645, -0.05609)
       E2 =    0.002488 ;% (0.00247, 0.002506)
       Q2 =     -0.4139 ;% (-0.4152, -0.4125)

%Goodness of fit:
 % SSE: 2.409e-05
 % R-square: 1
 % Adjusted R-square: 1
 % RMSE: 0.0007962
 
 %--------phi at 25??C
 %General model:
     %f(b) = ReducedMatrixFitOsmoticCoefficient(b, Q0, B0, C0, D0, E0)
%Coefficients (with 95% confidence bounds):
       B0 =      -22.51 ;% (-22.51, -22.5)
       C0 =     -0.7836 ;% (-0.7875, -0.7797)
       D0 =     0.09402 ;% (0.0931, 0.09494)
       E0 =   -0.004742 ;% (-0.004814, -0.004669)
       Q0 =       -82.9 ;% (-82.93, -82.88)

%Goodness of fit:
  %SSE: 3.271e-10
  %R-square: 1
  %Adjusted R-square: 1
  %RMSE: 2.934e-06
  
  
 %----Calculations with reduced matrixes
 
V_reduced = [A0, A1, A2, A3, A4; Q0, Q1, Q2, 0, 0; B0, B1, B2, 0, 0; C0, C1, C2, 0, 0; D0, D1, D2, 0, 0; E0, E1, E2, 0, 0];
T_gamma_reduced = [T0_gamma; T1_gamma; T2_gamma; T3_gamma; T4_gamma];
T_phi_reduced = [T0_phi; T1_phi; T2_phi; T3_phi; T4_phi];
T_L_reduced = [T0_L; T1_L; T2_L; T3_L; T4_L];
T_J_reduced = [T0_J; T1_J; T2_J; T3_J; T4_J];

phi_calc_reduced = double(1+B_phi*V_reduced*T_phi_reduced);
gamma_calc_reduced = double(exp(B_gamma*V_reduced*T_gamma_reduced));


%------ Plotting
%Phi as a function of b and T
figure(6)
[x1, y1, z1] = prepareSurfaceData( T_gamma_Clarke, b_Clarke, gamma_calc_reduced);
[x2, y2, z2] = prepareSurfaceData( T_gamma_Clarke, b_Clarke, gamma_Clarke);
h1 = scatter3(x1, y1, z1);
hold on
h2 = scatter3(x2, y2, z2);
legend([h1,h2],'Calculated','Clarke data', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'none' );
ylabel( 'b (mol/kg)', 'Interpreter', 'none' );
zlabel( 'phi', 'Interpreter', 'none' );
view( -53.2, 26.2 );
hold off
 
%{


%-------SECTION 2: FREEZING POINT

%-----Fit on Bodnar Data
%Linear model Poly4:
     %f(x) = t1*x^(1/2)+t2*x+t3*x^(3/2)+t4*x^2+ TfusH2O
%Coefficients (with 95% confidence bounds):
       t1 =      0.2874 ;% (0.2636, 0.3112)
       t2 =      -4.253 ;% (-4.304, -4.203)
       t3 =        1.04 ;% (1.006, 1.074)
       t4 =     -0.4541 ;% (-0.4615, -0.4468)
       TfusH2O = 273.15;
syms x
T_f(x) = t1*x^(1/2)+t2*x+t3*x^(3/2)+t4*x^2+ TfusH2O;

%Freezing Point of NaCl evaluated at Clarke molalities
T_fus = zeros(1,size(b_Clarke,1));
m = size(b_Clarke,1);
for i=1:m
    T_fus(1,i) = T_f(b_Clarke(i,1));
end


%estimation of properties at the freezing point: each property needs to be
%estimated at a given b and only the corresponding T_f

%-----Temperature matrix at freezing point for Lphi
T0_L_fus = zeros(size(T_fus));
T1_L_fus = -ones(size(T_fus));
T2_L_fus = -(T_fus-Theta);
T3_L_fus = -(1/2.*(T_fus-Theta).^2);
T4_L_fus = -(1/6.*(T_fus-Theta).^3);

T_L_fus = [T0_L_fus; T1_L_fus; T2_L_fus; T3_L_fus; T4_L_fus];

%-----Temperature matrix at freezing point for Jphi
T0_J_fus = zeros(size(T_fus));
T1_J_fus = zeros(size(T_fus));
T2_J_fus = -ones(size(T_fus));
T3_J_fus = -(T_fus-Theta);
T4_J_fus = -(1/2.*(T_fus-Theta).^2);

T_J_fus = [T0_J_fus; T1_J_fus; T2_J_fus; T3_J_fus; T4_J_fus];

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for gamma %---Problem with - signs
x_gamma_fus = (T_fus-Theta)./Theta;

T0_gamma_fus = -1/Theta.*ones(size(T_fus));
T1_gamma_fus = (1/Theta.*x_gamma_fus./(x_gamma_fus+1));
T2_gamma_fus = (x_gamma_fus.^2.*vpa(S2(x_gamma_fus)));
T3_gamma_fus = (1/2*Theta.*x_gamma_fus.^3.*vpa(S3(x_gamma_fus)));
T4_gamma_fus = (1/6*Theta^2.*x_gamma_fus.^4.*vpa(S4(x_gamma_fus)));


T_gamma_fus = [T0_gamma_fus; T1_gamma_fus; T2_gamma_fus; T3_gamma_fus; T4_gamma_fus];

%-----Temperature matrix for phi 
x_phi_fus = (T_fus-Theta)./Theta;

T0_phi_fus = -1/Theta.*ones(size(T_fus));
T1_phi_fus = (1/Theta.*x_phi_fus./(x_phi_fus+1));
T2_phi_fus = (x_phi_fus.^2.*vpa(S2(x_phi_fus)));
T3_phi_fus = (1/2*Theta.*x_phi_fus.^3.*vpa(S3(x_phi_fus)));
T4_phi_fus = (1/6*Theta^2.*x_phi_fus.^4.*vpa(S4(x_phi_fus)));

T_phi_fus = [T0_phi_fus; T1_phi_fus; T2_phi_fus; T3_phi_fus; T4_phi_fus];


%----43x43 matrixes where all parameters have been evaluated on all
%molalities at all temperatures
L_Tf = B_LCp*V_reduced*T_L_fus;
J_Tf = B_LCp*V_reduced*T_J_fus;
phi_Tf = 1+B_phi*V_reduced*T_phi_fus;
gamma_Tf = exp(B_gamma*V_reduced*T_gamma_fus);

%----43x1 matrixes where only the value at (b(Tf),Tf) is selected
L_fp = zeros(size(b_Clarke));
J_fp = zeros(size(b_Clarke));
phi_fp = zeros(size(b_Clarke));
gamma_fp = zeros(size(b_Clarke));
for i=1:m
    L_fp(i,1) = L_Tf(i,i);
    J_fp(i,1) = J_Tf(i,i);
    phi_fp(i,1) = phi_Tf(i,i);
    gamma_fp(i,1) = gamma_Tf(i,i);
end

L_NaCl_fp = horzcat(b_Clarke, L_fp);
J_NaCl_fp = horzcat(b_Clarke, J_fp);
writematrix(L_NaCl_fp ,'LNaClFP.csv')
writematrix(J_NaCl_fp ,'CpNaClFP.csv')
%}
%{
figure(7)
plot(b_Clarke, phi_fp, b_Clarke, phi_25_Clarke)

phi_NaCl_fp = horzcat(b_Clarke, phi_fp, phi_25_Clarke);
writematrix(phi_NaCl_fp ,'phiNaClFP.csv')
%}


gamma_NaCl_calc = horzcat(b_Clarke, gamma_calc_reduced);
writematrix(gamma_NaCl_calc,'gammaNaClCalc.csv')
gamma_NaCl_exp = horzcat(b_Clarke, gamma_Clarke);
writematrix(gamma_NaCl_exp,'gammaNaClExp.csv') 

phi_NaCl_calc = horzcat(b_Clarke, phi_calc_reduced);
writematrix(phi_NaCl_calc,'phiNaClCalc.csv')
phi_NaCl_exp = horzcat(b_Clarke, phi_Clarke);
writematrix(phi_NaCl_exp,'phiNaClExp.csv') 

%residuals
residuals_phi_NaCl = (phi_calc_reduced-phi_Clarke)./phi_Clarke*100;
residuals_gamma_NaCl = (gamma_calc_reduced-gamma_Clarke)./gamma_Clarke*100;
residuals_phi= horzcat(b_Clarke, residuals_phi_NaCl);
residuals_gamma= horzcat(b_Clarke, residuals_gamma_NaCl);
writematrix(residuals_phi,'phiNaClresiduals.csv') 
writematrix(residuals_gamma,'gammaNaClresiduals.csv')


