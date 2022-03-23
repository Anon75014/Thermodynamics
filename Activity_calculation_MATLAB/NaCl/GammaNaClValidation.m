%Calculation of reduced model at the average molality of Clarke, Partanen
%and Gibbard

reference_data_gamma = readmatrix('ReferenceActivityNaCl_clean.txt');
gamma_ref = reference_data_gamma(:,2:7);
b_averaged = reference_data_gamma(:,1);
T_gamma_ref = [273.15, 283.15, 288.15, 298.15, 313.15, 333.15];

 %---------- Molality Matrix for activity coefficient as in Clarke
b0_gamma = -(b_averaged.^(1/2)./(1+1.2*b_averaged.^(1/2))+2/1.2*log(1+1.2*b_averaged.^(1/2)));
b1_gamma = 1/2*(1-(1+2*b_averaged.^(1/2)-2*b_averaged).*exp(-2*b_averaged.^(1/2)));
b2_gamma = 2*b_averaged;
b3_gamma = 3/2*b_averaged.^2;
b4_gamma = 4/3*b_averaged.^3;
b5_gamma = 5/4*b_averaged.^4;

B_gamma = [b0_gamma, b1_gamma, b2_gamma, b3_gamma, b4_gamma, b5_gamma];

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

%-----Temperature matrix for gamma %---Problem with - signs
x_gamma = (T_gamma_ref-Theta)./Theta;

T0_gamma = -1/Theta.*ones(size(T_gamma_ref));
T1_gamma = (1/Theta.*x_gamma./(x_gamma+1));
T2_gamma = (x_gamma.^2.*vpa(S2(x_gamma)));
T3_gamma = (1/2*Theta.*x_gamma.^3.*vpa(S3(x_gamma)));
T4_gamma = (1/6*Theta^2.*x_gamma.^4.*vpa(S4(x_gamma)));

T_gamma = [T0_gamma; T1_gamma; T2_gamma; T3_gamma; T4_gamma];

%----------- COEFFICIENTS FROM CLARKE
%----------Ak coefficients from Clarke 1980 (nu_gk is not exactly u_k)
A0 = -0.391940*298.15;
A1 = 0.198653*298.15;
A2 = 0.772533;
A3 = 2/298.15*1.68848;
A4 = 6/(298.15)^2*1.99113;
%-----Reduced matrix fit
%------- Lphi at 25°C
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
 
  
  %------- Jphi at 25°C
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
 
 %--------phi at 25°C
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
gamma_calc = double(exp(B_gamma*V_reduced*T_gamma));


gamma_NaCl_calc = horzcat(b_averaged, gamma_calc);
writematrix(gamma_NaCl_calc,'gammaNaClCalc.csv')
gamma_NaCl_exp = horzcat(b_averaged, gamma_ref);
writematrix(gamma_NaCl_exp,'gammaNaClExp.csv') 

%residuals
residuals_gamma_NaCl = (gamma_calc-gamma_ref)./gamma_ref*100;
residuals_gamma= horzcat(b_averaged, residuals_gamma_NaCl);
writematrix(residuals_gamma,'gammaNaClresiduals.csv')