reference_data_gamma = readmatrix('ReferenceActivityNaCl_clean.txt');
gamma_ref = reference_data_gamma(:,2:7);
b_averaged = reference_data_gamma(:,1);
T_gamma_ref = [273.15, 283.15, 288.15, 298.15, 313.15, 333.15];

%Visualisation of phi as calculated from Clarke and agreement with Literature
%data
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
 
 %---------- Molality Matrix for gamma
 b1 = [0:0.02:0.1]';
 b2= [0.1:0.05:0.4]';
 b3= [0.4:0.2:6]';
 b=vertcat(b1,b2,b3);
 
 b0_gamma = -(b.^(1/2)./(1+1.2*b.^(1/2))+2/1.2*log(1+1.2*b.^(1/2)));
b1_gamma = 1/2*(1-(1+2*b.^(1/2)-2*b).*exp(-2*b.^(1/2)));
b2_gamma = 2*b;
b3_gamma = 3/2*b.^2;
b4_gamma = 4/3*b.^3;
b5_gamma = 5/4*b.^4;

B_gamma = [b0_gamma, b1_gamma, b2_gamma, b3_gamma, b4_gamma, b5_gamma];

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);
S5(x) = symsum((n/(4+n)*(-x)^(n-1)),n,0,Inf);
S6(x) = symsum((n/(5+n)*(-x)^(n-1)),n,0,Inf);
S7(x) = symsum((n/(6+n)*(-x)^(n-1)),n,0,Inf);

Theta = 298.15;
T = 273.15:1:338.15;
x_gamma = (T-Theta)./Theta;

T0_gamma = -1/Theta.*ones(size(T));
T1_gamma = (1/Theta.*x_gamma./(x_gamma+1));
T2_gamma = (x_gamma.^2.*vpa(S2(x_gamma)));
T3_gamma = (1/2*Theta.*x_gamma.^3.*vpa(S3(x_gamma)));
T4_gamma = (1/6*Theta^2.*x_gamma.^4.*vpa(S4(x_gamma)));
T5_gamma = (1/24*Theta^3.*x_gamma.^5.*vpa(S5(x_gamma)));
T6_gamma = (1/120*Theta^4.*x_gamma.^6.*vpa(S6(x_gamma)));

T_gamma = [T0_gamma; T1_gamma; T2_gamma; T3_gamma; T4_gamma; T5_gamma; T6_gamma];


gamma_calc = double(exp(B_gamma*V*T_gamma));


%Plot
figure(1)
h1 = surf(T, b, gamma_calc);
alpha 0.9
hold on
[x, y, z] = prepareSurfaceData( T_gamma_ref, b_averaged, gamma_ref);
h2 = scatter3(x, y, z, 'Marker','s',...
     'MarkerEdgeColor','k',...
     'MarkerFaceColor','y');
legend([h1,h2],'Calculated','Validation set', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'latex' );
ylabel( 'b (mol/kg)', 'Interpreter', 'latex' );
zlabel( '$\gamma_{NaCl}$', 'Interpreter', 'latex' );
ylim([0 5])
set(gca,'FontSize',14)
view( -112, 37 );
