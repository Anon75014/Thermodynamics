%Visualisation of Cpphi as calculated from Clarke and agreement with Literature
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
 
 %---------- Molality Matrix for Lphi, Cpphi and Jphi as in Clarke
 b1 = [0:0.02:0.1]';
 b2= [0.1:0.05:0.4]';
 b3= [0.4:0.2:6]';
 b=vertcat(b1,b2,b3);
b0_LCp = -4/1.2*log(1+1.2*b.^(1/2));
b1_LCp = 1-(1+2*b.^(1/2)).*exp(-2*b.^(1/2));
b2_LCp = 2*b;
b3_LCp = b.^2;
b4_LCp = 2/3*b.^3;
b5_LCp = 1/2*b.^4;

B_LCp = [b0_LCp, b1_LCp, b2_LCp, b3_LCp, b4_LCp, b5_LCp];

%-----Temperature matrix for Jphi
  Theta = 298.15;
  T = 273.15:1:338.15;
T0_J = zeros(size(T));
T1_J = zeros(size(T));
T2_J = ones(size(T));
T3_J = (T-Theta);
T4_J = (1/2.*(T-Theta).^2);
T5_J = (1/6.*(T-Theta).^3);
T6_J = (1/24.*(T-Theta).^4);

T_J = [T0_J; T1_J; T2_J; T3_J; T4_J; T5_J; T6_J];

R=8.314;
J_phi = R*B_LCp*V*T_J;

%Experimental data from Messikomer for comparison
Data = readmatrix('ApparentHeatCapacityCpphiNaCl.txt');
J_phi_Tanner = -Data(:,6:9)*4.184;
b_Tanner = Data(:,1);
T_Tanner = 273.15 + [5, 25, 45, 65];

%Plot
figure(1)
h1 = surf( T, b, J_phi);
alpha 0.7
hold on
[x, y, z] = prepareSurfaceData( T_Tanner, b_Tanner, J_phi_Tanner);
h2 = scatter3(x, y, z,  'Marker','s',...
     'MarkerEdgeColor','k',...
     'MarkerFaceColor','y');
legend([h1, h2],'Calculated', 'Messikomer', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'latex' );
ylabel( 'b (mol/kg)', 'Interpreter', 'latex' );
zlabel( '$C_p^{\phi}$ (J/mol/K)', 'Interpreter', 'latex' );
view( -143.2, 26.2 );
set(gca,'FontSize',14)
hold off


%residuals

b0_LCp_Tanner = -4/1.2*log(1+1.2*b_Tanner.^(1/2));
b1_LCp_Tanner = 1-(1+2*b_Tanner.^(1/2)).*exp(-2*b_Tanner.^(1/2));
b2_LCp_Tanner = 2*b_Tanner;
b3_LCp_Tanner = b_Tanner.^2;
b4_LCp_Tanner = 2/3*b_Tanner.^3;
b5_LCp_Tanner = 1/2*b_Tanner.^4;
B_LCp_Tanner = [b0_LCp_Tanner, b1_LCp_Tanner, b2_LCp_Tanner, b3_LCp_Tanner, b4_LCp_Tanner, b5_LCp_Tanner];

T0_J_Tanner = zeros(size(T_Tanner));
T1_J_Tanner = zeros(size(T_Tanner));
T2_J_Tanner = ones(size(T_Tanner));
T3_J_Tanner = (T_Tanner-Theta);
T4_J_Tanner = (1/2.*(T_Tanner-Theta).^2);
T5_J_Tanner = (1/6.*(T_Tanner-Theta).^3);
T6_J_Tanner = (1/24.*(T_Tanner-Theta).^4);
T_J_Tanner = [T0_J_Tanner; T1_J_Tanner; T2_J_Tanner; T3_J_Tanner; T4_J_Tanner; T5_J_Tanner; T6_J_Tanner];

J_phi_Tanner_calc = R*B_LCp_Tanner*V*T_J_Tanner;

%Plot residuals
%{
figure(2)
h3 = stem3( T_Tanner, b_Tanner, (J_phi_Tanner-J_phi_Tanner_calc)./J_phi_Tanner, 'Marker','s',...
     'MarkerEdgeColor','b',...
     'MarkerFaceColor','g');
legend([h3],'Residuals', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'none' );
ylabel( 'b (mol/kg)', 'Interpreter', 'none' );
zlabel( '%', 'Interpreter', 'none' );
view( -40, 40 );
set(gca,'FontSize',14)
%}