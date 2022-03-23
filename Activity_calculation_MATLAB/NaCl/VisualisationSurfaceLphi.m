%Visualisation of Lphi as calculated from Clarke and agreement with Literature
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

%-----Temperature matrix for Lphi
  Theta = 298.15;
  T = 298.15:1:333.15;
T0_L = zeros(size(T));
T1_L = ones(size(T));
T2_L = (T-Theta);
T3_L = (1/2.*(T-Theta).^2);
T4_L = (1/6.*(T-Theta).^3);
T5_L = (1/24.*(T-Theta).^4);
T6_L = (1/120.*(T-Theta).^5);

T_L = [T0_L; T1_L; T2_L; T3_L; T4_L; T5_L; T6_L];

R=8.314;
L_phi = -R*B_LCp*V*T_L;

%Experimental data from Messikomer for comparison
Data = readmatrix('ExcessEnthalpyLPhiNaClMessikomer.txt');
L_phi_Messikomer = Data(:,2:6)*4.184;
b_Messikomer = Data(:,1);
T_Messikomer = 273.15 + [25, 30, 40, 50, 60];

set(0,'defaulttextinterpreter','latex')

%Plot

%for k=1:1
 figure(1)
h1 = surf( T, b, L_phi);
 alpha 0.9
hold on
[x, y, z] = prepareSurfaceData( T_Messikomer, b_Messikomer, L_phi_Messikomer);
h2 = scatter3(x, y, z,  'Marker','s',...
     'MarkerEdgeColor','k',...
     'MarkerFaceColor','y');

legend([h1, h2],'Calculated', 'Messikomer', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( '$T$ (K)', 'Interpreter', 'latex' );
ylabel( '$b$ (mol/kg)', 'Interpreter', 'latex' );
zlabel( '$$\Delta H_{\phi}$$ (J/mol)', 'Interpreter', 'latex' );
xlim([290 340])
ylim([0 6])
zlim([-2000 2000])
view( -53.2, 26.2);
set(gca,'FontSize',14)
hold off
 %saveas(h1,sprintf('LphiNaCl-%d.png',k)); % will create FIG1, FIG2,...
%end



%Calculated values from Messikomer for comparison
T0_L_Messikomer = zeros(size(T_Messikomer));
T1_L_Messikomer = -ones(size(T_Messikomer));
T2_L_Messikomer = -(T_Messikomer-Theta);
T3_L_Messikomer = -(1/2.*(T_Messikomer-Theta).^2);
T4_L_Messikomer = -(1/6.*(T_Messikomer-Theta).^3);
T5_L_Messikomer = -(1/24.*(T_Messikomer-Theta).^4);
T6_L_Messikomer = -(1/120.*(T_Messikomer-Theta).^5);

T_L_Messikomer = [T0_L_Messikomer; T1_L_Messikomer; T2_L_Messikomer; T3_L_Messikomer; T4_L_Messikomer; T5_L_Messikomer; T6_L_Messikomer];

b0_LCp_Messikomer = -4/1.2*log(1+1.2*b_Messikomer.^(1/2));
b1_LCp_Messikomer = 1-(1+2*b_Messikomer.^(1/2)).*exp(-2*b_Messikomer.^(1/2));
b2_LCp_Messikomer = 2*b_Messikomer;
b3_LCp_Messikomer = b_Messikomer.^2;
b4_LCp_Messikomer = 2/3*b_Messikomer.^3;
b5_LCp_Messikomer = 1/2*b_Messikomer.^4;

B_LCp_Messikomer = [b0_LCp_Messikomer, b1_LCp_Messikomer, b2_LCp_Messikomer, b3_LCp_Messikomer, b4_LCp_Messikomer, b5_LCp_Messikomer];

L_phi_Messikomer_calc = R*B_LCp_Messikomer*V*T_L_Messikomer;

%Plot residuals
%{
%for k=1:360
 figure(1)
h3 = stem3( T_Messikomer, b_Messikomer, (L_phi_Messikomer-L_phi_Messikomer_calc)./L_phi_Messikomer, 'Marker','s',...
     'MarkerEdgeColor','k',...
     'MarkerFaceColor','y');
legend([h3],'Residuals', 'Location', 'NorthEast', 'Interpreter', 'latex' );
% Label axes
xlabel( 'T (K)', 'Interpreter', 'latex' );
ylabel( 'b (mol/kg)', 'Interpreter', 'latex' );
zlabel( '\%', 'Interpreter', 'latex' );
view( -90+k, 26.2);
set(gca,'FontSize',14)
hold off
 %saveas(h3,sprintf('residuals-%d.png',k)); % will create FIG1, FIG2,...
%end
%}