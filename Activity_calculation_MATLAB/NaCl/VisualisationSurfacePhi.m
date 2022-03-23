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
 
 %---------- Molality Matrix for phi
 b1 = [0:0.02:0.1]';
 b2= [0.1:0.05:0.4]';
 b3= [0.4:0.2:6]';
 b=vertcat(b1,b2,b3);
 b0_phi = -b.^(1/2)./(1+1.2*b.^(1/2));
b1_phi = b.*exp(-2*b.^(1/2));
b2_phi = b;
b3_phi = b.^2;
b4_phi = b.^3;
b5_phi = b.^4;

B_phi = [b0_phi, b1_phi, b2_phi, b3_phi, b4_phi, b5_phi];

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);
S5(x) = symsum((n/(4+n)*(-x)^(n-1)),n,0,Inf);
S6(x) = symsum((n/(5+n)*(-x)^(n-1)),n,0,Inf);
S7(x) = symsum((n/(6+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for phi
Theta = 298.15;
  T = 273.15:1:338.15;
x_phi = (T-Theta)./Theta;

T0_phi = -1/Theta.*ones(size(T));
T1_phi = (1/Theta.*x_phi./(x_phi+1));
T2_phi = (x_phi.^2.*vpa(S2(x_phi)));
T3_phi = (1/2*Theta.*x_phi.^3.*vpa(S3(x_phi)));
T4_phi = (1/6*Theta^2.*x_phi.^4.*vpa(S4(x_phi)));
T5_phi = (1/24*Theta^3.*x_phi.^5.*vpa(S5(x_phi)));
T6_phi = (1/120*Theta^4.*x_phi.^6.*vpa(S6(x_phi)));

T_phi = [T0_phi; T1_phi; T2_phi; T3_phi; T4_phi; T5_phi; T6_phi];

phi_calc = double(1+B_phi*V*T_phi); %double instead of sym for the surf() function

%Experimental data from Messikomer for comparison
Data = readmatrix('OsmoticCoefficientNaClGibbard.txt');
phi_Gibbard = Data(:,2:4);
b_Gibbard = Data(:,1);
T_Gibbard = 273.15 + [0, 25, 50];


%Plot

%for k=1:360
figure(1)
h1 = surf(T, b, phi_calc);
alpha 0.9
hold on
[x, y, z] = prepareSurfaceData( T_Gibbard, b_Gibbard, phi_Gibbard);
h2 = scatter3(x, y, z, 'Marker','s',...
     'MarkerEdgeColor','k',...
     'MarkerFaceColor','y');
legend([h1, h2],'Calculated', 'Gibbard', 'Location', 'NorthEast', 'Interpreter', 'latex');
% Label axes
xlabel( 'T (K)', 'Interpreter', 'latex' );
ylabel( 'b (mol/kg)', 'Interpreter', 'latex' );
zlabel( '$\phi$', 'Interpreter', 'latex' );
set(gca,'FontSize',14)
view( -53.2, 26.2)
hold off
 %saveas(h1,sprintf('phiNaCl-%d.png',k)); % will create FIG1, FIG2,...
%end


%residuals
b0_phi_Gibbard = -b_Gibbard.^(1/2)./(1+1.2*b_Gibbard.^(1/2));
b1_phi_Gibbard = b_Gibbard.*exp(-2*b_Gibbard.^(1/2));
b2_phi_Gibbard = b_Gibbard;
b3_phi_Gibbard = b_Gibbard.^2;
b4_phi_Gibbard = b_Gibbard.^3;
b5_phi_Gibbard = b_Gibbard.^4;

B_phi_Gibbard = [b0_phi_Gibbard, b1_phi_Gibbard, b2_phi_Gibbard, b3_phi_Gibbard, b4_phi_Gibbard, b5_phi_Gibbard];

 %----Taylor development for the calculation of Temperature Matrixes       
syms x n

S2(x) = symsum((n/(1+n)*(-x)^(n-1)),n,0,Inf);
S3(x) = symsum((n/(2+n)*(-x)^(n-1)),n,0,Inf);
S4(x) = symsum((n/(3+n)*(-x)^(n-1)),n,0,Inf);
S5(x) = symsum((n/(4+n)*(-x)^(n-1)),n,0,Inf);
S6(x) = symsum((n/(5+n)*(-x)^(n-1)),n,0,Inf);
S7(x) = symsum((n/(6+n)*(-x)^(n-1)),n,0,Inf);

%-----Temperature matrix for phi
x_phi_Gibbard = (T_Gibbard-Theta)./Theta;

T0_phi_Gibbard = -1/Theta.*ones(size(T_Gibbard));
T1_phi_Gibbard = (1/Theta.*x_phi_Gibbard./(x_phi_Gibbard+1));
T2_phi_Gibbard = (x_phi_Gibbard.^2.*vpa(S2(x_phi_Gibbard)));
T3_phi_Gibbard = (1/2*Theta.*x_phi_Gibbard.^3.*vpa(S3(x_phi_Gibbard)));
T4_phi_Gibbard = (1/6*Theta^2.*x_phi_Gibbard.^4.*vpa(S4(x_phi_Gibbard)));
T5_phi_Gibbard = (1/24*Theta^3.*x_phi_Gibbard.^5.*vpa(S5(x_phi_Gibbard)));
T6_phi_Gibbard = (1/120*Theta^4.*x_phi_Gibbard.^6.*vpa(S6(x_phi_Gibbard)));

T_phi_Gibbard = [T0_phi_Gibbard; T1_phi_Gibbard; T2_phi_Gibbard; T3_phi_Gibbard; T4_phi_Gibbard; T5_phi_Gibbard; T6_phi_Gibbard];

phi_calc_Gibbard = double(1+B_phi_Gibbard*V*T_phi_Gibbard);

%Plot residuals
%{
for k=1:360
figure(1)
h3 = stem3( T_Gibbard, b_Gibbard, (phi_Gibbard-phi_calc_Gibbard)./phi_Gibbard*100, 'Marker','s',...
     'MarkerEdgeColor','b',...
     'MarkerFaceColor','g');
legend([h3],'Residuals', 'Location', 'NorthEast', 'Interpreter', 'none');
% Label axes
xlabel( 'T (K)', 'Interpreter', 'latex','FontSize',14 );
ylabel( 'b (mol/kg)', 'Interpreter', 'latex','FontSize',14 );
zlabel( '\%', 'Interpreter', 'latex','FontSize',14 );
set(gca,'FontSize',14)
view( -90+k, 26.2)
 saveas(h3,sprintf('phiNaClresiduals-%d.png',k)); % will create FIG1, FIG2,...
end
%}