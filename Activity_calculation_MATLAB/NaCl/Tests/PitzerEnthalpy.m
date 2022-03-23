function z = PitzerEnthalpy(b,T, Q0, Q1, Q2, Q3, Q4, Q5, Q6, B0, B1, B2, B3, B4, B5, C0, C1, C2, C3, C4, C5, D0, D1, D2, D3, D4, E0, E1, E2)
%b: molality
%T: temperature
%---------- Ideal gas constant
R = 8.314; %J/K/mol
%---- Reference temperature
Theta = 298.15;

%----------Ak coefficients from Clarke 1980
A0 = 0.391940.*Theta;
A1 = 0.198653.*Theta;
A2 = 0.772533;
A3 = 2./Theta.*1.68848;
A4 = 6./(Theta).^2.*1.99113;

%---------- Molality Matrix as in Clarke
b0 = -4./1.2.*log(1+1.2.*b.^(1/2));
b1 = 1-(1+2.*b.^(1/2)).*exp(-2.*b.^(1/2));
b2 = 2.*b;
b3 = b.^2;
b4 = 2./3.*b.^3;
b5 = 1./2.*b.^4;

%B_L = [b0, b1, b2, b3, b4, b5];

%---- Temperature Matrix
T0 = 0;
T1 = 1;
T2 = T-Theta;
T3 = 1/2.*(T-Theta).^2;
T4 = 1/6.*(T-Theta).^3;
T5 = 1/24.*(T-Theta).^4;
T6 = 1/120.*(T-Theta).^5;

%T_L = [T0; T1; T2; T3; T4; T5; T6];

%Virial Matrix
%V = [A0, A1, A2, A3, A4, 0, 0; Q0, Q1, Q2, Q3, Q4, Q5, Q6; B0, B1, B2, B3, B4, B5, 0; C0, C1, C2, C3, C4, C5, 0; D0, D1, D2, D3, D4, 0, 0; E0, E1, E2, 0, 0, 0, 0];

%---- Function to fit (matrix multiplication doesn't seem to work in curve
%fitting tool)
%z= B_L.*V.*T_L;


%%%%Notes on curve fitting tool
%%%Apparent problem with matrix multiplication
%%%All multiplications and power require .* and .^ (error: dimenions don't
%%%agree)
%%% the notation +... isn't understood (error:statement not finished)
%%---- This snytax worked first in the curve fitting tool:
%z = 0*((-4/1.2*log(1+1.2*b^(1/2)))*0.391940+(1-(1+2*b^(1/2))*exp(-2*b^(1/2)))*Q0+2*b*B0+b^2*C0+2/3*b^3*D0+1/2*b^4*E0)+ 1*((-4/1.2*log(1+1.2*b^(1/2)))*(0.198653*298.15)+(1-(1+2*b^(1/2))*exp(-2*b^(1/2)))*Q1+2*b*B1+b^2*C1+2/3*b^3*D1+1/2*b^4*E1)+(T-298.15)*((-4/1.2*log(1+1.2*b^(1/2)))*0.772533+(1-(1+2*b^(1/2))*exp(-2*b^(1/2)))*Q2+2*b*B2+b^2*C2+2/3*b^3*D2+1/2*b^4*E2)+(1/2*(T-298.15)^2)*((-4/1.2*log(1+1.2*b^(1/2)))*2/298.15*1.68848+(1-(1+2*b^(1/2))*exp(-2*b^(1/2)))*Q3+2*b*B3+b^2*C3+2/3*b^3*D3)+(1/6*(T-298.15)^3)*((-4/1.2*log(1+1.2*b^(1/2)))*6/(298.15)^2*1.99113+(1-(1+2*b^(1/2))*exp(-2*b^(1/2)))*Q4+2*b*B4+b^2*C4+2/3*b^3*D4)+(1/24*(T-298.15)^4)*((1-(1+2*b^(1/2)))*exp(-2*b^(1/2))*Q5+2*b*B5+b^2*C5)+(1/120*(T-298.15)^5)*((1-(1+2*b^(1/2)))*exp(-2*b^(1/2))*Q6);
%%---- Than this one:
z = T0.*(b0.*A0+b1.*Q0+b2.*B0+b3.*C0+b4.*D0+b5.*E0)+T1.*(b0.*A1+b1.*Q1+b2.*B1+b3.*C1+b4.*D1+b5.*E1)+T2.*(b0.*A2+b1.*Q2+b2.*B2+b3.*C2+b4.*D2+b5.*E2)+T3.*(b0.*A3+b1.*Q3+b2.*B3+b3.*C3+b4.*D3)+T4.*(b0.*A4+b1.*Q4+b2.*B4+b3.*C4+b4.*D4)+T5.*(b1.*Q5+b2.*B5+b3.*C5)+T6.*(b1.*Q6);
%%----weird: results are different
%%---- But for now matrix multiplication deosn't work