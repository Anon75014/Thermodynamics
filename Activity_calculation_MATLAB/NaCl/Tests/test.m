function z=test(b,T,b1,c1,d1)
Theta = 298.15;

%----------Ak coefficients from Clarke 1980
A0 = 0.391940;
A1 = 0.198653*Theta;
A2 = 0.772533;
A3 = 2/Theta*1.68848;

%z=b1*A0.*b+c1*A1.*b.*T+d1*A2.*T.^2.*b+A3.*T*b.^2;

z = A0.*b1.*b.*T+A1.*c1.*b.^2.*A2.*T+d1.*b.*A3.*T.^2;