function vtb2_5(m,m0,e,z,rmin,rmax)

%VTB2_5 Rotating unbalance.
% VTB2_5(m,m0,e,zeta,rmin,rmax) plots the displacement and phase of 
% a system with rotating unbalance between the frequency 
% ratios rmin and rmax.  The mass of the system is m and the damping
% ratio is zeta.  The paramaters of the unbalance are the mass m0 and
% eccentricity e.

r=rmin:(rmax-rmin)/1000:rmax;

Xn=sqrt((r.^4)./((1-r.^2).^2+(2*z*r).^2));%(2.51)
X=Xn.*m0*e/m;%(2.51)
phi=atan2(2*z*r,1-r.^2);%(2.52)

plot(r,Xn)
grid on
xlabel('Frequency Ratio')
ylabel('Normalized Displacement Magnitude')
title('Normalized Displacement Magnitude versus Frequency Ratio')
pause
plot(r,X)
grid on
xlabel('Frequency Ratio')
ylabel('Displacement Magnitude')
title('Displacement Magnitude versus Frequency Ratio')
pause
plot(r,phi)
grid on
xlabel('Frequency Ratio')
ylabel('Phase')
title('Phase versus Frequency Ratio')

%Automatically check for updates
vtbchk
