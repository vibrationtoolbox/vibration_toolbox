function wn=vtb6_4(n,E,G,I,rho,A,L,K)
%VTB6_4 Comparison of natural frequencies for an Euler-Bernoulli
%beam, a beam including rotary inertia, and a beam including 
%shear deformation.
% w = VTB6_4(n,E,G,I,rho,A,L,K) will return the first n natural 
% frequencies for the three types of beams mentioned above. The input 
% parameters are:
% 
% E = Young's Modulus
% G = Shear Modulus
% I = moment of inertia about bending axis
% rho = density
% A = cross-sectional area
% L = length
% K = shear coefficient
%
% The output w is a n x 3 matrix.  The first column is the
% frequencies for an E-B beam, the second column corresponds to
% a beam with rotary inertia, and the third column is the beam
% with shear deformation.

alpha=sqrt(E*I/(rho*A));
r=sqrt(I/A);

for i=1:n

    wn(i,1)=sqrt((alpha^2*i^4*pi^4)/L^4);   %E_B beam Note (1)

    wn(i,2)=(alpha^2*i^4*pi^4)/L^4;
    wn(i,2)=sqrt(wn(i,2)*(1/(1+i^2*pi^2*r^2/L^2)));  %Rotary inertia

    wn(i,3)=(alpha^2*i^4*pi^4)/L^4;  %Shear deformation
    wn(i,3)=sqrt(wn(i,3)*(1/(1+(i^2*pi^2*r^2/L^2)*(E/(K*G)))));

end

%Automatically check for updates
vtbchk
