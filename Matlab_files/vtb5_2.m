function vtb5_2(k,Y,r,TR)
%VTB5_2 Base excitation force transmissibility for a SDOF system. 
% VTB5_2(k,Y,r,TR) returns the force transmitted through the base
% to a SDOF system.  The input variables are:
%               
%        k:  stiffness (N/m)
%        Y:  amplitude of base vibration (m)
%        r:  ratio of driving frequency to system frequency
%        TR: Transmissibility ratio

Ft=k*Y*r^2*TR;

disp(['The transmitted force is ',num2str(Ft),' newtons.'])


%Automatically check for updates
vtbchk
