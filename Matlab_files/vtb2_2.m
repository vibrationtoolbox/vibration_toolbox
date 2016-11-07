function vtb2_2(m,c,k,wdr,F0,tf)
%VTB2_2 Particular solution of an underdamped single degree of freedom system.
% VTB2_2(m,c,k,wdr,F0,tf) plots the response of an underdamped single
% degree of freedom system to a sinusoidal input with amplitude F0 and 
% frequency wdr.  The argument tf is the total time of the simulation.
% The system is a mass m, damping c, and stiffness k.
% VTB2_2(zeta,w,wdr,f0,tf) plots the response of an underdamped single
% degree of freedom system to a sinusoidal input with amplitude f0 and 
% frequency wdr.  The argument tf is the total time of the simulation.
% The system parameters are the damping ratio zeta and natural frequency w.


% This statement determines which type of input format you are using.
if nargin==5
  tf=F0;f0=wdr;wdr=k;w=c;z=m;
  m=1;c=2*z*w;k=w^2;F0=f0;
end
t=0:.005*tf:tf;
f0=F0/m;
w=sqrt(k/m);
z=c/2/w/m;
if z>=1
  disp('This system is NOT underdamped, sorry!')
else

A0=f0/sqrt((w^2-wdr^2)^2+(2*z*w*wdr)^2);
x=A0*cos(wdr*t-atan2(2*z*w*wdr,w^2-wdr^2));%(2.28)
aa=version;ll=length(aa);
plot(t,x)
grid on

xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')
end

%Automatically check for updates
vtbchk
