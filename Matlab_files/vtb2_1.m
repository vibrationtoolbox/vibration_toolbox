function vtb2_1(m,k,x0,v0,wdr,F0,tf)
%VTB2_1 Forced response of an undamped single degree of freedom system.
% VTB2_1(m,k,x0,v0,wdr,F0,tf) plots the response of an undamped single
% degree of freedom system to a sinusoidal input with amplitude F0 and
% frequency wdr.  The argument tf is the total time of the simulation.
% The initial displacement is x0 and initial velocity is v0.  The system
% is described by mass m and stiffness k.
% VTB2_1(w,x0,v0,wdr,f0,tf) plots the response of an undamped single
% degree of freedom system to a sinusoidal input with amplitude F0 and 
% frequency wdr.  The argument tf is the total time of the simulation.
% The natural frequency of the system, in rad/s, is w.


% This loop determines which type of input format you are using.
if nargin==6
  tf=F0;f0=wdr;wdr=v0;v0=x0;x0=k;w=m;k=m^2;m=1;F0=f0;
end
t=0:.005*tf:tf;
f0=F0/m;
w=sqrt(k/m);
x=v0/w*sin(w*t)+(x0-f0/(w^2-wdr^2))*cos(w*t)+f0/(w^2-wdr^2)*cos(wdr*t);%(2.11)
aa=version;ll=length(aa);
plot(t,x)
grid on
xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')

%Automatically check for updates
vtbchk
