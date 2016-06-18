function vtb3_1(m,c,k,Fo,tf)

%VTB3_1 Impulse response of a SDOF system.
% VTB3_1(m,c,k,Fo,tf) plots the response of the system to an
% impulse of magnitude Fo.  The input variables are the mass m,
% the stiffness k, and the damping c.  The total time of the 
% response is tf.
% VTB3_1(zeta,w,fo,tf) plots the response of the system 
% described by the damping ratio zeta and the 
% undamped natural frequency w to a normalized impulse (Fo/m) 
% of magnitude fo.
%
% Example
% vtb3_1(.005,100,1,20)
% Corrected by Richard Pappa, 4/2007
% Added settling time, 4/2007
   
if nargin==5
  w=sqrt(k/m);zeta=c/2/w/m;fo=Fo/m;
    else
  tf=Fo;zeta=m;w=c;fo=k;
end

tau=1/(zeta*w);
dt=2*pi/(10*w);
t=0:dt:tf;
wd=w*sqrt(1-zeta^2);
x=fo/wd.*exp(-zeta*w*t).*sin(wd*t);
plot(t,x)
grid on
%end
xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')
hold on
plot([3*tau 3*tau], [min(x) max(x)*0]/4,'k')
text(3*tau+.02*max(t),min(x)/4,'5% Settling time')
plot([4*tau 4*tau], [min(x)*0 max(x)]/4,'k')
text(4*tau+.02*max(t),-min(x)/4,'2% Settling time')
hold off

%Automatically check for updates
vtbchk
