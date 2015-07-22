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
   
if nargin==5
  w=sqrt(k/m);zeta=c/2/w/m;fo=Fo/m;
    else
  tf=Fo;zeta=m;w=c;fo=k;
end

t=0:tf/300:tf;
wd=w*sqrt(1-zeta^2);
x=fo/wd.*exp(-zeta*w*t).*sin(wd*t);
aa=version;ll=length(aa);
plot(t,x)
grid on
end
xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')
