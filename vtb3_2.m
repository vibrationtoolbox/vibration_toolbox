function vtb3_2(m,c,k,Fo,tf)

%VTB3_2 Step response of a SDOF system.
% VTB3_2(m,c,k,Fo,tf) plots the response of the system to an
% step of magnitude Fo.  The input variables are the mass m,
% the stiffness k, and the damping c.  The total time of the 
% response is tf.
% VTB3_2(zeta,w,fo,tf) plots the response of the system 
% described by the damping ratio zeta and the undamped 
% natural frequency w to a normalized step function (fo=Fo/m).    


if nargin==5
  w=sqrt(k/m);zeta=c/2/w/m;fo=Fo/m;
    else
  tf=Fo;zeta=m;w=c;fo=k;
end
t=0:tf/300:tf;
wd=w*sqrt(1-zeta^2);
if zeta~=1
  phi=atan(zeta/sqrt(1-zeta^2));
end
if (zeta<1 && zeta>=0)
  x=fo/w^2*(1-w/wd*exp(-zeta*w*t).*cos(wd*t-phi));
 elseif zeta==1
  lam=-w;
  A1=-fo/w^2;
  A2=-A1*lam;
  x=fo/w^2+A1*exp(lam*t)+A2*t.*exp(lam*t);
 elseif zeta>1
  lam1=-zeta*w-w*sqrt(zeta^2-1);
  lam2=-zeta*w+w*sqrt(zeta^2-1);
  A2=fo/w^2/(lam2/lam1-1);
  A1=-lam2/lam1*A2;
  x=fo/w^2+A1*exp(lam1*t)+A2*exp(lam2*t);
 else
  disp('Zeta should be greater than zero')
end
  


plot(t,x)
xlabel('Time')
ylabel('Displacement')
title('Displacement versus Time')
grid on

%Automatically check for updates
vtbchk
