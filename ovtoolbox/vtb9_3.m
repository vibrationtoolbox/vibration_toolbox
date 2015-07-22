function z = VTB9_3(rkf,u,t,x0)
%VTB9_3   VTB9_3(rkf,u,t,x0)
%       Runge-Kutta fourth order solution to a first order DE
%       t is a row vector from the initial time to the final time
%       step h.  x0 is the initial value of the function.
%       The force matrix u is ordered such that the nth column
%       of u is the force vector u evaluated at time n*dt.
%       rkf is a variable containing the name of the function 
%       file.  The equations in the function file mus be written
%       in first order state space form.
%       See example VTB9_3ex.m.
%       Example
%       t=0:.5:20;  % Creates time vector
%       u=[zeros(1,length(t));sin(t*1.1)];% Creates force matrix
%       x0=[1 ;0];  % Creates initial state vector.
%       x=vtb9_3('vtb9_3ex',u,t,x0); % Runs analysis.
%       plot(t,x(1,:)); % Plots displacement versus time.
%       plot(t,x(2,:)); % Plots velocity versus time.

disp('VTB9_3 has been grandfathered. Please use VTB1_3 in the future.')


n=length(t);
z=zeros(length(x0),length(t));
z(:,1)=x0;
h=t(2)-t(1);

for l1=1:(n-1);
   z1=z(:,l1);
   u1=u(:,l1);
   u2=u(:,l1+1);

   k1=h*feval(rkf,z1,u1,t(l1));
   k2=h*feval(rkf,z1+.5*k1,u1,t(l1)+.5*h);
   k3=h*feval(rkf,z1+.5*k2,u1,t(l1)+.5*h);
   k4=h*feval(rkf,z1+k3,u1,t(l1)+h);
   z(:,l1+1)=z(:,l1)+1/6*(k1+2*k2+2*k3+k4);
end
