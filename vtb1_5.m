function vtb1_5(m,k,dtype,dcoef,dt,tott,x0,v0)
%VTB1_5 Damping Simulations.
% VTB1_5(m,k,dtype,dcoef,dt,tott,x0,v0) plots the free decay of 
% single degree of freedom systems with different types of 
% damping. m is the mass, k is the stiffness, dtype is the 
% damping type number where 1 is linear viscous damping, 2 is 
% coulomb damping, and 3 is air damping. dt is the time step 
% size for the numerical simulation, and tott is the total 
% time of the simulation.  x0 is the initial displacement,
% v0 is the initial velocity.  The variable dcoef represents 
% the damping coefficient for the type of damping you have chosen.
% i.e. c for linear damping, mu N for coulomb damping, and
% alpha for air damping.  Note that the point at which 
% 'sticktion' occurs for coulomb damping is predicted, but
% this prediction is tenuous at best.
n=floor(tott/dt+1);
z(:,1)=[x0;v0];
h=dt;
t=0:dt:tott;
k1=[0;0];k2=[0;0];k3=[0;0];k4=[0;0];
if dtype==1
 for l1=1:(n-1);
   z1=z(:,l1);

   k1(1)=h*z1(2);
   k1(2)=h*(-k*z1(1)-dcoef*z1(2))/m;
   k2(1)=h*(z1(2)+.5*k1(2));
   k2(2)=h*(-k*(z1(1)+.5*k1(1))-dcoef*(z1(2)+.5*k1(2)))/m;
   k3(1)=h*(z1(2)+.5*k2(2));
   k3(2)=h*(-k*(z1(1)+.5*k2(1))-dcoef*(z1(2)+.5*k2(2)))/m;
   k4(1)=h*(z1(2)+k3(2));
   k4(2)=h*(-k*(z1(1)+k3(1))-dcoef*(z1(2)+k3(2)))/m;

   z(:,l1+1)=z(:,l1)+1/6*(k1+2*k2+2*k3+k4);
 end
coefval=['(c = ' num2str(dcoef) ' )'];
end

if dtype==2
 for l1=1:(n-1);
   z1=z(:,l1);

   k1(1)=h*(z1(2));
   k1(2)=h*(-k*z1(1)-dcoef*sign(z1(2)))/m;
   k2(1)=h*(z1(2)+.5*k1(2));
   k2(2)=h*(-k*(z1(1)+.5*k1(1))-dcoef*sign(z1(2)+.5*k1(2)))/m;
   k3(1)=h*(z1(2)+.5*k2(2));
   k3(2)=h*(-k*(z1(1)+.5*k2(1))-dcoef*sign(z1(2)+.5*k2(2)))/m;
   k4(1)=h*(z1(2)+k3(2));
   k4(2)=h*(-k*(z1(1)+k3(1))-dcoef*sign(z1(2)+k3(2)))/m;

   z(:,l1+1)=z(:,l1)+1/6*(k1+2*k2+2*k3+k4);
   if abs(z1(2))<abs(x0*sqrt(k/m)+v0)*.05 & dcoef>abs(k*z1(1))
     z(:,l1+1)=z(:,l1);
     z(2,l1+1)=0;
   end
 end
   coefval=['(\mu = ' num2str(dcoef) ' )'];
end

if dtype==3
 for l1=1:(n-1);
  z1=z(:,l1);

  k1(1)=h*(z1(2));
  k1(2)=h*(-k*z1(1)-dcoef*sign(z1(2))*(z1(2))^2)/m;
  k2(1)=h*(z1(2)+.5*k1(2));
  k2(2)=h*(-k*(z1(1)+.5*k1(1))-dcoef*sign(z1(2)+.5*k1(2))*(z1(2)+.5*k1(2))^2)/m;
  k3(1)=h*(z1(2)+.5*k2(2));
  k3(2)=h*(-k*(z1(1)+.5*k2(1))-dcoef*sign(z1(2)+.5*k2(2))*(z1(2)+.5*k2(2))^2)/m;
  k4(1)=h*(z1(2)+k3(2));
  k4(2)=h*(-k*(z1(1)+k3(1))-dcoef*sign(z1(2)+k3(2))*(z1(2)+k3(2))^2)/m;

  z(:,l1+1)=z(:,l1)+1/6*(k1+2*k2+2*k3+k4);
 end
coefval=['(\alpha = ' num2str(dcoef) ' )'];
end
aa=version;ll=length(aa);
plot(t,z(1,:))
xlabel('Time')
ylabel('Displacement')

title([ 'Displacement versus Time ' coefval])
grid on
zoom on

%Automatically check for updates
vtbchk
