function vtb2_6(m,k,dtype,dcoef,dt,tott,x0,v0)
%VTB2_6 Damping Simulations.
% VTB2_6(m,k,dtype,dcoef,dt,tott,x0,v0) plots the free decay of 
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
% 

disp('VTB2_6 has been grandfathered. Please use VTB1_5 in the future.')

vtb1_5(m,k,dtype,dcoef,dt,tott,x0,v0)

%Automatically check for updates
vtbchk
