function [zd]=vtb1_3ex(z,u,t)
%    function for    2
%                  dx      dx
%                m --  + c -- +k x = f(t)
%                    2     dt
%                  dt                     dx
%    where m=2,k=1,c=.1, and z(1)=x, z(2)=--
%                                         dt
% z(1) is the displacement
% z(2) is the velocity
% u(2) is the mass normalized force
x=z(1);
v=z(2);
m=2;
k=3;
c=.1;
f=u(2);
zd=[z(2);
    -c/m*v-k/m*x+f];%This is the part solved for \ddot{x}
