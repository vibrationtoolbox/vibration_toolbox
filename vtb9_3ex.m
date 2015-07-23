function [zd]=rfunc(z,u,t)
%    function for    2
%                  dx      dx
%                m --  + c -- +k x = f(t)
%                    2     dt
%                  dt                     dx
%    where m=1,k=1,c=.1, and z(1)=x, z(2)=--
%                                         dt

disp('VTB9_3ex has been grandfathered. Please use VTB1_3ex in the future.')

zd=[z(2);
    -.1*z(2)-z(1)+u(2)];
