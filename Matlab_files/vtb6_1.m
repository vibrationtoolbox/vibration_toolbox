function [w,U]=vtb6_1(E,rho,L,n,bctype,plotpar)

%VTB6_1 Natural frequencies and mode shapes of a uniform bar.
% [w,U]= VTB6_1(E,rho,L,n,bctype) will return the natural
% frequencies (w) and mode shapes (U) for the the first n modes
% of a uniform bar.  The boundary condition is specified according
% to the value of bctype:
% 
% bctype = 1: free-free
% bctype = 2: fixed-free
% bctype = 3: fixed-fixed
%
% The material properties of the bar are Young's Modulus (E) and
% the density (rho).  The length of the beam is denoted L.
% [w,U]=VTB6_1(E,rho,L,n,bctype,1) will also plot the mode shapes.
% The modal amplitude is normalized to be equal to one at x = L for
% fixed-free and free-free and x=L/2 for fixed-fixed.


if nargin==5
   plotpar=0;
end

x=[0:.01*L:L]';
c=sqrt(E/rho);

%Calculates natural frequencies from equation 6.63

if bctype ==1
   for i=1:n
      w(i,1)=(i*pi*c)/L;
      U(:,i)=cos(i*pi*x/L);
      U(:,i)=U(:,i)/U(101,i);
   end
elseif bctype==2
   for i=1:n
      w(i,1)=(2*i-1)*pi*c/(2*L);
      U(:,i)=sin((2*i-1)*pi*x/(2*L));
      U(:,i)=U(:,i)/U(101,i);
   end
elseif bctype==3
   for i=1:n
      w(i,1)=i*pi*c/L;
      U(:,i)=sin(i*pi*x/L);
      U(:,i)=U(:,i)/U(52,i);
   end
end

%Plotting routine if so chosen.
if plotpar==1
   for i=1:n
      plot(x,U(:,i))
      title(['Mode ',int2str(i),'     Natural Frequency = ',num2str(w(i)),' rad/s'])
      ylabel('Modal Amplitude')
      xlabel('Length along bar - x')
      grid on
      pause
   end
end

%Automatically check for updates
vtbchk
