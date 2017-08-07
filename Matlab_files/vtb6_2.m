function [w,U]=vtb6_2(n,bctype,bmpar,cstype,cspar,plotpar)
%VTB6_2 Torsional natural frequencies and mode shapes of a 
% uniform bar.
% [w,U]=VTB6_2(n,bctype,bmpar,cstype,cspar) will return the first
% n natural frequencies (w) and mode shapes (U) for a uniform bar 
% undergoing torsion.  The parameters are input into the function
% in the following manner:
% 
% bctype defines the boundary condition
% 1 = free-free bar, 2 = fixed-free bar, 3 = fixed-fixed
% 
% bmpar is a vector of beam parameters defined as:
% bmpar = [G J rho L] 
%
% cstype and cmpar defined the cross-section type and parameters.
% cstype = 1 is a circular shaft and cspar = R.
% cstype = 2 is a hollow circular shaft and cspar = [R1 R2]
% cstype = 3 is a square shaft and cspar = a
% cstype = 4 is a hollow rectangular shaft and cspar = [a b A B]
%
% All variables names are consistent with Section 6.4 of the book
% and Tables 6.2 and 6.3.
%
% [w,U]=VTB6_2(n,bctype,bmpar,cstype,cspar,1) will also plot the
% mode shapes.

if nargin==5
  plotpar=0;
end

%First,calculate gamma from Table 6.2

if cstype==1    %Circular shaft
      R=cspar(1);
      g=(pi*R^4)/2;
elseif cstype==2				%hollow circular shaft
      R1=cspar(1);
      R2=cspar(2);
      g=(pi/2)*(R2^4-R1^4);
elseif cstype==3				%square shaft
      a=cspar(1);
      g=.1406*a^4;
elseif cstype==4    %hollow square shaft
      if length(cspar)~=4
          error('Not enough parameters for cstype 4')
      end
      a=cspar(1);
      b=cspar(2);
      A=cspar(3);
      B=cspar(4);
      g=(2*A*B*(a-A)^2*(b-B)^2)/(a*A+b*B-A^2-B^2);
end

if g<= 0
    error('The constant gamma is less than or equal to zero.')
end

G=bmpar(1);
J=bmpar(2);
rho=bmpar(3);
L=bmpar(4);
c=sqrt(G*g/(rho*J));

x=0:.01*L:L;
x=x';

%calculate the natural frequencies and mode shapes
%for the various boundary conditions

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
