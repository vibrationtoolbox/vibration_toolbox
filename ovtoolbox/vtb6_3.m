function [w,U]=vtb6_3(n,bctype,bmpar,plotpar)

%VTB6_3 Natural frequencies and mode shapes for an Euler-
% Bernoulli beam with a given boundary condition.
% [w,U]=VTB6_3(n,bctype,bmpar) will return the first n natural 
% frequencies (w) and mode shapes (U) of an Euler-Bernoulli beam.
% The boundary condition is defined as follows:
%
% bctype = 1 free-free
% bctype = 2 clamped-free
% bctype = 3 clamped-pinned
% bctype = 4 clamped-sliding
% bctype = 5 clamped-clamped
% bctype = 6 pinned-pinned
%
% The beam parameters are input through the vector bmpar:
% bmpar = [E I rho A L];
% where the variable names are consistent with Section 6.5 of the 
% text.
% [w,U]=vtb6_3(n,bctype,bmpar,1) will also plot the mode shapes.


if nargin==3
   plotpar=0;
end

E=bmpar(1);
I=bmpar(2);
rho=bmpar(3);
A=bmpar(4);
L=bmpar(5);

len=[0:.01:1]';  %Normalized length of the beam

%Determine natural frequencies and mode shapes depending on the
%boundary condition.

if bctype==1
   desc=['Free-Free '];
   Bnl=[4.73 7.85 11 14.14 17.28];
   if n>5
      for i=6:n
         Bnl(i)=(2*i+1)*pi/2;
      end
   end
   for i=1:n
      sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      b=Bnl(i)*len;
      U(:,i)=cosh(b)+cos(b)-sig*(sinh(b)+sin(b));
      U(:,i)=U(:,i)/U(101,i);
   end

elseif bctype==2
   desc=['Clamped-Free '];
   Bnl=[1.88 4.69 7.85 10.99 14.14];
   if n>5
      for i=6:n
         Bnl(i)=(2*i-1)*pi/2;
      end
   end
   for i=1:n
      sig=(sinh(Bnl(i))-sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      b=Bnl(i)*len;
      U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
      U(:,i)=U(:,i)/U(101,i);
   end

elseif bctype==3
   desc=['Clamped-Pinned '];
   Bnl=[3.93 7.07 10.21 13.35 16.49];
   if n>5
      for i=6:n
         Bnl(i)=(4*i+1)*pi/4;
      end
   end
   for i=1:n
      sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      b=Bnl(i)*len;
      U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
      U(:,i)=U(:,i)/U(52,i);
   end

elseif bctype==4
   desc=['Clamped-Sliding '];
   Bnl=[2.37 5.50 8.64 11.78 14.92];
   if n>5
      for i=6:n
         Bnl(i)=(4*i-1)*pi/4;
      end
   end
   for i=1:n
      sig=(sinh(Bnl(i))+sin(Bnl(i)))/(cosh(Bnl(i))-cos(Bnl(i)));
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      b=Bnl(i)*len;
      U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
      U(:,i)=U(:,i)/U(101,i);
   end

elseif bctype==5
   desc=['Clamped-Clamped']
   Bnl=[4.73 7.85 11 14.14 17.28];
   if n>5
      for i=6:n
         Bnl(i)=(2*i+1)*pi/2;
      end
   end
   for i=1:n
      sig=(cosh(Bnl(i))-cos(Bnl(i)))/(sinh(Bnl(i))-sin(Bnl(i)));
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      b=Bnl(i)*len;
      U(:,i)=cosh(b)-cos(b)-sig*(sinh(b)-sin(b));
      U(:,i)=U(:,i)/U(52,i);
   end
elseif bctype==6
   desc=['Pinned-Pinned'];
   for i=1:n
      Bnl(i)=i*pi;
      w(i,1)=(Bnl(i)^2)*sqrt(E*I/(rho*A*L^4));
      U(:,i)=sin(i*pi*len);
   end
end


x=len*L;
%Plotting routine if so chosen.
if plotpar==1
   for i=1:n
      plot(x,U(:,i))
      title([desc,'  ','Mode ',int2str(i),'     Natural Frequency = ',num2str(w(i)),' rad/s'])
      ylabel('Modal Amplitude')
      xlabel('Length along bar - x')
      grid on
      pause
   end
end