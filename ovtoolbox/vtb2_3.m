function VTB2_3(z,rmin,rmax,opt)
%VTB2_3 Steady state magnitude and phase of a single degree of freedom 
% damped system.  
% VTB2_3(zeta,rmin,rmax) plots the response of a single degree of 
% freedom system with damping ratio zeta between the frequency 
% ratios rmin and rmax.
% zeta can be a list of damping ratios.
% VTB2_3(zeta,rmin,rmax,opt) allows plotting of magnitude
% (opt=1), phase (opt=2), or magnitude and phase (opt=3, default).
%
% Example:
% vtb2_3([0:.2:1],0,2,3)
%graw('reset\n')
clg
if nargin==1
	opt=z;
	z=[.0:.05:.25];
	rmin=0;
	rmax=2;
	disp('Demo mode. Type ''help vtb2_3'' to learn how to enter values.')
end
if nargin==0 
	z=[.0:.2:2];
	rmin=0;
	rmax=2;
	disp('Demo mode. Type ''help vtb2_3'' to learn how to enter values.')
	opt=3;
end

if nargin==3, opt=3;end

r=rmin:(rmax-rmin)/1000:rmax;

z=z+eps;z;
for i=1:length(z)
  A0(i,:)=1./(1-r.^2+2*j*r*z(i));
end

grid('on')
gset nokey

if opt==3 
  subplot(2,1,1)
end

if opt ~=2
  xlabel('Frequency Ratio')
  ylabel('Normalized Amplitude (dB)')
  title('Normalized Amplitude versus Frequency Ratio')
  plot(r,20*log10(abs(A0)))
end

if opt==3
  subplot(2,1,2)
end

if opt ~=1
  xlabel('Frequency Ratio')
  ylabel('Phase lag (deg)')
  title('Phase versus Frequency Ratio')
  plot(r,-angle(A0)/pi*180)
end

