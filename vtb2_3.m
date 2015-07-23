function vtb2_3(z,rmin,rmax,opt)
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
% vtb2_3([0.001:.04:.2],0,2,3)

if ishold==0
	clf
end
if nargin==1
	opt=z;
	z=.0:.05:.5;
	rmin=0;
	rmax=4;
	disp('Demo mode. Type ''help vtb2_3'' to learn how to enter values.')
end
if nargin==0 
	z=.0:.2:2;
	rmin=0;
	rmax=4;
	disp('Demo mode. Type ''help vtb2_3'' to learn how to enter values.')
	opt=3;
end

if nargin==3, opt=3;end

r=rmin:(rmax-rmin)/1000:rmax;
z=z+eps;
for i=1:length(z)
	A0(i,:)=1./(1-r.^2+2*j*r*z(i));
	a{i}=['\zeta = ' num2str(z(i))];
end

if opt==3 
	subplot(2,1,1)
end
if opt ~=2
%semilogy(r,log(abs(A0)))
plot(r,20*log10(abs(A0)))
%plot(r,log10(abs(A0)))
%labels=str2num(get(gca,'yticklabel'));
%for i=1:length(labels)
%	labels(i,:)=10^labels(i);
%end
%labels
%labels=num2str(labels)
%set(gca,'yticklabel',labels,'yticklabelmode','auto')
legend(char(a),1)
xlabel('Frequency Ratio')
ylabel('Normalized Amplitude (dB)')
title('Normalized Amplitude versus Frequency Ratio')
%set(gca,'YTick',0:22.5:180)
grid on
zoom on
end

%pause
if opt==3
subplot(2,1,2)
end
if opt ~=1
plot(r,-angle(A0)/pi*180)
legend(char(a),4)
xlabel('Frequency Ratio')
ylabel('Phase lag (\circ)')
title('Phase versus Frequency Ratio')
set(gca,'YTick',0:22.5:180)
grid on
zoom on
axis(axis+[0 0 -10 10])
%set(gca,'Buttondownfcn',['[x,y]=ginput(1);xlabel([''Frequency Ratio '' num2str(x)]);ylabel([''Phase lag (\circ) '' num2str(y)])'])
end



% break
% % Old code
% for i=1:length(z)
% 	z=zz(i);
% 	A0=(1)./sqrt((1-r.^2).^2+(2*z*r).^2);%(2.30)
% 	c
% 	hold on
% end
% xlabel('Frequency Ratio')
% ylabel('Normalized Amplitude')
% title('Normalized Amplitude versus Frequency Ratio')
% grid on
% zoom on
% hold off
% pause
% for i=1:length(zz)
% 	z=zz(i);
% 	phi=atan2(2*z*r,1-r.^2);%(2.30)
% 	plot(r,phi/pi*180)
% 	hold on	
% end
% 
% xlabel('Frequency Ratio')
% ylabel('Phase lag (\circ)')
% title('Phase versus Frequency Ratio')
% set(gca,'YTick',0:22.5:180)
% grid on
% zoom on
% hold off

%legend(char(a),4)

%Automatically check for updates
vtbchk
