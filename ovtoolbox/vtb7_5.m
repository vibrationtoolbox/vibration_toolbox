function [freqout,recep,mobil,inert]=vtb7_5(M,D,K,numin,numout,freq)
%VTB7_5 Transfer Function from second order system matrices.
% [Freq,Recep,Mobil,Inert] = VTB7_5(M,D,K,NUMIN,NUMOUT,Freq) 
% returns the Compliance, Mobility, and Inertance Transfer 
% Functions (FRF) between a force at degree of freedom 
% NUMIN and a response at degree of freedom NUMOUT.
% M, D, and K are the mass, damping, and stiffness matrices
% repectively. Freq is a vector of frequencies over which
% the evaluated transfer function is desired (in Hz).
%
% VTB7_5(M,D,K,NUMIN,NUMOUT,Freq) plots the Transfer Functions 
% if there are no output arguments.  Click in the region of 
% interest to zoom in.  Each click will double the size of 
% the plot. Double click to return to full scale.
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K; 
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,1,linspace(0,.5,1024));
% vtb7_5(M,C,K,1,1,linspace(0,.5,1024))
%
% See also VTB7_2.

% Copyright Joseph C. Slater, 1995
% Modifies from sostf in professional Vibration Toolbox
% Simplified for student understanding. Capabilities
% for large matrices removed. 10/7/99
% Renamed from tf to sostf 9/23/98 to avoid conflict 
% with control toolbox
% Switched call from f_sspace to ssit 9/23/98

n=256;
if exist('freq')~=1
%  nf=(sqrt(eig(K,M)))/2/pi
%  [nf2,shape]=f_sspace(K,M,1,ones(max(size(M)),1));
  [nf2,shape]=ssit(M,K,1);
  nf2=nf2*2*pi*2*pi;
  nfmax=1/nf2*1.3;
  [nf2,shape]=ssit(M,K,1);
  nfmin=nf2/4;
  freq=(nfmin:(nfmax-nfmin)/(n-1):nfmax)';
 elseif length(freq)==2
  fmin=freq(1);fmax=freq(2);
  freq=(fmin:(fmax-fmin)/(n-1):fmax)';
 elseif size(freq,1)==1
  freq=freq';
end

omega=freq*2*pi;


%adsign=(-1)^(numin+numout);
tfunc1=omega;


for i=1:length(omega)
  MDK=K+j*omega(i)*D-omega(i)^2*M;
  MDKi=inv(MDK);  
  tfunc1(i)=MDKi(numin,numout);
end




tfunc2=tfunc1.*omega*j;
tfunc3=-tfunc1.*omega.^2;

% If no left hand arguments then plot results
if nargout==0
 subplot(211)
 plot(freq,20*log10(abs(tfunc1)))
 title('Compliance Transfer Function')
 xlabel('Frequency (Hz)')
 ylabel('Mag (dB)')
 grid on
 zoom on
 subplot(212)
 phase=[angle(tfunc1(1)) ; unwrap(angle(tfunc1(2:length(tfunc1))))]*180/pi;
 plot(freq,phase)
 xlabel('Frequency (Hz)')
 ylabel('Phase (deg)')
 grid on
 sphase=sort(phase);
 numnan=sum(isnan(sphase));size(numnan);
 sphase=sphase(1:length(sphase)-numnan);
 phmin_max=[floor(min(sphase)/45)*45-5 ceil(max(sphase)/45)*45+5];
 set(gca,'YLim',phmin_max)
 gridmin_max=round(phmin_max/90)*90;
 set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
% set(gca,'GridLineStyle','--')
% gridmin_max=round(phmin_max/45)*45;
% set(gca,'YTick',gridmin_max(1):45:gridmin_max(2))
 set(gca,'GridLineStyle',':')
 set(gca,'YTickLabel',gridmin_max(1):90:gridmin_max(2))
 zoom on
% uicontrol('style','pushbutton','units','normal','position',[.91 .95 .075 .05],'string','Print','callback','print')
 pause

 subplot(211)
 plot(freq,20*log10(abs(tfunc2)))
 title('Mobility Transfer Function')
 xlabel('Frequency (Hz)')
 ylabel('Mag (dB)')
 grid on
 zoom on
 subplot(212)
 if isnan(angle(tfunc2(1)))==1
   tfunc2(1)=0;
 end
 angle(tfunc2(1:10));
 
 phase=[angle(tfunc2(1)) ; unwrap(angle(tfunc2(2:length(tfunc2))))]*180/pi;
 plot(freq,phase)
 xlabel('Frequency (Hz)')
 ylabel('Phase (deg)')
 grid on
 sphase=sort(phase);
 numnan=sum(isnan(sphase));
 sphase=sphase(1:length(sphase)-numnan);
 phmin_max=[floor(min(sphase)/45)*45-5 ceil(max(sphase)/45)*45+5];
 set(gca,'YLim',phmin_max)
 gridmin_max=round(phmin_max/90)*90;
 set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
% set(gca,'GridLineStyle','--')
% gridmin_max=round(phmin_max/45)*45;
% set(gca,'YTick',gridmin_max(1):45:gridmin_max(2))
 set(gca,'GridLineStyle',':')
 set(gca,'YTickLabel',gridmin_max(1):90:gridmin_max(2))
 zoom on
 pause

 subplot(211)
 plot(freq,20*log10(abs(tfunc3)))
 title('Inertance Transfer Function')
 xlabel('Frequency (Hz)')
 ylabel('Mag (dB)')
 grid on
 zoom on
 subplot(212)
 if isnan(angle(tfunc3(1)))==1
   tfunc3(1)=0;
 end
 phase=[angle(tfunc3(1)) ; unwrap(angle(tfunc1(2:length(tfunc3))))]*180/pi;
 plot(freq,phase)
 xlabel('Frequency (Hz)')
 ylabel('Phase (deg)')
 grid on
 sphase=sort(phase);
 numnan=sum(isnan(sphase));
 sphase=sphase(1:length(sphase)-numnan);
 phmin_max=[floor(min(sphase)/45)*45-5 ceil(max(sphase)/45)*45+5];
 if phmin_max(1)==phmin_max(2)
  phmin_max(1)=-.000000000001+phmin_max(1);
  phmin_max(2)=.000000000001+phmin_max(2);
 end
 set(gca,'YLim',phmin_max)
 gridmin_max=round(phmin_max/90)*90;
 set(gca,'YTick',gridmin_max(1):90:gridmin_max(2))
% set(gca,'GridLineStyle','--')
% gridmin_max=round(phmin_max/45)*45;
% set(gca,'YTick',gridmin_max(1):45:gridmin_max(2))
 set(gca,'GridLineStyle',':')
 set(gca,'YTickLabel',gridmin_max(1):90:gridmin_max(2))
 zoom on

 return
end

freqout=freq;
recep=tfunc1;
mobil=tfunc2;
inert=tfunc3;

