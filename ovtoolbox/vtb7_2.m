function [f,TXF]=VTB7_2(x,f,dt,n)

%VTB7_2 
%  [W,TF] = VTB7_2(x,f,dt,n)
%  generates frequency response function (H(iw))
%  dt is the time step of the sampled data
%  n is the number of points in the fft.
%  recommended value is 256 unless zooming
%  use higher n for zooming.  Load the data
%  in vtb7_2ex.mat for example.

lf = length(f);
if (nargin == 3), n=256; end
if n < 128
  n=256;
end

w=(sin(pi*(0:lf-1)/(lf-1)).^2)';

x=x(:);
f=f(:);
xw=x.*w;
fw=f.*w;
FX=fft(xw);
FF=fft(fw);
SXF=FF.*conj(FX);
SXX=FX.*conj(FX);
SFF=FF.*conj(FF);
SFX=FX.*conj(FF);
TXF=SXX./SXF;TXF2=SFX./SFF;
lTXF2=length(TXF)/2;
freq=[0 ((1:lTXF2-1)/lTXF2/2/dt)]';
TXF=TXF(1:lTXF2);
mag=abs(TXF);
ang=angle(TXF)*180/pi;
aa=version;ll=length(aa);
semilogy(freq,mag)
title('Txf - Transfer function magnitude'), ...
xlabel('Frequency (Hz)')
grid on
pause
plot(freq,ang)
title('Txf - Phase'), ...
xlabel('Frequency (Hz)')
grid on
pause
Coh=(abs(SXF).^2)./SXX./SFF;
plot(freq,Coh(1:lTXF2))
title('Txf - Coherence'), ...
xlabel('Frequency (Hz)')
grid on
pause



