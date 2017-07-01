function shakerfrf(x,fmin,fmax,beamparams)
%SHAKERFRF(X,FMIN,FMAX,BEAMPARAMS) plots the frequency response of a base excited shker
%between an accelerometer at the shaker head and an accelerometer at a
%location X on the beam.
%[F,H]=SHAKERFRF(X,FMIN,FMAX,BEAMPARAMS) returns the frequency response of a base excited shker
%between an accelerometer at the shaker head and an accelerometer at a
%location X on the beam.
% BEAMPARAMS=[E I rho A L]; 
%
% Example: 30 cm long aluminum beam with h=1.5 cm, b=3 cm, tip response
% E=7.31e10;
% I=1/12*.03*.015^3;
% rho=2747;
% A=.015*.03;
% L=0.3;
% shakerfrf(.2,100,[E I rho A L])
%
wn=0;
i=0;
w=(fmin:(fmax-fmin)/3000:fmax)'*2*pi;

E=beamparams(1);
I=beamparams(2);
rho=beamparams(3);
A=beamparams(4);
L=beamparams(5);
zeta=0;
while wn<5*(fmax*2*pi)%Including contributions of mode with frequencies 5x the max
    i=i+1;
    [wn,xx,U]=vtb6_3(i,2,beamparams,200);
    modeint=sum(U)*L;
    a(:,i)=rho *A*modeint.*w.^2./(wn^2-w.^2+2*zeta*wn*w*sqrt(-1))*spline(xx,U,x);
    f(i)=wn/2/pi;
end
plot(w/2/pi,20*log10(abs(sum(a,2))))
f;
grid on
xlabel('Frequency (Hz)')
ylabel('FRF (dB)')
