function VTB7_1
%VTB7_1 is an example of taking a power spectral density
% of data. The power spectral density of a signal containing
% .5 Hz and 1.59 Hz signals is plotted.

t=0:.1:100;t=t';
y=sin(pi*t)+sin(10*t);
figure
aa=version;ll=length(aa);
plot(t,y);
title('Sin(pi*t)+Sin(10*t) versus Time')
xlabel('Time')
ylabel('F(t)')
x=fft(y,512);
grid on
pause
Pyy=x.*conj(x)/512;
f=10*(0:255)/512;
plot(f,Pyy(1:256))
grid on
title('Power Spectral Density of F(t)')
xlabel('Frequency (Hz)')
ylabel('Power Spectral Density')
pause
