function [z,nf,a,com]=VTB7_4(f,TF)
%[z,nf,a,com]=VTB7_4(f,TF) Curve fit to SDOF FRF.
% f is the frequency vector in Hz. It does not have to 
%    start at 0 Hz.
% TF is the complex transfer function.
% z and nf are the damping ratio and natural frequency (Hz)
% a is the product of the residues of the coordinates the 
% transfer function is between. (For example, in the example 
% below, a1 times a2 is returned. See equation 7.42)
% If com is returned as a real number, then it is the 
% compliance between the two coordinates. 
% Only one peak may exist in the segment of the FRF passed to 
% VTB7_4. No zeros may exist withing this segment. Otherwise, 
% curve fitting becomes unreliable.
%
% EXAMPLE:
% M=eye(2);
% K=[2 -1;-1 2];
% C=.01*K;
% [Freq,Recep,Mobil,Inert]=vtb7_5(M,C,K,1,2,linspace(0,.5,1024));
% figure(1)
% n=250;
% f2=Freq((1:n)+450);
% R2=Recep((1:n)+450);
% R2=R2+.1*randn(n,1)+.1*randn(n,1)*i;% Poorly Simulated Noise
% [z,nf,a,com]=vtb7_4(f2,R2)
%
% Note that by changing the parts of Freq and Recep used
% We can curve fit to other modes.

% Copyright Joseph C. Slater, 10/8/99
% Updated 11/8/99 to improve robustness
% Updated 05/9/00 to improve robustness
%    vtb7_4 now includes proprietary software from the 
%    professional edition of the vibration toolbox and has 
%    been encoded in the form of a p-file. 
