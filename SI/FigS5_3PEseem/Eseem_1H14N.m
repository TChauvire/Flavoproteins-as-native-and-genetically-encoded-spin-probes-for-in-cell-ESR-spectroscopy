%%%
%This script is used for fitting 3PEseem spectrum in frequency domain (or %alternatively in time %domain with a Background subtraction step).
%Two homemade functions are used for the fit: 'FrequencyFitEseem' or %'TimeFitEseem'.
%
%Timothee Chauvire November, 11th, 2024
%%%
clear all
Filename1 = 'A112118_N14Aer.dat';
Filename2 = 'A112118_N14Aer_Time.dat';
freq = importdata(Filename1);
time = importdata(Filename2);
 
tau = 0.21;
Exp.Sequence = '3pESEEM';
Exp.Field = 1192.8;
Exp.dt = 0.016;
Exp.tau = tau;
Exp.T = 0.1;
Exp.mwFreq = 33.67;
Exp.ExciteWidth =5e8;
 
Sys.g = [2.00436 2.00402 2.00228];
Sys.Nucs = '14N,14N';
Sys.A = [1.3 -0.8 58.6; 3.6 4.1 21.3];
Sys.T1 = [4000];
Sys.T2 = [200];
 
Opt.GridSize = 480; % A lot of orientation are needed as the hfccs of the 14N/15N are strongly anisotropic
Opt.TimeDomain = 0 ; 
Opt.ProductRule = 0;
Opt.ZeroFillFactor = 4;
Opt.Method = 'matrix';
Opt.Verbosity = 2;
Vary.A = [2 2 20;4 4 20];
 
Exp.nPoints = size(freq(:,2),1)/2; % A fit on half of the frequency domain was achieved
spcFreq = freq(1:length(freq)/2,2)/max(freq(1:length(freq)/2,2));
esfit(spcFreq,@FrequencyFitEseem,{Sys,Exp,Opt},{Vary});
 
%Either a fit in frequency or a fit in the time domain can be used:
%Exp.nPoints = size(time(:,2),1)/2;
%spcTime = time(1:length(time)/2,2)/max(time(1:length(time)/2,2)) ;
%esfit('TimeFitEseem',spcTime,Sys,Vary,Exp,Opt); 
 
%%%
%This function is used as a global fit procedure for simulating 3PEseem spectrum in 
%frequency domain.
% 
%Timothee Chauvire November, 11th, 2024
%%%
function fout = FrequencyFitEseem(Sys,Exp,Opt)
[x,y,out] = saffron(Sys,Exp,Opt); %Simulate the 3P-ESEEM spectrum in the time domain but extract the frequency domain via the 'out' variable
% Get the limit of half of the positive quadrant in the frequency domain
N1 = length(out.fd)/2 ;
N2 = length(out.fd)/4 ; 
fout = abs(out.fd(:,N1+1:N1+N2))/max(abs(out.fd(:,N1+1:N1+N2))); % extract the data in the frequency domain and normalized it
end
%%%
%This function is used as a global fit procedure for simulating 3PEseem spectrum in 
%time domain with a background subtraction (Exponential of a polynomial function).
%An Homemade function called 'bcgd_eseem' is used for the background subtraction.
% 
%Timothee Chauvire November, 11th, 2024
%%%
function y_new = TimeFitEseem(Sys,Exp,Opt)
Exp.nPoints = 512;
[x,y,out] = saffron(Sys,Exp,Opt);
xnew = 1000.*x; % Convert the ns time scale in us
[y2,background]= bcgd_eseem(xnew,real(y),3); % Homemade function to subtract an exponential polynome in the time domain
y_new = y2/max(y2);
end
 
%%%
%This function enables to fit the background via an exponential of a nth degree 
%polynomial.
%Three inputs are needed : x (the x absciss), y (the y absciss), and nPoly the 
%order of polynomial needed for the fit.
%Two outputs are generated : the new values of y with background subtraction y_new, 
%and the 'background' variable used for the subtraction.
%
%Timothee Chauvire November, 07th, 2019
%%%
function [y_new,background]= bcgd_eseem(x,y,nPoly)
ylog = log(y); % convert the time domain in logarithmic scale
p = polyfit(x,ylog,nPoly); %polynomial fit in the logarithmic domain
f = polyval(p,x); %reconstruct the polynomial function in the logarithmic domain
y_new = exp(ylog)./exp(f)-1; % subtract the background
background = exp(f); % export the background
end
