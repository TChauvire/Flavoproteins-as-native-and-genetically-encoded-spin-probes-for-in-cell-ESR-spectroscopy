%%%
%This script is used for globally fitting isotopically modified (14N/15N) cw-ESR %spectrum.
%The parameters employed are fitted via a variable 'Param' inserted in the %classically defined %spin system 'Sys' and used in a 'user-defined' function (here %named 'Global14N15N_M4G').
%
%Timothee Chauvire November, 11th, 2024
%%%
clear all;
Filename1= 'A040619_AerD2O_M=4G.dat'; %load experimental spectrum Aer(2H/14N) 
Filename2 = 'A041019_AerN15D2O_M=4G.dat'; %load experimental spectrum Aer(2H/15N) 
temp1 = importdata(Filename1);
temp2 = importdata(Filename2);
 
B1 = temp1(:,1)/10; % Convert Magnetic Field in mT
spc14N = temp1(:,2)/max(temp1(:,2)); % Normalize the spectra to simulate
B2 = temp2(:,1)/10; % Convert Magnetic Field in mT
spc15N = temp2(:,2)/max(temp2(:,2)); % Normalize the spectra to simulate
 
% Definition of the 2 set of spin systems to simulate the 14N and 15N hyperfine interaction
Sys14N = struct('g',[2.00436 2.00402 2.00228],'S',0.5,'lwpp',[0.22 0]);
Sys15N = struct('g',[2.00436 2.00402 2.00228],'S',0.5,'lwpp',[0.3 0]);
Sys14N.Nucs = '14N, 14N';
Sys14N.n = [1 1]; 
Sys15N.Nucs = '15N, 15N';
Sys15N.n = [1 1]; 
Sys14N.A = [1.666300205	1.217980273	56.34843717;	4.517784059	4.44814188	23.12313496]; % define the hfccs of the two 14N used for the simulation
 
% Adjust the different spin system value by using a set of coefficient parameters
Param = [1.0 1.0 1.0 1.0 1.0 1.0 1.0 1.0 0.84];
Sys14N.Param = Param; % Create a Variable that will be incorporate in the spin system structure used by EasySpin for the simulation
gyr14N15N = 1.4027; % Gyromagnetic ration for the conversion between 14N and 15N
Sys15N.A = Sys14N.A.*gyr14N15N; % copy the hfccs of the two nitrogens 
Sys15N.Param = Param;
% define the experimental parameters:
Exp14N = struct('mwFreq',9.1965,'Harmonic',1,'Temperature',295,...
     'Range',[B1(1) B1(end)],'ModAmp',0.4,'nPoints',1024);
Exp15N = struct('mwFreq',9.1925,'Harmonic',1,'Temperature',295,...
     'Range',[B2(1) B2(end)],'ModAmp',0.4,'nPoints',1024);
Exp = Exp14N;
 
% define the optimization parameters used by the function pepper
Opt.Verbosity = 2;
Opt.nKnots = 20;
Opt.Method = 'hybrid';         
Vary.Param = [0.4 0.4 0.4 0.4 0.4 0.4 0.2 0.2 0.2];
spc = [spc14N;spc15N] ; % concatenate the 2 experimental spectrum Aer(2H/14N) and Aer (2H/15N) in one variable 'spc'
Exp.nPoints = 2048; % Use the size of the variable 'spc' for the fit
 
% Execute the fitting procedure by employing an homemade function 'Global14N15N_M4G'
esfit('Global14N15N_M4G',spc,Sys14N,Vary,Exp,Opt)
 
%%%
%This function is used as a global fit procedure for simulating isotopically 
%modified (14N/15N) cw-ESR spectrum.
%The parameters employed are stored in a variable called 'Param'.
%The trick is then to duplicate the 14N spin system to a 15N spin system and adjust %the fit via %the variable 'Param'.
%The two spectra simulated are then concatenated via an output variable y.
%
%Timothee Chauvire November, 11th, 2024
%%%
function y = Global14N15N_M4G(Sys,Exp,Opt)
fullSys14N = Sys; %extract the spin system used for the fitting
fullSys15N = fullSys14N; %duplication of the spin system for Aer(2H/15N) spectrum
Param = Sys.Param; %extract the fitting parameters 
%Adjust the hfccs tensor and linewidth of the simulated spectrum with the fitting parameters
fullSys14N.A(1,1) = Sys.A(1,1)*Param(1); 
fullSys14N.A(1,2) = Sys.A(1,2)*Param(2); 
fullSys14N.A(1,3) = Sys.A(1,3)*Param(3); 
fullSys14N.A(2,1) = Sys.A(2,1)*Param(4); 
fullSys14N.A(2,2) = Sys.A(2,2)*Param(5); 
fullSys14N.A(2,3) = Sys.A(2,3)*Param(6); 
fullSys14N.lwpp = Sys.lwpp*Param(7);
Exp.nPoints = 1024; 
gyr14N15N = 1.4027; % Gyromagnetic ration for the conversion between 14N and 15N
 
fullSys15N.A = fullSys14N.A.*gyr14N15N;
fullSys15N.Nucs = '15N,15N';
fullSys15N.lwpp = Sys.lwpp*Param(8);
 
[x14N,y14N] = pepper(fullSys14N,Exp,Opt);
Exp.mwFreq = 9.166; %Compensate the magnetic field shift via a frequency shift
[x15N,y15N] = pepper(fullSys15N,Exp,Opt);
y = [y14N'; Param(9).*y15N']; %Use of a parameter for a relative scaling of the spectra
end
