clc; clear all; close all;

[A] = readmatrix('ENDOR_N14_Aer_A_W.txt')
f = A(:,1);
I = A(:,2);


figure; 
plot(f,I/max(I),'k','Linewidth',2)
box off
set(gca,'FontSize',14)
set(gca,'fontname','Arial') 

xlim([30 75])
ylim([-0.15 1.1])
xlabel('Frequency (MHz)')
ylabel('Normalized Intensity')