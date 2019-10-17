%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
%%%%  Evaluate coherence of all sampling patterns for spherical harmonics matrix
%%%%  including proposed sampling.
%%%%  For well-known sampling points,i.e., spiral,fibonacci,hammersley are
%%%%  given in "A Comparison of Popular Point Configurations on S2"
%%%%  https://github.com/gradywright/spherepts
%%%%  Created by Arya Bangun at TI RWTH Aachen 2018 31.08.2018
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Check Coherence SH
clear all
close all
 
load SH_N100.mat 
 
N_SH=B^2;
%% Product
 

%%
Legendre_bound1  = zeros(1,length(m));
Welchbound_SH = Legendre_bound1;
Coh_equi = Legendre_bound1;
Coh_spi = Legendre_bound1;
Coh_fib = Legendre_bound1;
Coh_Hamm = Legendre_bound1;
Coh_proposed = Legendre_bound1;
for ii=1:length(m);
    %% Welchbound
    Welchbound_SH(ii)=sqrt((N_SH-m(ii))/(m(ii)*(N_SH-1)));
    %% Create matrix from equiangular
    ang_equi=total_angles.equi{ii};
    [SH_equi]=SH_matrix(ang_equi,B);
    Coh_equi(ii)=Tes_Coherence(SH_equi);
    %% Create matrix from spiral
 
    ang_spiral=total_angles.spiral{ii};
    [SH_spiral]=SH_matrix(ang_spiral,B);
    Coh_spi(ii)=Tes_Coherence(SH_spiral);
    %% Create matrix from Fibonacci
     ang_fibo=  total_angles.fibo{ii}; 
    [SH_fibo ]=SH_matrix(ang_fibo,B);
    Coh_fib(ii)=Tes_Coherence(SH_fibo);
    %% Create matrix from Hammersley
    
    ang_Hammersley=  total_angles.hammersley{ii}; 
    [SH_hammersley]=SH_matrix(ang_Hammersley,B);
    Coh_Hamm(ii)=Tes_Coherence(SH_hammersley);
 
    %% Create matrix from proposed  
    ang_proposed=total_angles.proposed{ii};
    [SH_proposed]=SH_matrix(ang_proposed,B);
    Coh_proposed(ii)=Tes_Coherence(SH_proposed);
  
    %% Legendre bound
    l2=B-1;
    l1=l2-2;
    Pl1=legendreP(l1,cos(ang_proposed(:,1)));
    Pl2=legendreP(l2,cos(ang_proposed(:,1)));
    Legendre_bound1(ii)=abs(Pl1'*Pl2)/(norm(Pl1)*norm(Pl2));
      ii
 end
FS=60;
MS=30;
figure;plot(m,Coh_equi,'-ob','LineWidth',6,'MarkerSize',MS);
hold on;
grid on
plot(m,Coh_spi,'-or','LineWidth',6,'MarkerSize',MS);
plot(m,Coh_fib,'-ok','LineWidth',6,'MarkerSize',MS);
plot(m,Coh_Hamm,'-om','LineWidth',6,'MarkerSize',MS);
 plot(m,Coh_proposed,'-oc','LineWidth',6,'MarkerSize',MS);
 plot(m,Legendre_bound1,'--sr','LineWidth',6,'MarkerSize',MS);
plot(m,Welchbound_SH,'--sg','LineWidth',6,'MarkerSize',MS);
title(['Coherence of sampling patterns on the sphere (N =', num2str(N),')'],'Interpreter','latex','FontSize',FS);
ylabel('Coherence','Interpreter','latex','FontSize',FS);
xlabel('Samples (m)','Interpreter','latex','FontSize',FS);
lgd=legend('Equiangular sampling','Spiral sampling','Fibonacci sampling','Hammersley sampling'...
    ,'Proposed sampling','Proposition 2','Welch bound');
lgd.FontSize = 60;
set(gca,'fontsize',50);
xlim([m(1) m(end)])
ylim([0 1])

 
