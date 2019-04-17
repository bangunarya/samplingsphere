%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This algorithms try to search pair angles on S2(sphere) and SO(3) rotation group
%%% that have lower coherence property in a matrix from spherical harmonics
%%% and Wigner-D functions . Given
%%% deterministic angles on elevation (theta) which is equispaced sampling
%%% patterns i.e cos(theta_p) = (2p-m-1)/(m-1) for p in [m].
%%%
%%% The goal now is to find angles on azimuth (phi_p) for p in [m] for
%%% spherical harmonics and azimuth (phi_p) and polarization (chi_p) for p in [m]
%%% for Wigner-D functions
%%% N(size column) for spherical harmonics matrix N_SH = B^2
%%% N(size column) for Wigner-D functions matrix N_W = B*(2*B-1)*(2*B+1)/3;
%%% Created by Arya Bangun at TI RWTH Aachen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1
clear all;
close all

clc;
%% Generate degree and number of samples and iteration
B=10; % Bandlimited degree

[lmn,lm]=degree_order(B); %Generate combination of degree and orders
%% Preallocation
MC=10; % Number of iteration
m=17:4:97; % Number of measurement
Coh_all=zeros(1,length(m));
ang_all=cell(1,length(m));
Coh_val_temp=1; %Initial coherence
type = 'SH'; % SH or Wigner
for jj=1:length(m);
    %% Initialization
   
    %% Fix theta
    initial(:,1)=acos(linspace(-1,1,m(jj)))';
    
    %% Lower bound
    PB1=legendreP(B-1,cos(initial(:,1)));
    PB3=legendreP(B-3,cos(initial(:,1)));
    norm_PB1=norm(PB1);
    norm_PB3=norm(PB3);
    Best_coh_PS = abs(PB1'*PB3)/(norm_PB1*norm_PB3);
    %%%%%% Start Iteration %%%%%%%%%
    for ii=1:MC;
        if isequal('SH',type) == 1;
            %%%%%%% Search algorithm Spherical Harmonics
            [y_ps,Coh_val_ps] = azimuth_search(initial,m(jj),lm);
        else
            %%%%%%% Search algorithm Wigner-D functions
            [y_ps,Coh_val_ps] = az_pol_search(initial,m(jj),lmn);
            %%%%%%%%%%%%%% Check how good the coherence is
        end
        dist=abs(Best_coh_PS - Coh_val_ps);
        if dist < 1e-2;
            newx_ps=y_ps;
            
            disp(['M (measurement size) = ', num2str(m(jj)),', Lower bound = ',num2str(Best_coh_PS),' Best Achievable Coherence = ',num2str(Coh_val_ps),', Distance  = ',num2str(dist), ', MC (Monte Carlo) = ',num2str(ii)]);
            break
        else
            
            if Coh_val_ps < Coh_val_temp;
                newx_ps=y_ps;
                Coh_val_temp = Coh_val_ps;
            end
            
            disp(['M (measurement size) = ', num2str(m(jj)),', Lower bound = ',num2str(Best_coh_PS),', Actual Coherence = ',num2str(Coh_val_ps),', Best Achievable Coherence = ',num2str(Coh_val_temp),', Distance  = ',num2str(dist), ', MC (Monte Carlo) = ',num2str(ii), ' selesai']);
        end
        %% Try to change initialization
        x = newx_ps + rand(m(jj),size(newx_ps,2))*1e-2;
    end
    %%%%%%%%%%%%%% Save the best Coherence
    Coh_all(jj)=Best_coh_PS;
    ang_all{jj}=[initial(:,1) newx_ps];
     clear y_ps ang_new_ps newx_ps initial
    %disp(['Simulasi m (measurement size) = ', num2str(m(jj)),' Best Coherence = ',num2str(Best_coh_PS),'Actual Coherence = ',num2str(Coh_val_temp),', Jarak coherence  = ',num2str(dist), ', MC (Monte Carlo) = ',num2str(ii), ' selesai']);
end

