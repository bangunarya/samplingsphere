%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% This algorithm tries to search pair angles on S2(sphere) and SO(3) rotation group
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
type_matrix = "SH";
if type_matrix == "SH"
    B = 10; % Bandlimited degree
    N = B^2; % Dimension of column matrix
    m = 17:4:N;
    dim = 1;
else
    B = 4;
    [lmn,lm]=degree_order(B); %Generate combination of degree and orders
    N = length(lmn);
    m=17:4:N; % Number of measurement
    dim = 2;
end
Coh_all=zeros(1,length(m));
ang_all=cell(1,length(m));

max_iter = 10;
init_sample = @(m_sample,dim) rand(m_sample,dim)*2*pi;
for jj=1:length(m);
    %% Initialization
    Coh_val_temp=1; %Initial coherence
    x = init_sample(m(jj),dim); 
    %% Fix theta
    initial(:,1)=acos(linspace(-1,1,m(jj)))';
   
    %% Lower bound
    PB1=legendreP(B-1,cos(initial(:,1)));
    PB3=legendreP(B-3,cos(initial(:,1)));
    norm_PB1=norm(PB1);
    norm_PB3=norm(PB3);
    Best_coh_PS = abs(PB1'*PB3)/(norm_PB1*norm_PB3);
    
    
   
    
    k = 1;
    while k < max_iter && Coh_val_temp - Best_coh_PS > 1e-2
	    if type_matrix ==  "SH"
		%%%%%%%% Search Algorithm Spherical Harmonics %%%%%%%
	   
		[y_ps,Coh_val_ps] = azimuth_search(initial,x,B);
		
		%%%%%%%%%% Check whether we have better coherence
            if Coh_val_ps < Coh_val_temp
                Coh_val_temp = Coh_val_ps;
                newx_ps = mod(y_ps,2*pi);	
                x = mod(y_ps,2*pi);
            else
                new_ps = mod(y_ps,2*pi);
                x = init_sample(m(jj),dim);
            end	
	    else
		
		%%%%%%% Search algorithm Wigner-D functions %%%%%%%%%
    	
        [y_ps,Coh_val_ps] = az_pol_search(initial,x,lmn);
        
        %%%%%%% Check whether we have better coherence
            if Coh_val_ps < Coh_val_temp
                Coh_val_temp = Coh_val_ps;
                newx_ps = y_ps;
                x = mod(y_ps,2*pi);
            else
                new_ps = mod(y_ps,2*pi);
                x = init_sample(m(jj),dim);
            end
	    end
        disp(['M (measurement size) = ', num2str(m(jj)),', Lower bound = ',num2str(Best_coh_PS),', Actual Coherence = ',num2str(Coh_val_ps),', Best Achievable Coherence = ',num2str(Coh_val_temp),', Distance  = ',num2str(abs(Coh_val_temp - Best_coh_PS)), ', Iteration  = ',num2str(k)]);
 
        k = k +1;
        
        
    end
    %%%%%%%%%%%%%% Save the best Coherence
    Coh_all(jj)=Best_coh_PS;
    ang_all{jj}=[initial(:,1) newx_ps];
    clear y_ps ang_new_ps newx_ps initial
end





