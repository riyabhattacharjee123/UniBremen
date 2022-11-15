%%%%%%%%%%%% steeringVectors.m %%%%%%%%%%%%
%% calculate the steering vectors for each channel established between %%
%% satellites and GS %%
gs_sat_distance = [] ;
gs_sat_polar_theta = [];
gs_sat_azimuth_phi = [];
sat_AoD_THETA=[];% radians, angle of departure of satellite antenna

for fa = 1:Ns
    % satellite distance from GS
    gs_sat_distance(1,fa) = r_start_sat_gs_spherical(1,(3*fa-2)); % meters
    % satellite polar/inclination angle with y-axis of Ground stn plane
    gs_sat_polar_theta(1,fa) = r_start_sat_gs_spherical(1,(3*fa-1));%radians
    % satellite elevation angle with x-axis of Ground stn x-y plane
    gs_sat_elevation_theta(1,fa) = pi/2 - gs_sat_polar_theta(1,fa);%radians    
    % satellite azimuth angle with x-axis of Ground stn x-z plane
    gs_sat_azimuth_phi(1,fa) = r_start_sat_gs_spherical(1,(3*fa-0));%radians
    
    % satellite angle of departure to Ground stn from y-axis
    gs_sat_AoD_THETA_y(1,fa) = gs_sat_elevation_theta(1,fa) - pi/2;%radians
    % satellite angle of departure to Ground stn from x-axis
    gs_sat_AoD_THETA_x(1,fa) = gs_sat_azimuth_phi(1,fa);%radians   
end

%% Calculate the Channel gain matrix (NrXNs)
sum_channel_gain = zeros();
for m_gs = 1:(size(rcvr_pos_eci,1))
    for l_ns = 1:(size(r_start_sat,2)/3)
        %% Calculate the distance between each satellite and each GS antenna
        del_x = r_start_sat(1,(3*l_ns-2))- rcvr_pos_eci(m_gs,1); % meters
        del_y = r_start_sat(1,(3*l_ns-1))- rcvr_pos_eci(m_gs,2); % meters
        del_z = r_start_sat(1,(3*l_ns-0))- rcvr_pos_eci(m_gs,3); % meters
        
        distance(m_gs,l_ns)= sqrt(del_x^2+del_y^2+del_z^2); % meters  
        
        %% Calculate the pathloss in dB and linear scale
        pathloss_dB(m_gs,l_ns) = -20*log10(lambda/(4*pi ...
                                   * 600000)) ... %distance(m_gs,l_ns)))...
                                   - g_tx_dB - g_rx_dB;
                            
        pathloss_linear(m_gs,l_ns) = 10^(pathloss_dB(m_gs,l_ns)/10);
        
        %% Calculate the Complex valued Channel Gain
        channel_gain(m_gs,l_ns) = (1/sqrt(pathloss_linear(m_gs,l_ns))) ...
                                   * exp(-1i*nu* 600000 );
                               %distance(m_gs,l_ns));                    
        
        
    end    
    
end
channel_gain_avg = mean(channel_gain,'all');
sum_channel_gain = sum(channel_gain,1);


% Calculate the steering vectors for each NtXNr path
% Each satellite(1 Nt antenna) with GS (many Nr antennas)
al_1 = [];
al_2= [];
al_3 = [];
for cns = 1:(length(gs_sat_polar_theta)) 
    %% calculate steering vectors from GS to Sat. uplink %%    
    % steering vector in x-direction
    ul_x = 1*cos(gs_sat_elevation_theta(1,cns))...
            *cos(gs_sat_azimuth_phi(1,cns)); % (lambda/(Nrx*D_Ant))
    
    % steering vector in y-direction
    ul_y = 1*sin(gs_sat_elevation_theta(1,cns)); % (lambda/(Nrx*D_Ant))
    
    % steering vector in z-direction
    ul_z = 1*cos(gs_sat_elevation_theta(1,cns))...
            *sin(gs_sat_azimuth_phi(1,cns));
    
    a_ul_x=[];
    a_ul_y=[];
    for cnrx=1:Nrx
        a_ul_x= (1/sqrt(Nrx))*[a_ul_x (exp(-1i*(cnrx-1)*ul_x))];
        a_ul_y = (1/sqrt(Nrx))*[a_ul_y (exp(-1i*(cnrx-1)*ul_y))];        
    end
    
    if cns==1
        al_1 =  kron((a_ul_x.'),(a_ul_y.'));
    elseif cns==2
        al_2 =  kron((a_ul_x.'),(a_ul_y.'));
    else
        al_3 =  kron((a_ul_x.'),(a_ul_y.'));
    end
    
    %% Calculate steering vectors from Sat. to GS downlink %%

    % steering vector in x-direction
    vl_x = 1*sin(gs_sat_AoD_THETA_y(1,cns))...
            *cos(gs_sat_AoD_THETA_x(1,cns)) ;
    
    % steering vector in y-direction
    vl_y = 1*cos(gs_sat_AoD_THETA_y(1,cns));
    
    % steering vector in z-direction
    vl_z = 1*sin(gs_sat_AoD_THETA_y(1,cns))...
            *sin(gs_sat_AoD_THETA_x(1,cns)) ;
    
    b_vl_x=[];
    b_vl_y=[];
    
    for cntx=1:N_sat_ant
        b_vl_x= (1/sqrt(N_sat_ant))*[b_vl_x (exp(-1i*(cntx-1)*vl_x))];
        b_vl_y = (1/sqrt(N_sat_ant))*[b_vl_y (exp(-1i*(cntx-1)*vl_y))];        
    end
    
    if cns==1
        bl_1 =  kron((b_vl_x.'),(b_vl_y.'));
    elseif cns==2
        bl_2 =  kron((b_vl_x.'),(b_vl_y.'));
    else
        bl_3 =  kron((b_vl_x.'),(b_vl_y.'));
    end
    
end

% Steering matrix by concatenating all steering vectors columnwise
steering_matrix_A = ([al_1 al_2 al_3]) ;

%for cnr=1:Nr
    %steering_matrix_B(:,cnr) = ([bl_1 bl_2 bl_3])' ; %;    
%end

steering_matrix_B = eye(Ns);

% Calculate the Estimated Channel matrix %
Channel_matrix_H = (channel_gain_avg)*(steering_matrix_A) *(steering_matrix_B);
%Channel_matrix_H = (channel_gain)*(steering_matrix_B)';
    




% calculating discreet fourier transform matrix
%dft_1 = conj((al_2.').')*(al_1.');
%dft_1 = conj((al_1.').')*(al_2.');
%sum_dft_1 = sum(dft_1(:));

%dft_2 = conj((al_2.').')*(al_3.');
%sum_dft_2 = sum(dft_2(:));

%dft_3 = conj((al_3.').')*(al_1.');
%sum_dft_3 = sum(dft_3(:));

%% Calculate the channel precoder matrix

G_geo = sqrt(P_tx_gs/N_sat_ant).* steering_matrix_B ;

% Calculate the channel response vectors for each satellite to GS antenna path

%G_ul_1 = (al_1)* 10^(gain_channel_ul/10) ;
%G_ul_2 = (al_2)* 10^(gain_channel_ul/10) ;
%G_ul_3 = (al_3)* 10^(gain_channel_ul/10) ;

%G_ul = [G_ul_1 G_ul_2 G_ul_3] ;
G_dl = (Channel_matrix_H)*(G_geo) ;
G_dl_hermitian = (G_dl)' ;


