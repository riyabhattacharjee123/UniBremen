%%%%%%%%%%%% steeringVectors.m %%%%%%%%%%%%
%% calculate the steering vectors for each channel established between %%
%% satellites and GS %%
gs_sat_distance = [] ;
gs_sat_polar_theta = [];
gs_sat_azimuth_phi = [];
sat_AoD_THETA=[];% radians, angle of departure of satellite antenna

for fa = 1:Nt
    % satellite distance from GS
    gs_sat_distance(1,fa) = r_start_sat_gs_spherical(1,(3*fa-2)); % meters
    % satellite polar/inclination angle with y-axis of Ground stn plane
    gs_sat_polar_theta(1,fa) = r_start_sat_gs_spherical(1,(3*fa-1));%radians
    % satellite elevation angle with x-axis of Ground stn x-y plane
    gs_sat_elevation_theta(1,fa) = pi/2 - gs_sat_polar_theta(1,fa);%radians
    
    % satellite angle of departure to Ground stn 
    gs_sat_AoD_THETA(1,fa) = gs_sat_elevation_theta(1,fa) - pi/2;%radians
    
    % satellite azimuth angle with x-axis of Ground stn x-z plane
    gs_sat_azimuth_phi(1,fa) = r_start_sat_gs_spherical(1,(3*fa-0));%radians
end

% Calculate the steering vectors for each NtXNr path
% Each satellite(1 Nt antenna) with GS (many Nr antennas)
al_1 = [];
al_2= [];
al_3 = [];
for cns = 1:(length(gs_sat_polar_theta)) 
    
    % steering vector in x-direction
    ul = (lambda/(Nrx*D_Ant))*cos(gs_sat_elevation_theta(1,cns))...
            *cos(gs_sat_azimuth_phi(1,cns));
    
    % steering vector in y-direction
    vl = (lambda/(Nrx*D_Ant))*sin(gs_sat_elevation_theta(1,cns));
    
    a_ul=[];
    a_vl=[];
    for cnrx=1:Nrx
        a_ul= (1/sqrt(Nrx))*[a_ul (exp(-1i*(cnrx-1)*ul))];
        a_vl = (1/sqrt(Nrx))*[a_vl (exp(-1i*(cnrx-1)*vl))];        
    end
    
    if cns==1
        al_1 =  kron((a_ul.'),(a_vl.'));
    elseif cns==2
        al_2 =  kron((a_ul.'),(a_vl.'));
    else
        al_3 =  kron((a_ul.'),(a_vl.'));
    end   
end

% Steering matrix by concatenating all steering vectors columnwise
steering_matrix_A = [al_1 al_2 al_3] ;

% calculating discreet fourier transform matrix
%dft_1 = conj((al_2.').')*(al_1.');
dft_1 = conj((al_1.').')*(al_2.');
sum_dft_1 = sum(dft_1(:));

%dft_2 = conj((al_2.').')*(al_3.');
%sum_dft_2 = sum(dft_2(:));

%dft_3 = conj((al_3.').')*(al_1.');
%sum_dft_3 = sum(dft_3(:));

% Calculate the channel response vectors for each satellite to GS antenna path

G_ul_1 = (al_1)* 10^(gain_channel_ul/10) ;
G_ul_2 = (al_2)* 10^(gain_channel_ul/10) ;
G_ul_3 = (al_3)* 10^(gain_channel_ul/10) ;

G_ul = [G_ul_1 G_ul_2 G_ul_3] ;
G_ul_hermitian = transpose(conj(G_ul));


