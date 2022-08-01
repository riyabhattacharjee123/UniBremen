Nr = length(rcvr_pos_eci(:,1)); % number of receiver antenna
Nt = (size(r_start_sat,2)/3); % number of transmitting antenna
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Generate Channel matrix for each channel
H = [];
for n=1:size(rcvr_pos_eci,1)
    for m=1:size(sat_trx_pos_eci,1)
               
        del_x = sat_trx_pos_eci(m,1)- rcvr_pos_eci(n,1); % meters
        del_y = sat_trx_pos_eci(m,2)- rcvr_pos_eci(n,2); % meters
        del_z = sat_trx_pos_eci(m,3)- rcvr_pos_eci(n,3); % meters        
        distance(n,m) = sqrt(del_x^2+del_y^2+del_z^2); % meters        
       
        % calculate the pathloss in dB
        pathloss_dB(n,m) = -20*log10(lambda/(4*pi*distance(n,m)))- g_tx_dB - g_rx_dB;        
        pathloss_linear(n,m) = 10^(pathloss_dB(n,m)/10);
        
        % calculate the elements of channel matrix
        H(n,m) = (1/sqrt(pathloss_linear(n,m)))* exp(-1i*nu*distance(n,m));     
        
    end    
end

% Channel Rate
H_hermitian = transpose(conj(H));
I = eye(Nr,Nr); % Identity matrix

for s = 1:size(SNR_dB,2)  
    SNR = 10^(SNR_dB(s)/10); % converting SNR from dB to linear scale
    R =  log2(real(det(I+((SNR/Nt)*H*H_hermitian)))); % bits/sec/Hz   
    Y(s)=R;
end   

