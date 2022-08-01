%% Estimate received signal and Linear receiver at Satellite %%

% Normalize the x_signal sent along each row.
% each row is sent by each GS antenna, and has 8 bits of information
x_signal_normalized = normalize(x_signal,2);
%for ngs = 1:Ngs     
    %x_signal_normalized(:,m) = (x_signal(:,m)-mean(x_signal(:,m)))...
     %                            /std(x_signal(:,m));        
    %x_mean(ngs,:) = mean(x_signal_normalized(ngs,:));
     %x_var(ngs,:) = var(x_signal_normalized(ngs,:));
      %x_std(ngs,:) = std(x_signal_normalized(ngs,:));
%end    

% generate AWGN noise
% convert Noise power from dBW to linear scale
P_noise_linear = 10^(P_noise_dbW / 10); % Watts , sigma^2
for ns=1:Ns
    for m=1:M
        noise_signal(ns,m) = P_noise_linear;
    end
end


numSamples = Ns;
noiseVar   = 1;

mA = sqrt(noiseVar / 2) * (randn(numSamples,M) + (1i * randn(numSamples,M)));

% Normalize noise signal 
noise_signal_normalized = wgn(Ns,M,0);%normalize(noise_signal,2);%normalize(noise_signal,2);ones(Ns,M)
for ns=1:Ns
    %n_mean(ns,:)=mean(noise_signal_normalized(ns,:));
    %n_var(ns,:)=var(noise_signal_normalized(ns,:));
    ma_var = var(mA(ns,:))
    %n_std(ns,:)=std(noise_signal_normalized(ns,:));
end

% received signal y
y_received_signal_sat = transpose(G_ul) ...
                         *sqrt(P_tx_gs)*x_signal_normalized...
                         + mA ; %noise_signal_normalized ;
                     
rho_ul_inv = (1/P_tx_gs)* ma_var;% n_var; 
                     
% linear receiver equalizer matrix for maximum SINR
W = inv(G_ul_hermitian*G_ul + rho_ul_inv) * G_ul_hermitian ;

% Recovered signal at satellite
x_signal_recovered = transpose(W)*y_received_signal_sat;
                     
                     
