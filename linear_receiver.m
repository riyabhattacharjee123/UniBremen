%% Estimate received signal and Linear receiver at Satellite %%

% Normalize the x_signal sent along each row.
% each row is sent by each GS antenna, and has 8 bits of information
%run('input_signal.m')
plot_lims = [-3 3];
%run('input_signal_qpsk.m');x_signal=qpsk_sig;
run('QPSK_symbols_complex.m');x_signal=tx_sig;
s_signal=tx_sig;
x_signal_normalized = x_signal %normalize(x_signal,1);

%% Signal transmitted by satellites
x_tx_signal = G_geo * s_signal;

%x_signal_normalized = normalize(tx_sig,2);
%for ngs = 1:Ngs     
    %x_signal_normalized(:,m) = (x_signal(:,m)-mean(x_signal(:,m)))...
     %                            /std(x_signal(:,m));        
    %x_mean(ngs,:) = mean(x_signal_normalized(ngs,:));
     %x_var(ngs,:) = var(x_signal_normalized(ngs,:));
      %x_std(ngs,:) = std(x_signal_normalized(ngs,:));
%end    

% generate AWGN noise
% convert Noise power from dBW to linear scale
%P_noise_linear = 10^(P_noise_dbW / 10); % Watts , sigma^2
%for ns=1:Ns
 %   for m=1:M
  %      noise_signal(ns,m) = P_noise_linear;
   % end
%end
% actual power in signal vector
L=1;
Power_vector=L*sum(sum(abs(s_signal).^2))/length(s_signal); 
N0=Power_vector/SNR_linear; %Find the noise spectral density
noiseVar = N0/2; %P_noise_var; %sigma_sq ; %1;

numSamples = Ngs;
sample_factor = size(x_signal_normalized,2);
t = [0:sample_factor-1];
% AWGN Complex Circular Noise
awgn_cn = sqrt(N0/2) * ...
    (randn(numSamples,sample_factor) + ...
    (1i * randn(numSamples,sample_factor)));

%awgn_cn = awgn(s_signal,SNRdB_user_input,'measured');

% Normalize noise signal 
%noise_signal_normalized = wgn(Ns,sample_factor,0);%normalize(noise_signal,2);%normalize(noise_signal,2);ones(Ns,M)
%for ns=1:Ns
    %n_mean(ns,:)=mean(noise_signal_normalized(ns,:));
    %n_var(ns,:)=var(noise_signal_normalized(ns,:));
 %   ma_var(ns) = var(awgn_cn(ns))
    %n_std(ns,:)=std(noise_signal_normalized(ns,:));
%end
%mult1 = (G_ul)*sqrt(P_tx_gs)*(x_signal_normalized);

%% received signal y at GS 
y_received_signal_sat = (Channel_matrix_H)*(G_geo)*(s_signal)+awgn_cn; 
                     
rho_ul_inv = (1/P_tx_gs)* noiseVar;% n_var; 
                     

% linear receiver equalizer matrix for maximum SINR
Sum_G_dl_k=zeros();
for ut_k = 1:Ngs
    G_dl_k(ut_k,:) = G_dl(ut_k,:);
    G_dl_k_Hermition(:,ut_k)=transpose(conj(G_dl_k(ut_k,:)));
    
    for ut_i = 1:Ngs
        if ut_i ~= ut_k
            G_dl_i = G_dl(ut_i,:);
            G_dl_i_H = transpose(conj(G_dl_i));
            mult2 = G_dl_i * G_dl_i_H ;
            Sum_G_dl_k = Sum_G_dl_k + mult2;            
        end        
    end    
    sum_element = Sum_G_dl_k + rho_ul_inv* eye(1);
    sum_element_inv = inv(sum_element) ;    
    W_k(:,ut_k)= G_dl_k_Hermition(:,ut_k)*sum_element_inv ;
end
%W = G_ul_hermitian * sum_element_inv ;
%W = inv(G_ul_hermitian*G_ul + rho_ul_inv*eye) * G_ul_hermitian ;
term1 = steering_matrix_A*(conj(transpose(steering_matrix_A)))...
        + rho_ul_inv*eye(Nr)

%W_lin = [w_l1 w_l2 w_l3] ;

%% calculate W_lin %% Used this code

%for ns_l = 1: Ns

W_lin =[];

for ns_l = 1:Ns
    g_l_geo = G_geo(ns_l,:);
    g_l_geo_H = conj(g_l_geo.');
    
    H_l = Channel_matrix_H(:,ns_l);
    H_l_H = conj(H_l.');
    sum_H_g = zeros(64,64);
    
    for ns_i = 1:Ns
        if ns_i ~= ns_l
            H_i = Channel_matrix_H(:,ns_i);
            H_i_H = conj(H_i.');
            
            %g_i_geo = G_geo(ns_i,:);
            %g_i_geo_H = conj(g_i_geo.');
            
            mult_term1 = H_i * H_i_H ;
            sum_H_g = mult_term1 + sum_H_g ;            
        end      
    end
    
    mult_term2 = (sum_H_g + noiseVar * eye(Nr))^(-1); % nooisevar/Ptx     
    mult_term3 = H_l_H ;
   
    W_lin(ns_l,:) = mult_term3 * mult_term2;     
end




%% calculate the estimated signal at GS
s_est_signal = (W_lin) * (y_received_signal_sat) ;


% Recovered signal at satellite
x_signal_recovered = s_est_signal; %(W_k)*y_received_signal_sat;


% plot the signals
%figure();
%plot(t,x_signal(1,:));
%title('Signal transmitted from Satellite 1 / Desired Signal');
%xlabel('number of bits per symbol');
%ylabel(' binary signal sent');


%figure();
%subplot(2,2,1)
%plot(real(y_received_signal_sat(1,:)),imag(y_received_signal_sat(1,:)),'.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 1 with AWGN');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,2)
%plot(real(y_received_signal_sat(2,:)),imag(y_received_signal_sat(2,:)),'.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 2 with AWGN');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,3)
%plot(real(y_received_signal_sat(3,:)),imag(y_received_signal_sat(3,:)),'.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 3 with AWGN');
%xlabel('real part');
%ylabel('imaginary part');


%figure();
%subplot(2,2,1)
%plot(real(s_est_signal(1,:)), imag(s_est_signal(1,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 1 estimated');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,2)
%plot(real(s_est_signal(2,:)), imag(s_est_signal(2,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 2 estimated');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,3)
%plot(real(s_est_signal(3,:)), imag(s_est_signal(3,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 3 estimated');
%xlabel('real part');
%ylabel('imaginary part');
                     
                     
