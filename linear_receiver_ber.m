%%%% linear_receiver_ber.m %%%
% calculate the BER of the uncoded transmitted and received signals %


%% Demodulate the received QPSK signals%%
% QPSK demodulator at the Receiver
uncoded_bits_rx = zeros(Ns,2*size(s_est_signal,2));
for ns= 1:Ns
        B4(ns,:) = (real(s_est_signal(ns,:))<0);
        B3(ns,:) = (imag(s_est_signal(ns,:))<0);  
        
        uncoded_bits_rx(ns,(1:2:end)) = B3(ns,:);
        uncoded_bits_rx(ns,(2:2:end)) = B4(ns,:);

end

qpsk_sig_demodulated=uncoded_bits_rx;

%% Calculate the BER of the recived signal from satellite 1 %%
 T_Errors = zeros();
 T_bits = zeros();
 for ns = 1:Ns
     while T_Errors < 500
         % Calculate Bit Errors
         diff(ns,:) = uncoded_bits(ns,:) - uncoded_bits_rx(ns,:);
         T_Errors = T_Errors + sum(abs(diff(ns,:)));
         T_bits = T_bits + length(x_signal(2));
     end     
     % Calculate Bit Error Rate
     BER(ns,:) = T_Errors / T_bits;     
 end
 
 
 figure();
plot_lims = [-3 3];
subplot(2,2,1)
plot(real(qpsk_sig_demodulated(1,:)), imag(qpsk_sig_demodulated(1,:)), '.');
xlim(plot_lims);
ylim(plot_lims);
title('QPSK constellation from satelite 1 demodulated');
xlabel('real part');
ylabel('imaginary part');

subplot(2,2,2)
plot(real(qpsk_sig_demodulated(2,:)), imag(qpsk_sig_demodulated(2,:)), '.');
xlim(plot_lims);
ylim(plot_lims);
title('QPSK constellation from satelite 2 demodulated');
xlabel('real part');
ylabel('imaginary part');

subplot(2,2,3)
plot(real(qpsk_sig_demodulated(3,:)), imag(qpsk_sig_demodulated(3,:)), '.');
xlim(plot_lims);
ylim(plot_lims);
title('QPSK constellation from satelite 3 demodulated');
xlabel('real part');
ylabel('imaginary part');

 
 
 

figure();
semilogy(SNR_dB,BER(1,:),'or');
hold on;
xlabel('SNR (dB)');
ylabel('BER');
title('SNR Vs BER plot for QPSK Modualtion in Rayleigh Channel for satellite 1 data');