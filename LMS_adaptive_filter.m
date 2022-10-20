%%lms_adaptive_filter.m%%
%% This program estimates the channel matrix based on LMS method%%

n_seq = size(x_signal,2); % length of sequence
m_taps = 1;  % number of taps
mu_step = 0.0010 ; % step size mu
t=[0:n_seq-1];
Weight = zeros(m_taps,n_seq); % weight vector

% Desired signal
D_signal = s_signal ;

% Corrupted received_signal
A_signal = (Channel_matrix_H)*(G_geo)*(D_signal) + awgn_cn ;
E = []; % error signal
for ix=(m_taps+1):n_seq
    a_array = A_signal((ix-(m_taps)+1):ix) % 1 row 
    E(ix) = D_signal(ix) - a_array*Weight(:,ix); % 1 row - weight * array
    Weight(:,ix+1) = Weight(:,ix) - mu_step * E(ix)*(a_array)';    
end
for ix=(m_taps+1):n_seq
    a_array = A_signal((ix-(m_taps)+1):ix)
    %Received signal from satellite.
    yd(:,ix) = a_array*Weight(:,ix);  
end

%figure();plot(t,D_signal(1,:));title('Desired Signal');
figure();
plot(real(A_signal(1,:)),imag(A_signal(1,:)),'.');
xlim(plot_lims);
ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');
title('Signal received at GS 1 of LMS Adaptive Filter');


figure();
plot(real(E), imag(E),'.');
title('Error signal in LMS Adaptive filter');
xlim(plot_lims);
ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');

figure();
plot(real(yd),imag(yd),'.');
title('Adaptive Desired output / Estimate signal');
xlim(plot_lims);
ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');

%% Demodulate the received QPSK signals%%
% QPSK demodulator at the Receiver
uncoded_bits_rx_AF = zeros(1,2*size(yd,2)); % later change 1 to Ns
for ns= 1:1
    B4_AF(ns,:) = (real(x_signal_recovered(ns,:))<0);
    B3_AF(ns,:) = (imag(x_signal_recovered(ns,:))<0);          
    uncoded_bits_rx_AF(ns,(1:2:end)) = B4_AF(ns,:);
    uncoded_bits_rx_AF(ns,(2:2:end)) = B3_AF(ns,:);
end


%% Calculate the BER of the recived signal from satellite 1 %%
 T_Errors_AF = zeros();
 T_bits_AF = zeros();
 for ns = 1:1
     while T_Errors_AF < 500
         % Calculate Bit Errors
         diff(ns,:) = uncoded_bits(ns,:) - uncoded_bits_rx_AF(ns,:);
         T_Errors_AF = T_Errors_AF + sum(abs(diff(ns,:)));
         T_bits_AF = T_bits_AF + length(x_signal(2));
     end     
     % Calculate Bit Error Rate
     BER_AF(ns,:) = T_Errors_AF / T_bits_AF;     
 end
