%%lms_adaptive_filter.m%%
%% This program estimates the channel matrix based on LMS method%%

n_seq = size(qpsk_sig,2); % length of sequence
m_taps = 1;  % number of taps
mu_step = 0.0005; %3.6896e+27; %0.09 ; % step size mu
t=[0:n_seq-1];
n_symb = [1:num_symbols]
Weight_real = zeros(Ns,n_seq); % weight vector Real part
Weight_imag = zeros(Ns,n_seq); % weight vector Imaginary part
E_real=[]; % error vector real part
E_imag=[]; % error vector Imaginary part
%threshold = 0.05;

% Desired signal
D_signal = qpsk_sig ;
D_signal_real = real(qpsk_sig);
D_signal_imag = imag(qpsk_sig);

% Corrupted received_signal
A_signal = (Channel_matrix_H)*(G_geo)*(D_signal) + awgn_cn ; % input vector
A_signal_amplified = (1e3)*A_signal;
U_signal_real = real(A_signal_amplified); % input vector real part
U_signal_imag = imag(A_signal_amplified); % input vector imaginary part



%% Find Optimum W. Used this code%%
%W_af = zeros(Ngs,Ns); % Ngs=64
%error=zeros(Ns,n_seq); %Ns = 3
%WA = [];
%x_estimated=zeros(Ns,n_seq);

%% LMS adaptive filter for complex signals %%
for ins = 1:Ns   
    for ix = 1:n_seq
        for ings = 1:Ngs
            % Calculate output Real part
            y_real(ins,ix) = sum(Weight_real(ins,ix)' ...
                            * U_signal_real(ings,ix)) ...
                            - sum(Weight_imag(ins,ix)' ...
                            * U_signal_imag(ings,ix)) ; 
            
            % Calculate output Imag part            
            y_imag(ins,ix) = sum(Weight_real(ins,ix)' ...
                            * U_signal_imag(ings,ix)) ...
                            + sum(Weight_imag(ins,ix)' ...
                            * U_signal_real(ings,ix)) ;
                        
            % Calculate error vector Real and imaginary parts            
            E_real(ins,ix) = D_signal_real(ins,ix) - y_real(ins,ix);
            E_imag(ins,ix) = D_signal_imag(ins,ix) - y_imag(ins,ix);
            E(ins,ix) = complex(E_real(ins,ix),E_imag(ins,ix));
            
            % Calculate Weight vector Real and imaginary parts
            Weight_real(ins,ix+1) = Weight_real(ins,ix) + mu_step ...
                                    *( E_real(ins,ix) ...
                                    * U_signal_real(ings,ix) ...
                                    - E_imag(ins,ix) ...
                                    * U_signal_imag(ings,ix));
                                
            Weight_imag(ins,ix+1) = Weight_imag(ins,ix) + mu_step ...
                                    *( E_real(ins,ix) ...
                                    * U_signal_imag(ings,ix) ...
                                    - E_imag(ins,ix) ...
                                    * U_signal_real(ings,ix));
        end
    end
end

W_AF_real = zeros(Ns,Ngs);
W_AF_imag = zeros(Ns,Ngs);
for ins = 1:Ns  
    for ings = 1:Ngs
        W_AF_real(ins,ings) = Weight_real(ins,n_seq+1); 
        W_AF_imag(ins,ings) = Weight_imag(ins,n_seq+1); 
    end
end

s_est_signal_af_real = zeros(Ns,n_seq);
s_est_signal_af_imag = zeros(Ns,n_seq);

%% Calculate the estimated signal in Adaptive Filter %%
for ins = 1:Ns   
    for ix = 1:n_seq
        for ings = 1:Ngs
            %s_est_signal_af() = (W_af)' * (A_signal) ;
            s_est_signal_af_real(ins,ix) = sum(W_AF_real(ins,ings) ...
                            * U_signal_real(ings,ix)) ...
                            - sum(W_AF_imag(ins,ings) ...
                            * U_signal_imag(ings,ix)) ; 
            
            % Calculate output Imag part            
            s_est_signal_af_imag(ins,ix) = sum(W_AF_real(ins,ings) ...
                            * U_signal_imag(ings,ix)) ...
                            + sum(W_AF_imag(ins,ings) ...
                            * U_signal_real(ings,ix)) ;
        end
    end
end

figure();plot(t,real(D_signal(1,:)));title('Desired Signal');
figure();
%plot(real(A_signal(1,:)),imag(A_signal(1,:)),'.');
plot(t,real(A_signal(1,:)));
%xlim(plot_lims);
%ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');
title('Signal received at GS 1 of LMS Adaptive Filter');

figure();
%plot(real(A_signal(1,:)),imag(A_signal(1,:)),'.');
plot(t,real(A_signal_amplified(1,:)));
%xlim(plot_lims);
%ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');
title('Amplified signal received at GS 1 of LMS Adaptive Filter');


figure();
%plot(real(E(1,:)), imag(E(1,:)),'.');
plot(t,real(E_real(1,:)));
title('Error signal in LMS Adaptive filter');
%xlim(plot_lims);
%ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');

figure();
%plot(real(D_signal_est(1,:)),imag(D_signal_est(1,:)),'.');
plot(t,real(y_real(1,:)));
title('Adaptive desired output signal');
%xlim(plot_lims);
%ylim(plot_lims);
xlabel('real part');
ylabel('                                                                                                                                                     imaginary part');

figure();
%plot(real(D_signal_est(1,:)),imag(D_signal_est(1,:)),'.');
plot(t,s_est_signal_af_real(1,:));
title('Estimated output signal Adaptive Filter');
%xlim(plot_lims);
%ylim(plot_lims);
xlabel('real part');
ylabel('imaginary part');


%% Demodulate the received QPSK signals%%
% QPSK demodulator at the Receiver
uncoded_bits_rx_AF = zeros(1,2*size(A_signal,2)); % later change 1 to Ns
for ins= 1:Ns
    B4_AF(ins,:) = (s_est_signal_af_real(ins,:)<0);%(y_real(ins,:)<0);
    B3_AF(ins,:) = (s_est_signal_af_imag(ins,:)<0);%(y_imag(ins,:)<0);          
    uncoded_bits_rx_AF(ins,(1:2:end)) = B4_AF(ins,:);
    uncoded_bits_rx_AF(ins,(2:2:end)) = B3_AF(ins,:);
end

qpsk_sig_demodulated_AF = uncoded_bits_rx_AF;

%% Calculate the BER of the recived signal from satellite 1 %%
 T_Errors_AF = zeros(Ns,num_symbols);
 T_bits_AF = zeros(Ns,num_symbols);
 for ins = 1:Ns
     while T_Errors_AF(ins,:) < 500
         % Calculate Bit Errors
         diff(ins,:) = uncoded_bits(ins,:) - uncoded_bits_rx_AF(ins,:);
         T_Errors_AF(ins,:) = T_Errors_AF(ins,:) + sum(abs(diff(ins,:)));
         T_bits_AF(ins,:) = T_bits_AF(ins,:) + length(x_signal(2));
     end     
     % Calculate Bit Error Rate
     BER_AF(ins,:) = T_Errors_AF(ins,:) / T_bits_AF(ins,:);     
 end

figure();
plot_lims = [-3 3];
subplot(2,2,1)
plot(y_real(1,:), y_imag(1,:), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 1 demodulated AF');
xlabel('real part');
ylabel('imaginary part');% figure();

subplot(2,2,2)
plot(y_real(2,:), y_imag(2,:), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 2 demodulated AF');
xlabel('real part');
ylabel('imaginary part');

subplot(2,2,3)
plot(y_real(3,:), y_imag(3,:), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 3 demodulated AF');
xlabel('real part');
ylabel('imaginary part');

figure();
subplot(2,2,1)
plot(real(qpsk_sig_demodulated_AF(1,:)), imag(qpsk_sig_demodulated_AF(1,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 1 demodulated');
xlabel('real part');
ylabel('imaginary part');

subplot(2,2,2)
plot(real(qpsk_sig_demodulated_AF(2,:)), imag(qpsk_sig_demodulated_AF(2,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 2 demodulated');
xlabel('real part');
ylabel('imaginary part');

subplot(2,2,3)
plot(real(qpsk_sig_demodulated_AF(3,:)), imag(qpsk_sig_demodulated_AF(3,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
title('QPSK constellation from satelite 3 demodulated');
xlabel('real part');
ylabel('imaginary part');

figure();
plot(t,abs(E(1,:)) );
xlabel("Number of adaptations");
ylabel("Instantaneous Error");
title("Instantaneous Error signal with iterations");

figure();
plot((n_symb), abs(diff(1,:)) );
xlabel("Number of symbols");
ylabel("Instantaneous Bit Error");
title("Instantaneous BER with iterations");

