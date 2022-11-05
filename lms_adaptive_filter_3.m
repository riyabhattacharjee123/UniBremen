%%lms_adaptive_filter_3.m%%
%% This program estimates the channel matrix based on LMS method%%

n_seq = size(x_signal,2); % length of sequence
m_taps = 1;  % number of taps
mu_step = 0.1 ; % step size mu
t=[0:n_seq-1];
Weight = zeros(Ns,n_seq); % weight vector
E=[]
threshold = 0.05;

% Desired signal
D_signal = s_signal ;

% Corrupted received_signal
A_signal = (Channel_matrix_H)*(G_geo)*(D_signal) + awgn_cn ;

%% Find Optimum W. Used this code%%
W_af = zeros(Ngs,Ns); % Ngs=64
error=zeros(Ns,n_seq); %Ns = 3
WA = [];
WB=[];
x_estimated=zeros(Ns,n_seq);

%% AF code
for ins = 1:Ns
    for ix = 1:n_seq
        for ings = 1:Ngs
            E(ins,ix) = D_signal(ins,ix) - Weight(ins,ix)' ...
                                           * A_signal(ings,ix);
            
            Weight(ins,ix+1) = Weight(ins,ix) + mu_step ...
                                                * E(ins,ix) ...
                                                * A_signal(ings,ix);  
        end
    end
    for ix = 1:n_seq
        for ings = 1:Ngs
            D_signal_est(ins,ix) = sum(Weight(ins,ix+1)' ...
                                       * A_signal(ings,ix)) ;
        end
    end
end

%% Calculate the estimated signal in Adaptive Filter %%
%s_est_signal_af = (W_af)' * (A_signal) ;


%figure();plot(t,D_signal(1,:));title('Desired Signal');
%figure();
%plot(real(A_signal(1,:)),imag(A_signal(1,:)),'.');
%xlim(plot_lims);
%ylim(plot_lims);
%xlabel('real part');
%ylabel('imaginary part');
%title('Signal received at GS 1 of LMS Adaptive Filter');

%figure();
%plot(real(E(1,:)), imag(E(1,:)),'.');
%title('Error signal in LMS Adaptive filter');
%xlim(plot_lims);
%ylim(plot_lims);
%xlabel('real part');
%ylabel('imaginary part');

%figure();
%plot(real(D_signal_est(1,:)),imag(D_signal_est(1,:)),'.');
%title('Estimate signal');
%xlim(plot_lims);
%ylim(plot_lims);
%xlabel('real part');
%ylabel('imaginary part');



%% Demodulate the received QPSK signals%%
% QPSK demodulator at the Receiver
uncoded_bits_rx_AF = zeros(1,2*size(A_signal,2)); % later change 1 to Ns
for ns= 1:1
    B4_AF(ns,:) = (real(D_signal_est(ns,:))<0);
    B3_AF(ns,:) = (imag(D_signal_est(ns,:))<0);          
    uncoded_bits_rx_AF(ns,(1:2:end)) = B4_AF(ns,:);
    uncoded_bits_rx_AF(ns,(2:2:end)) = B3_AF(ns,:);
end


%% Calculate the BER of the recived signal from satellite 1 %%
 T_Errors_AF = zeros();
 T_bits_AF = zeros();
 for ns = 1:3
     while T_Errors_AF < 500
         % Calculate Bit Errors
         diff(ns,:) = uncoded_bits(ns,:) - uncoded_bits_rx_AF(ns,:);
         T_Errors_AF = T_Errors_AF + sum(abs(diff(ns,:)));
         T_bits_AF = T_bits_AF + length(x_signal(2));
     end     
     % Calculate Bit Error Rate
     BER_AF(ns,:) = T_Errors_AF / T_bits_AF;     
 end
