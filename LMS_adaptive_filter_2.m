%%lms_adaptive_filter_2.m%%
n_seq = size(x_signal_normalized,2); % length of sequence
l_path = size(x_signal_normalized,1); % number of paths
m_taps = 1; % number of taps
mu_step = 0.0005 ; % step size mu
t=[0:n_seq-1];
Weight = zeros(m_taps,n_seq); % weight vector

% Desired signal
D_signal = x_signal ;

% noise
awgn = sqrt(noiseVar) * (randn(Nr,n_seq) + (1i * randn(Nr,n_seq)));

% Corrupted received_signal
A_signal =D_signal + awgn ; %D_signal + awgn ; %y_received_signal_sat ;x_signal_recovered
E = []; % error signal

for iy=1:l_path
    for ix=(m_taps+1):n_seq        
        a_array(iy,:) = A_signal((ix-(m_taps)+1):ix)
    end
end
for iy=1:l_path
    for ix=(m_taps+1):n_seq
        E(iy,ix) = D_signal(iy,ix) - ...
                                     transpose(conj(Weight(iy,ix)))...
                                                .*a_array(iy,:) ;
    end
end
for iy=1:l_path
    for ix=(m_taps+1):n_seq
        Weight(iy,ix+1) = Weight(iy,ix) + mu_step * ...
                                        E(iy,ix).* ...
                                        (a_array(iy,:))';
    end
end
    
for iy=1:l_path
    for ix=(m_taps+1):n_seq
        a_array(iy,:) = A_signal((ix-(m_taps)+1):ix)
        yd(iy,ix) = a_array(iy,ix).*Weight(iy,ix);  
    end       
end


%figure();plot(t,D_signal(1,:));title('Desired Signal');
figure();plot(t,A_signal(1,:));title('Input Signal+Noise');
figure();plot(t,E(1,:));title('Error');
figure();plot(t,yd(1,:));title('Adaptive Desired output / Estimate signal');

