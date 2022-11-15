%%% QPSK_symbols_complex.m %%%

num_symbols = 1e3; %1e4; %bit count, frame rate
M = num_symbols
%amplitude = 1
%for ns = 1:3
 %   int_symbols(ns,:) = randi([1, 4], 1, num_symbols);
  %  A = sqrt(amplitude);
   % qpsk_symbols(ns,:) = zeros(size(int_symbols(ns,:)));
    %qpsk_symbols(int_symbols == 1) =   A + 1i*A;
    %qpsk_symbols(int_symbols == 2) =   A - 1i*A;
    %qpsk_symbols(int_symbols == 3) = - A + 1i*A;
    %qpsk_symbols(int_symbols == 4) = - A - 1i*A;
    %tx_sig(ns,:) = qpsk_symbols(ns,:);    
%end

% Generate some information bits
uncoded_bits = round(rand(Ns,num_symbols));
for ns = 1:Ns
    % Split the stream into two streams, for Quadrature Carriers
    B1(ns,:) = uncoded_bits(ns,(1:2:end));
    B2(ns,:) = uncoded_bits(ns,(2:2:end));
    
    % QPSK modulator set to pi/4 radians constellation
     % If you want to change the constellation angles
     % just change the angles. (Gray Coding)
    qpsk_sig(ns,:) = ((B1(ns,:)==0).*(B2(ns,:)==0)*(exp(i*pi/4)) ...
                        +(B1(ns,:)==0).*(B2(ns,:)==1)*(exp(3*i*pi/4))...
                        +(B1(ns,:)==1).*(B2(ns,:)==1)*(exp(5*i*pi/4))...
                        +(B1(ns,:)==1).*(B2(ns,:)==0)*(exp(7*i*pi/4))); 

    tx_sig(ns,:) = qpsk_sig(ns,:);   % transmitted signal

end


%figure();
%plot_lims = [-3 3];
%subplot(2,2,1)
%plot(real(qpsk_sig(1,:)), imag(qpsk_sig(1,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 1 without noise');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,2)
%plot(real(qpsk_sig(2,:)), imag(qpsk_sig(2,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 2 without noise');
%xlabel('real part');
%ylabel('imaginary part');

%subplot(2,2,3)
%plot(real(qpsk_sig(3,:)), imag(qpsk_sig(3,:)), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title('QPSK constellation from satelite 3 without noise');
%xlabel('real part');
%ylabel('imaginary part');

%snr = 15; %// in dB
%rx_sig = awgn(tx_sig, snr, 'measured');

%fh2 = figure;
%plot(real(rx_sig), imag(rx_sig), '.');
%xlim(plot_lims);
%ylim(plot_lims);
%title(['QPSK constellation at an SNR of ' num2str(snr) ' dB']);
%xlabel('real part');
%ylabel('imaginary part');