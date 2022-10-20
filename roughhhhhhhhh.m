
close all
clear all
clc
Ns = 3
num_symbols = 1e3
%uncoded_bits = zeros()
uncoded_bits = round(rand(Ns,num_symbols));
for ns = 1:3
     % Split the stream into two streams, for Quadrature Carriers
     B1 = uncoded_bits(ns,(1:2:end));
     B2 = uncoded_bits(ns,(2:2:end));
        
     % QPSK modulator set to pi/4 radians constellation
     % If you want to change the constellation angles
     % just change the angles. (Gray Coding)
     qpsk_sig(ns,:) = ((B1==0).*(B2==0)*(exp(i*pi/4))+(B1==0).*(B2==1)...
                *(exp(3*i*pi/4))+(B1==1).*(B2==1)*(exp(5*i*pi/4))...
                +(B1==1).*(B2==0)*(exp(7*i*pi/4))); 

     tx_sig(ns,:) = qpsk_sig(ns,:);   

end