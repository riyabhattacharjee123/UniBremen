close all;
clear all;
clc;

L=41; %filter Length
R=2.5e9; %data Rate = 2.5Gbps
Fs=8*R; %oversampling by 8
T=1/R; %pulse duration
Ts=1/Fs; %sampling duration
alpha =1; % Design Factor for Raised Cosing Filter
upsampling = Fs/R;

%----------------------------------------------------------
%Raised Cosing Filter Design
%----------------------------------------------------------
if mod(L,2)==0 
 M=L/2 ; % for even value of L
else
 M=(L-1)/2; % for odd value of L 
end 
g=zeros(1,L); %Place holder for RC filter's transfer function 
for n=-M:M 
 num=sin(pi*n*Ts/T)*cos(alpha*pi*n*Ts/T);
 den=(pi*n*Ts/T)*(1-(2*alpha*n*Ts/T)^2);
 g(n+M+1)=num/den; 
 if (1-(2*alpha*n*Ts/T)^2)==0 
 g(n+M+1)=pi/4*sin(pi*n*Ts/T)/(pi*n*Ts/T);
 end 
 if (pi*n*Ts/T)==0 
 g(n+M+1)=cos(alpha*pi*n*Ts/T)/(1-(2*alpha*n*Ts/T)^2); 
 end 
end 
%----------------------------------------------------------
%Generate data of random 1s and 0s
data1 = round(rand(1,2^7)); % Generate a random bit stream
%data2 = (-2*data1) + 1; %generate [-1,1] for BPSK
output1=upsample(data1,Fs/R);
%output2=upsample(data2,Fs/R);
%y=conv(g,output); %Convolving the data signal with the Raised Cosine Filter
y1=filter(g,1,output1); %you can use either Conv function or filter function to obtain the output
%y2=filter(g,1,output2);
%----------------------------------------------------------
%Plot data and RC filtered Output
figure;
subplot(2,1,1);
stem(data1);
title('Input data to the Raised Cosine Filter');
xlabel('Samples');
ylabel('Amplitude');
subplot(2,1,2);
plot(y1);
title('Response of the Raised Cosine Filter for the given Input');
xlabel('Samples');
ylabel('Amplitude');
%Plot data and RC filtered Output
%figure;
%subplot(2,1,1);
%stem(data2);
%title('Input data to the Raised Cosine Filter');
%xlabel('Samples');
%ylabel('Amplitude');
%subplot(2,1,2);
%plot(y2);
%title('Response of the Raised Cosine Filter for the given Input');
%xlabel('Samples');
%ylabel('Amplitude');
%Obtain FFT of the output to plot its frequency response.
Fn=Fs/2;NFFY=2.^(ceil(log(length(y1))/log(2)));
FFTY=fft(y1,NFFY);%pad with zeros
NumUniquePts=ceil((NFFY+1)/2); 
FFTY=FFTY(1:NumUniquePts);
MY=abs(FFTY);
MY=MY*2;
MY(1)=MY(1)/2;
MY(length(MY))=MY(length(MY))/2;
MY=MY/length(y1);
f1=(0:NumUniquePts-1)*2*Fn/NFFY;
%Plot Frequency spectrum
figure;
plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(Hz)');ylabel('DB');
grid
title('Frequency domain plots') 
%Obtain FFT of the output to plot its frequency response.
%Fn=Fs/2;NFFY=2.^(ceil(log(length(y2))/log(2)));
%FFTY=fft(y2,NFFY);%pad with zeros
%NumUniquePts=ceil((NFFY+1)/2); 
%FFTY=FFTY(1:NumUniquePts);
%MY=abs(FFTY);
%MY=MY*2;
%MY(1)=MY(1)/2;
%MY(length(MY))=MY(length(MY))/2;
%MY=MY/length(y2);
%f1=(0:NumUniquePts-1)*2*Fn/NFFY;
%Plot Frequency spectrum
%figure;
%plot(f1,20*log10(abs(MY).^2));xlabel('FREQUENCY(Hz)');ylabel('DB');
%grid
%title('Frequency domain plots') 
% Plot Eye Diagram
%ploteye(y1,Fs/R);
%ploteye(y2,Fs/R);