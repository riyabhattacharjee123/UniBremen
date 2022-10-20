%%%% input_signal_qpsk.m%%%

%x = randi([0 1],1,8) ; % input sequence

nr_data_bits=50;
M = nr_data_bits
nr_symbols=nr_data_bits/2;
b_data=(randn(1,nr_data_bits)>.5);
b=[b_data];
x=[b_data];

% Bits to polar mapping
for pidx = 1:length(x)
    if x(pidx) ==0
        p(pidx)=-1;
    else
        p(pidx)=+1;
    end    
end

figure();
stem(x(1,:));
title('Input data binary');
xlabel('Samples');
ylabel('Amplitude');

% Separation of even and odd sequences
even_seq = p(1:2:length(x));
odd_seq = p(2:2:length(x));

% NRZ polar line coded signal generation
pidx = 1;
mx = 2:2:length(x);
t=0:0.01:length(x);
for jx = 1:length(t)
    if t(jx)<=mx(pidx)
        even_ps(jx) = even_seq(pidx);
    else
        even_ps(jx) = even_seq(pidx);
        pidx=pidx+1;
    end   
end


pidx = 1;
mx = 2:2:length(x);
for jx = 1:length(t)
    if t(jx)<=mx(pidx)
        odd_ps(jx) = odd_seq(pidx);
    else
        odd_ps(jx) = odd_seq(pidx);
        pidx=pidx+1;
    end   
end

%figure(1);
%subplot(2,1,1);
%plot(t,even_ps,'r');
%subplot(2,1,2);
%plot(t,odd_ps,'k');

% quadrature carrier signal generator
c1 = cos(2*pi*1*t);
c2 = sin(2*pi*1*t);

%figure(2);
%subplot(2,1,1);
%plot(t,c1,'r');
%subplot(2,1,2);
%plot(t,c2,'k');

% QPSK waveform generator
r1 = even_ps.*c1;
r2 = odd_ps.*c2;
qpsk_sig = r1-r2;

figure(2);
%subplot(3,1,1);
%plot(t,r1,'r');
%subplot(3,1,2);
%plot(t,r2,'k');
%subplot(3,1,3);
plot(t,qpsk_sig,'g');
title('QPSK Waveform');
xlabel('Samples');
ylabel('Amplitude');

% QPSK constellation diagram generator
d=zeros(1,length(b)/2);

for n=1:length(b)/2;
    p=b(2*n);
    imp=b(2*n-1);
    if (imp==0)&(p==0)
        d(n)=exp(j*pi/4);
    end
    if(imp==1)&(p==0)
        d(n)=exp(j*3*pi/4);
    end
    if (imp==1)&(p==1)
        d(n)=exp(j*5*pi/4);
    end
    if (imp==0)&(p==1)
        d(n)=exp(j*7*pi/4);
    end
end
qpsk=d;
figure(3);
plot(d,'o');
axis([-2 2 -2 2]);
grid on;
xlabel('real');
ylabel('imag');
title('QPSK constellation')

