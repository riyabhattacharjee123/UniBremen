clear all;
 close all;
 l=10000;
 snrdb=1:1:10;
 snrlin=10.^(snrdb/10);
 for snrdb=1:1:10
     si=2*(round(rand(1,l))-0.5);
     sq=2*(round(rand(1,l))-0.5);
     s=si+j*sq;
     w=awgn(s,snrdb,'measured');
     r=w;
     si_=sign(real(r));
     sq_=sign(imag(r));
     ber1=(l-sum(si==si_))/l;
     ber2=(l-sum(sq==sq_))/l;
     ber(snrdb)=mean([ber1 ber2]);
 end
 
figure();
plot_lims = [-3 3];
subplot(2,2,1)
plot(real(s), imag(s), '.');
xlim(plot_lims);
ylim(plot_lims);
title('QPSK constellation from satelite 1 without noise');
xlabel('real part');
ylabel('imaginary part');
 
figure();
plot(real(w), imag(w), '.');



figure();
 %semilogy(snrdb, ber,'o-')
 snrdb=1:1:10;
 snrlin=10.^(snrdb./10);
 tber=0.5.*erfc(sqrt(snrlin));
 semilogy(snrdb,ber,'-bo',snrdb,tber,'-mh')
 title('QPSK with awgn');
 xlabel('Signal to noise ratio');
 ylabel('Bit error rate');
 grid on;
 
 
 
 Mod = 4;
N = 256;
x_n = randi([0 Mod-1],N,1);

P = 2; % Set desired constellation power
s_n = sqrt(P) * pskmod(x_n,Mod,pi/Mod);
figure()
plot(real(s_n), imag(s_n));