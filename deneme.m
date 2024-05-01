%% Theroical simulation of the binary signal noise ratio
Eb = 1;
% make a vector of SNR values from 0 to 11 dB
SNRdb = 0:0.01:11;
%make Pb vector using SNRdb
Pb =qfunc(sqrt(2*10.^(SNRdb/10)));
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb);
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;
%you can see that the Pb decreases as the SNR increases as example at 10 db the Pb is 10^-6
% its mean there is 1 error in 1 million bits
% To test that we have to send 10 million bits
% and check the number of errors bu sending 100 million bits is better

% make random bits 100 million
bits = randi([0 1],1,10^8);

% plot the distribution of the bits
figure;
histogram(bits,100);
xlabel('bits');
ylabel('count');
title('bits distribution');
grid on;

% make the signal
signal = 2*bits - 1;
figure;
plot(signal(1:10));
xlabel('time');
ylabel('signal');
title('signal');
grid on;

% make the noise
noise = randn(1,10^8);
figure;
histogram(noise,100);
xlabel('noise');
ylabel('count');
title('noise distribution');
grid on;


% make the received signal
SNRdb = 10;
SNR = 10^(SNRdb/10);
received_signal = signal + noise/sqrt(SNR);
figure;
plot(received_signal(1:10));
xlabel('time');
ylabel('received signal');
title('received signal');
grid on;




