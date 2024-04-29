%% Theroical simulation of the binary signal noise ratio
Eb = 1;
% make a vector of SNR values from 0 to 11 dB
SNRdb = 0:11;
%make Pb vector using SNRdb
Pb =qfunc(sqrt(2*10.^(SNRdb/10)));
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb);
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;

