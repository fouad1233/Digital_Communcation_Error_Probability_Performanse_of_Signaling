
A = 1;
T = 2*10^-3;

s1 = @(t) A.*sin((2*pi*t)/T) .* (0<=t & t<=T/2);
s2 = @(t) -s1(t- T/2) .* (T/2<=t & t<=T);
Ts = 0.00001;
t = 0:Ts:T-Ts;
fs = 1/Ts;
figure;
% make subplot for the 2 signals
subplot(2,1,1);
plot(t,s1(t));
xlabel('time');
ylabel('s1(t)');
title('s1(t)');
grid on;

subplot(2,1,2);
plot(t,s2(t));
xlabel('time');
ylabel('s2(t)');
title('s2(t)');
grid on;

Eb = A^2*T/4;
% make a vector of SNR values from 0 to 11 dB
SNRdb = 0:0.01:11;
Pb =qfunc(sqrt(10.^(SNRdb/10)));

%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb);
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;

% make random bits 100 thousand
bit_num = 10^6;
bits = randi([0 1],1,bit_num);
% if the bit is 1 then the signal is s1(t) if 0 the signal is s2(t)
signal = zeros(1,bit_num * length(t));
for i = 1:bit_num
    if bits(i) == 1
        signal((i-1)*length(t)+1:i*length(t)) = s1(t);
    else
        signal((i-1)*length(t)+1:i*length(t)) = s2(t);
    end
end
t_signal = 0:Ts:(bit_num*T);
%plot the signal for 5 T time
figure;
plot(t_signal(1:5*length(t)),signal(1:5*length(t)));
xlabel('time');
ylabel('signal');
title('signal');
grid on;

% add additive white gaussion noise to the signal for all SNR db values
SNRdb = 0:0.01:11;
% make the noise using awgn function
SNR = 10.^(SNRdb/10);
received_signal = awgn(signal,1/1000000,'measured');
% plot the received signal for 5 T time
figure;
plot(t_signal(1:5*length(t)),received_signal(1:5*length(t)));
xlabel('time');
ylabel('received signal');
title('received signal');
grid on;

% make the decision for the received signal and obtain received bits
received_bits = zeros(1,bit_num);
% we will use correlation to decide the received bits
% r(t) * (s1(t) - s2(t))  integral from 0 to T will be positive if the bit is 1
% and negative if the bit is 0

comporator_gama = 0;
for i = 1:bit_num
    % get the received signal for the i th bit
    received_signal_i = received_signal((i-1)*length(t)+1:i*length(t));
    % calculate the correlation
    correlation = sum(received_signal_i .* (s1(t) - s2(t)) .* Ts ) ;
    if correlation > comporator_gama
        received_bits(i) = 1;
    else
        received_bits(i) = 0;
    end
end
% plot received bits and original bits
figure;
plot(bits(1:100));
hold on;
plot(received_bits(1:100));
xlabel('bit index');
ylabel('bit value');
title('received bits vs original bits');
grid on;
legend('original bits','received bits');


% calculate the bit error rate
bit_error_rate = sum(abs(received_bits - bits))/bit_num;
% print sum of received bits and original bits
disp(['sum of (received bits - original bits) = ',num2str(sum(received_bits - bits))]);


disp(['bit error rate = ',num2str(bit_error_rate)]);


