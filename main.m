
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
bit_num = 10^3;
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
bit_error_rate_sim = zeros(1,length(SNR));
for j = 1:length(SNR)
    received_signal = awgn(signal,SNR(j));
    % plot the received signal for 5 T time
    
    %figure;
    %plot(t_signal(1:5*length(t)),received_signal(1:5*length(t)));
    %xlabel('time');
    %ylabel('received signal');
    %title('received signal');
    %grid on;

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
        correlation = sum(received_signal_i .* (s1(t) - s2(t)) * Ts);
        if correlation > comporator_gama
            received_bits(i) = 1;
        else
            received_bits(i) = 0;
        end
    end

    % calculate the bit error rate
    bit_error_rate = sum(abs(received_bits - bits))/bit_num;
    %display bit error for each SNR write SNR and bit error rate
    display(['SNR = ',num2str(SNRdb(j)),' bit error rate = ',num2str(bit_error_rate)]);
    bit_error_rate_sim(j) = bit_error_rate;



end

% plot the bit error rate simulation
figure;
semilogy(SNRdb,bit_error_rate_sim);
xlabel('SNR (dB)');
ylabel('bit error rate');
title('SNR vs bit error rate');
grid on;

%% print the bit error from theorical and simulation
% the theorical bit error rate is Pb
% the simulation bit error rate is bit_error_rate_sim

figure;
semilogy(SNRdb,Pb);
hold on;
semilogy(SNRdb,bit_error_rate_sim);
xlabel('SNR (dB)');
ylabel('bit error rate');
title('SNR vs bit error rate');
legend('theorical','simulation');
grid on;

