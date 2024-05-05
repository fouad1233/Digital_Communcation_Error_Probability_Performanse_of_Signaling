A = 1;
T = 2*10^-3;

% Define functions properly
s1 = @(t) A.*sin((2*pi*t)/T) .* (0<=t & t<=T/2);
s2 = @(t) -s1(t- T/2) .* (T/2<=t & t<=T);

Ts = 0.00001;
t = 0:Ts:T-Ts;
fs = 1/Ts;

% Plot original signals
figure;
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

% Calculate SNR and Pb
Eb = A^2*T/4;
SNRdb = 0:0.01:11;
Pb = qfunc(sqrt(10.^(SNRdb/10)));

% Plot SNR vs Pb
figure;
semilogy(SNRdb,Pb);
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;

% Generate random bits
bit_num = 10^6;
bits = randi([0 1],1,bit_num);

% Generate signal
signal = zeros(1,bit_num * length(t));
for i = 1:bit_num
    if bits(i) == 1
        signal((i-1)*length(t)+1:i*length(t)) = s1(t);
    else
        signal((i-1)*length(t)+1:i*length(t)) = s2(t);
    end
end

% Add AWGN to signal
received_signal = awgn(signal,1);

%print 10 log10(Eb)
disp(['10 log10(Eb) = ',num2str(10*log10(Eb))]);
% Plot received signal
t_signal = 0:Ts:(bit_num*T);
figure;
plot(t_signal(1:5*length(t)),received_signal(1:5*length(t)));
xlabel('time');
ylabel('signal');
title('received signal');
grid on;


% Make the decision for the received signal and obtain received bits
received_bits = zeros(1,bit_num);
s1_minus_s2 = s1(t) - s2(t);
gamma = 0;
error_bits = 0 ;
for i = 1:bit_num
    received_signal_i = received_signal((i-1)*length(t)+1:i*length(t));
    correlation = sum(received_signal_i .* s1_minus_s2 .* Ts );
    error_bits = error_bits + (correlation>gamma) .* (1- bits(i)) + (correlation<=gamma) .* bits(i);
    if correlation > gamma
        received_bits(i) = 1;
    else
        received_bits(i) = 0;
    end

end
% print error_bits
disp(['error_bits = ',num2str(error_bits)]);

% Calculate bit error rate
bit_error_rate = sum(abs(received_bits - bits))/bit_num;

% Plot received bits and original bits
figure;
plot(bits(1:100));
hold on;
plot(received_bits(1:100));
xlabel('bit index');
ylabel('bit value');
title('received bits vs original bits');
grid on;
legend('original bits','received bits');

% Print bit error rate
disp(['Bit error rate = ',num2str(bit_error_rate)]);
