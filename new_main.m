
A = 2;
T = 1;

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
SNRdb = 0:1:14;
Pb =qfunc(sqrt(10.^(SNRdb/10)));
%%
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb);

ylim([10^-6,5*10^-1])
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb teorical for P(1)=P(0)=0.5');
grid on;
%%
% make random bits 100 thousand
bit_num = 10^8;
bits = randi([0 1],1,bit_num);
%calculate a1 and a2
a1 = A^2*T/4;
a2 = -a1;
%if the bit is 1 then z = a1 + n if 0 then z = a2 + n
% if bit in bits is 1 then ai(T) = a1 else ai(T) = a2
ai = zeros(1,bit_num);
ai(bits == 1) = a1;
ai(bits == 0) = a2;

%let's calculate z = ai + n by adding gaussian noise for different SNR
%values
z = zeros(1,bit_num);
shat = zeros(1,bit_num);
Pb_sim = zeros(1,length(SNRdb));
comparator_gama = 0;
comparator_gamas = comparator_gama * ones(1,bit_num);
%find N0's value for given SNR
N0 = Eb.*(10.^(-SNRdb/10));
for i = 1:length(N0)

    z = ai + (sqrt(N0(i)*(A^2)*T/4) ).*randn(1,bit_num);
    %find shat by using comparator
    
    shat = double(z>comparator_gamas);
    %calculate Pb_sim
    Pb_sim(i) = sum(abs(bits-shat))/bit_num;


end
%%
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb_sim);

ylim([10^-6,5*10^-1])
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb from simulation for P(1)=0.5 P(0)=0.5');
grid on;

%plot the two graphs
figure;
semilogy(SNRdb,Pb);
hold on;
semilogy(SNRdb,Pb_sim);

ylim([10^-6,5*10^-1])
legend('Pb from theory','Pb simulation')
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;
%%
%% For different Probobilities


A = 2;
T = 1;

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
SNRdb = 0:1:14;
%calculate a1 and a2
a1 = A^2*T/4;
a2 = -a1;
N0 = Eb.*(10.^(-SNRdb/10));

gama0 = 0.549.*N0;
sigma0 = A.*sqrt(N0.*T)/2;


Pb =0.25*(1-qfunc((gama0 - a1) ./sigma0))+...
    0.75*qfunc( (gama0 - a2) ./sigma0 ) ;
%%
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb);
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb from theoric for P(1)=0.25 P(0)=0.75');
grid on;
%%

% make random bits 100 thousand
bit_num = 10^8;
random_numbers = rand(1,bit_num);
P1 = 1/4;
bits = random_numbers< P1;
%calculate a1 and a2
a1 = A^2*T/4;
a2 = -a1;
%if the bit is 1 then z = a1 + n if 0 then z = a2 + n
% if bit in bits is 1 then ai(T) = a1 else ai(T) = a2
ai = zeros(1,bit_num);
ai(bits == 1) = a1;
ai(bits == 0) = a2;

%let's calculate z = ai + n by adding gaussian noise for different SNR
%values
z = zeros(1,bit_num);
shat = zeros(1,bit_num);
Pb_sim = zeros(1,length(SNRdb));
%find N0's value for given SNR
N0 = Eb.*(10.^(-SNRdb/10));

for i = 1:length(N0)
    comparator_gamas = 0.549*N0(i) * ones(1,bit_num);


    z = ai + (sqrt(N0(i)*(A^2)*T/4) ).*randn(1,bit_num);
    %find shat by using comparator
    
    shat = double(z>comparator_gamas);
    %calculate Pb_sim
    Pb_sim(i) = sum(abs(bits-shat))/bit_num;


end
%%
%plot SNRdb vs Pb
figure;
semilogy(SNRdb,Pb_sim);
ylim([10^-6,5*10^-1])
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb from simulation for P(1)=0.25 P(0)=0.75');
grid on;

%plot the two graphs
figure;
semilogy(SNRdb,Pb);
hold on;
semilogy(SNRdb,Pb_sim);
ylim([10^-6,5*10^-1])
legend('Pb from theory','Pb simulation')
xlabel('SNR (dB)');
ylabel('Pb');
title('SNR vs Pb');
grid on;



