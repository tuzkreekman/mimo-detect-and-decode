global LEN ITERS snr K R N qamBitSize qamSize normAnt normConst precode

LEN = 2; % how many K-bit length messages we will send (per tx/rx)
ITERS = 10000;
snr = 0:10;
n = 2; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;

addpath('./samples/polar');
addpath('./samples/polar/functions');

initPC(N,K,'AWGN',0); % changd snr
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

subplot(2,1,1);
hold on;
perfectKnowledge = 1;

n=2;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.r','DisplayName','2x2 Linear');
plot(snr,polarBER,':r','DisplayName','2x2 Polar Only');

n=4;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.g','DisplayName','4x4 Linear');
plot(snr,polarBER,':g','DisplayName','4x4 Polar Only');

n=8;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.b','DisplayName','8x8 Linear');
plot(snr,polarBER,':b','DisplayName','8x8 Polar Only');

title('BER vs SNR for MIMO Systems, K=16, Perfect Channel Knowledge');
ylabel('BER');
xlabel('SNR');
set(gca,'YScale','log');
legend(gca,'show');

subplot(2,1,2);
perfectKnowledge = 0;
hold on;

n=2;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.r','DisplayName','2x2 Linear');
plot(snr,polarBER,':r','DisplayName','2x2 Polar Only');

n=4;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.g','DisplayName','4x4 Linear');
plot(snr,polarBER,':g','DisplayName','4x4 Polar Only');

n=8;
[qamBER, linearBER, polarBER] = GenerateBER(n, perfectKnowledge);
plot(snr,linearBER,'-.b','DisplayName','8x8 Linear');
plot(snr,polarBER,':b','DisplayName','8x8 Polar Only');

title('BER vs SNR for MIMO Systems, K=16, Estimated Channel');
ylabel('BER');
xlabel('SNR');
set(gca,'YScale','log');
legend(gca,'show');
hold off;




saveas(gca,'baselineBER.png');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [qamBER,combinedBER,polarBER] = GenerateBER(n,perfectKnowledge)
global LEN ITERS snr K R N qamBitSize qamSize normAnt normConst precode

qamBER       = zeros(11,ITERS);
combinedBER  = zeros(11,ITERS);
polarBER     = zeros(11,ITERS);

for (SNR=snr)
  for (i=1:ITERS)
    [ber0,ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,perfectKnowledge,precode);
    qamBER(SNR+1,i) = ber0;
    combinedBER(SNR+1,i) = ber1;
    polarBER(SNR+1,i) = ber2;
  end
end

qamBER = mean(qamBER,2);
combinedBER = mean(combinedBER,2);
polarBER = mean(polarBER,2);

end
