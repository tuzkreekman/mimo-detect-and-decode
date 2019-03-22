LEN = 4; % how many K-bit length messages we will send (per tx/rx)
ITERS = 1000;
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

linearBER    = zeros(11,ITERS);
iterativeBER = zeros(11,ITERS);
n=2;

for (SNR=snr)
  for (i=1:ITERS)
    [ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,0,precode);
    linearBER(SNR+1,i) = ber1;
    iterativeBER(SNR+1,i) = ber2;
  end
end

linearBER = mean(linearBER,2);
iterativeBER = mean(iterativeBER,2);

fig = figure();
hold on;
plot(snr,linearBER,'-.r','DisplayName','2x2 Linear');
plot(snr,iterativeBER,':r','DisplayName','2x2 Polar Only');

% 4x4 mimo

linearBER    = zeros(11,ITERS);
iterativeBER = zeros(11,ITERS);
n=4;

for (SNR=snr)
  for (i=1:ITERS)
    [ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,0,precode);
    linearBER(SNR+1,i) = ber1;
    iterativeBER(SNR+1,i) = ber2;
  end
end

linearBER = mean(linearBER,2);
iterativeBER = mean(iterativeBER,2);

plot(snr,linearBER,'-.g','DisplayName','4x4 Linear');
plot(snr,iterativeBER,':g','DisplayName','4x4 Polar Only');

% 8x8 mimo

linearBER    = zeros(11,ITERS);
iterativeBER = zeros(11,ITERS);
n=8;

for (SNR=snr)
  for (i=1:ITERS)
    [ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,0,precode);
    linearBER(SNR+1,i) = ber1;
    iterativeBER(SNR+1,i) = ber2;
  end
end

linearBER = mean(linearBER,2);
iterativeBER = mean(iterativeBER,2);

plot(snr,linearBER,'-.b','DisplayName','8x8 Linear');
plot(snr,iterativeBER,':b','DisplayName','8x8 Polar Only');
title('BER vs SNR for MIMO Systems, K=16, Estimated Channel');
ylabel('BER');
xlabel('SNR');
set(gca,'YScale','log');
legend(gca,'show');
hold off;


saveas(fig,'estimateBER.png');


