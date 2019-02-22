LEN = 10; % how many K-bit length messages we will send (per tx/rx)
ITERS = 2000;
snr = 0:10;
n = 2; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 0;
normConst = 0;
precode = 0;

addpath('./samples/polar');
addpath('./samples/polar/functions');

linearBER    = zeros(11,ITERS);
iterativeBER = zeros(11,ITERS);

for (SNR=snr)
  for (i=1:ITERS)
    [ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,precode);
    linearBER(SNR+1,i) = ber1;
    iterativeBER(SNR+1,i) = ber2;
  end
end

linearBER = mean(linearBER,2);
iterativeBER = mean(iterativeBER,2);

fig = figure();
hold on;
plot(snr,linearBER,'DisplayName','Linear');
plot(snr,iterativeBER,'DisplayName','Iterative');
title('BER vs SNR');
ylabel('BER');
xlabel('SNR');
set(gca,'YScale','log');
legend();
hold off;

saveas(fig,'BER.png');

function [linearBER, iterativeBER] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,precode)
initPC(N,K,'AWGN',0); % changd snr
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

% Create constallation table
qamTab = ConstellationTable(qamSize, normConst);

% Channel matrix - Gaussian
H = randn(n).*exp(-1i*2*pi*rand(n,n));

% Can do precoding with full or statistical channel information for point-to-point MIMO systems
% For now, assume TX has perfect initial channel knowledge
% Can reach channel capacity with perfect channel knowledge and SVD precoding
% We aren't actually precoding right now though
H_known = H;

% Create MIMO data
B = MIMOGenerator(n, LEN, K);

% Polar encode and modulate
[X, newLen, enc, enc_old] = ApplyPolarQAM(B, n, LEN, N, K, R, qamBitSize, qamTab, precode, H_known);


% Receive antenna noise - AWGN
noiseVal = 10^(-SNR/10);% CHANGED *K/N;
noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
% Apply channel
Y = H*X + noiseVec; % Nonfading gaussian channel

%Hest = ChannelEstimate(rxPilots, txPilots);

Hest=H; % perfect CSI

% MIMO Detect
[Yhat,wzf,zf] = LinearMIMODecoder(n, newLen, N, Y, qamTab, Hest, normAnt);
%Yhat - original polar bits it guessed
%wzf - equalized qam
%zf - recentered qam

% Linear Polar Decode - uses centered qams
Bhat1 = PolarDecoder(n, LEN, K, N, SNR,  zf);
% Iterative Polar Decode - uses equalized qams
Bhat2 = PolarDecoder(n, LEN, K, N, SNR, wzf);


linearBER = sum(sum(sum(abs(B-Bhat1))))/(K*LEN*n);
iterativeBER = sum(sum(sum(abs(B-Bhat2))))/(K*LEN*n);

end

