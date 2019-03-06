LEN = 4; % how many K-bit length messages we will send (per tx/rx)
ITERS = 100;
snr = 10;
n = 8; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;
perfectKnowledge = 1;

addpath('./samples/polar');
addpath('./samples/polar/functions');

initPC(N,K,'AWGN',0); % changd snr
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

qamBER    = zeros(ITERS,1);
linearBER = zeros(ITERS,1);
polarBER  = zeros(ITERS,1);

[ber0,ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,perfectKnowledge,precode);

parfor (i=1:ITERS)
    [ber0,ber1,ber2] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,perfectKnowledge,precode);
    qamBER(i)       = ber0;
    linearBER(i)    = ber1;
    polarBER(i)     = ber2;
end

qamBER    = mean(qamBER);
linearBER = mean(linearBER);
polarBER  = mean(polarBER);

disp('MIMO Decode only BER');
disp(qamBER);
disp('BER - Linear (QAM + Polar)');
disp(linearBER);
disp('BER - Polar Only');
disp(polarBER);
