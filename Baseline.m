LEN = 20; % how many K-bit length messages we will send (per tx/rx)
SNR = 300;
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


initPC(N,K,'AWGN',SNR);
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
noiseVal = 10^(-SNR/10)*K/N;
noiseVec = sqrt(noiseVal)*randn(n,1);
        
% Apply channel
Y = H*X + noiseVec;

%Hest = ChannelEstimate(rxPilots, txPilots);

Hest=H; % perfect CSI

% MIMO Detect
[Yhat,zf] = LinearMIMODecoder(n, newLen, N, Y, qamTab, Hest, normAnt);
%Yhat - original polar bits it guessed
%zf - recentered qam

% Polar Decode
Bhat = PolarDecoder(n, LEN, K, N, SNR, zf);


disp('How off the post-encoded bits and pre-decoded bits are');
disp(sum(sum(sum(abs(enc-Yhat)))));
disp('How off the pre coded bits and post decoded bits are');
disp(sum(sum(sum(abs(B-Bhat)))));
