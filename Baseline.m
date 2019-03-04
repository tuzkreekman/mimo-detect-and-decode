LEN = 200; % how many K-bit length messages we will send (per tx/rx)
SNR = 10;
n = 8; % number of tx and rx antennas
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
noiseVec = sqrt(noiseVal)*randn(n,newLen+n); % Each symbol is received noisily

antennaNorm = 1;
if (normAnt)
    antennaNorm = 1/sqrt(n);
end

pilotData = hadamard(n);
X = [pilotData, X]; %Xpilots,Xdata

% Apply channel
Y = antennaNorm*H*X + noiseVec; % Nonfading gaussian channel

p = 1:n;

%Hest = ChannelEstimate(rxPilots, txPilots);
%Hest = (H*X + noiseVec)/X;
%Hest = sqrt(n)*((Y(:,p))'.*(X(:,p)))*(inv((X(:,p))'*(X(:,p))));
Hest = sqrt(n)*Y(:,p)*X(:,p)'*inv(X(:,p)*X(:,p)');

Y =  Y(:,n+1:end); % Ydata

% MIMO Detect
[Yhat,wzf,zf] = LinearMIMODecoder(n, newLen, N, Y, qamTab, Hest, normAnt);
%Yhat - original polar bits it guessed
%wzf - equalized qam
%zf - recentered qam

% Linear Polar Decode - uses centered qams
Bhat1 = PolarDecoder(n, LEN, K, N, SNR,  zf);
% Iterative Polar Decode - uses equalized qams
Bhat2 = PolarDecoder(n, LEN, K, N, SNR, wzf);

disp('MIMO Decode only BER');
disp(sum(sum(sum(abs(enc-Yhat))))/(newLen*n));
disp('BER - Linear');
disp(sum(sum(sum(abs(B-Bhat1))))/(K*LEN*n));
disp('BER - Iterative');
disp(sum(sum(sum(abs(B-Bhat2))))/(K*LEN*n));




