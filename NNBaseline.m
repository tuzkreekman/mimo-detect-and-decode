LEN = 1; % how many K-bit length messages we will send (per tx/rx)
SNR = 1;
n = 2; % number of tx and rx antennas
K = 2; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 0;
normConst = 0;
precode = 0;

TRAIN_SIZE = 200; 
TEST_SIZE = 20;

training.B = zeros(n, K, TRAIN_SIZE);
training.Y = zeros(n, N/qamBitSize, 1, TRAIN_SIZE);

testing.B = zeros(n, K, TEST_SIZE);
testing.Y = zeros(n, N/qamBitSize, 1, TEST_SIZE);

addpath('./samples/polar');
addpath('./samples/polar/functions');

layers = [
    imageInputLayer([n N/qamBitSize 1])
    
    fullyConnectedLayer(512)
    reluLayer
    fullyConnectedLayer(256)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(128)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(64)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(32)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(16)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(K)
    softmaxLayer %should be sigmoid? 
    pixelClassificationLayer  %maybe??
];

options =  trainingOptions('adam', ...
    'InitialLearnRate',3e-4, ...
    'SquaredGradientDecayFactor',0.99, ...
    'MaxEpochs',20, ...
    'MiniBatchSize',64, ...
    'Plots','training-progress');

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

% Channel estimate
B = MIMOGenerator(n, LEN, K);
[X, newLen, enc, enc_old] = ApplyPolarQAM(B, n, LEN, N, K, R, qamBitSize, qamTab, precode, H_known);
noiseVal = 10^(-SNR/10)*K/N;
noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
Hest = (H*X + noiseVec)/X; 

% Create Training data
for (i=1:TRAIN_SIZE)
  % Create MIMO data
  B = MIMOGenerator(n, LEN, K);

  % Polar encode and modulate
  [X, newLen, enc, enc_old] = ApplyPolarQAM(B, n, LEN, N, K, R, qamBitSize, qamTab, precode, H_known);


  % Receive antenna noise - AWGN
  noiseVal = 10^(-SNR/10)*K/N;
  noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
  % Apply channel
  Y = Hest*X + noiseVec; % Nonfading gaussian channel

  training.Y(:,:,1,i) = Y;
  training.B(:,:,i) = B; 
end

% Create Test data
for (i=1:TEST_SIZE)
  % Create MIMO data
  B = MIMOGenerator(n, LEN, K);

  % Polar encode and modulate
  [X, newLen, enc, enc_old] = ApplyPolarQAM(B, n, LEN, N, K, R, qamBitSize, qamTab, precode, H_known);


  % Receive antenna noise - AWGN
  noiseVal = 10^(-SNR/10)*K/N;
  noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
  % Apply channel
  Y = H*X + noiseVec; % Nonfading gaussian channel

  testing.Y(:,:,1,i) = Y;
  testing.B(:,:,i) = B;
end

%net = trainNetwork(training.Y,training.B,layers,options);
%Bhat = classify(net,testing.Y);

%ber = sum(sum(sum(abs(Bhat-testing.B))))/(n*TEST_SIZE*K);

%disp('BER - NN');
%disp(mean(ber));

data.training = training;
data.testing = testing;

save('data.mat','data');


