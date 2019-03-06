LEN = 1; % how many K-bit length messages we will send (per tx/rx)
SNR = 10;
n = 2; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;

TRAIN_SIZE = 1024; %2^K; 
TEST_SIZE = 1024; %2^K;
INPUT_SIZE = 2*n*(N/qamBitSize + n);
OUTPUT_SIZE = K;

training.B = zeros(n*K, TRAIN_SIZE);
training.Y = zeros(n*(2*N/qamBitSize + 2*n), TRAIN_SIZE);

size(training.B)
size(training.Y)

testing.B = zeros(n*K, TEST_SIZE);
testing.Y = zeros(n*(2*N/qamBitSize + 2*n), TEST_SIZE);

addpath('./samples/polar');
addpath('./samples/polar/functions');

initPC(N,K,'AWGN',0);
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

% Create constallation table
qamTab = ConstellationTable(qamSize, normConst);

% Channel matrix - Gaussian
H = eye(n);
%H = randn(n).*exp(-1i*2*pi*rand(n,n));

% Can do precoding with full or statistical channel information for point-to-point MIMO systems
% For now, assume TX has perfect initial channel knowledge
% Can reach channel capacity with perfect channel knowledge and SVD precoding
% We aren't actually precoding right now though
H_known = H;

% Channel estimate
noiseVal = 10^(-SNR/10)*K/N;
noiseVec = sqrt(noiseVal)*randn(n,n); % Each symbol is received noisily

antennaNorm = 1;
if (normAnt)
    antennaNorm = 1/sqrt(n);
end

pilotData = hadamard(n);
% Apply channel
Y = antennaNorm*H*pilotData + noiseVec; % Nonfading gaussian channel

%Hest = 1/antennaNorm*Y*pilotData'*inv(pilotData*pilotData');
Hest = H;
flatHest = [real(Hest); imag(Hest)];
flatHest = reshape(flatHest,[],1);

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
  Y = antennaNorm*Hest*X + noiseVec; % Nonfading gaussian channel

  Y = [real(Y); imag(Y)];

  training.Y(:,i) = [reshape(Y,[],1); flatHest];
  training.B(:,i) = reshape(B, [], 1);
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
  Y = antennaNorm*H*X + noiseVec; % Nonfading gaussian channel

  Y = [real(Y); imag(Y)];

  testing.Y(:,i) = [reshape(Y,[],1); flatHest];
  testing.B(:,i) = reshape(B, [], 1);
end

data.training = training;
data.testing = testing;

layers = [
    imageInputLayer([INPUT_SIZE 1 1])

    fullyConnectedLayer(300)
    reluLayer
    fullyConnectedLayer(300)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(300)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(300)
    reluLayer
    fullyConnectedLayer(300)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(300)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(OUTPUT_SIZE)
    sigmoidLayer
    bitClassificationLayer
];


options =  trainingOptions('adam', ...
    'InitialLearnRate',3e-4, ...
    'MaxEpochs',200, ...
    'MiniBatchSize',64, ...
    'Plots','training-progress');

training.B = training.B(1:K,:);
testing.B = testing.B(1:K,:);
training.Y(1:K,:) = training.B;
testing.Y(1:K,:) = testing.B;

training.B = reshape(training.B, 1, 1, K, TRAIN_SIZE);
testing.B = reshape(testing.B, 1, 1, K, TEST_SIZE);
training.Y = reshape(training.Y, INPUT_SIZE, 1, 1, TRAIN_SIZE);
testing.Y = reshape(testing.Y, INPUT_SIZE, 1, 1, TEST_SIZE);

net = trainNetwork(training.Y,training.B,layers,options);
Bhat = predict(net,testing.Y);

disp(mean(mean(abs(squeeze(Bhat>.5) - squeeze(testing.B)'))));

