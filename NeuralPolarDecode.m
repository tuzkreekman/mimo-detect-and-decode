global H Hest flatHest TRAIN_SIZE TEST_SIZE SNR n K R N INPUT_SIZE qamTab antennaNorm options

SNR = 10;
n = 4; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;

TRAIN_SIZE = 2^K;
TEST_SIZE = 1024; 

INPUT_SIZE = 2*n*(N/qamBitSize + n);
OUTPUT_SIZE = K;

addpath('./samples/polar');
addpath('./samples/polar/functions');

initPC(N,K,'AWGN',0);
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

% Create constallation table
qamTab = ConstellationTable(qamSize, normConst);

% Channel matrix - Gaussian
%H = eye(n);
H = randn(n).*exp(-1i*2*pi*rand(n,n));

% Channel estimate setup
noiseVal = 10^(-SNR/10)*K/N;
noiseVec = sqrt(noiseVal)*randn(n,n); % Each symbol is received noisily

antennaNorm = 1;
if (normAnt)
    antennaNorm = 1/sqrt(n);
end

% Channel estimation
pilotData = hadamard(n);
% Apply channel
Y = antennaNorm*H*pilotData + noiseVec; % Nonfading gaussian channel
%Hest = 1/antennaNorm*Y*pilotData'*inv(pilotData*pilotData');
Hest = H; %perfect channel knowledge
flatHest = [real(Hest); imag(Hest)];
flatHest = reshape(flatHest,[],1);

layers = [
    imageInputLayer([INPUT_SIZE 1 1])

    fullyConnectedLayer(512)
    reluLayer
    fullyConnectedLayer(256)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(128)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(64)
    reluLayer
    fullyConnectedLayer(32)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(16)
    batchNormalizationLayer
    reluLayer
    fullyConnectedLayer(OUTPUT_SIZE)
    sigmoidLayer
    sigmoidClassificationLayer
];


options =  trainingOptions('sgdm', ...
    'MaxEpochs',1, ...
    'Verbose',false, ...
    'Plots','none',...
    'MiniBatchSize',64);

SNR = 30;

net = TrainEpoch(layers);

for i=1:(10^4)
    layers = net.Layers;
    net = TrainEpoch(layers);
end
disp('first 10000');

SNR = 8;

for i=1:(2*10^4)
    layers = net.Layers;
    net = TrainEpoch(layers);
end
disp('next 20000');
SNR = 0;

for i=1:(2*10^4)
    layers = net.Layers;
    net = TrainEpoch(layers);
end


save ('polar_decoder.mat','net')


function [net] = TrainEpoch(layers)
global H Hest flatHest TRAIN_SIZE TEST_SIZE SNR n K R N INPUT_SIZE qamTab antennaNorm options

training = GenerateData(Hest,Hest,flatHest,TRAIN_SIZE,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm);
testing = GenerateData(H,Hest,flatHest,TEST_SIZE,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm);

net = trainNetwork(training.Y,training.B,layers,options);
Bhat = predict(net,testing.Y);
%Bhat = 1./(1+exp(-Bhat));

disp(mean(mean(abs(squeeze(Bhat>.5) - squeeze(testing.B)'))));
end

function [data] = GenerateData(H,Hest,flatHest,LEN,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm)
% Create MIMO data
B = MIMOGenerator(n, LEN, K);

% Polar encode and modulate
[X, newLen, enc, enc_old] = ApplyPolarQAM(B, n, LEN, N, K, R, qamTab.qamBitSize, qamTab, 0, 0);

% Receive antenna noise - AWGN
noiseVal = 10^(-SNR/10)/n;    % scaled by 1/n so that the noise power per receiver is relative to unit power
noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
% Apply channel
Y = antennaNorm*H*X + noiseVec; % Nonfading gaussian channel

% CHannel equalize
Y = pinv(Hest)*Y;

Y = [real(Y); imag(Y)];

Y = [reshape(Y,[],LEN); repmat(flatHest,1,LEN)];
B = reshape(B,[],LEN);

data.Y = reshape(Y, INPUT_SIZE, 1, 1, LEN);
data.B = reshape(B(1:K,:), 1, 1, K, LEN);
end




