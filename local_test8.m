global H Hest flatHest TRAIN_SIZE TEST_SIZE SNR n K R N INPUT_SIZE qamTab antennaNorm options BER

SNR = 30;
n = 2; % number of tx and rx antennas
K = 8
%K = 16; % bits per msg
R = .5; % polar rate
%N = 32
disp('N is')
BER = [];
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;

TRAIN_SIZE = 2^K; 
TEST_SIZE = 40; 

disp('Input size is')
INPUT_SIZE = 2*n*(N/(n*qamBitSize) + n)

OUTPUT_SIZE = K;

addpath('./samples/polar');
addpath('./samples/polar/functions');

initPC(N,K,'AWGN',0);
%SNR: Default: 0dB;  := Eb/N0,  where (K*Eb/N) is the energy used during BPSK modulation of coded-bits)

% Create constallation table
qamTab = ConstellationTable(qamSize, normConst);

% Channel matrix - Gaussian
H = eye(n);
%H = randn(n).*exp(-1i*2*pi*rand(n,n));

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
Hest = H;
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
    %sigmoidLayer
    sigmoidClassificationLayer
];


options =  trainingOptions('adam', ...
    'InitialLearnRate',3e-3, ...
    'GradientDecayFactor',0.99,...
    'MaxEpochs',1, ...
    'MiniBatchSize',64);


%disp('Lr is')
lr = 3e-3;
SNR = 0;
snr = 0
net = TrainEpoch(layers);
%BER = [];
for i=1:1000
    lr = lr*0.99;
    disp('iter is')
    disp(i)
    %disp('b is')
    b = mod(i,100);
    if b == 0
        snr = snr+1;
        disp('new snr')
        SNR = snr
    end
    %options.InitialLearnRate = lr
    options =  trainingOptions('adam', ...
    'InitialLearnRate',lr, ...
    'MaxEpochs',1, ...
    'MiniBatchSize',64);
    for j = 1:1
    %disp('iter is')
    %disp(i*j)
    options =  trainingOptions('adam', ...
    'InitialLearnRate',lr, ...
    'MaxEpochs',100, ...
    'MiniBatchSize',64);
   % disp('i is')
   % disp(i)
        layers = net.Layers;
        net = TrainEpoch(layers);
    end
end

disp('BER is')
disp(BER)
save("BER8.mat","BER")
%disp('BER is')
%disp(BER)
%{
SNR = 8;

for i=1:1
    layers = net.Layers;
    net = TrainEpoch(layers);
end

SNR = 6;

for i=1:1
    layers = net.Layers;
    net = TrainEpoch(layers);
end

%}
function [net] = TrainEpoch(layers)
global H Hest flatHest TRAIN_SIZE TEST_SIZE SNR n K R N INPUT_SIZE qamTab antennaNorm options BER

disp('Training begins')

training = GenerateData(Hest,flatHest,TRAIN_SIZE,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm);
disp('Testing begins')
testing = GenerateData(H,flatHest,TEST_SIZE,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm);

net = trainNetwork(training.Y,training.B,layers,options);
%disp('Bhat is')
Bhat = predict(net,testing.Y);
Bhat = 1./(1+exp(-Bhat));
disp('We print here')
disp(squeeze(Bhat>.5))
A = (mean(mean(abs(squeeze(Bhat>.5) - squeeze(testing.B)'))))
%disp('BER is')
BER = [BER ; A];

end

function [data] = GenerateData(H,flatHest,LEN,SNR,n,K,R,N,INPUT_SIZE,qamTab,antennaNorm)
% Create MIMO data
B = parMIMOGenerator(n, LEN, K);

% Polar encode and modulate
[X, newLen, enc, enc_old] = ParPolarQAM(B, n, LEN, N, K, R, qamTab.qamBitSize, qamTab, 0, 0);

% Receive antenna noise - AWGN
noiseVal = 10^(-SNR/10)/n;    % scaled by 1/n so that the noise power per receiver is relative to unit power
noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
% Apply channel
Y = antennaNorm*H*X + noiseVec; % Nonfading gaussian channel

Y = [real(Y); imag(Y)];

Y = [reshape(Y,[],LEN); repmat(flatHest,1,LEN)];
B = reshape(B,[],LEN);
%disp('size of Y is')
%disp(size(Y))
%disp('LEN is')
%disp(LEN)
data.Y = reshape(Y, INPUT_SIZE, 1, 1, LEN);
%disp(size(Y))
%disp('size of B is')
data.B = reshape(B(1:K,:), 1, 1, K, LEN);
%disp(size(B))

end

%save("trainedNet.mat","net")



