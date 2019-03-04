LEN = 1; % how many K-bit length messages we will send (per tx/rx)
SNR = 1;
n = 2; % number of tx and rx antennas
K = 16; % bits per msg
R = .5; % polar rate
N = (2^nextpow2(K))/R; % bits per coded symbol
qamBitSize = 1;
qamSize = 2^qamBitSize;
normAnt = 1;
normConst = 1;
precode = 0;

TRAIN_SIZE = 200; 
TEST_SIZE = 20;

training.B = zeros(n*K, TRAIN_SIZE);
training.Y = zeros(n*(N/qamBitSize + n), TRAIN_SIZE);

size(training.B)
size(training.Y)

testing.B = zeros(n*K, TEST_SIZE);
testing.Y = zeros(n*(N/qamBitSize + n), TEST_SIZE);

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

Hest = Y*pilotData'*inv(pilotData*pilotData');

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

  %Y = [real(Y); imag(Y)];

  training.Y(:,i) = [reshape(Y,[],1); reshape(Hest,[],1)];
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

  %Y = [real(Y); imag(Y)];

  testing.Y(:,i) = [reshape(Y,[],1); reshape(Hest,[],1)];
  testing.B(:,i) = reshape(B, [], 1);
end

data.training = training;
data.testing = testing;

save('data.mat','data');
save('train.mat','training');
save('test.mat','testing');


