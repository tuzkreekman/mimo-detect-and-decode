function [qamBER, linearBER, polarBER] = GetBER(LEN,SNR,n,K,R,N,qamBitSize,qamSize,normAnt,normConst,perfectKnowledge,precode)

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
noiseVal = 10^(-SNR/10)/n;    % scaled by 1/n so that the noise power per receiver is relative to unit power
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
Hest = 1/antennaNorm*Y(:,p)*X(:,p)'*inv(X(:,p)*X(:,p)');

if (perfectKnowledge)
    Hest = H;
end

Y = Y(:,n+1:end); % Ydata


% MIMO Detect
[Yhat,wzf,zf] = LinearMIMODecoder(n, newLen, N, Y, qamTab, Hest, antennaNorm);
%Yhat - original polar bits it guessed
%wzf - equalized qam
%zf - recentered qam

% Linear Polar Decode - uses centered qams
Bhat1 = PolarDecoder(n, LEN, K, N, SNR,  zf);
% Simple Polar Decode - uses equalized qams with noise
Bhat2 = PolarDecoder(n, LEN, K, N, SNR, wzf);

qamBER = sum(sum(sum(abs(enc-Yhat))))/(newLen*n);
linearBER = sum(sum(sum(abs(B-Bhat1))))/(K*LEN*n);
polarBER = sum(sum(sum(abs(B-Bhat2))))/(K*LEN*n);

end


