addpath('./samples/polar');
addpath('./samples/polar/functions');

n=4;
SNR=10;
LEN=100;

B=MIMOGenerator(n,LEN,16);
initPC(32,16,'AWGN',0);

qamTab = ConstellationTable(4,0);
[X,newLen,enc,enc_old] = ApplyPolarQAM(B,n,LEN,32,16,.5,2,qamTab,0,0);

% Channel matrix - Gaussian
H = randn(n).*exp(-1i*2*pi*rand(n,n));

% Receive antenna noise - AWGN
noiseVal = 10^(-SNR/10);% CHANGED *K/N;
noiseVec = sqrt(noiseVal)*randn(n,newLen); % Each symbol is received noisily
        
% Apply channel
Y = H*X + noiseVec; % Nonfading gaussian channel

k=1;

for i=1:n
    for j=1:newLen
        x(k) = real(Y(i,j));
        y(k) = imag(Y(i,j));
        z(k) = i;
        k=k+1;
    end
end

scatter3(x,y,z)


