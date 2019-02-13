function SimZFMLComplex(n, SNRvec, qamSize, normAnt, normConst)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate ZF and ML decoding in MIMO systems.
% Includes both full ML decoding (CSI at RX only) and ML decoding with 
% SVD precoding (CSI at TX and RX). ZF and MMSE assume no precoding. 
%
% Inputs:
%   n - number of antennas.
%   SNRcec - vector of the inverses of noise variances, in dB.
%            Note: if the average Tx power is unity, then this corresponds
%            to the SNR values.
%   qamSize - Size of QAM constellation to be transnmitted in each antenna.
%             2 corresponds to {+-1}, 4 corresponds to {+-1, +-i},...
%   normAnt - a binary flag for normalizing the Tx power by 1/sqrt(n).
%   normConst - a binary flag for normalizing the Tx constellation to have
%               an average unit power.
%
%   Output - performance plot.
%
%   Example: SimZFMLComplex(4,[-5:5:25],4,1,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LEN = 2000;
NoiseVarLin = 10.^(SNRvec/10);
PC_ML_BER = zeros(LEN,length(NoiseVarLin));
ZF_BER = zeros(LEN,length(NoiseVarLin));
MMSE_BER = zeros(LEN,length(NoiseVarLin));
ML_BER = zeros(LEN,length(NoiseVarLin));

% Create constallation table
[tab,normedEnergy, avgEnergy] = genConstalTab(qamSize, normConst);
antennaNorm = 1;
if (normAnt)
    antennaNorm = 1/sqrt(n);
end

sigX = (avgEnergy/normedEnergy^2)/antennaNorm^2;

for (jj = 1 : length(NoiseVarLin))
    noiseVal = 1./NoiseVarLin(jj);
    snrVal = sigX/noiseVal;
     
    for (kk = 1 : LEN)
        % Channel matrix
        H = randn(n).*exp(-1i*2*pi*rand(n,n));
        [U,D,V] = svd(H);
        
        noiseVec = sqrt(noiseVal)*randn(n,1);
        idxVec = randi(qamSize, n, 1);
        x = tab(idxVec).'/antennaNorm;

        %%% We consider two Tx schemes - with pre-coding and without
        % Precoding
        Ypc = H*V*x + noiseVec;
        % Clean
        Y = H*x + noiseVec;
        
        %%% Decoding based on pre-coding transmission
        rml = U'*Ypc; % ML
                
        %%% Decoding with MMSE
        wzf = pinv(H)*Y; % ZF
        wmmse = ((1/snrVal)*eye(n) + H'*H)\H'*Y; % MMSE
        if (2 <= n) && (n <= 4)
            wml = MLDec(Y,H,n,tab.'/antennaNorm);
            ML_BER(kk,jj) = length(find( wml ~= x )) / n;
        end
        
        PC_ML_BER(kk,jj) = length(find(qamDetector(rml./diag(D), qamSize, antennaNorm, normedEnergy) ~= x))/n;
        mmseTmp = length(find(qamDetector(wmmse, qamSize, antennaNorm, normedEnergy) ~= x))/n;
        MMSE_BER(kk,jj) = mmseTmp;
        zfTmp = length(find(qamDetector(wzf, qamSize, antennaNorm, normedEnergy) ~= x))/n;
        ZF_BER(kk,jj) = zfTmp;
    end
end

PC_ML_BERAvg = mean(PC_ML_BER);
MMSE_BERAvg = mean(MMSE_BER);
ZF_BERAvg = mean(ZF_BER);
ML_BERAvg = mean(ML_BER);

% plotting
figure;
if (2 <= n) && (n <= 4)
    semilogy(SNRvec,PC_ML_BERAvg,'b',SNRvec,MMSE_BERAvg,'r',SNRvec, ZF_BERAvg,'k',SNRvec, ML_BERAvg,'m');
    grid on;
    legend('ML PC','MMSE', 'ZF', 'ML');
else
    semilogy(SNRvec,PC_ML_BERAvg,'b',SNRvec,MMSE_BERAvg,'r',SNRvec, ZF_BERAvg,'k');
    grid on;
    legend('ML PC','MMSE', 'ZF');
end
xlabel('SNR [dB]');
ylabel('Pe');
title(['n = ' num2str(n) ', normAnt = ' num2str(normAnt) ', qamSize = ' num2str(qamSize) ', normConst = ' num2str(normConst)]);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function x = MLDec(Y,H,n,consTab)
% Exhaustively check all possible symbols vectors and return the single
% vector that is closest in L2 norm to Y.
if n == 2
    xcand = cartprod( consTab, consTab ).';
elseif n == 3
    xcand = cartprod( consTab, consTab, consTab ).';
elseif n == 4
    xcand = cartprod( consTab, consTab, consTab, consTab ).';
end
        
errMat = repmat(Y,1,length(xcand)) - H*xcand;
errVec = sum(errMat.*conj(errMat)).';
[~,idx] = min(errVec);
x = xcand(:,idx);
end

%-----------------------------------------------------------------%
function [tab, normFactor, avgEnergy] = genConstalTab(qamSize, normConst)
% Generate M-QAM constellation and return the average energy of the
% constellation
m = sqrt(qamSize);
assert( m == round(m) );

pam = 2*[1:m]-m-1;
prod = cartprod(pam, pam);
twoSided = (prod(:,1) + prod(:,2)*1i).';
avgEnergy = twoSided*twoSided'/length(twoSided);
normFactor = 1;
if (normConst)
    normFactor = sqrt(avgEnergy);
end

tab = twoSided/normFactor;

end

function r = qamDetector(y, qamSize, antennaNorm, constNorm)
    pamSize = sqrt(qamSize);
    
    re = pamDetector( real(y), pamSize, antennaNorm, constNorm );
    im = pamDetector( imag(y), pamSize, antennaNorm, constNorm );
    
    r = re + 1i*im;
end

%-----------------------------------------------------------------%
function r = pamDetector(y, pamSize, antennaNorm, constNorm)

r = zeros(size(y));

% calc constallation norm factor andnormalize back
oneSided = 2*[1:pamSize/2]-1;

ytmp = abs(y)*constNorm*antennaNorm;

for (kk = 1 : length(ytmp))
    
    if (ytmp(kk) > oneSided(end))
        r(kk) = oneSided(end);
        continue;
    end
    
    r(kk) = 2*floor(ytmp(kk)/2) + 1;    
    
end

r = r.*sign(y)/constNorm/antennaNorm;

end
