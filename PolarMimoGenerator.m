function PolarMIMOGenerator(n, K, R, qamSize, normAnt, normConst, precode, H)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate polar coding of bits in MIMO systems.
% Optionally can apply SVD precoding based on known channel matrix
%
% Inputs:
%   n - number of antennas.
%   K - number of bits per symbol per antenna
%   R - polar coding rate
%       e.g. K=16, R=.5 -> 32 bits output of polar codes
%   qamSize - Size of QAM constellation to be transnmitted in each antenna.
%             2 corresponds to {+-1}, 4 corresponds to {+-1, +-i},...
%   normAnt - a binary flag for normalizing the Tx power by 1/sqrt(n).
%   normConst - a binary flag for normalizing the Tx constellation to have
%               an average unit power.
%   precode - a binary flag for SVD precoding
%   H - the known channel, if available, used for precoding
%
%   Output - MIMO signal about to be transmitted
%
%   Example: PolarMIMOGenerator(4,4,1,1,1,ones(4))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

LEN = 20;

% Create constallation table
[qamTable, normedEnergy, avgEnergy] = genConstalTab(qamSize, normConst);
antennaNorm = 1;
if (normAnt)
    antennaNorm = 1/sqrt(n);
end

sigX = (avgEnergy/normedEnergy^2)/antennaNorm^2;


data = zeros(n, K, LEN); % antennas, bits, nMessages
enc = zeros(n, nextpow2(K/R), LEN); %antennas, polarized bits, nMessages

for (kk = 1 : LEN)
    if (precode)
        [U,D,V] = svd(H);
    end

    data(:,:,kk) = randi([0, 1], n, K); %generate bits
    for (i = 1 : n)
        enc(i,:,kk) = nrPolarEncode(data(i,:,kk), K/R);	
    end

end

enc = reshape(enc, n, LEN*polarSize/qamSize, qamSize);

collapsed = bi2de(enc);

output = qamTable(collapsed);

if (precode)
    output = V.*output;
end

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


%-----------------------------------------------------------------%
