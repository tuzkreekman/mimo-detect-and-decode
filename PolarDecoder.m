function dec = PolarDecoder(n, LEN, K, N, SNR, Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate polar decoding of bits in MIMO systems.
%
% Inputs:
%   n - number of antennas.
%   LEN - number of symbols to transmit
%   K - number of bits per symbol per antenna
%   R - polar coding rate
%       e.g. K=16, R=.5 -> 32 bits output of polar codes
%   Y - data to decode
%
%   Output - MIMO signal after QAM detection linearly
%
%   Example: PolarDecoder(4,20,16,.5,Y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = reshape(Y, n, LEN, N);

dec = zeros(n, LEN, K);

for (i=1:n)
    dec(i,:,:) = pdecode(squeeze(data(i,:,:)), 'AWGN', SNR);
end

dec = transpose(dec, [1,3,2]);

end

