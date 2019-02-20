function h_hat = estimate_channel( H, snr_dB, k )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Perform least-squares based channel estimation
% Inputs:
%   H - complex channel gain matrix, must be square
%   snr_dB - ratio of pilot power to channel noise in dB
%   k - (Optional) number of pilots, default is the size of n
%
% Outputs:
%   h_hat - least-squares estimate of H
%
% Example: 
%   h = randn(4,4).*exp(-1i*2*pi*rand(4,4));
%   h_hat = estimate_channel( h, 10 );
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Make sure H is square and k is set 
    [n,m] = size(H);
    assert(n == m );
    if nargin < 3 || k < n
        k = n;
    end
    
    %Pilots are typically a DFT matrix
    P = dftmtx(n);
    %P = hadamard(n); %This is what the paper says but not sure it matters
    % much since abs(det(detmtx(n))) == abs(det(hadamard(n))
    if k > n
        P = [P exp(-1i*2*pi*randi(n, n, k-n)/n)];
    end
    
    %Compute noise power
    P_power = sum(abs(P(:,1).^2)) / n;
    snr_lin = 10.^(snr_dB/10);
    noise_var = P_power / snr_lin;
    
    %noise is circularly symmetric complex Gaussian
    E = noise_var*randn(n,k).*exp(-1i*2*pi*rand(n,k));
    
    S = H*P+E;
    
    h_hat = S/P;
end

function d = dftmtx(n)
% n-dimensional DFT matrix
  assert( round(n) == n );
  d = fft(eye(n));
end

function h = hadamard( n )
% n-dimensional Hadamard matrix (i.e. \pm 1 valued matrix with maximal
% determinant.  Simple construction that requires n is a power of 2.
    assert( round(log2(n)) == log2(n) );
    
    h2 = [1 1; 1 -1];
    h  = h2;
    for ii = 2:log2(n)
        h = kron(h, h2);
    end
end
