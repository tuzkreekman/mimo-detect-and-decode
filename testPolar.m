function [data,enc,dec] = testPolar(n,LEN,K,N)

data = randi([0, 1], n, LEN, K);

enc = zeros(n, LEN, N); %antennas, nMessages, polarized bits

for (kk = 1 : LEN)
    for (i = 1 : n)
        % polar encode each set of K bits at a time
        enc(i,kk,:) = pencode(squeeze(data(i,kk,:)));
    end

end

dec = zeros(n, LEN, K);

SNR=0;

noiseVal = 10^(-SNR/10);
 
for (kk = 1 : LEN)
    enc(:,kk,:) = enc(:,kk,:) + sqrt(noiseVal)*rand(1,n,N);
    for (i = 1 : n)
        dec(i,kk,:) = pdecode(squeeze(enc(i,kk,:)), 'AWGN', SNR);
    end
end


disp('# bits with errors');
disp(sum(sum(sum(abs(dec-data)))));

end

