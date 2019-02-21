function data = MIMOGenerator(n, LEN, K)
% antennas, nMessages, bits
data = randi([0, 1], n, LEN, K);
end

