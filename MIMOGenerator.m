function data = MIMOGenerator(n, LEN, K)
% antennas, bits, nMessages
data = randi([0, 1], n, LEN, K);
end

