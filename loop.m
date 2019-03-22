final_BER = [];
fig = figure();
SNR = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
sum = 0;
for i = 1:10
    sum = 0;
    disp('sum reset')
    for j = 0:99
        sum = sum + BER(i+j);
    end
    sum = sum/100;
    final_BER = [final_BER ; sum];
end
disp(final_BER)
plot(SNR, final_BER);

savefig(avgplot);
