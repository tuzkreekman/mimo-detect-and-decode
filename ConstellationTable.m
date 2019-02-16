function qam = ConstellationTable(qamSize, normConst)
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

qam.table = tab;
qam.normFactor = normFactor;
qam.avgEnergy = avgEnergy;

end


