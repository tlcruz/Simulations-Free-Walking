function [per] = GetPerturbation(rNs, params)

pert = zeros(rNs,1);
perti = randi(rNs, [floor(rNs/(3*params.pertD)) 1]);
perti = sort(perti);
indsEx = find(diff(perti)<2.5*params.pertD);
perti(indsEx) = [];

for i = 1 : length(perti)
    pert(perti(i):(perti(i)+params.pertD)) = 1;
end

pert = pert * params.SizePerturbation;

if params.PertOn == 1
    per.on = 1;
else
    per.on = 0;
end

switch params.PerP
    case 'FTBL'
        per.FTBL = pert;
        per.FTBR = zeros(length(pert),1);
        per.BTFL = zeros(length(pert),1);
        per.BTFR = zeros(length(pert),1);
    case 'FTBR'
        per.FTBL = zeros(length(pert),1);
        per.FTBR = pert;
        per.BTFL = zeros(length(pert),1);
        per.BTFR = zeros(length(pert),1);
    case 'BTFL'
        per.FTBL = zeros(length(pert),1);
        per.FTBR = zeros(length(pert),1);
        per.BTFL = pert;
        per.BTFR = zeros(length(pert),1);
    case 'BTFR'
        per.FTBL = zeros(length(pert),1);
        per.FTBR = zeros(length(pert),1);
        per.BTFL = zeros(length(pert),1);
        per.BTFR = pert;
    otherwise
        per.FTBL = zeros(length(pert),1);
        per.FTBR = zeros(length(pert),1);
        per.BTFL = zeros(length(pert),1);
        per.BTFR = zeros(length(pert),1);  
end

per.pert = pert/(max(pert)-min(pert));

end