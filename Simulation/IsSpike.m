function [bl] = IsSpike(XX, x, y, params)
if length(XX) > 25 && rand < (params.spkPa + params.spkPb*exp(params.spkPl*sqrt(x^2 + y^2)))
    bl = 1;
else
    bl = 0;
end
end