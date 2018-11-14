function [angAdj, vvf] = VMInt(vF, mF, params, per, n)
vvf = 0;
if per.on == 1
    if per.FTBL(n) == 0
        adjFTBL = params.kVisualSize*vF.vFFTBL + 0.5*params.WeightFTB*min(mF,0);
        vvf = vvf + params.kVisualSize*vF.vFFTBL;
    else
        adjFTBL = per.FTBL(n);
    end
    if per.FTBR(n) == 0
        adjFTBR = params.kVisualSize*vF.vFFTBR + 0.5*params.WeightFTB*max(mF,0);
        vvf = vvf + params.kVisualSize*vF.vFFTBR;
    else
        adjFTBR = -per.FTBR(n);
    end
    if per.BTFL(n) == 0
        adjBTFL = params.kVisualSize*vF.vFBTFL + 0.5*params.WeightBTF*max(mF,0);
        vvf = vvf + params.kVisualSize*vF.vFBTFL;
    else
        adjBTFL = -per.BTFL(n);
    end
    if per.BTFR(n) == 0
        adjBTFR = params.kVisualSize*vF.vFBTFR + 0.5*params.WeightBTF*min(mF,0);
        vvf = vvf + params.kVisualSize*vF.vFBTFR;
    else
        adjBTFR = per.BTFR(n);
    end
else
    adjFTBL = params.kVisualSize*vF.vFFTBL + 0.5*params.WeightFTB*min(mF,0);
    adjFTBR = params.kVisualSize*vF.vFFTBR + 0.5*params.WeightFTB*max(mF,0);
    adjBTFL = params.kVisualSize*vF.vFBTFL + 0.5*params.WeightBTF*max(mF,0);
    adjBTFR = params.kVisualSize*vF.vFBTFR + 0.5*params.WeightBTF*min(mF,0);
    vvf = params.kVisualSize*vF.vFFTBL + params.kVisualSize*vF.vFFTBR + ...
        params.kVisualSize*vF.vFBTFL + params.kVisualSize*vF.vFBTFR;
end
angAdj = adjFTBL+adjFTBR+adjBTFL+adjBTFR;
end