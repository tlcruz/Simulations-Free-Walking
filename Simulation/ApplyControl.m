function [ang, mF, vF] = ApplyControl(vr, ang, params, per, i)
mF = 0;
vF = 0;
if length(vr) > 32
    addpath('HR Model');
    if params.kVisualSize ~= 0
        if isnan(params.WeightR)
            [vF] = GetVisualFeedbackHR(vr((end-30):end), params);
            vF = params.kVisualSize*vF;
            mF = -0.1*params.kernelMotorSize*params.kernelMotor*vr((end-30):end);
            ang = ang + mF + vF;
        else
%             [vF, ~, ~, params] = GetVisualFeedbackHRLR(vr((end-30):end), params);
            [vF, ~, ~, params] = GetVisualFeedbackHRLR(vr((end-10):end), params);
            mF = -0.1*params.kernelMotorSize*params.kernelMotor*vr((end-30):end);
            [angAdj, vF] = VMInt(vF, mF, params, per, i);
            ang = ang + angAdj;
        end
    else
        mF = -0.1*params.kernelMotorSize*params.kernelMotor*vr((end-30):end);
        ang = ang + mF;
    end 
end
end