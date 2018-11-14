function [Vf, Vs, Vr, Vt, Angle, nJ] = ClearJumps(Vf, Vs, Vr, Angle)
if (max(abs(diff(Vr))) > 250)
    [~,locs] = findpeaks(abs(diff(Vr)),'MinPeakHeight',250, 'MinPeakDistance',6);
    locs = vertcat(1,locs, length(Vf));
    nJ = length(locs-1);
    for i = 2 : length(locs)
        if (locs(i)+5 < length(Vf))
            if (locs(i) > 6)
                Vf(locs(i)-2:locs(i)+2) = mean(Vf(locs(i)-5:locs(i)-2));
                Vs(locs(i)-2:locs(i)+2) = mean(Vs(locs(i)-5:locs(i)-2));
                Vr(locs(i)-2:locs(i)+2) = mean(Vr(locs(i)-5:locs(i)-2));
            elseif (locs(i) > 2)
                Vf(locs(i)-2:locs(i)+2) = Vf(locs(i)-2);
                Vs(locs(i)-2:locs(i)+2) = Vs(locs(i)-2);
                Vr(locs(i)-2:locs(i)+2) = Vr(locs(i)-2);
            else
                Vf(locs(i):locs(i)+2) = Vf(locs(i)+4);
                Vs(locs(i):locs(i)+2) = Vs(locs(i)+4);
                Vr(locs(i):locs(i)+2) = Vr(locs(i)+4);
            end
        else
            Vf(locs(i)) = Vf(locs(i)-4);
            Vs(locs(i)) = Vs(locs(i)-4);
            Vr(locs(i)) = Vr(locs(i)-4);
        end
    end
end
Vt = sqrt(Vf.*Vf + Vs.*Vs);
nJ = 0;
end

