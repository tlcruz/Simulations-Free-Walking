function [sig] = SpikeDirection(x0, y0, x, y)
if x0 > 0
    if y0 > 0
        if y0 > x0
            if x > x0
                sig = -1;
            else
                sig = 1;
            end
        else
            if y > y0
                sig = 1;
            else
                sig = -1;
            end
        end
    else
        if abs(y0) > x0
            if y0 > y
                sig = -1;
            else
                sig = 1;
            end
        else
            if x > x0
                sig = 1;
            else
                sig = -1;
            end
        end
    end
else
    if y0 > 0
        if y0 > abs(x0)
            if x0 > x
                sig = 1;
            else
                sig = -1;
            end
        else
            if y > y0
                sig = -1;
            else
                sig = 1;
            end
        end
    else
        if abs(y0) > abs(x0)
            if x0 > x
                sig = -1;
            else
                sig = 1;
            end
        else
            if y0 > y
                sig = 1;
            else
                sig = -1;
            end
        end
    end
end
end