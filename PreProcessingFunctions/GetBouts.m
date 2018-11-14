function [actState, ActBouts] = GetBouts(Vr, Vf, Vs, dts, dtw, tvr, tvf)

if nargin  < 4
    dts = 20;
    dtw = 10;
    tvr = 50;
    tvf = 0.5;
end

%% Find Movement Bouts Treadmilll
actState = zeros(length(Vr),1);
ind = find(abs(Vr) > tvr | abs(Vf) > tvf | abs(Vs)>tvf);
a = 1;
for i = 1 : length(actState)
    if a < length(ind) && i == ind(a)
        actState(i) = 1;
        a = a + 1;
    end
end
a = 1;
thisBout = [];
for i = 1 : length(actState) - 1
    if actState(i) == 1 && actState(i + 1) == 1
        thisBout = vertcat(thisBout, i);
    elseif actState(i) == 1 && actState(i + 1) == 0
        thisBout = vertcat(thisBout, i);
        if length(thisBout) < dtw
            actState(thisBout) = 0;
        end
        thisBout = [];
    end
end
thisBout = [];
for i = 1 : length(actState) - 1
    if actState(i) == 0 && actState(i + 1) == 0
        thisBout = vertcat(thisBout, i);
    elseif actState(i) == 0 && actState(i + 1) == 1
        thisBout = vertcat(thisBout, i);
        if length(thisBout) < dts
            actState(thisBout) = 1;
        end
        thisBout = [];
    end
end

%% Get Bouts
ind1 = find(actState == 1);
ind0 = find(actState == 0);
b = 1;
for j = 1 : (length(ind1) - 1)
    if (ind1(j+1) - ind1(j) == 1)
    else
        b = b + 1;
    end
end
c = 1;
a = 1;
ActBouts = cell(b, 1);
for j = 1 : (length(ind1) - 1)
    if (ind1(j+1) - ind1(j) == 1)
        inds(a) = ind1(j);
        a = a + 1;
    else
        inds(a) = ind1(j);
        if length(inds) > dtw
            ActBouts{c} = inds;
            c = c + 1;
        end
        a = 1;
        clearvars inds;
    end
end
if length(ind1)~= 0 && length(ind0)~= 0
    if exist('inds')
        if length(inds) > dtw
            ActBouts{c} = inds;
        end
    end
elseif length(ind1)~= 0 && length(ind0)== 0
    if exist('inds')
        if length(inds) > dtw
            ActBouts{c} = inds;
        end
    end
end
clearvars inds ind1 ind2;
while(isempty(ActBouts{end}))
    ActBoutsT = ActBouts;
    clearvars ActBouts1
    ActBouts = cell(length(ActBoutsT)-1,1);
    for i = 1 : length(ActBoutsT)-1
        ActBouts(i) = ActBoutsT(i);
    end
    clearvars ActBoutsT
    if length(ActBouts) < 1
        break;
    end
end

end

