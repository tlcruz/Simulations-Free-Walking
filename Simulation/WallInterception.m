function [pi] = WallInterception(p1, p2, r)
dx = p2(1)-p1(1);
dy = p2(2)-p1(2);
dr = sqrt(dx*dx + dy*dy);
D = p1(1)*p2(2)-p1(2)*p2(1);
Dt = sqrt(r*r*dr*dr-D*D);

x1 = (D*dy + dx*Dt)/(dr*dr);
y1 = (-D*dx + abs(dy)*Dt)/(dr*dr);
x2 = (D*dy - dx*Dt)/(dr*dr);
y2 = (-D*dx - abs(dy)*Dt)/(dr*dr);

bx = abs(x2-x1);
bxa = abs(p1(1)-x1);
bxb = abs(p2(1)-x1);
if(bxa > bx || bxb > bx)    
    pi(1) = x2;
else
    pi(1) = x1;
end

by = abs(y2-y1);
bya = abs(p1(2)-y1);
byb = abs(p2(2)-y1);
if(bya > by || byb > by)
    pi(2) = y2;
else
    pi(2) = y1;
end
end