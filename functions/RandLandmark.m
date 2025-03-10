function [l, sDev] = RandLandmark(p,N)

wp      = p.WayPoints;
length  = p.Length;

sDev    = max(ceil(length/15),3);

idx     = randi([1, size(wp,2)],[1, N]);
lPx     = wp(1,idx);
lPy     = wp(2,idx);

r       = sqrt(rand(1, N)) * sDev;
theta   = rand(1, N) * 2 * pi;

l(1,:)  = r .* cos(theta) + lPx;
l(2,:)  = r .* sin(theta) + lPy;

end