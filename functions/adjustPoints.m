function [xq, yq] = adjustPoints(x, y, n)

t = 1:length(x);

tq = linspace(1, length(x), n);

if n > length(x)
    xq = interp1(t, x, tq, 'linear');
    yq = interp1(t, y, tq, 'linear');
elseif n < length(x)
    step = round(length(x) / n);
    xq = x(1:step:end);
    yq = y(1:step:end);
else
    xq = x;
    yq = y;
end
end
