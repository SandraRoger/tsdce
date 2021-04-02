function [w1,w2,line1,line2] = wlsfreqest(u1,v1,wi)
% Simple function to estimate 2D frequencies of complex sinusoid
%
% Input:
% u1 = complex sinusoid in direction 1
% v1 = complex sinusoid in direction 2
% wi = weighting values for weighted LS.
%
% Output:
% w1, w2: Estimated frequencies


% Number of points
nt = length(u1);

% ===============================================
% Frequency of u1
% ===============================================

% Force 0 angle at first element before wrapping
Yd1 = angle(u1(2:end).*conj(u1(1:end-1)));

% Analyze variance for each wrapping alternative
% wraptopi works best for frequencies close to 0
% wrapto2pi works best for frequencies close to pi
vwrap_Yd1(1) = var(Yd1);
vwrap_Yd1(2) = var(wrapTo2Pi(Yd1));

if vwrap_Yd1(2)<vwrap_Yd1(1)
    Yd1 = wrapTo2Pi(Yd1);
end

line1 = [0; cumsum(Yd1)];

% Weighted least squares
wi = wi(:);
xi = (0:nt-1);
xi = xi(:);
barxw = sum(wi.*xi)/sum(wi);
yi = line1;
baryw = sum(wi.*yi)/sum(wi);
w1 = sum(wi.*(xi-barxw).*(yi-baryw));
w1 = w1./(sum(wi.*((xi-barxw).^2)));

w1 = wrapToPi(w1);

% ===============================================
% Frequency of u2
% ===============================================

Yd2 = angle(v1(2:end).*conj(v1(1:end-1)));

vwrap_Yd2(1) = var(Yd2);
vwrap_Yd2(2) = var(wrapTo2Pi(Yd2));

if vwrap_Yd2(2)<vwrap_Yd2(1)
    Yd2 = wrapTo2Pi(Yd2);
end

line2 = [0; cumsum(Yd2)];

wi = wi(:);
xi = (0:nt-1);
xi = xi(:);
barxw = sum(wi.*xi)/sum(wi);
yi = line2;
baryw = sum(wi.*yi)/sum(wi);
w2 = sum(wi.*(xi-barxw).*(yi-baryw));
w2 = w2./(sum(wi.*((xi-barxw).^2)));

w2 = wrapToPi(w2);

