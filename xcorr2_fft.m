function c = xcorr2_fft(a,b)
%XCORR2_FFT Two-dimensional cross-correlation evaluated with FFT algorithm.

if nargin == 1
	b = a;
end

% Matrix dimensions
adim = size(a);
bdim = size(b);

% Cross-correlation dimension
cdim = adim+bdim-1;

bpad = zeros(cdim);
apad = zeros(cdim);

apad(1:adim(1),1:adim(2)) = a;
bpad(1:bdim(1),1:bdim(2)) = b;
ffta = fft2(apad);
fftb = fft2(bpad);
c = ifft2(ffta.*conj(fftb));