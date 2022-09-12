% A self-contained script that demonstrates smoothing optimization for the
% measurement of multiple noisy Gaussian peaks. Compares sliding-average,
% triangular, pseudo-Gaussian, and Savitzky-Golay smooths. 
% See http://terpconnect.umd.edu/~toh/spectrum/Smoothing.html and 

% Tom O'Haver (toh@umd.edu). Version 4, January 2021
% You can change the signal characteristics in lines 14-19
format short g
format compact
clear
increment=1;
x=[1:increment:18000]; % Change to 180000 to get 1000 peaks

% For each simulated peak, compute the amplitude, position, and width
Height=1; % For simplicity, the heights of all the peaks are equal to 1.00.
separation=150.01; % Separation between peaks, in x units
pos=[200:separation:17800];   % Positions of the peaks (Change to 178000 for 1000 peaks)
amp=Height.*ones(size(pos));  % Amplitudes of the peaks  (Change if desired)
wid=20.*ones(size(pos));   % Half-width (FWHM) of the peaks (Change if desired)
Noise=.01; % Amount of random noise added to the signal. (Change if desired) 

% A = matrix containing one of the unit-amplidude peak in each of its rows
A = zeros(length(pos),length(x));
ActualPeaks=[0 0 0 0 0];
p=1;
for k=1:length(pos)
    % Create a series of peaks of different x-positions
    A(k,:)=exp(-((x-pos(k))./(0.6005615.*wid(k))).^2); % Gaussian peaks, wid=halfwidth
    % A(k,:)=ones(size(x))./(1+((x-pos(k))./(0.5.*wid(k))).^2);  % Lorentzian peaks
    % Assembles actual parameters into ActualPeaks matrix: each row = 1
    % peak; columns are Peak #, Position, Height, Width, Area
    ActualPeaks(p,:) = [p pos(k) amp(k) wid(k) 1.0646.*amp(k)*wid(k)];
    p=p+1;
end
z=amp*A;  % Multiplies each row by the corresponding amplitude and adds them up
NoiseArray=Noise.*randn(size(z)); % Creates random white noise array
y=z+NoiseArray; % Adds constant random noise
% y=z+Noise.*sqrtnoise(z);  % Adds signal-dependent random noise
% y=y+5.*gaussian(x,0,4000); % Optionally adds a broad background signal
% y=y+gaussian(x,0,4000); % Optionally adds a broad background signal
% datamatrix=[x' y']; % Assembles x and y vectors into data matrix

d2y=-derivxy(x,derivxy(x,y));
NoiseArray=-derivxy(x,derivxy(x,NoiseArray));

figure(1);
clf
plot(x(100:2000),d2y(100:2000),'.')
title('Second derivative of raw data: The first 12 peaks of raw data')
xlabel('x');ylabel('Second derivative, d(d/dx)/dx')
grid
drawnow
    
% Smooth widths vector
SmoothWidth=[1 3 5 7 9 11 13 15 17 21 29 41 49]; 
lsw=length(SmoothWidth);
% Note: If you change the length of this vector, search and change "(2:9)" to match.

% A: sliding average 
tic
for trial=1:length(SmoothWidth)
    if SmoothWidth(trial)>1
        sy=fastsmooth(d2y,SmoothWidth(trial),1);
        SmoothedNoise=fastsmooth(NoiseArray,SmoothWidth(trial),1);
        stdnoiseA(trial)=std(SmoothedNoise);
    else
        sy=d2y;
    end
    stdnoiseA(1)=std(NoiseArray);
    % Determine half-width of the impulse response function
    for k=1:length(SmoothWidth)
        df=[zeros(1,100) 1 zeros(1,100)]; % Unit-height impulse
        ir=fastsmooth(df,SmoothWidth(k),1); % Impulse response function
        tswA(k)=halfwidth(1:length(df),ir,100); % half-width of impulse response function
        NoiseReductionA(k)=1./sum(ir.^2);
        NoiseReductionA(k)=sqrt(NoiseReductionA(k));
    end

    for n=1:length(pos)
        range=val2ind(x,pos(n)-wid(n)):val2ind(x,pos(n)+wid(n));
        mamp(n)=max(sy(range));
        pindex=range(1)+val2ind(sy(range),mamp(n));
        mpos(n)=pindex(1);
        mwid(n)=halfwidth(x,sy,round(pos(n)));
        peakrange=val2ind(x,pos(n)-2*wid(n)):val2ind(x,pos(n)+2*wid(n)); % range spanning each entire peak
        mareaA(n)=trapz(x(peakrange),sy(peakrange));
    end
    PercentPositionErrorA(trial)=100.*mean(abs(pos-mpos)./separation);
    meanmampA(trial)=mean(mamp);
    rsdmampA(trial)=rsd(mamp);
    PercentHeightErrorA(trial)=100.*mean((meanmampA-Height)./Height);
    PercentRSDamplitudeA(trial)=100.*rsd(meanmampA);
    PercentWidthErrorA(trial)=100.*mean((mwid-wid)/wid);
    PercentRSDwidthA(trial)=100.*rsd(mwid);
    AverageArea=1.0646.*mean(amp).*mean(wid);
    MeanAreaErrorA=mean(abs((1.0646.*amp.*wid)-mareaA));
    RelativeAreaErrorA(trial)=100.*(MeanAreaErrorA./AverageArea);
    PercentAreaRSDA(trial)=100.*rsd(mareaA);
end
elapsedtimeA=toc;

% B: Triangular (two passes of sliding average) 
tic
for trial=1:length(SmoothWidth)
    if SmoothWidth(trial)>1
        sy=fastsmooth(d2y,SmoothWidth(trial),2);
        SmoothedNoise=fastsmooth(NoiseArray,SmoothWidth(trial),2);
        stdnoiseB(trial)=std(SmoothedNoise);
    else
        sy=d2y;
    end
    stdnoiseB(1)=std(NoiseArray);
    % Determine half-width of the impulse response function
    for k=1:length(SmoothWidth)
        df=[zeros(1,100) 1 zeros(1,100)]; % Unit-height impulse
        ir=fastsmooth(df,SmoothWidth(k),2); % Impulse response function
        tswB(k)=halfwidth(1:length(df),ir,100); % half-width of impulse response function
        NoiseReductionB(k)=1./sum(ir.^2);
        NoiseReductionB(k)=sqrt(NoiseReductionB(k));
    end
    for n=1:length(pos)
        range=val2ind(x,pos(n)-wid(n)):val2ind(x,pos(n)+wid(n));
        mamp(n)=max(sy(range));
        pindex=range(1)+val2ind(sy(range),mamp(n));
        mpos(n)=pindex(1);
        mwid(n)=halfwidth(x,sy,round(pos(n)));
        peakrange=val2ind(x,pos(n)-2*wid(n)):val2ind(x,pos(n)+2*wid(n)); % range spanning each entire peak
        mareaB(n)=trapz(x(peakrange),sy(peakrange));
    end
    PercentPositionErrorB(trial)=std(pos-mpos);
    meanmampB(trial)=mean(mamp);
    rsdmampB(trial)=rsd(mamp);
    PercentHeightErrorB(trial)=100.*mean((meanmampB-Height)./Height);
    PercentRSDamplitudeB(trial)=100.*rsd(meanmampB);
    PercentWidthErrorB(trial)=100.*mean((mwid-wid)/wid);
    PercentRSDwidthB(trial)=100.*rsd(mwid);
    AverageArea=1.0646.*mean(amp).*mean(wid);
    MeanAreaErrorB=mean(abs((1.0646.*amp.*wid)-mareaB));
    RelativeAreaErrorB(trial)=100.*(MeanAreaErrorA./AverageArea);
    PercentAreaRSDB(trial)=100.*rsd(mareaB);
end
elapsedtimeB=toc;

% C: Gaussian (three passes of sliding average) 
tic
for trial=1:length(SmoothWidth)
    if SmoothWidth(trial)>1
        sy=fastsmooth(d2y,SmoothWidth(trial),3);
        SmoothedNoise=fastsmooth(NoiseArray,SmoothWidth(trial),3);
        stdnoiseC(trial)=std(SmoothedNoise);
    else
        sy=d2y;
    end
    stdnoiseC(1)=std(NoiseArray);
    % Determine half-width of the impulse response function
     for k=1:length(SmoothWidth)
        df=[zeros(1,100) 1 zeros(1,100)]; % Unit-height impulse
        ir=fastsmooth(df,SmoothWidth(k),3); % Impulse response function
        tswC(k)=halfwidth(1:length(df),ir,100); % half-width of impulse response function
        NoiseReductionC(k)=1./sum(ir.^2);
        NoiseReductionC(k)=sqrt(NoiseReductionC(k));
    end
    for n=1:length(pos)
        range=val2ind(x,pos(n)-wid(n)):val2ind(x,pos(n)+wid(n));
        mamp(n)=max(sy(range));
        pindex=range(1)+val2ind(sy(range),mamp(n));
        mpos(n)=pindex(1);
        mwid(n)=halfwidth(x,sy,round(pos(n)));
        peakrange=val2ind(x,pos(n)-2*wid(n)):val2ind(x,pos(n)+2*wid(n)); % range spanning each entire peak
        mareaC(n)=trapz(x(peakrange),sy(peakrange));
    end
    PercentPositionErrorC(trial)=100.*mean(abs(pos-mpos)./separation);
    meanmampC(trial)=mean(mamp);
    rsdmampC(trial)=rsd(mamp);
    PercentHeightErrorC(trial)=100.*mean((meanmampC-Height)./Height);
    PercentRSDamplitudeC(trial)=100.*rsd(meanmampC);
    PercentWidthErrorC(trial)=100.*mean((mwid-wid)/wid);
    PercentRSDwidthC(trial)=100.*rsd(mwid);
    AverageArea=1.0646.*mean(amp).*mean(wid);
    MeanAreaErrorC=mean(abs((1.0646.*amp.*wid)-mareaC));
    RelativeAreaErrorC(trial)=100.*(MeanAreaErrorC./AverageArea);
    PercentAreaRSDC(trial)=100.*rsd(mareaC);
end
elapsedtimeC=toc;

% D. Savitzky-Golay filter
tic
for trial=2:length(SmoothWidth)
    if SmoothWidth(trial)>1
        FrameSize=makeodd(2.*SmoothWidth);
        PolynomialOrder=2;
        DerivativeOrder=2;
        sy=-savitzkyGolayFilt(y,PolynomialOrder,DerivativeOrder,FrameSize(trial));
        SmoothedNoise=savitzkyGolayFilt(NoiseArray,PolynomialOrder,DerivativeOrder,FrameSize(trial));
        stdnoiseD(trial)=std(SmoothedNoise);
        % Determine half-width of the impulse response function
        for k=2:length(SmoothWidth)
            df=[zeros(1,100) 1 zeros(1,100)]; % Unit-height impulse
            ir=-savitzkyGolayFilt(df,PolynomialOrder,DerivativeOrder,FrameSize(k)); % Impulse response function
            tswD(k)=halfwidth(1:length(df),ir,101); % half-width of impulse response function
            NoiseReductionD(k)=1./sum((ir.^2));
            NoiseReductionD(k)=sqrt(NoiseReductionD(k));
        end
    else
        sy=y;
    end
    stdnoiseD(1)=std(NoiseArray);
    for n=1:length(pos)
        range=val2ind(x,pos(n)-wid(n)):val2ind(x,pos(n)+wid(n));
        mamp(n)=max(sy(range));
        pindex=range(1)+val2ind(sy( range),mamp(n));
        mpos(n)=pindex(1);
        mwid(n)=halfwidth(x,sy,round(pos(n)));
        peakrange=val2ind(x,pos(n)-2*wid(n)):val2ind(x,pos(n)+2*wid(n)); % range spanning each entire peak
        mareaD(n)=trapz(x(peakrange),sy(peakrange));
    end
    PercentPositionErrorD(trial)=100.*mean(abs(pos-mpos)./separation);
    meanmampD(trial)=mean(mamp);
    rsdmampD(trial)=rsd(mamp);
    PercentHeightErrorD(trial)=100.*mean((meanmampD-Height)./Height);
    PercentRSDamplitudeD(trial)=100.*rsd(meanmampD);
    PercentWidthErrorD(trial)=100.*mean((mwid-wid)/wid);
    PercentRSDwidthD(trial)=100.*rsd(mwid);
    AverageArea=1.0646.*mean(amp).*mean(wid);
    MeanAreaErrorD=mean(abs((1.0646.*amp.*wid)-mareaD));
    RelativeAreaErrorD(trial)=100.*(MeanAreaErrorA./AverageArea);
    PercentAreaRSDD(trial)=100.*rsd(mareaD);
end
elapsedtimeD=toc;


figure(4);
clf
plot(x(100:2000),sy(100:2000),'.')
title('Smoothed second derivative of raw data: The first 12 peaks of raw data')
xlabel('x');ylabel('y')
grid
drawnow

figure(2)
clf
subplot(2,2,1)
 plot(tswA,PercentPositionErrorA,tswB,PercentPositionErrorB,tswC,PercentPositionErrorC,tswD(2:lsw),PercentPositionErrorD(2:lsw),'linewidth',2)
 title('Peak position error')
 xlabel(['Total smooth halfwidth.    Peak width: ' num2str(mean(wid)) ])
 ylabel('Percent error')
 legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
subplot(2,2,2)
 plot(tswA,meanmampA,tswB,meanmampB,tswC,meanmampC,tswD(2:lsw),meanmampD(2:lsw),'linewidth',2)
 title('Negative second derivative peak height')
  xlabel(['Total smooth halfwidth.    Peak width: ' num2str(mean(wid)) ])
  legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
 subplot(2,2,3)
 plot(tswA,stdnoiseA,tswB,stdnoiseB,tswC,stdnoiseC,tswD(2:lsw),stdnoiseD(2:lsw),'linewidth',2)
 title('Standard deviation of the noise')
  xlabel(['Total smooth halfwidth.    Peak width: ' num2str(mean(wid)) ])
  legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
 subplot(2,2,4)
 plot(tswA,rsdmampA,tswB,rsdmampB,tswC,rsdmampC,tswD(2:lsw),rsdmampD(2:lsw),'linewidth',2)
title('Relative standard deviation of derivative peak height')
 xlabel(['Total smooth halfwidth.    Peak width: ' num2str(mean(wid)) ])
 legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
grid

% Plot noise reduction on the abscissa
figure(3)
clf
subplot(2,2,1)
 plot(NoiseReductionA,PercentPositionErrorA,NoiseReductionB,PercentPositionErrorB,NoiseReductionC,PercentPositionErrorC,NoiseReductionD(2:lsw),PercentPositionErrorD(2:lsw),'linewidth',2)
 title('Peak position error')
 xlabel('Noise reduction')
 ylabel('Percent error')
 legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
subplot(2,2,2)
 plot(NoiseReductionA,meanmampA,NoiseReductionB,meanmampB,NoiseReductionC,meanmampC,NoiseReductionD(2:lsw),meanmampD(2:lsw),'linewidth',2)
 title('Negative second derivative peak height')
 xlabel('Noise reduction')
  legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
 subplot(2,2,3)
 plot(NoiseReductionA,stdnoiseA,NoiseReductionB,stdnoiseB,NoiseReductionC,stdnoiseC,NoiseReductionD(2:lsw),stdnoiseD(2:lsw),'linewidth',2)
 title('Standard deviation of the noise')
 xlabel('Noise reduction')
  legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
 grid
 subplot(2,2,4)
 plot(NoiseReductionA,rsdmampA,NoiseReductionB,rsdmampB,NoiseReductionC,rsdmampC,NoiseReductionD(2:lsw),rsdmampD(2:lsw),'linewidth',2)
title('Relative standard deviation of derivative peak height')
xlabel('Noise reduction')
 legend('Moving average','Triangular','Gaussian','Savitsky-Golay')
grid

% Internal functions
function [FWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,y,xo)
% function [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y,xo) computes the full
% width at half maximum of the maximum of any shape peak that peaks near xo
% and has a zero baseline. If xo is omitted, it computes the halfwidth from
% the maximum y value. Not highly accurate if the function is too sparsely
% sampled. Also optionally returns the leading and trailing edge slopes
% (slope1 and slope 2, respectively) and the leading and trailing edge
% half-width-at-half-maximum, if those output arguments are incuded. 
% Tom O'Haver (toh@umd.edu) Version 2, July 2016
%
% Example 1:
% x=-5:.01:5;
% y=sinc(x).^2;
% FWHM=halfwidth(x,y,0) % FWHM of Central peak
% FWHM=halfwidth(x,y,1.5) % FWHM of Positive side peak
%
% Example 2:
% x=[0:.1:10];
% W=3; % W is the true half-width
% y=gaussian(x,5,W);
% FWHM=halfwidth(x,y,5)
% FWHM=halfwidth(x,y) % gives the same result since there is only one peak
%
% Example 3: % Returns NaN because signal does not contain the halfwidth
% x=[0:.1:1];
% y=exp(-(x).^2)
% halfwidth(x,y) 
% 
% Example 4: Also returns the leading and trailing edge slopes (+1 and -1)
% and the leading and trailing edge half-width-at-half-maximum
% x=[0:.1:10];
% y=triangle(x,5,1);
% [FWHM,slope1,slope2,x1,x2]=halfwidth(x,y)
%
% Example 5. Comparison of Gaussian (g) and BiGaussian peaks (bg)
% x=[1:100];
% g=gaussian(x,50,20);
% bg=bigaussian(x,50,20,3);
% plot(x,g,x,bg)
% [GaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,g)
% [BiGaussianFWHM,slope1,slope2,hwhm1,hwhm2]=halfwidth(x,bg)

if nargin==2,xo=x(val2ind(y,max(y)));end
try   
    indmax=val2ind(x,xo);
    maxy=y(indmax);
    oy=y-maxy/2;
    
    n=indmax;
    while oy(n)>0
        n=n-1;
    end
    x1=interp1([y(n) y(n+1)],[x(n) x(n+1)],maxy/2);
    slope1=(y(n+1)-y(n))./(x(n+1)-x(n));
    % polyfit((x(n+1)-x(n)),(y(n+1)-y(n)),1)
    
    n=indmax;
    while oy(n)>0
        n=n+1;
    end
    x2= interp1([y(n-1) y(n)],[x(n-1) x(n)],maxy/2);
    slope2=(y(n-1)-y(n))./(x(n-1)-x(n));
    % polyfit((x(n-1)-x(n)),(y(n-1)-y(n)),1)
    FWHM=x2-x1;
    hwhm1=xo-x1;
    hwhm2=x2-xo;
catch
    FWHM=NaN;
end
end

function [index,closestval]=val2ind(x,val)
% Returns the index and the value of the element of vector x that is closest to val
% If more than one element is equally close, returns vectors of indicies and values
% Tom O'Haver (toh@umd.edu) October 2006
% Examples: If x=[1 2 4 3 5 9 6 4 5 3 1], then val2ind(x,6)=7 and val2ind(x,5.1)=[5 9]
% [indices values]=val2ind(x,3.3) returns indices = [4 10] and values = [3 3]
dif=abs(x-val);
index=find((dif-min(dif))==0);
closestval=x(index);
end

function y=savitzkyGolayFilt(x,N,DN,F,W,DIM)
%savitzkyGolayFilt Savitzky-Golay Filtering.
%   savitzkyGolayFilt(X,N,DN,F) filters the signal X using a Savitzky-Golay 
%   (polynomial) filter.  The polynomial order, N, must be less than the
%   frame size, F, and F must be odd.  DN specifies the differentiation
%   order (DN=0 is smoothing). For a DN higher than zero, you'll have to
%   scale the output by 1/T^DN to acquire the DNth smoothed derivative of
%   input X, where T is the sampling interval. The length of the input X
%   must be >= F.  If X is a matrix, the filtering is done on the columns
%   of X.
%
%   Note that if the polynomial order N equals F-1, no smoothing
%   will occur.
%
%   savitzkyGolayFilt(X,N,DN,F,W) specifies a weighting vector W with
%   length F containing real, positive valued weights employed during the
%   least-squares minimization. If not specified, or if specified as
%   empty, W defaults to an identity matrix.
%
%   savitzkyGolayFilt(X,N,DN,F,[],DIM) or savitzkyGolayFilt(X,N,DN,F,W,DIM)
%   operates along the dimension DIM.
%
%   See also savitzkyGolay, FILTER, sgolayfilt

%   References:
%     [1] Sophocles J. Orfanidis, INTRODUCTION TO SIGNAL PROCESSING,
%              Prentice-Hall, 1995, Chapter 8.

%   Author(s): R. Losada
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.11.4.4 $  $Date: 2009/08/11 15:47:54 $

error(nargchk(4,6,nargin,'struct'));

% Check if the input arguments are valid
if round(F) ~= F, error(generatemsgid('MustBeInteger'),'Frame length must be an integer.'), end
if rem(F,2) ~= 1, error(generatemsgid('SignalErr'),'Frame length must be odd.'), end
if round(N) ~= N, error(generatemsgid('MustBeInteger'),'Polynomial order must be an integer.'), end
if N > F-1, error(generatemsgid('InvalidRange'),'The Polynomial order must be less than the frame length.'), end
if DN > N, error(generatemsgid('InvalidRange'),'The Differentiation order must be less than or equal to the Polynomial order.'), end

if nargin < 5 || isempty(W)
   % No weighting matrix, make W an identity
   W = ones(F,1);
else
   % Check for right length of W
   if length(W) ~= F, error(generatemsgid('InvalidDimensions'),'The weight vector must be of the same length as the frame length.'),end
   % Check to see if all elements are positive
   if min(W) <= 0, error(generatemsgid('InvalidRange'),'All the elements of the weight vector must be greater than zero.'), end
end

if nargin < 6, DIM = []; end

% Compute the projection matrix B
pp = fix(-F./2):fix(F./2);
B = savitzkyGolay(pp,N,DN,pp,W);

if ~isempty(DIM) && DIM > ndims(x)
	error(generatemsgid('InvalidDimensions'),'Dimension specified exceeds the dimensions of X.')
end

% Reshape X into the right dimension.
if isempty(DIM)
	% Work along the first non-singleton dimension
	[x, nshifts] = shiftdim(x);
else
	% Put DIM in the first dimension (this matches the order 
	% that the built-in filter function uses)
	perm = [DIM,1:DIM-1,DIM+1:ndims(x)];
	x = permute(x,perm);
end

if size(x,1) < F, error(generatemsgid('InvalidDimensions'),'The length of the input must be >= frame length.'), end

% Preallocate output
y = zeros(size(x));

% Compute the transient on (note, this is different than in sgolayfilt,
% they had an optimization leaving out some transposes that is only valid
% for DN==0)
y(1:(F+1)/2-1,:) = fliplr(B(:,(F-1)/2+2:end)).'*flipud(x(1:F,:));

% Compute the steady state output
ytemp = filter(B(:,(F-1)./2+1),1,x);
y((F+1)/2:end-(F+1)/2+1,:) = ytemp(F:end,:);

% Compute the transient off
y(end-(F+1)/2+2:end,:) = fliplr(B(:,1:(F-1)/2)).'*flipud(x(end-(F-1):end,:));

% Convert Y to the original shape of X
if isempty(DIM)
	y = shiftdim(y, -nshifts);
else
	y = ipermute(y,perm);
end
end