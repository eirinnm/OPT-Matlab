function [ fitted, offset ] = fit_sine( y )
%FIT_SINE sine fitting function from 
% https://uk.mathworks.com/matlabcentral/answers/121579-curve-fitting-to-a-sinusoidal-function
original_y = y;
y = sgolayfilt(original_y, 3, 9);
x = [1:size(y)]';
yu = max(y);
yl = min(y);
yr = (yu-yl)/2;                               % Range of ‘y’
yz = y-yu+(yr);  % y centred around the origin
z = yz .* circshift(yz,[0 1]);
zx = x(z <= 10);     % Find zero-crossings
zdiff = diff(zx);
per = 360; %2*mean(zdiff(zdiff>1));                     % Estimate period
ym = nanmean(y);                               % Estimate offset
%estimate the phase - this helps a lot with the fitting function
[~, max_ix] = max(y);
if max_ix < 90
    pha = 90 - max_ix;
else
    pha = 450 - max_ix;
end
fit = @(b,x)  b(1).*(sin(deg2rad(360*x./b(2) + b(3)))) + b(4);    % Function to fit
fcn = @(b) nansum((fit(b,x) - y).^2);                              % Least-Squares cost function
options = optimset('MaxFunEvals',1000);
[sinfit,fval,exitflag,output] = fminsearch(fcn, [yr;  per;  pha;  ym], options);                    % Minimise Least-Squares
 %fit to a sine curve, returns amplitude, period, phase, offset
% phase = sinfit(2)/(2*sinfit(3)); 
% phase=0;
if(exitflag==1)
    fprintf('Sine fitting complete after %d iterations\n',output.iterations);
else
    error('Sine fitting failed.');
end
phase=sinfit(3);
fitted = sin(deg2rad(phase + ([1:size(y,1)])*360/sinfit(2)))*sinfit(1);
vertical_offset = sinfit(4);
fitted = fitted + vertical_offset;
fitted = fitted';
offset = fitted-original_y; %how much do we need to warp the sinogram?
offset = fillmissing(offset, 'linear');
fprintf('Sine fitted with max divergence %0.2f pixels and std %0.4f\n',max(abs(offset)),std(offset));

end

