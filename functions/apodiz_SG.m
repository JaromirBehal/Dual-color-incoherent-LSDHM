function aperture=apodiz_SG(x,y,radius,value,exponent)
% Super-Gaussian apodization filter 'apodiz_SG' with values in range [0,1]
% 'x','y','radius' must be in the same units

% 'x','y' ...  coordinate vectors
% 'radius' ... sets area of radius 'radius' with borders depending on 'value'
% 'value'  ... determines decay where 'apodiz_SG'=='value'
% 'exponent' ... order of SG function

[X,Y]=meshgrid(x,y);
sigma = sqrt(-(radius^2)^exponent/(2*log(value)));
aperture = exp(-(X.^2+Y.^2).^exponent/(2*sigma^2));

end
