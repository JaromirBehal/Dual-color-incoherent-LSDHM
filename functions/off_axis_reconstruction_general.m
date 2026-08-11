function [Urec0,F0,aper_four]=off_axis_reconstruction_general(I0,fx,fy,radius_freq,threshold,order)
% I0 ... recorded hologram (CENTERED VALUABLE DIFFRACTION ORDER !!! )


aper_four=apodiz_SG(fx,fy,radius_freq,threshold,order);

% Spectra
F0=ifftshift(fft2(I0)).*aper_four;
Urec0=(ifft2(ifftshift(F0)));

end