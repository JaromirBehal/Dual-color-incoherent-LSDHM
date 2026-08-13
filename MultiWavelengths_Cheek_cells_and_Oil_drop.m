clear all; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%         OPTICAL PARAMETERS        %%%%%%%%%%%%%%%%%%%%%%
%addpath('C:\Users\....\functions'); % FOLDER WITH USED FUNCTIONS

    %columns                % rows
    M=3096;                 N=2080;         % camera pixels
    m=-M/2:(M/2-1);         n=-N/2:(N/2-1);
    dx=2.4*10^(-6);         dy=dx;              % camera pixel size
    x=m.*dx;                y=n.*dy;        
    [X,Y]=meshgrid(x,y);
    lambdaG=532*10^(-9);    % central wavelength 532 nm
    lambdaB=450*10^(-9);    % central wavelength 450 nm

    fx=m/(M*dx);            fy=n/(N*dy);       % spatial frequencies (in image space)
    [FX,FY]=meshgrid(fx,fy);
    NA=0.25;                % numerical aperture of microscope objective (MO) 
    magMO=10;               % nominal magnification of the MO
    fTL_ideal=0.180;        % nominal focal length of ideal Tube lens (for which the MO is designed)
    fTL=0.300;              % focal length of used Tube lens
    beta=magMO*fTL/fTL_ideal*(150/200);     % overall lateral magnification of the sample-imaging optical path   

    xobj=x/beta*1000000;    % x dimensions in object space in micrometers
    yobj=y/beta*1000000; % y dimensions in object space in micrometers

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%      INTERFERENCE PATTERNS - CHEEK CELLS       %%%%%%%%%%%%%%%%

directory_images='C:\Users\....\data'; % folder with camera records

    % interference pattern with CHEEK cells
    I=apodization_for_propag(fitsread([directory_images '\2026-07-09-0742_9-CapObj_0000.FIT']));

    % interference pattern without sample
    Iref=apodization_for_propag(fitsread([directory_images '\2026-07-09-0741_5-CapObj_0000.FIT']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%      CPLX amplitude reconstruction      %%%%%%%%%%%%%%%%%%% 

    coef = 0.6; % to avoid overlap between neighbour holographic orders

    row_R1=1393;          col_R1=2141;
    U1=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R1).*X+fy(row_R1).*Y)),fx,fy,coef*(NA/lambdaB)/beta,0.5,24); % cplx. ampl for lambdaB
    U1ref=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R1).*X+fy(row_R1).*Y)),fx,fy,coef*(NA/lambdaB)/beta,0.5,24); % ref cplx. ampl for lambdaB
    U1corpom=U1./U1ref; % phase background correction
    U1cor=U1corpom.*exp(-1i*mean2(angle(U1corpom(1340:1500,1100:1250)))); % getting rid of residual offset

    row_R2=506+2;         col_R2=648;
    U2=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R2).*X+fy(row_R2).*Y)),fx,fy,coef*(NA/lambdaG)/beta,0.5,24); % cplx. ampl for lambdaG
    U2ref=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R2).*X+fy(row_R2).*Y)),fx,fy,coef*(NA/lambdaG)/beta,0.5,24); % ref cplx. ampl for lambdaG
    U2corpom=U2./U2ref; % phase background correction
    U2cor=U2corpom.*exp(-1i*mean2(angle(U2corpom(1340:1500,1100:1250)))); % getting rid of residual offset


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%              FIGURES              %%%%%%%%%%%%%%%%%%%%%%

figure(1); imagesc(xobj,yobj,I); axis image; colormap gray; colorbar; clim([3000 65000])
figure(2); imagesc(I(1600+(1:52),900+(1:52))); axis image; colormap gray; colorbar; 
figure(3); imagesc(fx*beta/1000,fy*beta/1000,log1p(abs(fftshift(fft2(I))))); axis image; colormap gray; clim([9 22]); colorbar
figure(4); imagesc(angle(U1cor)); axis image; colormap parula; colorbar; clim([-pi pi])
figure(5); imagesc(xobj,yobj,angle(U2cor)); axis image; colormap parula; colorbar; clim([-pi pi])








%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       INTERFERENCE PATTERNS - OIL      %%%%%%%%%%%%%%%%%%%%%%

    % interference pattern with drops of OIL
    I=apodization_for_propag(fitsread([directory_images '\2026-07-03-1309_7-CapObj_0000.FIT']));

    % interference pattern without sample
    Iref=apodization_for_propag(fitsread([directory_images '\2026-07-03-1310_5-CapObj_0000.FIT']));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%      CPLX amplitude reconstruction      %%%%%%%%%%%%%%%%%%%    
    
    row_R1=1393;          col_R1=2141;
    U1=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R1).*X+fy(row_R1).*Y)),fx,fy,coef*(NA/lambdaB)/beta,0.5,24); % cplx. ampl for lambdaB
    U1ref=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R1).*X+fy(row_R1).*Y)),fx,fy,coef*(NA/lambdaB)/beta,0.5,24); % ref cplx. ampl for lambdaB
    U1corpom=U1./U1ref; % phase background correction
    U1cor=U1corpom.*exp(-1i*mean2(angle(U1corpom(750:810,1950:2060)))); % getting rid of residual offset


    row_R2=506+2;         col_R2=648;
    U2=off_axis_reconstruction_general(I.*exp(-1i*2*pi*(fx(col_R2).*X+fy(row_R2).*Y)),fx,fy,coef*(NA/lambdaG)/beta,0.5,24); % cplx. ampl for lambdaG
    U2ref=off_axis_reconstruction_general(Iref.*exp(-1i*2*pi*(fx(col_R2).*X+fy(row_R2).*Y)),fx,fy,coef*(NA/lambdaG)/beta,0.5,24); % ref cplx. ampl for lambdaG
    U2corpom=U2./U2ref; % phase background correction
    U2cor=U2corpom.*exp(-1i*mean2(angle(U2corpom(750:810,1950:2060)))); % getting rid of residual offset

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%       1D cross-sections           %%%%%%%%%%%%%%%%%%%%%%

    coli=[590, 2200];
    rowi=[1725, 976];
    
    c1 = improfile(angle(U1cor),coli,rowi);
    c2 = improfile(angle(U2cor),coli,rowi);
    c3 = improfile(angle((U1cor)./(U2cor)),coli,rowi);
    c4 = improfile(angle((U1cor).*(U2cor)),coli,rowi);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%              FIGURES              %%%%%%%%%%%%%%%%%%%%%%

figure(11); 
imagesc(angle(U1cor)); axis image; colormap parula; colorbar; clim([-pi pi])
hold on
line(coli,rowi,'Color','k','LineStyle','-','Linewidth',3);

figure(12); 
imagesc(angle(U2cor)); axis image; colormap parula; colorbar; clim([-pi pi])
hold on
line(coli,rowi,'Color','k','LineStyle','-','Linewidth',3);

figure(13); 
imagesc(angle((U1cor)./(U2cor))); axis image; colormap parula; colorbar; clim([-pi pi])
hold on
line(coli,rowi,'Color','k','LineStyle','-','Linewidth',3);

figure(14); 
imagesc(angle((U1cor).*(U2cor))); axis image; colormap parula; colorbar; clim([-pi pi])
hold on
line(coli,rowi,'Color','k','LineStyle','-','Linewidth',3);

figure(15); 
subplot(4,1,1);plot(c1,'b-','linewidth',2); axis tight
hold on
subplot(4,1,2);plot(c2,'r-','linewidth',2); axis tight
hold on
subplot(4,1,3);plot(c3,'k-','linewidth',3); axis tight
hold on
subplot(4,1,4);plot(c4,'k-','linewidth',3); axis tight
