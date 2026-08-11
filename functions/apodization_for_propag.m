function [ Im_out ] = apodization_for_propag( Im_in )

% apodization often used for propagation
% input: non apodized image
% output: apodized image

dim=size(Im_in);
Nx=dim(1);Ny=dim(2);
a=0;
y1=0; y2=Nx+1; x1=0; x2=Ny+1;
for j=1:25;
    a=a+0.04;
    y1=y1+1; y2=y2-1; x1=x1+1; x2=x2-1;
    Im_in(y1,:)=Im_in(y1,:).*a;
    Im_in(y2,:)=Im_in(y2,:).*a;
    Im_in(:,x1)=Im_in(:,x1).*a;
    Im_in(:,x2)=Im_in(:,x2).*a;
end
Im_out=Im_in;
end
