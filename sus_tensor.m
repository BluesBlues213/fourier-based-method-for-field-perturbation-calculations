function [A]= sus_tensor(PHI,THETA,B0,sus, AA,BB,CC)

% INPUT PHI between 0 and pi/2
%       THETA between 0 and pi/2
%       B0 in Tesla
%       sus in ppb
%       AA,BB,CC correspond to diagonal elements of 3x3 matrix.
gamma=42.576*10^6; %Hz/Tesla
sus=sus*10^-9;
[IX IY IZ]=meshgrid(-120:1:120);
cx=0;
cy=0;
cz=0;

[xx yy]=meshgrid(-120:1:120);
R = sqrt( (xx-cx).^2 + (yy-cy).^2);
R(R>10)=inf;
R(R~=inf)=sus;
R(R==inf)=0;
N=241;
for i=1:N;
    RM(:,:,i)=R;
end
MASK1=logical(RM);

PHI=pi-PHI; 
phi=ones(size(MASK1)).*(PHI);
theta=ones(size(MASK1)).*(THETA); 

%% field generation from tensor components AA,BB,CC

TH=sin(theta);
PA=((AA*(cos(phi).^2))+(BB*(sin(phi).^2)))*sus; 
PB=((-AA*cos(phi).*sin(phi)) + (BB*cos(phi).*sin(phi)) )*sus;
PC=phi~=0;
PC=PC*sus*CC;

PA=PA.*MASK1;
PB=PB.*MASK1;
PC=PC.*MASK1;

partONE=fftshift(fftn(fftshift( PA - ((cos(theta).^2).*(PA-PC)) ))).*(1/3);
JR=ifftshift(ifftn(ifftshift(partONE)));

FTB=fftshift(fftn(fftshift(PB.*TH)));
FTA=fftshift(fftn(fftshift(PA.*TH)));
FTC=fftshift(fftn(fftshift(PC.*cos(theta))));

kx=(IY-cy).*FTB;
ky=(IX-cx).*FTA;
kz=(IZ-cz).*FTC;
tot=kx+ky+kz;

partTWO= ((TH.*(IX-cx) + cos(theta).*IZ)./((IX-cx).^2 + (IY-cy).^2 +(IZ-cz).^2)) .* tot;
partTWO(isnan(partTWO))=0;
iftpt= ifftshift(ifftn(ifftshift(partTWO)));

IN=partONE-partTWO;
ISO=ifftshift(ifftn(ifftshift(IN)))*B0*gamma;
A=real(ISO);


%% figures
close all
figure;imagesc(squeeze(A(:,:,120)));axis square;colorbar;axis off;title('field')
