function [k]=covar(x, x0, model, c, vnu)
% [k]=Covar(x,x0,model,c)
% Fonction pour calculer les covariances avec des mod`eles spécifiées comme avec cokri
% La fonction calcule pour toutes les positions de x (n1) avec toutes les positions de x0 (n2)  K est donc n1xn2
% auteur D. Marcotte avril 2002

% here we define the equations for the various covariograms. Any new model
% can be added here.
nugget = @(h)      (h == 0);  %1 - Nugget
expon  = @(h)      exp(-h);   %2 - Exponential
gaus   = @(h)      exp(-(h).^2); %3 - Gaussian
spher  = @(h)      1 - (1.5 * min(h, 1) / 1 - 0.5 * (min(h, 1) / 1).^3); %4 - Spherical
cubic  = @(h)      1 - (7 * min(h, 1).^2 - 8.75 * min(h, 1).^3 + 3.5 * min(h, 1).^5 - 0.75 * min(h, 1).^7); %5 - Cubique
penta  = @(h)      1 - 22 / 3 * min(h, 1).^2 + 33 * min(h, 1).^4 - 77 / 2 * min(h, 1).^5 + 33 / 2 * min(h, 1).^7 - 11 / 2 * min(h, 1).^9 + 5 / 6 * min(h, 1).^11; %6 - Penta
cauchy = @(h, nu)  (h.^2 + 1).^(-nu); %7 - Cauchy with shape parameter nu
KBess  = @(h, nu)  1 / (gamma(nu) * 2^(nu - 1)) * max(h, eps).^nu .* besselk(nu, max(h, eps)); %8 - Matern
linear = @(h)      1 - h; %9 - linear
spline = @(h)      (h.^2).*log(max(h,eps)); %10 - Spline plaque mince
HoleS  = @(h)      sin(max(eps,h*2*pi))./max(eps,h*2*pi);  %11 - Hole Effect (Sinus)
HoleC  = @(h)      cos(h*2*pi); %12 - Hole Effect (cosinus)
chris  = @(h)      1-h.^2./(1+h.^2); %13 - Christakos
wand0  = @(h, nu)  (1-min(h,1)).^nu; % 14 Wendland 0 with shape parameter nu
wand1  = @(h, nu)  (1+h*(nu+1))*(1-min(h,1)).^(nu+1); % 15 Wendland 0 with shape parameter nu
wand2  = @(h, nu)  (1+h*(nu+2)/3+h.^2*(nu^2+4*nu+3))*(1-min(h,1)).^(nu+2); % 16 Wendland 0 with shape parameter nu
Bohman = @(h)      (h==0)+(h<1).*((1-h).*sin(2*pi*h)./(2*pi*max(h,eps))+(1-cos(2*pi*h))./(2*pi^2*max(h,eps))); % 17 Bohman covariance

Gam = {nugget, expon, gaus, spher, cubic, penta,  cauchy, KBess,...
    linear, spline, HoleS, HoleC, chris, wand0, wand1, wand2, Bohman};

% Definition of some constants
[n1,~]=size(x); % d dimension de l'espace
[n2,d]=size(x0);
[rp,p]=size(c);
r=rp/p;  % nombre de structures
nm=size(model,2);

% ne pas permettre des portées de 0 en input pour éviter divisions par 0
if nm>2
    model(:,2:1+d)=max(model(:,2:1+d),100*eps);
else
    model(:,2)=max(model(:,2),100*eps);
end


% calculer les covariances

k=zeros(n1*p,n2*p);
for i=1:r

    % calculation of matrix of reduced rotated distances H

    [t1]=trans(x(:,1:d),model,i);
    [t2]=trans(x0,model,i);
    h=0;
    for id=1:d
        h=h+(t1(:,id)*ones(1,n2)-ones(n1,1)*t2(:,id)').^2;
    end
    h=sqrt(h);
    ji=(i-1)*p+1; js=i*p ;

    % evaluation of the current basic structure
    gh = Gam{model(i, 1)};
    % Apply the covariance model

    % If Wendland, Cauchy or Matern
    if model(i, 1) == 7 || model(i, 1) == 8 || model(i, 1) == 14 || model(i, 1) == 15 || model(i, 1) == 16
        nu = vnu(i);
        g = gh(h, nu);
    else % Others
        g = gh(h);
    end

    k=k+kron(g,c(ji:js,:));
end
end

function [cx,rot]=trans(cx,model,im)
% function [cx,rot]=trans(cx,model,im);
%
% TRANS is called from COKRI2. It takes as input original coordinates and
%       return the rotated and reduced coordinates following specifications
%       described in model(im,:)
% Rotations are all performed anticlockwise with the observer located on the positive side of
% the axis and looking toward the origin. In 3D, rotations are performed first along z,
% then along rotated y and then along twice rotated x.
% Author: D. Marcotte
% Version 2.1  97/aug/18

% some constants are defined

[~,d]=size(cx);
[~,p]=size(model);

% check for 1-D or isotropic model

if p-1>d

    % perform rotation counterclockwise

    if d==2
        ang=model(im,4); cang=cos(ang/180*pi); sang=sin(ang/180*pi);
        rot=[cang,-sang;sang,cang];
    else

        % rotation matrix in 3-D is computed around z, y and x in that order

        angz=model(im,7); cangz=cos(angz/180*pi); sangz=sin(angz/180*pi);
        angy=model(im,6); cangy=cos(angy/180*pi); sangy=sin(angy/180*pi);
        angx=model(im,5); cangx=cos(angx/180*pi); sangx=sin(angx/180*pi);
        rotz=[cangz,-sangz,0;sangz,cangz,0;0 0 1];
        roty=[cangy,0,sangy;0 1 0;-sangy,0,cangy];
        rotx=[1 0 0;0 cangx -sangx;0 sangx cangx];
        rot=rotz*roty*rotx;
    end

    % rotation is performed around z, y and x in that order, the other coordinates are left unchanged.

    dm=min(3,d);
    cx(:,1:dm)=cx(:,1:dm)*rot;
    t=[model(im,2:1+dm),ones(d-dm,1)];
    t=diag(t);
else
    t=eye(d)*model(im,2);
end

% perform contractions or dilatations (reduced h)

cx=cx/t;
end