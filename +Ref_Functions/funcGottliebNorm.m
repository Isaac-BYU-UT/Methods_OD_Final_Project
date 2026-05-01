%% Function to Calculate the Acceleration due to Earth's Spherical Harmonics
% March 2026
% Jack Li, jack.li@utexas.edu
%%

% Calculates the acceleration due to Earth's gravity using NASA's normalised Gottlieb algorithm
%% Inputs
% mu: the standard gravitational parameter of the Earth in km^3/s^2
% re: the radius of the Earth in km
% xin: the state vector of the spacecraft in position and velocity [6x1] in km
% C: the C matrix for the spherical harmonics model
% S: the S matrix for the spherical harmonics model
% nax: the maximum degree for the spherical harmonics model in n
% max: the maximum degree for the spherical harmonics model in m
% rnp: rotation matrix from ECEF to ECI, set to identity since rotation is calculated externally
%% Outputs
% aGravity: the acceleration vector due to two-body and spherical harmonics in x, y, and z
function aGravity = funcGottliebNorm(mu, re, xin, C, S, nax, max, rnp)
for n = 2:nax+1
    norm1(n) = sqrt((2*n+1)/(2*n-1));
    norm2(n) = sqrt((2*n+1)/(2*n-3));
    norm11(n) = sqrt((2*n+1)/(2*n))/(2*n-1);
    normn10(n) = sqrt((n+1)*n/2);
    for m = 1:n %RAE
        norm1m(n,m) = sqrt((n-m)*(2*n+1)/((n+m)*(2*n-1))); %RAE
        norm2m(n,m) = sqrt((n-m)*(n-m-1)*(2*n+1)/((n+m)*(n+m-1)*(2*n-3))); %RAE
        normn1(n,m) = sqrt((n+m+1)*(n-m)); %RAE
    end %RAE
end %RAE
x=rnp*xin; %RAE
r = sqrt(x(1)^2+x(2)^2+x(3)^2);
ri=1/r;
xor=x(1)*ri;
yor=x(2)*ri;
zor=x(3)*ri;
ep=zor;
reor=re*ri;
reorn=reor;
muor2=mu*ri*ri;
p(1,1) = 1; %RAE
p(1,2) = 0; %RAE
p(1,3) = 0; %RAE
p(2,2) = sqrt(3); %RAE %norm
p(2,3) = 0; %RAE
p(2,4) = 0; %RAE
for n = 2:nax %RAE
    ni = n+1; %RAE
    p(ni,ni) = norm11(n)*p(n,n)*(2*n-1); %RAE %norm
    p(ni,ni+1) = 0; %RAE
    p(ni,ni+2) = 0; %RAE
end
ctil(1)=1; %RAE
stil(1)=0; %RAE
ctil(2)=xor; %RAE
stil(2)=yor; %RAE
sumh=0;
sumgm=1;
sumj=0;
sumk=0;
p(2,1) = sqrt(3)*ep; %RAE %norm
for n=2:nax
    ni=n+1; %RAE
    reorn=reorn*reor;
    n2m1=n+n-1;
    nm1=n-1;
    np1=n+1;
    p(ni,n) = normn1(n,n-1)*ep*p(ni,ni); %RAE %norm
    p(ni,1) = (n2m1*ep*norm1(n)*p(n,1)-nm1*norm2(n)*p(nm1,1))/n; %RAE %norm
    p(ni,2) = (n2m1*ep*norm1m(n,1)*p(n,2)-n*norm2m(n,1)*p(nm1,2))/(nm1); %RAE %norm
    sumhn=normn10(n)*p(ni,2)*C(ni,1); %norm %RAE
    sumgmn=p(ni,1)*C(ni,1)*np1; %RAE
if (max>0)
    for m = 2:n-2
        mi = m+1; %RAE
        p(ni,mi) = (n2m1*ep*norm1m(n,m)*p(n,mi)-...
        (nm1+m)*norm2m(n,m)*p(nm1,mi))/(n-m); %RAE %norm
    end
    sumjn=0;
    sumkn=0;
    ctil(ni)=ctil(2)*ctil(ni-1)-stil(2)*stil(ni-1); %RAE
    stil(ni)=stil(2)*ctil(ni-1)+ctil(2)*stil(ni-1); %RAE
    if(n<max)
        lim=n;
    else
        lim=max;
    end
    for m=1:lim
        mi=m+1; %RAE
        mm1=mi-1; %RAE
        mp1=mi+1; %RAE
        mxpnm=m*p(ni,mi); %RAE
        bnmtil=C(ni,mi)*ctil(mi)+S(ni,mi)*stil(mi); %RAE
        sumhn=sumhn+normn1(n,m)*p(ni,mp1)*bnmtil; %RAE %norm
        sumgmn=sumgmn+(n+m+1)*p(ni,mi)*bnmtil; %RAE
        bnmtm1=C(ni,mi)*ctil(mm1)+S(ni,mi)*stil(mm1); %RAE
        anmtm1=C(ni,mi)*stil(mm1)-S(ni,mi)*ctil(mm1); %RAE
        sumjn=sumjn+mxpnm*bnmtm1;
        sumkn=sumkn-mxpnm*anmtm1;
    end
    sumj=sumj+reorn*sumjn;
    sumk=sumk+reorn*sumkn;
end
    sumh = sumh+reorn*sumhn;
    sumgm = sumgm+reorn*sumgmn;
end
lambda=sumgm+ep*sumh;
g(1,1)=-muor2*(lambda*xor-sumj);
g(2,1)=-muor2*(lambda*yor-sumk);
g(3,1)=-muor2*(lambda*zor-sumh);
aGravity=rnp*g; %RAE
end
