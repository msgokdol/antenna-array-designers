function data = reflectarrayPattern3D(A,phaseDeg,targetThetaDeg,targetPhiDeg)
%REFLECTARRAYPATTERN3D Sample the normalized upper-hemisphere field.
%   The regular display mesh is refined around the commanded direction
%   and includes that direction exactly, even for narrow off-grid beams.

if ~any(A.active(:))
    data = struct('x',[],'y',[],'z',[],'db',[], ...
        'thetaDeg',[],'phiDeg',[],'fieldMagnitude',[]);
    return;
end
apertureSpanLambda = max([max(A.X(:))-min(A.X(:)) ...
    max(A.Y(:))-min(A.Y(:))]) / ...
    A.operatingLambdaMM;
thetaHalfWidth = min(4,rad2deg(2/max(apertureSpanLambda,1)));
phiHalfWidth = min(8,thetaHalfWidth / ...
    max(sind(targetThetaDeg),0.25));
thetaLocal = max(0,min(90, ...
    targetThetaDeg+linspace(-thetaHalfWidth,thetaHalfWidth,21)));
phiLocal = mod(targetPhiDeg+ ...
    linspace(-phiHalfWidth,phiHalfWidth,21),360);
theta = unique([0:2:90,thetaLocal]);
% An elevation cut includes the commanded phi plane and its opposite
% half-plane. Keep both in the 3D mesh even when the regular phi grid
% misses one of them.
phi = unique([0:4:360,phiLocal,mod(targetPhiDeg+180,360)]);
[TH,PH] = meshgrid(theta,phi);
fieldMagnitude = abs(reflectarrayField(A,phaseDeg,TH,PH));
db = 20*log10(max(fieldMagnitude,1e-8));
db = db-max(db(:));
radius = max(0.03,1+max(db,-40)/40);
data = struct('x',radius.*sind(TH).*cosd(PH), ...
    'y',radius.*sind(TH).*sind(PH), ...
    'z',radius.*cosd(TH),'db',db, ...
    'thetaDeg',theta,'phiDeg',phi, ...
    'fieldMagnitude',fieldMagnitude);
end
