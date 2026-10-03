function maps = reflectarrayPhaseMaps(A,phaseDeg,S,samples)
%REFLECTARRAYPHASEMAPS Phase-law terms and a continuous reference.
%   All phases are wrapped to [0,360) degrees at the frequency at which
%   the cell phases were set. The continuous map uses the nominal surface
%   outline, so it remains a reference when individual cells are edited.

if nargin < 4, samples = 181; end
if ~isnumeric(samples) || ~isreal(samples) || ...
        ~isscalar(samples) || ~isfinite(samples) || ...
        samples < 2 || mod(samples,1) ~= 0
    error('Reflectarray:PhaseMap','Sample count must be an integer >= 2.');
end
if ~isnumeric(phaseDeg) || ~isreal(phaseDeg) || ...
        ~isequal(size(phaseDeg),size(A.active)) || ...
        any(~isfinite(phaseDeg(A.active)))
    error('Reflectarray:PhaseMap','Cell phase map size must match the surface.');
end
k = 2*pi*A.phaseReferenceGHz/299.792458;
ux = sind(S.beamThetaDeg)*cosd(S.beamPhiDeg);
uy = sind(S.beamThetaDeg)*sind(S.beamPhiDeg);

% Incident field phase, desired outgoing phase gradient, and
% actual cell reflection phase (which may have been edited by the user).
maps.spatialDeg = mod(-rad2deg(k*A.incidentPathMM),360);
maps.progressiveDeg = mod(-rad2deg(k*(A.X*ux+A.Y*uy)),360);
maps.reflectionDeg = mod(phaseDeg,360);

x = linspace(min(A.baseX(:))-A.dxMM/2, ...
    max(A.baseX(:))+A.dxMM/2,samples);
y = linspace(min(A.baseY(:))-A.dyMM/2, ...
    max(A.baseY(:))+A.dyMM/2,samples);
[X,Y] = meshgrid(x,y);
if strcmp(A.illuminationMode,'Horn')
    F = A.feedMM;
    incidentPathMM = sqrt((X-F(1)).^2+(Y-F(2)).^2+F(3)^2);
else
    uiX = sind(A.incidentThetaDeg)*cosd(A.incidentPhiDeg);
    uiY = sind(A.incidentThetaDeg)*sind(A.incidentPhiDeg);
    incidentPathMM = -(X*uiX+Y*uiY);
end
maps.continuousDeg = mod(rad2deg(k*(incidentPathMM-X*ux-Y*uy)),360);
maps.continuousMask = reflectarrayShapeMask(X,Y,S.shape);
if ~any(A.active(:)), maps.continuousMask(:) = false; end
maps.continuousXMM = x;
maps.continuousYMM = y;
end
