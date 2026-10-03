function A = reflectarrayAperture(P)
%REFLECTARRAYAPERTURE Geometry, incident illumination, and cell phases.
%   The surface is in the xy plane. All positions are in millimetres and
%   theta is measured from +z. This is an ideal scalar, phase-only model;
%   no physical unit-cell reflection coefficient or absolute gain is used.

required = {'M','N','dxLambda','dyLambda','designFreqGHz', ...
    'opFreqGHz','retunePhase','feedMM','feedQ', ...
    'beamThetaDeg','beamPhiDeg','shape','taper'};
for j = 1:numel(required)
    if ~isfield(P,required{j})
        error('Reflectarray:MissingSetting','Missing setting: %s.',required{j});
    end
end
if ~isnumeric(P.M) || ~isnumeric(P.N) || ...
        ~isreal(P.M) || ~isreal(P.N) || ...
        ~isscalar(P.M) || ~isscalar(P.N) || ...
        any(~isfinite([P.M P.N])) || ...
        any([P.M P.N] < 1) || any(mod([P.M P.N],1) ~= 0)
    error('Reflectarray:Grid','M and N must be positive integers.');
end
spacingAndFrequency = [P.dxLambda P.dyLambda ...
    P.designFreqGHz P.opFreqGHz];
if ~isnumeric(spacingAndFrequency) || ~isreal(spacingAndFrequency) || ...
        numel(spacingAndFrequency) ~= 4 || ...
        any(~isfinite(spacingAndFrequency)) || ...
        any(spacingAndFrequency <= 0)
    error('Reflectarray:Dimensions','Spacing and frequency must be positive.');
end
if ~isscalar(P.retunePhase) || ...
        ~ismember(P.retunePhase,[false true])
    error('Reflectarray:PhaseMode','Phase mode must be on or off.');
end
if ~isnumeric(P.feedMM) || ~isreal(P.feedMM) || ...
        numel(P.feedMM) ~= 3 || any(~isfinite(P.feedMM)) || P.feedMM(3) <= 0
    error('Reflectarray:Feed','The feed must have finite x/y and z > 0.');
end
if ~isnumeric(P.feedQ) || ~isreal(P.feedQ) || ...
        ~isscalar(P.feedQ) || ~isfinite(P.feedQ) || P.feedQ < 0
    error('Reflectarray:FeedPattern','Feed exponent q must be nonnegative.');
end
illuminationMode = settingOr(P,'illuminationMode','Horn');
incidentThetaDeg = settingOr(P,'incidentThetaDeg',0);
incidentPhiDeg = settingOr(P,'incidentPhiDeg',0);
if ~any(strcmp(illuminationMode,{'Horn','Plane wave'})) || ...
        ~isnumeric(incidentThetaDeg) || ~isreal(incidentThetaDeg) || ...
        ~isscalar(incidentThetaDeg) || ~isfinite(incidentThetaDeg) || ...
        incidentThetaDeg < 0 || incidentThetaDeg > 85 || ...
        ~isnumeric(incidentPhiDeg) || ~isreal(incidentPhiDeg) || ...
        ~isscalar(incidentPhiDeg) || ~isfinite(incidentPhiDeg) || ...
        incidentPhiDeg < 0 || incidentPhiDeg > 360
    error('Reflectarray:Illumination','Invalid illumination mode or incidence angle.');
end
if ~isnumeric(P.beamThetaDeg) || ~isnumeric(P.beamPhiDeg) || ...
        ~isreal(P.beamThetaDeg) || ~isreal(P.beamPhiDeg) || ...
        ~isscalar(P.beamThetaDeg) || ~isscalar(P.beamPhiDeg) || ...
        any(~isfinite([P.beamThetaDeg P.beamPhiDeg])) || ...
        P.beamThetaDeg < 0 || P.beamThetaDeg > 90
    error('Reflectarray:Beam','Beam theta must be between 0 and 90 degrees.');
end

designLambdaMM = 299.792458 / P.designFreqGHz;
operatingLambdaMM = 299.792458 / P.opFreqGHz;
dxMM = P.dxLambda*designLambdaMM;
dyMM = P.dyLambda*designLambdaMM;
k = 2*pi/operatingLambdaMM;
if P.retunePhase
    phaseReferenceGHz = P.opFreqGHz;
else
    phaseReferenceGHz = P.designFreqGHz;
end
kPhase = 2*pi*phaseReferenceGHz/299.792458;
x = ((1:P.N) - (P.N+1)/2)*dxMM;
y = ((1:P.M) - (P.M+1)/2)*dyMM;
[baseX,baseY] = meshgrid(x,y);
active = reflectarrayShapeMask(baseX,baseY,P.shape);
if isfield(P,'activeMask') && ~isempty(P.activeMask)
    if ~isequal(size(P.activeMask),[P.M P.N]) || ...
            ~(islogical(P.activeMask) || isnumeric(P.activeMask)) || ...
            any(~isfinite(double(P.activeMask(:)))) || ...
            any(~ismember(P.activeMask(:),[0 1]))
        error('Reflectarray:ActiveMask','Cell mask must be an M-by-N logical grid.');
    end
    active = logical(P.activeMask);
end
gridAngle = settingOr(P,'gridAngleDeg',90);
rowStagger = settingOr(P,'rowStaggerLambda',0);
staggerAngle = settingOr(P,'staggerAngleDeg',0);
offsetX = settingOr(P,'xOffsetLambda',zeros(P.M,P.N));
offsetY = settingOr(P,'yOffsetLambda',zeros(P.M,P.N));
if ~isnumeric(offsetX) || ~isnumeric(offsetY) || ...
        ~isreal(offsetX) || ~isreal(offsetY) || ...
        ~isequal(size(offsetX),[P.M P.N]) || ...
        ~isequal(size(offsetY),[P.M P.N]) || ...
        any(~isfinite([offsetX(:);offsetY(:)])) || ...
        any(~isfinite([gridAngle rowStagger staggerAngle])) || ...
        gridAngle <= 0 || gridAngle >= 180
    error('Reflectarray:Lattice','Invalid lattice angle, stagger, or offset.');
end
% Brick lattice: every second row shifts along x. The displayed stagger
% angle is atan2(stagger,dy); it is not another geometric rotation.
rowShiftMM = repmat(mod((1:P.M).',2) == 0,1,P.N) * ...
    rowStagger*designLambdaMM;
X = baseX + baseY*cosd(gridAngle) + ...
    rowShiftMM + offsetX*designLambdaMM;
Y = baseY*sind(gridAngle) + ...
    offsetY*designLambdaMM;

F = reshape(P.feedMM,1,3);
R = sqrt((X-F(1)).^2 + (Y-F(2)).^2 + F(3)^2);
theta = P.beamThetaDeg;
phi = P.beamPhiDeg;
ux = sind(theta)*cosd(phi);
uy = sind(theta)*sind(phi);
if strcmp(illuminationMode,'Horn')
    incidentPathMM = R;
    % The horn is aimed at the aperture centre. cos^q is a FIELD-pattern
    % model, not a power-pattern exponent; free-space field falls as 1/R.
    toCellX = X-F(1); toCellY = Y-F(2); toCellZ = -F(3);
    toCentre = -F/norm(F);
    cosFeedAngle = (toCellX*toCentre(1) + toCellY*toCentre(2) + ...
        toCellZ*toCentre(3))./R;
    illumination = max(cosFeedAngle,0).^P.feedQ ./ R;
    illumination(cosFeedAngle <= 0) = 0;
else
    % Incidence angles specify the direction TO the source, measured from
    % +z. The incoming propagation direction is its negative. A plane
    % wave has uniform field amplitude and phase k*r_transverse*u_inc.
    uiX = sind(incidentThetaDeg)*cosd(incidentPhiDeg);
    uiY = sind(incidentThetaDeg)*sind(incidentPhiDeg);
    incidentPathMM = -(X*uiX+Y*uiY);
    illumination = ones(P.M,P.N);
end
pathDifference = incidentPathMM-X*ux-Y*uy;
idealPhaseDeg = mod(rad2deg(k*pathDifference),360);
referencePhaseDeg = mod(rad2deg(kPhase*pathDifference),360);

switch char(P.taper)
    case 'None'
        taper = ones(P.M,P.N);
    case 'Hann (ideal amplitude)'
        % Bin-centred raised cosine avoids forced zero edge cells.
        wx = 0.5 - 0.5*cos(2*pi*((1:P.N)-0.5)/P.N);
        wy = 0.5 - 0.5*cos(2*pi*((1:P.M)-0.5)/P.M);
        taper = wy(:)*wx(:).';
        taper = taper/max(taper(:));
    otherwise
        error('Reflectarray:Taper','Unknown amplitude taper.');
end
amplitude = illumination.*taper;
amplitude(~active) = 0;
peak = max(amplitude(:));
if peak == 0 && any(active(:))
    error('Reflectarray:Unlit','The incident field does not illuminate the surface.');
end
if peak > 0, amplitude = amplitude/peak; end
illuminationPeak = max(illumination(:));
if illuminationPeak > 0
    normalizedIllumination = illumination/illuminationPeak;
else
    normalizedIllumination = zeros(size(illumination));
end

A = struct('X',X,'Y',Y,'x',x,'y',y,'active',active,'R',R, ...
    'incidentPathMM',incidentPathMM, ...
    'idealPhaseDeg',idealPhaseDeg, ...
    'referencePhaseDeg',referencePhaseDeg,'amplitude',amplitude, ...
    'illumination',normalizedIllumination, ...
    'lambdaMM',designLambdaMM,'designLambdaMM',designLambdaMM, ...
    'operatingLambdaMM',operatingLambdaMM, ...
    'dxMM',dxMM,'dyMM',dyMM,'k',k,'feedMM',F, ...
    'illuminationMode',illuminationMode, ...
    'incidentThetaDeg',incidentThetaDeg, ...
    'incidentPhiDeg',incidentPhiDeg, ...
    'phaseReferenceGHz',phaseReferenceGHz, ...
    'beamThetaDeg',theta,'beamPhiDeg',phi, ...
    'baseX',baseX,'baseY',baseY);
end

function value = settingOr(P,name,fallback)
if isfield(P,name), value = P.(name); else, value = fallback; end
end
