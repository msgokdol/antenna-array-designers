function E = reflectarrayField(A,phaseDeg,thetaDeg,phiDeg)
%REFLECTARRAYFIELD Scalar upper-hemisphere pattern of an illuminated surface.
%   E is normalized to the sum of active reflected-field amplitudes. It is
%   not absolute gain. thetaDeg and phiDeg must be equally sized arrays.

if ~isequal(size(phaseDeg),size(A.X))
    error('Reflectarray:PhaseSize','Phase map must match the cell grid.');
end
if ~isequal(size(thetaDeg),size(phiDeg))
    error('Reflectarray:DirectionSize','Theta and phi grids must match.');
end
if ~isnumeric(phaseDeg) || ~isreal(phaseDeg) || ...
        ~isnumeric(thetaDeg) || ~isreal(thetaDeg) || ...
        ~isnumeric(phiDeg) || ~isreal(phiDeg)
    error('Reflectarray:RealInput', ...
        'Cell phases and directions must be real numeric arrays.');
end
if any(~isfinite(phaseDeg(A.active))) || ...
        any(~isfinite(thetaDeg(:))) || any(~isfinite(phiDeg(:)))
    error('Reflectarray:Finite','Active phases and directions must be finite.');
end

U = sind(thetaDeg).*cosd(phiDeg);
V = sind(thetaDeg).*sind(phiDeg);
E = complex(zeros(size(U)));
idx = find(A.active & A.amplitude > 0);
if isempty(idx), return; end
for j = 1:numel(idx)
    i = idx(j);
    arg = -A.k*A.incidentPathMM(i) + deg2rad(phaseDeg(i)) + ...
        A.k*(A.X(i)*U + A.Y(i)*V);
    E = E + A.amplitude(i)*exp(1i*arg);
end
E = E/sum(A.amplitude(idx));
E(thetaDeg < 0 | thetaDeg > 90) = 0;
end
