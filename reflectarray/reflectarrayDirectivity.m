function D = reflectarrayDirectivity(A,phaseDeg)
%REFLECTARRAYDIRECTIVITY Exact normalization for this scalar aperture model.
%   The existing field is zero below the xy plane. For two cells in that
%   plane, the upper-hemisphere integral of their interference term is
%   2*pi*sin(k*rho)/(k*rho), with the limit 2*pi at rho = 0. Summing those
%   pair terms gives the radiated-power integral without angular sampling.
%   D.factor converts |reflectarrayField|^2 to linear directivity. This
%   model does not include unit-cell patterns, spillover, or losses.

if ~isnumeric(phaseDeg) || ~isreal(phaseDeg) || ...
        ~isequal(size(phaseDeg),size(A.X)) || ...
        any(~isfinite(phaseDeg(A.active)))
    error('Reflectarray:Directivity','Cell phases must be a finite real grid.');
end
idx = find(A.active & A.amplitude > 0);
if isempty(idx)
    D = struct('factor',NaN,'pairPower',0);
    return;
end
x = A.X(idx);
y = A.Y(idx);
weights = A.amplitude(idx);
excitation = weights .* exp(1i*(deg2rad(phaseDeg(idx)) - ...
    A.k*A.incidentPathMM(idx)));
pairPower = 0;
blockSize = 256;
for first = 1:blockSize:numel(idx)
    rows = first:min(first+blockSize-1,numel(idx));
    kr = A.k*hypot(x(rows)-x.',y(rows)-y.');
    kernel = ones(size(kr));
    nonzero = kr > 1e-10;
    kernel(nonzero) = sin(kr(nonzero))./kr(nonzero);
    pairPower = pairPower + real(excitation(rows).' * ...
        (kernel*conj(excitation)));
end
if ~isfinite(pairPower)
    error('Reflectarray:Directivity','Radiated-power integral is invalid.');
end
if pairPower <= 1e-12*sum(weights)^2
    % Coincident cells with opposing phases can cancel every direction.
    D = struct('factor',NaN,'pairPower',max(pairPower,0));
    return;
end
D = struct('factor',2*sum(weights)^2/pairPower, ...
    'pairPower',pairPower);
end
