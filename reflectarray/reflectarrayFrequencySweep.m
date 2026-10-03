function result = reflectarrayFrequencySweep(P,phaseDeg,frequenciesGHz)
%REFLECTARRAYFREQUENCYSWEEP Evaluate one fixed phase map across frequency.
%   The physical lattice, source, and assigned cell phases are held fixed.
%   Each sample recomputes propagation and scalar modeled directivity; it
%   does not model frequency-dependent cell reflection or realized gain.

if ~isnumeric(frequenciesGHz) || ~isreal(frequenciesGHz) || ...
        ~isvector(frequenciesGHz) || isempty(frequenciesGHz) || ...
        any(~isfinite(frequenciesGHz(:))) || ...
        any(frequenciesGHz(:) <= 0)
    error('Reflectarray:SweepFrequency', ...
        'Sweep frequencies must be a finite positive real vector.');
end
frequenciesGHz = reshape(frequenciesGHz,1,[]);
base = reflectarrayAperture(P);
if ~isnumeric(phaseDeg) || ~isreal(phaseDeg) || ...
        ~isequal(size(phaseDeg),size(base.X)) || ...
        any(~isfinite(phaseDeg(base.active)))
    error('Reflectarray:SweepPhase', ...
        'Fixed cell phases must match the surface and be finite.');
end

count = numel(frequenciesGHz);
[thetaGrid,phiGrid] = meshgrid(0:2:90,0:4:356);
result = struct('frequencyGHz',frequenciesGHz, ...
    'targetDbi',nan(1,count),'peakDbi',nan(1,count), ...
    'peakThetaDeg',nan(1,count),'peakPhiDeg',nan(1,count), ...
    'pointingErrorDeg',nan(1,count));
for n = 1:count
    settings = P;
    settings.opFreqGHz = frequenciesGHz(n);
    aperture = reflectarrayAperture(settings);
    D = reflectarrayDirectivity(aperture,phaseDeg);
    if ~isfinite(D.factor), continue; end

    targetField = abs(reflectarrayField(aperture,phaseDeg, ...
        P.beamThetaDeg,P.beamPhiDeg));
    result.targetDbi(n) = 10*log10(max( ...
        D.factor*targetField^2,realmin));

    magnitude = abs(reflectarrayField(aperture,phaseDeg, ...
        thetaGrid,phiGrid));
    [peakField,index] = max(magnitude(:));
    peakTheta = thetaGrid(index);
    peakPhi = phiGrid(index);
    % Include the exact target so a coarse global grid cannot report a
    % sampled peak below the field at the commanded direction.
    if targetField > peakField
        peakField = targetField;
        peakTheta = P.beamThetaDeg;
        peakPhi = P.beamPhiDeg;
    end
    thetaLocal = max(0,peakTheta-2):0.25:min(90,peakTheta+2);
    phiLocal = mod(peakPhi+(-4:0.5:4),360);
    [thetaFine,phiFine] = meshgrid(thetaLocal,phiLocal);
    fineMagnitude = abs(reflectarrayField(aperture,phaseDeg, ...
        thetaFine,phiFine));
    [finePeak,fineIndex] = max(fineMagnitude(:));
    if finePeak > peakField
        peakField = finePeak;
        peakTheta = thetaFine(fineIndex);
        peakPhi = phiFine(fineIndex);
    end

    result.peakDbi(n) = 10*log10(max( ...
        D.factor*peakField^2,realmin));
    result.peakThetaDeg(n) = peakTheta;
    if peakTheta < 0.25
        % Phi is not defined for a beam along the surface normal.
        result.peakPhiDeg(n) = NaN;
    else
        result.peakPhiDeg(n) = mod(peakPhi,360);
    end
    cosine = cosd(peakTheta)*cosd(P.beamThetaDeg) + ...
        sind(peakTheta)*sind(P.beamThetaDeg)* ...
        cosd(peakPhi-P.beamPhiDeg);
    result.pointingErrorDeg(n) = acosd(max(-1,min(1,cosine)));
end
end
