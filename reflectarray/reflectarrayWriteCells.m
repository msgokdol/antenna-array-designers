function count = reflectarrayWriteCells(A,phaseDeg,S,destination,kind)
%REFLECTARRAYWRITECELLS Export active ideal cells, not a physical reflector.
%   CSV keeps the reflection and incident phase separate. CST-shaped TSV
%   is an equivalent driven aperture: Phase = reflection + incident phase.

if ~ismember(kind,{'csv','tsv'})
    error('Reflectarray:Export','Unknown export format.');
end
if ~isnumeric(phaseDeg) || ~isreal(phaseDeg) || ...
        ~isequal(size(phaseDeg),size(A.active)) || ...
        any(~isfinite(phaseDeg(A.active)))
    error('Reflectarray:Export','Cell phases must be finite and match the surface.');
end
[cols,rows] = find(A.active.');
idx = sub2ind(size(A.active),rows,cols);
count = numel(idx);
if count == 0
    error('Reflectarray:Export','There are no active cells to export.');
end
reflection = mod(phaseDeg(idx),360);
incident = mod(-rad2deg(A.k*A.incidentPathMM(idx)),360);
equivalent = mod(reflection+incident,360);
phaseReferenceGHz = A.phaseReferenceGHz;
if isfield(S,'manualEnabled') && S.manualEnabled
    % Manually assigned phases have no common synthesis frequency.
    phaseReferenceGHz = NaN;
end
if strcmp(kind,'csv')
    T = table(rows,cols,A.X(idx),A.Y(idx),zeros(count,1), ...
        A.amplitude(idx),reflection,incident,equivalent, ...
        repmat(phaseReferenceGHz,count,1), ...
        repmat(S.opFreqGHz,count,1), ...
        'VariableNames',{'row','column','x_mm','y_mm','z_mm', ...
        'relative_amplitude','reflection_phase_deg', ...
        'incident_phase_deg','equivalent_excitation_phase_deg', ...
        'phase_reference_GHz', ...
        'operating_frequency_GHz'});
    writetable(T,destination);
else
    fid = fopen(destination,'w');
    if fid < 0, error('Reflectarray:Export','Cannot open export file.'); end
    cleanupFile = onCleanup(@()fclose(fid));
    fprintf(fid,'# Equivalent aperture from Reflectarray Designer; NOT a feed or physical reflector model\n');
    fprintf(fid,'# unit: meters\n');
    fprintf(fid,'# design frequency: %.10g GHz\n',S.designFreqGHz);
    if isnan(phaseReferenceGHz)
        fprintf(fid,'# phase reference frequency: manual / unspecified\n');
    else
        fprintf(fid,'# phase reference frequency: %.10g GHz\n',phaseReferenceGHz);
    end
    fprintf(fid,'# operating frequency: %.10g GHz\n',S.opFreqGHz);
    fprintf(fid,'# illumination: %s\n',A.illuminationMode);
    fprintf(fid,'# Phase = reflection phase + incident field phase at operating frequency\n');
    fprintf(fid,'# Element\tX\tY\tZ\tMagnitude\tPhase\tPhi\tTheta\tGamma\n');
    for t = 1:count
        fprintf(fid,'%d,%d\t%.12g\t%.12g\t0\t%.12g\t%.12g\t0\t0\t0\n', ...
            cols(t),rows(t),A.X(idx(t))*1e-3, ...
            A.Y(idx(t))*1e-3,A.amplitude(idx(t)), ...
            equivalent(t));
    end
    clear cleanupFile;
end
end
