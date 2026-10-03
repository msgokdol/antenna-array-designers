function testPhasedArrayDesignerElementReference
% Compare diagonal element cuts with R2025b Phased Array Toolbox responses.
% The values below are peak-relative dB at 65-degree az/el beamwidths;
% they are fixed references, so this test runs without that toolbox.
here = fileparts(mfilename('fullpath'));
addpath(here,'-begin');
cleanup = onCleanup(@()close(findall(0,'Type','figure'),'force')); %#ok<NASGU>
phasedArrayDesigner;
f = findall(0,'Type','figure','Name','Planar Phased Array Designer');
assert(isscalar(f),'Designer window did not open.');
auto = one(f,'cbAuto'); auto.Value = false;
for tag = {'spM','spN'}
    sp = one(f,tag{1}); sp.Value = 1; sp.ValueChangedFcn(sp,[]);
end
follow = one(f,'cbCutFollow');
follow.Value = false; follow.ValueChangedFcn(follow,[]);
cutPhi = one(f,'spCutAzimuth');
cutPhi.Value = 45; cutPhi.ValueChangedFcn(cutPhi,[]);
dd = one(f,'ddEF');
compute = one(f,'btnCompute');
models = {'Gaussian','Sinc','3GPP TR 38.901 shape'};
referenceDb = [-11.4095 -22.0124; ...
               -10.9047 -17.9753; ...
               -11.3705 -21.9371];
referenceDirDb = [9.8569 9.6478 9.8257];
for k = 1:numel(models)
    dd.Value = models{k}; dd.ValueChangedFcn(dd,[]);
    compute.ButtonPushedFcn(compute,[]);
    line = one(f,'cutTotal');
    [~,i0] = min(abs(line.XData));
    for j = 1:2
        angle = [60 80];
        [~,ix] = min(abs(line.XData-angle(j)));
        actualDb = line.YData(ix)-line.YData(i0);
        assert(abs(actualDb-referenceDb(k,j))<0.12, ...
            '%s diagonal cut at %d deg: %.3f dB, expected %.3f dB.', ...
            models{k},angle(j),actualDb,referenceDb(k,j));
    end
    actualDirDb = str2double(one(f,'cardDirectivity').Text);
    assert(abs(actualDirDb-referenceDirDb(k))<0.12, ...
        '%s directivity: %.3f dBi, expected %.3f dBi.', ...
        models{k},actualDirDb,referenceDirDb(k));
end

dd.Value = 'Cardioid'; dd.ValueChangedFcn(dd,[]);
compute.ButtonPushedFcn(compute,[]);
line = one(f,'cutTotal');
[~,i0] = min(abs(line.XData));
[~,i90] = min(abs(line.XData-90));
horizonDb = line.YData(i90)-line.YData(i0);
assert(abs(horizonDb+3.0103)<0.12, ...
    'Cardioid horizon should be -3 dB, not the narrower textbook shape.');
assert(abs(str2double(one(f,'cardDirectivity').Text)-3.0103)<0.12, ...
    'Cardioid should have the gallery reference model''s 3.01 dBi directivity.');
fprintf('Off-axis element patterns and gallery Cardioid match references.\n');
end

function h = one(parent,tag)
h = findall(parent,'Tag',tag);
assert(isscalar(h),'Missing or duplicated %s.',tag);
end
