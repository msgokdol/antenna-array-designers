function testPhasedArrayDesigner2D
% Focused smoke test for 2D cuts and selectable pattern factors.
addpath(fileparts(mfilename('fullpath')),'-begin');
close(findall(0,'Type','figure'),'force');
cleanup = onCleanup(@()close(findall(0,'Type','figure'),'force')); %#ok<NASGU>
phasedArrayDesigner;
main = findall(0,'Type','figure','Name','Planar Phased Array Designer');
button = findall(main,'Tag','btnPattern2D');
assert(isscalar(button),'2D pattern ribbon button missing.');
button.ButtonPushedFcn(button,[]);
win = findall(main,'Tag','pattern2DWin');
assert(isscalar(win),'Embedded 2D pattern pane did not open.');
assert(numel(findall(0,'Type','figure'))==1, ...
    '2D pattern created a separate window.');
viewRibbon = findall(main,'Tag','btnPinRef');
viewSettings = findall(main,'Tag','spCutAzimuth');
assert(strcmp(viewRibbon.Visible,'on') && strcmp(viewSettings.Visible,'on'), ...
    '2D view hid the View ribbon or its settings.');
viewTab = findall(main,'Tag','tabView');
viewTab.ButtonPushedFcn(viewTab,[]);
workspace = findall(main,'Tag','resultWorkspace');
assert(strcmp(workspace.Visible,'off'), ...
    'View tab did not restore the array and 3D workspace.');
patternTab = findall(main,'Tag','tab2D Pattern');
patternTab.ButtonPushedFcn(patternTab,[]);
assert(strcmp(win.Visible,'on') && strcmp(viewRibbon.Visible,'on'), ...
    'Returning to the 2D tab lost the plot or View ribbon.');

mode = findall(win,'Tag','ddPattern2DMode');
fixed = findall(win,'Tag','spPattern2DFixed');
line = findall(win,'Tag','pattern2DTotal');
assert(strcmp(mode.Value,'Elevation pattern') && ...
    strcmp(line.Parent.ThetaZeroLocation,'top'), ...
    'The 2D ribbon button should open the elevation cut by default.');
follow = findall(main,'Tag','cbCutFollow');
steerTheta = findall(main,'Tag','spSteerFirst');
steerPhi = findall(main,'Tag','spSteerSecond');
assert(follow.Parent == steerPhi.Parent.Parent && ...
    follow.Layout.Row == steerPhi.Parent.Layout.Row+1 && ...
    follow.Value && strcmp(fixed.Enable,'off'), ...
    'The 2D cut follow control must sit below the steering angles.');
steerTheta.Value = 30; steerTheta.ValueChangedFcn(steerTheta,[]);
steerPhi.Value = 37; steerPhi.ValueChangedFcn(steerPhi,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(fixed.Value == 37 && strcmp(fixed.Enable,'off') && ...
    contains(line.Parent.Title.String,'phi = 37.0°'), ...
    'The open 2D elevation view must follow changes to steering phi.');
follow.Value = false; follow.ValueChangedFcn(follow,[]);
fixed.Value = 15; fixed.ValueChangedFcn(fixed,[]);
steerPhi.Value = 90; steerPhi.ValueChangedFcn(steerPhi,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(strcmp(fixed.Enable,'on') && fixed.Value == 15 && ...
    contains(line.Parent.Title.String,'phi = 15.0°'), ...
    'Unlinked 2D elevation cuts must keep their selected plane.');
follow.Value = true; follow.ValueChangedFcn(follow,[]);
assert(fixed.Value == 90 && strcmp(fixed.Enable,'off'), ...
    'Re-linking the 2D cut must use the current steering phi.');
steerTheta.Value = 0; steerTheta.ValueChangedFcn(steerTheta,[]);
steerPhi.Value = 0; steerPhi.ValueChangedFcn(steerPhi,[]);
mode.Value = 'Azimuth pattern';
mode.ValueChangedFcn(mode,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(isscalar(line) && numel(line.RData)==721, ...
    'Azimuth polar cut was not sampled around the full circle.');
polar = line.Parent;
assert(strcmp(polar.ThetaZeroLocation,'right') && ...
    strcmp(polar.ThetaDir,'counterclockwise'), ...
    'Azimuth zero must remain on the +x/right axis.');
az0 = line.RData;
fixed.Value = 30;
fixed.ValueChangedFcn(fixed,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(any(abs(line.RData-az0)>1e-4), ...
    'Changing fixed elevation did not change the azimuth pattern.');
button.ButtonPushedFcn(button,[]);
assert(fixed.Value==30 && strcmp(mode.Value,'Azimuth pattern'), ...
    'Reopening 2D pattern reset the chosen cut.');

mode.Value = 'Elevation pattern';
mode.ValueChangedFcn(mode,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(isscalar(line) && min(line.ThetaData)<0 && max(line.ThetaData)>0, ...
    'Elevation polar cut did not span both sides of broadside.');
polar = line.Parent;
assert(strcmp(polar.ThetaZeroLocation,'top') && ...
    strcmp(polar.ThetaDir,'clockwise') && ...
    strcmp(polar.ThetaTickLabel{1},'0'), ...
    'Theta zero must point up and positive theta must run to the right.');

mode.Value = 'U pattern';
mode.ValueChangedFcn(mode,[]);
line = findall(win,'Tag','pattern2DTotal');
assert(isscalar(line) && abs(line.XData(1)+1)<1e-12 && ...
    abs(line.XData(end)-1)<1e-12 && all(isfinite(line.YData)), ...
    'U cut must cover finite total-field levels from -1 to +1.');
absolute = findall(win,'Tag','cbPattern2DAbs');
absolute.Value = false;
absolute.ValueChangedFcn(absolute,[]);
source = findall(win,'Tag','ddPattern2DSource');
assert(strcmp(source.Value,'Total (EF x AF)'), ...
    '2D pattern should default to the total field.');
line = findall(win,'Tag','pattern2DTotal');
mainCut = findall(main,'Tag','cutTotal');
assert(isscalar(mainCut),'Main total cut missing.');
[~,iU] = min(abs(line.XData));
[~,iCut] = min(abs(mainCut.XData));
assert(abs(line.YData(iU)-mainCut.YData(iCut))<1e-8, ...
    'U=0 and the main broadside cut disagree.');
for pair = {'Array factor only','cutAF'; 'Element factor only','cutEF'}'
    source.Value = pair{1};
    source.ValueChangedFcn(source,[]);
    line = findall(win,'Tag','pattern2DTotal');
    mainCut = findall(main,'Tag',pair{2});
    [~,iU] = min(abs(line.XData));
    [~,iCut] = min(abs(mainCut.XData));
    assert(abs(line.YData(iU)-mainCut.YData(iCut))<1e-8, ...
        '%s disagrees with the main cut at broadside.',pair{1});
end
polChoice = findall(win,'Tag','ddPattern2DPol');
assert(strcmp(polChoice.Enable,'on'), ...
    'Element-factor polarization control should remain available.');
source.Value = 'Array factor only'; source.ValueChangedFcn(source,[]);
assert(strcmp(polChoice.Enable,'off'), ...
    'Array factor has no polarization to select.');
source.Value = 'Total (EF x AF)'; source.ValueChangedFcn(source,[]);

mainAbsolute = findall(main,'Tag','cbAbs');
mainAbsolute.Value = true; mainAbsolute.ValueChangedFcn(mainAbsolute,[]);
absolute.Value = true; absolute.ValueChangedFcn(absolute,[]);
for pair = {'Array factor only','cutAF'; 'Element factor only','cutEF'; ...
        'Total (EF x AF)','cutTotal'}'
    source.Value = pair{1};
    source.ValueChangedFcn(source,[]);
    line = findall(win,'Tag','pattern2DTotal');
    mainCut = findall(main,'Tag',pair{2});
    [~,iU] = min(abs(line.XData));
    [~,iCut] = min(abs(mainCut.XData));
    assert(abs(line.YData(iU)-mainCut.YData(iCut))<1e-8, ...
        'Absolute %s disagrees with the main cut.',pair{1});
end
mainAbsolute.Value = false; mainAbsolute.ValueChangedFcn(mainAbsolute,[]);
absolute.Value = false; absolute.ValueChangedFcn(absolute,[]);

convention = findall(main,'Tag','ddAngleConvention');
convention.Value = 'Azimuth / elevation';
convention.ValueChangedFcn(convention,[]);
mode.Value = 'Azimuth pattern';
mode.ValueChangedFcn(mode,[]);
label = findall(win,'Tag','lblPattern2DFixed');
assert(strcmp(label.Text,'Elevation (°)') && fixed.Value==0, ...
    '2D fixed angle did not convert to elevation.');
menuU = findall(main,'Tag','menuPattern2DU');
menuU.MenuSelectedFcn(menuU,[]);
assert(strcmp(mode.Value,'U pattern'),'View menu did not select the U cut.');

theme(main,'dark'); drawnow;
darkLine = findall(win,'Tag','pattern2DTotal');
darkColor = darkLine.Color;
assert(strcmp(main.Theme.BaseColorStyle,'dark'), ...
    'Main app did not use the dark theme.');
theme(main,'light'); drawnow;
lightLine = findall(win,'Tag','pattern2DTotal');
assert(strcmp(main.Theme.BaseColorStyle,'light') && ...
    ~isequal(lightLine.Color,darkColor), ...
    '2D pattern trace did not recolor with the light theme.');

% The crossed-dipole guide depicts the two circular-polarization lobes,
% whereas its combined power has a shallow (3 dB) equatorial waist. A
% 40 dB radius made that combined 2D cut look almost circular.
closeTab = findall(main,'Tag','close2D Pattern');
assert(isprop(closeTab,'ImageClickedFcn') && ...
    isequal(closeTab.Parent,patternTab.Parent), ...
    'The close glyph should sit inside the result tab.');
closeTab.ImageClickedFcn(closeTab,[]);
assert(isempty(findall(main,'Tag','pattern2DWin')) && ...
    strcmp(closeTab.Visible,'off') && strcmp(workspace.Visible,'off'), ...
    'Closing 2D Pattern did not restore and hide its result tab.');
auto = findall(main,'Tag','cbAuto'); auto.Value = false;
for tag = {'spM','spN'}
    spArray = findall(main,'Tag',tag{1});
    spArray.Value = 1; spArray.ValueChangedFcn(spArray,[]);
end
element = findall(main,'Tag','ddEF');
element.Value = 'Crossed dipole (RHCP)';
element.ValueChangedFcn(element,[]);
full = findall(main,'Tag','cbFullSphere');
full.Value = true; full.ValueChangedFcn(full,[]);
compute = findall(main,'Tag','btnCompute');
compute.ButtonPushedFcn(compute,[]);
button.ButtonPushedFcn(button,[]);
win = findall(main,'Tag','pattern2DWin');
mode = findall(win,'Tag','ddPattern2DMode');
mode.Value = 'Elevation pattern'; mode.ValueChangedFcn(mode,[]);
polChoice = findall(win,'Tag','ddPattern2DPol');
scale = findall(win,'Tag','ddPattern2DScale');
assert(strcmp(polChoice.Value,'Total (any pol)') && ...
    strcmp(scale.Value,'dB (dynamic range)'), ...
    'The 2D pattern should open with total power on a dB radius.');
scale.Value = 'Linear power'; scale.ValueChangedFcn(scale,[]);
combined = findall(win,'Tag','pattern2DTotal');
[~,iFront] = min(abs(combined.ThetaData));
[~,iSide] = min(abs(combined.ThetaData-pi/2));
[~,iBack] = min(abs(combined.ThetaData-pi));
assert(abs(combined.RData(iFront)-1)<0.02 && ...
    abs(combined.RData(iSide)-0.5)<0.02 && ...
    abs(combined.RData(iBack)-1)<0.02, ...
    'Combined crossed-dipole power should have a 3 dB equatorial waist.');
polChoice.Value = 'RHCP + LHCP'; polChoice.ValueChangedFcn(polChoice,[]);
right = findall(win,'Tag','pattern2DRHCP');
left = findall(win,'Tag','pattern2DLHCP');
assert(isscalar(right) && isscalar(left),'Both CP lobes were not plotted.');
assert(abs(right.RData(iFront)-1)<0.02 && right.RData(iBack)<0.01 && ...
    left.RData(iFront)<0.01 && abs(left.RData(iBack)-1)<0.02 && ...
    abs(right.RData(iSide)-0.25)<0.02 && ...
    abs(left.RData(iSide)-0.25)<0.02, ...
    'Crossed-dipole RHCP/LHCP lobes do not match the guide.');
scale.Value = 'dB (dynamic range)'; scale.ValueChangedFcn(scale,[]);
assert(strcmp(findall(win,'Tag','cbPattern2DAbs').Visible,'on'), ...
    'Absolute dBi control should be available in the dB-radius view.');
fprintf('2D azimuth, elevation, U, and crossed-dipole CP views passed.\n');
end
