function reflectarrayDesigner
%REFLECTARRAYDESIGNER Standalone ideal phase-only reflectarray prototype.
%   Run from this folder, or add this folder to the MATLAB path. The
%   phased-array application is deliberately not modified by this tool.

BG = [0.12 0.12 0.12]; PANEL = [0.15 0.15 0.15];
PLOT = [0.065 0.065 0.07]; WHITE = [0.91 0.93 0.96];
MUTED = [0.69 0.72 0.77]; BLUE = [0.20 0.58 0.91];
NAVY = [0.015 0.24 0.42];
iconRoot = fullfile(fileparts(mfilename('fullpath')),'icons');
defaults = struct('M',20,'N',20,'dxLambda',0.5,'dyLambda',0.5, ...
    'designFreqGHz',10,'opFreqGHz',10,'freqUnit','GHz', ...
    'retunePhase',true, ...
    'shape','Circle','feedMM',[-40 0 120],'feedQ',4, ...
    'illuminationMode','Horn','incidentThetaDeg',0,'incidentPhiDeg',0, ...
    'beamThetaDeg',0,'beamPhiDeg',0,'angleConvention','Theta / phi', ...
    'taper','None','manualEnabled',false,'manualPhaseDeg',[], ...
    'gridAngleDeg',90,'rowStaggerLambda',0,'staggerAngleDeg',0, ...
    'activeMask',[],'xOffsetLambda',zeros(20), ...
    'yOffsetLambda',zeros(20), ...
    'selected',zeros(0,2),'selectedCells',zeros(0,2));
S = defaults;
A = reflectarrayAperture(S);
phaseDeg = A.referencePhaseDeg;
patternData = struct();
directivity = struct('factor',NaN,'pairPower',0);
idealDirectivity = [];
selectedGeometryHighlight = gobjects(1);
editorControls = struct();
shapePreviewHeight = 122;
shapeRecent = {'Diamond','Custom','Circle'};
shapePreviewButtons = gobjects(1,3);
shapeGallery = gobjects(1);
projectFile = '';
dirty = false;
lastDeletedCells = zeros(0,2);
lastDeletedSelection = zeros(0,2);
geometryViewInitialized = false;
patternViewInitialized = false;
sweepStartGHz = 0.9*S.opFreqGHz;
sweepStopGHz = 1.1*S.opFreqGHz;
sweepResult = [];

fig = uifigure('Name','Reflectarray Designer', ...
    'Position',[90 65 1480 920],'WindowState','maximized', ...
    'Color',BG,'Tag','reflectarrayDesigner');
fig.CloseRequestFcn = @(~,~)closeDesigner();
root = uigridlayout(fig,[3 1]);
root.RowHeight = {76,'1x',26};
root.Padding = [8 5 8 5]; root.RowSpacing = 6;
root.BackgroundColor = BG;

bar = uigridlayout(root,[1 11]);
bar.ColumnWidth = {88,88,88,88,128,'1x',68,68,68,68,68};
bar.RowHeight = {'1x'}; bar.Padding = [3 3 3 3];
bar.ColumnSpacing = 4; bar.BackgroundColor = NAVY;
viewNames = {'Geometry','Phase Map','3D Pattern','2D Pattern', ...
    'Frequency Sweep'};
viewButtons = gobjects(1,numel(viewNames));
for j = 1:numel(viewNames)
    name = viewNames{j};
    viewButtons(j) = uibutton(bar,'Text',name,'Tag',['tab' strrep(name,' ','')], ...
        'FontSize',13,'FontWeight','bold','FontColor',WHITE, ...
        'BackgroundColor',NAVY,'ButtonPushedFcn',@(~,~)showView(name));
    viewButtons(j).Layout.Column = j;
end
hdr = uilabel(bar,'Text','REFLECTARRAY','FontSize',13, ...
    'FontWeight','bold','FontColor',WHITE,'HorizontalAlignment','right');
hdr.Layout.Column = 6;
bNew = uibutton(bar,'Text','New','FontColor',WHITE, ...
    'Icon',fullfile(iconRoot,'toolbar_new.png'),'IconAlignment','top', ...
    'BackgroundColor',NAVY,'Tag','reflectarrayNew', ...
    'ButtonPushedFcn',@(~,~)newProject());
bNew.Layout.Column = 7;
bSave = uibutton(bar,'Text','Save','FontColor',WHITE, ...
    'Icon',fullfile(iconRoot,'toolbar_save.png'),'IconAlignment','top', ...
    'BackgroundColor',NAVY,'Tag','reflectarraySave', ...
    'ButtonPushedFcn',@(~,~)saveProject());
bSave.Layout.Column = 8;
bOpen = uibutton(bar,'Text','Open','FontColor',WHITE, ...
    'Icon',fullfile(iconRoot,'toolbar_open.png'),'IconAlignment','top', ...
    'BackgroundColor',NAVY,'Tag','reflectarrayOpen', ...
    'ButtonPushedFcn',@(~,~)openProject());
bOpen.Layout.Column = 9;
bCsv = uibutton(bar,'Text','CSV','FontColor',WHITE, ...
    'Icon',fullfile(iconRoot,'toolbar_csv.png'),'IconAlignment','top', ...
    'BackgroundColor',NAVY,'Tag','reflectarrayCSV', ...
    'ButtonPushedFcn',@(~,~)exportCells('csv'));
bCsv.Layout.Column = 10;
bTsv = uibutton(bar,'Text','CST .tsv','FontColor',WHITE, ...
    'Icon',fullfile(iconRoot,'toolbar_cst.png'),'IconAlignment','top', ...
    'BackgroundColor',NAVY,'Tag','reflectarrayCST', ...
    'Tooltip','Export cell positions, amplitude and phase as an equivalent aperture. Not a CST feed-port model.', ...
    'ButtonPushedFcn',@(~,~)exportCells('tsv'));
bTsv.Layout.Column = 11;

body = uigridlayout(root,[1 2]);
body.ColumnWidth = {350,'1x'}; body.RowHeight = {'1x'};
body.Padding = [0 0 0 0]; body.ColumnSpacing = 8;
body.BackgroundColor = BG;

leftShell = uigridlayout(body,[2 1]);
leftShell.RowHeight = {36,'1x'};
leftShell.Padding = [0 0 0 0]; leftShell.RowSpacing = 5;
leftShell.BackgroundColor = BG;
leftTabs = uigridlayout(leftShell,[1 3]);
leftTabs.Layout.Row = 1;
leftTabs.ColumnWidth = {'1x','1x','1x'};
leftTabs.RowHeight = {'1x'};
leftTabs.Padding = [0 0 0 0]; leftTabs.ColumnSpacing = 4;
leftTabs.BackgroundColor = NAVY;
bDesignSection = uibutton(leftTabs,'Text','DESIGN', ...
    'FontColor',WHITE,'FontWeight','bold','FontSize',10, ...
    'BackgroundColor',BG, ...
    'Tag','reflectarrayDesignSection', ...
    'ButtonPushedFcn',@(~,~)showControlSection('Design'));
bDesignSection.Layout.Column = 1;
bGeometrySection = uibutton(leftTabs,'Text','ARRAY GEOMETRY', ...
    'FontColor',WHITE,'FontWeight','bold','FontSize',10, ...
    'BackgroundColor',NAVY, ...
    'Tag','reflectarrayArrayGeometrySection', ...
    'ButtonPushedFcn',@(~,~)showControlSection('Array Geometry'));
bGeometrySection.Layout.Column = 2;
bBeamSection = uibutton(leftTabs,'Text','BEAM', ...
    'FontColor',WHITE,'FontWeight','bold','FontSize',10, ...
    'BackgroundColor',NAVY, ...
    'Tag','reflectarrayBeamSection', ...
    'ButtonPushedFcn',@(~,~)showControlSection('Beam'));
bBeamSection.Layout.Column = 3;
leftWorkspace = uigridlayout(leftShell,[1 1]);
leftWorkspace.Layout.Row = 2;
leftWorkspace.Padding = [0 0 0 0]; leftWorkspace.BackgroundColor = BG;
designPanel = uipanel(leftWorkspace,'BorderType','none', ...
    'BackgroundColor',BG,'Tag','reflectarrayDesignPanel');
designPanel.Layout.Row = 1; designPanel.Layout.Column = 1;
geometryPanel = uipanel(leftWorkspace,'BorderType','none', ...
    'BackgroundColor',BG,'Visible','off', ...
    'Tag','reflectarrayArrayGeometryPanel');
geometryPanel.Layout.Row = 1; geometryPanel.Layout.Column = 1;
beamPanel = uipanel(leftWorkspace,'BorderType','none', ...
    'BackgroundColor',BG,'Visible','off', ...
    'Tag','reflectarrayBeamPanel');
beamPanel.Layout.Row = 1; beamPanel.Layout.Column = 1;
left = uigridlayout(designPanel,[4 1]);
left.RowHeight = {132,82,165,'1x'};
left.Padding = [0 0 0 0]; left.RowSpacing = 6;
left.BackgroundColor = BG;
beamLeft = uigridlayout(beamPanel,[4 1]);
beamLeft.RowHeight = {132,98,90,'1x'};
beamLeft.Padding = [0 0 0 0]; beamLeft.RowSpacing = 6;
beamLeft.BackgroundColor = BG;

% Design frequency fixes physical cell pitch; operating frequency changes
% propagation without moving the already designed surface.
pCfg = section(left,'CONFIGURATION'); pCfg.Layout.Row = 1;
gCfg = uigridlayout(pCfg,[4 3]);
gCfg.ColumnWidth = {130,'1x',65};
gCfg.RowHeight = {25,25,18,22};
gCfg.Padding = [7 2 7 2]; gCfg.RowSpacing = 2;
gCfg.ColumnSpacing = 4; gCfg.BackgroundColor = PANEL;
cfgLabel = uilabel(gCfg,'Text','Design frequency', ...
    'FontColor',WHITE,'FontSize',11);
cfgLabel.Layout.Row = 1; cfgLabel.Layout.Column = 1;
spDesignFreq = number(gCfg,S.designFreqGHz,[0.01 300],0.5,1,2, ...
    @(src,~)setDesignFrequency(src.Value));
spDesignFreq.Tag = 'reflectarrayDesignFrequency';
ddFreqUnit = dropdown(gCfg,{'GHz','MHz'},S.freqUnit,1,3, ...
    @(~,~)setFrequencyUnit());
ddFreqUnit.Tag = 'reflectarrayFrequencyUnit';
cfgLabel = uilabel(gCfg,'Text','Operating frequency', ...
    'FontColor',WHITE,'FontSize',11);
cfgLabel.Layout.Row = 2; cfgLabel.Layout.Column = 1;
spOperatingFreq = number(gCfg,S.opFreqGHz,[0.01 300],0.5,2,2, ...
    @(src,~)setOperatingFrequency(src.Value));
spOperatingFreq.Layout.Column = [2 3];
spOperatingFreq.Tag = 'reflectarrayOperatingFrequency';
lblDesignLambda = uilabel(gCfg,'Text','', ...
    'FontColor',MUTED,'FontSize',10, ...
    'Tag','reflectarrayDesignWavelength');
lblDesignLambda.Layout.Row = 3; lblDesignLambda.Layout.Column = [1 3];
cbFixedPhase = uicheckbox(gCfg, ...
    'Text','Keep phases fixed at design frequency', ...
    'Value',~S.retunePhase,'FontColor',WHITE,'FontSize',11, ...
    'Tag','reflectarrayFixedPhase', ...
    'ValueChangedFcn',@(~,~)setPhaseMode());
cbFixedPhase.Layout.Row = 4; cbFixedPhase.Layout.Column = [1 3];
cbFixedPhase.Tooltip = ['Checked: use phases synthesized at the design ' ...
    'frequency while evaluating the pattern at the operating frequency.'];

pGeom = section(left,'SURFACE'); pGeom.Layout.Row = 2;
g = sectionGrid(pGeom,2);
label(g,'Rows M',1,1); spM = number(g,S.M,[2 60],1,1,2,@(~,~)setGridSize());
label(g,'Cols N',1,3); spN = number(g,S.N,[2 60],1,1,4,@(~,~)setGridSize());
spM.Tag = 'reflectarrayRows'; spN.Tag = 'reflectarrayColumns';
label(g,'dx (λ)',2,1); spDx = number(g,S.dxLambda,[0.01 10],0.05,2,2, ...
    @(src,~)setSetting('dxLambda',src.Value));
spDx.Tag = 'reflectarrayDxLambda';
spDx.Tooltip = 'Physical dx equals this value times the design wavelength.';
label(g,'dy (λ)',2,3); spDy = number(g,S.dyLambda,[0.01 10],0.05,2,4, ...
    @(src,~)setSetting('dyLambda',src.Value));
spDy.Tag = 'reflectarrayDyLambda';
spDy.Tooltip = 'Physical dy equals this value times the design wavelength.';

pFeed = section(beamLeft,'ILLUMINATION'); pFeed.Layout.Row = 1;
g = sectionGrid(pFeed,3);
label(g,'Source',1,1);
ddIllumination = dropdown(g,{'Horn','Plane wave'},S.illuminationMode, ...
    1,2,@(~,~)setIllumination());
ddIllumination.Layout.Column = [2 4];
ddIllumination.Tag = 'reflectarrayIlluminationMode';
lblFeedX = label(g,'X (mm)',2,1); spFx = number(g,S.feedMM(1),[-10000 10000],5,2,2, ...
    @(src,~)setFeed(1,src.Value));
spFx.Tag = 'reflectarrayFeedX';
lblFeedY = label(g,'Y (mm)',2,3); spFy = number(g,S.feedMM(2),[-10000 10000],5,2,4, ...
    @(src,~)setFeed(2,src.Value));
lblFeedZ = label(g,'Z (mm)',3,1); spFz = number(g,S.feedMM(3),[1 10000],5,3,2, ...
    @(src,~)setFeed(3,src.Value));
lblFeedQ = label(g,'Feed q',3,3); spQ = number(g,S.feedQ,[0 30],0.5,3,4, ...
    @(src,~)setSetting('feedQ',src.Value));
spQ.Tooltip = ['The horn aims at the surface centre. Field illumination ' ...
    'varies as cos^q(angle) divided by distance. Phase correction does not flatten it.'];
lblIncidentTheta = label(g,'Inc. θ (°)',2,1);
spIncidentTheta = number(g,S.incidentThetaDeg,[0 85],1,2,2, ...
    @(src,~)setIncident('incidentThetaDeg',src.Value));
spIncidentTheta.Tag = 'reflectarrayIncidentTheta';
spIncidentTheta.Tooltip = 'Angle toward the plane-wave source, measured from +z. 0° is normal incidence.';
lblIncidentPhi = label(g,'Inc. φ (°)',2,3);
spIncidentPhi = number(g,S.incidentPhiDeg,[0 360],1,2,4, ...
    @(src,~)setIncident('incidentPhiDeg',src.Value));
spIncidentPhi.Tag = 'reflectarrayIncidentPhi';
spIncidentPhi.Tooltip = 'Azimuth toward the plane-wave source, measured from +x.';
planeNote = uilabel(g,'Text','Uniform incident amplitude · propagation toward the surface', ...
    'FontColor',MUTED,'FontSize',10,'WordWrap','on');
planeNote.Layout.Row = 3; planeNote.Layout.Column = [1 4];

pBeam = section(beamLeft,'INTENDED BEAM'); pBeam.Layout.Row = 2;
g = sectionGrid(pBeam,2);
label(g,'Angles',1,1);
ddConvention = dropdown(g,{'Theta / phi','Elevation / azimuth'}, ...
    S.angleConvention,1,2,@(~,~)setConvention());
ddConvention.Tag = 'reflectarrayAngleConvention';
ddConvention.Layout.Column = [2 4];
lblTheta = label(g,'Theta (°)',2,1);
spTheta = number(g,S.beamThetaDeg,[0 90],1,2,2,@(~,~)setBeam());
spTheta.Tag = 'reflectarrayBeamTheta';
lblPhi = label(g,'Phi (°)',2,3);
spPhi = number(g,S.beamPhiDeg,[0 360],1,2,4,@(~,~)setBeam());
spPhi.Tag = 'reflectarrayBeamPhi';
spPhi.Tooltip = 'Azimuth angle measured from +x in the xy plane.';

pTaper = section(beamLeft,'AMPLITUDE'); pTaper.Layout.Row = 3;
gt = uigridlayout(pTaper,[2 2]); gt.ColumnWidth = {115,'1x'};
gt.RowHeight = {27,'1x'}; gt.Padding = [7 5 7 5];
gt.BackgroundColor = PANEL;
label(gt,'Extra taper',1,1);
ddTaper = uidropdown(gt,'Items',{'None','Hann taper'}, ...
    'ItemsData',{'None','Hann (ideal amplitude)'},'Value',S.taper, ...
    'ValueChangedFcn',@(src,~)setSetting('taper',src.Value),'FontSize',11);
ddTaper.Tooltip = 'Hann applies an additional mathematical amplitude taper. Phase-only cells cannot impose it independently.';
ddTaper.Layout.Row = 1; ddTaper.Layout.Column = 2;
taperNote = uilabel(gt,'Text','Phase correction aligns cells but does not flatten horn illumination. Hann adds extra amplitude taper.', ...
    'FontColor',MUTED,'FontSize',10,'WordWrap','on');
taperNote.Layout.Row = 2; taperNote.Layout.Column = [1 2];
syncIlluminationControls();

pCell = section(left,'SELECTED CELL & MODEL'); pCell.Layout.Row = 3;
gc = uigridlayout(pCell,[5 2]); gc.ColumnWidth = {125,'1x'};
gc.RowHeight = {24,27,32,23,'1x'}; gc.Padding = [7 3 7 3];
gc.BackgroundColor = PANEL;
selectedText = uilabel(gc,'Text','No cell selected','FontColor',WHITE, ...
    'FontWeight','bold','Tag','reflectarraySelectedCellLabel');
selectedText.Layout.Row = 1;
selectedText.Layout.Column = [1 2];
label(gc,'Reflection (°)',2,1);
spCell = number(gc,0,[-36000 36000],5,2,2, ...
    @(src,~)editSelected(src.Value));
spCell.Tag = 'reflectarrayCellPhase';
spCell.Enable = 'off';
btnIdeal = uibutton(gc,'Text','Restore synthesized phases', ...
    'BackgroundColor',NAVY,'FontColor',WHITE, ...
    'ButtonPushedFcn',@(~,~)restoreIdeal());
btnIdeal.Tag = 'reflectarrayRestoreSynthesizedPhases';
btnIdeal.Layout.Row = 3; btnIdeal.Layout.Column = [1 2];
modeText = uilabel(gc,'Text','Synthesized phases','FontColor',WHITE, ...
    'FontWeight','bold');
modeText.Layout.Row = 4; modeText.Layout.Column = [1 2];
infoText = uilabel(gc,'Text', ...
    'Normalized scalar model; no absolute gain or physical cell losses.', ...
    'FontColor',MUTED,'WordWrap','on','VerticalAlignment','top');
infoText.Layout.Row = 5; infoText.Layout.Column = [1 2];

buildInlineGeometryControls();

% Views share one workspace; result plots never create extra windows.
stack = uigridlayout(body,[1 1]); stack.Padding = [0 0 0 0];
stack.BackgroundColor = BG;
views = gobjects(1,numel(viewNames));
for j = 1:numel(viewNames)
    views(j) = uipanel(stack,'BorderType','line','BackgroundColor',PANEL, ...
        'Visible','off','Tag',['view' strrep(viewNames{j},' ','')]);
    views(j).Layout.Row = 1; views(j).Layout.Column = 1;
end
geometryGrid = uigridlayout(views(1),[2 2]);
geometryGrid.RowHeight = {31,'1x'};
geometryGrid.ColumnWidth = {170,'1x'};
geometryGrid.Padding = [10 10 10 10];
geometryGrid.BackgroundColor = PANEL;
showGeometryBeam = uicheckbox(geometryGrid,'Text','Show 3D beam', ...
    'Value',false,'FontColor',WHITE,'FontSize',12, ...
    'Tooltip','Beam radius uses a display scale, not physical millimetres.', ...
    'Tag','reflectarrayShowGeometryBeam', ...
    'ValueChangedFcn',@(~,~)drawGeometry());
showGeometryBeam.Layout.Row = 1;
showGeometryBeam.Layout.Column = 1;
geometryHint = uilabel(geometryGrid, ...
    'Text','Click a blue cell to select it.', ...
    'FontColor',MUTED,'FontSize',11,'HorizontalAlignment','right');
geometryHint.Layout.Row = 1;
geometryHint.Layout.Column = 2;
axGeom = uiaxes(geometryGrid,'Tag','reflectarrayGeometryAxes');
axGeom.Layout.Row = 2;
axGeom.Layout.Column = [1 2];
styleAxes(axGeom);

mapView = uigridlayout(views(2),[2 1]);
mapView.RowHeight = {35,'1x'};
mapView.Padding = [7 7 7 7]; mapView.RowSpacing = 5;
mapView.BackgroundColor = PANEL;
mapTools = uigridlayout(mapView,[1 3]);
mapTools.ColumnWidth = {45,155,'1x'};
mapTools.Padding = [0 0 0 0]; mapTools.ColumnSpacing = 5;
mapTools.BackgroundColor = PANEL;
label(mapTools,'Show',1,1);
ddMapMode = dropdown(mapTools,{'Phase map','Illumination map'}, ...
    'Phase map',1,2,@(~,~)showMapMode());
ddMapMode.Tag = 'reflectarrayMapMode';
mapNote = uilabel(mapTools,'Text','', ...
    'FontColor',MUTED,'FontSize',11);
mapNote.Layout.Row = 1; mapNote.Layout.Column = 3;
phaseGrid = uigridlayout(mapView,[1 2]);
phaseGrid.Layout.Row = 2;
phaseGrid.ColumnWidth = {'1x',270};
phaseGrid.Padding = [0 0 0 0]; phaseGrid.ColumnSpacing = 7;
phaseGrid.BackgroundColor = PANEL;
phasePlots = uigridlayout(phaseGrid,[2 2]);
phasePlots.Layout.Column = 1;
phasePlots.RowHeight = {'1x','1x'};
phasePlots.ColumnWidth = {'1x','1x'};
phasePlots.Padding = [0 0 0 0]; phasePlots.RowSpacing = 2;
phasePlots.ColumnSpacing = 2; phasePlots.BackgroundColor = PANEL;
axSpatial = uiaxes(phasePlots,'Tag','reflectarraySpatialDelayAxes');
axSpatial.Layout.Row = 1; axSpatial.Layout.Column = 1;
styleAxes(axSpatial);
axProgressive = uiaxes(phasePlots,'Tag','reflectarrayProgressivePhaseAxes');
axProgressive.Layout.Row = 1; axProgressive.Layout.Column = 2;
styleAxes(axProgressive);
axPhase = uiaxes(phasePlots,'Tag','reflectarrayPhaseAxes');
axPhase.Layout.Row = 2; axPhase.Layout.Column = 1;
styleAxes(axPhase);
axContinuous = uiaxes(phasePlots,'Tag','reflectarrayContinuousPhaseAxes');
axContinuous.Layout.Row = 2; axContinuous.Layout.Column = 2;
styleAxes(axContinuous);
phaseTableBox = uigridlayout(phaseGrid,[1 1]);
phaseTableBox.Layout.Column = 2;
phaseTableBox.RowHeight = {'1x'};
phaseTableBox.Padding = [0 0 0 0];
phaseTableBox.BackgroundColor = PANEL;
phaseTable = uitable(phaseTableBox,'Tag','reflectarrayPhaseTable', ...
    'CellEditCallback',@onTableEdit, ...
    'CellSelectionCallback',@onTableSelect);
phaseTable.Layout.Row = 1;
phaseTable.BackgroundColor = PLOT;
phaseTable.ForegroundColor = WHITE;
illuminationPlots = uigridlayout(mapView,[1 2]);
illuminationPlots.Layout.Row = 2;
illuminationPlots.ColumnWidth = {'1x','1x'};
illuminationPlots.Padding = [0 0 0 0];
illuminationPlots.ColumnSpacing = 7;
illuminationPlots.BackgroundColor = PANEL;
illuminationPlots.Visible = 'off';
axIllumination = uiaxes(illuminationPlots, ...
    'Tag','reflectarrayIlluminationAxes');
axIllumination.Layout.Column = 1;
styleAxes(axIllumination);
axApertureAmplitude = uiaxes(illuminationPlots, ...
    'Tag','reflectarrayApertureAmplitudeAxes');
axApertureAmplitude.Layout.Column = 2;
styleAxes(axApertureAmplitude);

patternGrid = uigridlayout(views(3),[2 1]);
patternGrid.RowHeight = {32,'1x'};
patternGrid.Padding = [10 10 10 10];
patternGrid.BackgroundColor = PANEL;
patternTools = uigridlayout(patternGrid,[1 2]);
patternTools.ColumnWidth = {175,'1x'};
patternTools.Padding = [0 0 0 0]; patternTools.BackgroundColor = PANEL;
cbShow3DCut = uicheckbox(patternTools,'Text','Show selected 2D cut', ...
    'Value',false,'FontColor',WHITE,'FontSize',11, ...
    'Tag','reflectarrayShow3DCut', ...
    'ValueChangedFcn',@(~,~)drawPattern3D());
cbShow3DCut.Layout.Column = 1;
cutOverlayInfo = uilabel(patternTools,'Text','', ...
    'FontColor',MUTED,'FontSize',11);
cutOverlayInfo.Layout.Column = 2;
axPattern = uiaxes(patternGrid,'Tag','reflectarray3DAxes');
axPattern.Layout.Row = 2;
styleAxes(axPattern);

cutGrid = uigridlayout(views(4),[2 1]);
cutGrid.RowHeight = {43,'1x'}; cutGrid.Padding = [8 8 8 8];
cutGrid.BackgroundColor = PANEL;
cutTools = uigridlayout(cutGrid,[1 8]);
cutTools.ColumnWidth = {27,95,43,95,95,42,70,'1x'};
cutTools.Padding = [0 0 0 0]; cutTools.ColumnSpacing = 4;
cutTools.BackgroundColor = PANEL;
label(cutTools,'Cut',1,1);
ddCut = dropdown(cutTools,{'Theta cut','Phi cut'},'Theta cut',1,2, ...
    @(~,~)refreshCutViews());
ddCut.Tag = 'reflectarrayCut';
ddCut.Tooltip = ['Theta cut sweeps theta in a selected phi plane; ' ...
    'Phi cut sweeps phi from 0° to 360° at a fixed theta.'];
label(cutTools,'Display',1,3);
ddPlotStyle = dropdown(cutTools, ...
    {'Polar','Cartesian'},'Polar',1,4, ...
    @(~,~)refresh2D());
ddPlotStyle.Tag = 'reflectarrayPlotStyle';
cbAbsoluteDbi = uicheckbox(cutTools,'Text','Absolute dBi', ...
    'Value',false,'FontColor',WHITE,'FontSize',11, ...
    'Tooltip','Show modeled directivity in dBi instead of peak-relative field level.', ...
    'Tag','reflectarrayAbsoluteDbi', ...
    'ValueChangedFcn',@(~,~)refresh2D());
cbAbsoluteDbi.Layout.Column = 5;
lblCutTheta = label(cutTools,'Cut θ',1,6);
spCutTheta = number(cutTools,S.beamThetaDeg,[0 90],1,1,7, ...
    @(~,~)refreshCutViews());
spCutTheta.Tag = 'reflectarrayCutTheta';
lblCutPhi = label(cutTools,'Cut φ',1,6);
spCutPhi = number(cutTools,S.beamPhiDeg,[0 360],1,1,7, ...
    @(~,~)refreshCutViews());
spCutPhi.Tag = 'reflectarrayCutPhi';
cbFollowBeamCut = uicheckbox(cutTools,'Text','Follow beam', ...
    'Value',true,'FontColor',WHITE,'FontSize',11, ...
    'Tooltip','Use the commanded angle for the selected cut plane.', ...
    'Tag','reflectarrayFollowBeamCut', ...
    'ValueChangedFcn',@(~,~)refreshCutViews());
cbFollowBeamCut.Layout.Column = 8;
plotHost = uigridlayout(cutGrid,[1 1]);
plotHost.Layout.Row = 2;
plotHost.Padding = [8 8 8 8];
plotHost.BackgroundColor = PANEL;
axCartesian = uiaxes(plotHost,'Tag','reflectarray2DCartesian');
axCartesian.Layout.Row = 1; axCartesian.Layout.Column = 1;
styleAxes(axCartesian);
axPolar = polaraxes(plotHost);
axPolar.Tag = 'reflectarray2DPolar';
axPolar.Layout.Row = 1; axPolar.Layout.Column = 1;
axPolar.Color = PLOT; axPolar.ThetaColor = WHITE;
axPolar.RColor = WHITE; axPolar.GridColor = MUTED;
axPolar.Visible = 'off';

sweepGrid = uigridlayout(views(5),[2 1]);
sweepGrid.RowHeight = {58,'1x'};
sweepGrid.Padding = [10 10 10 10]; sweepGrid.RowSpacing = 7;
sweepGrid.BackgroundColor = PANEL;
sweepTools = uigridlayout(sweepGrid,[2 9]);
sweepTools.RowHeight = {27,20};
sweepTools.ColumnWidth = {45,85,45,85,42,55,90,'1x',85};
sweepTools.Padding = [0 0 0 0]; sweepTools.ColumnSpacing = 5;
sweepTools.BackgroundColor = PANEL;
label(sweepTools,'Start',1,1);
spSweepStart = number(sweepTools,sweepStartGHz,[0.01 300],0.1,1,2, ...
    @(src,~)setSweepLimit('start',src.Value));
spSweepStart.Tag = 'reflectarraySweepStart';
label(sweepTools,'Stop',1,3);
spSweepStop = number(sweepTools,sweepStopGHz,[0.01 300],0.1,1,4, ...
    @(src,~)setSweepLimit('stop',src.Value));
spSweepStop.Tag = 'reflectarraySweepStop';
label(sweepTools,'Points',1,5);
spSweepPoints = number(sweepTools,21,[3 41],2,1,6,@(~,~)invalidateSweep());
spSweepPoints.RoundFractionalValues = 'on';
spSweepPoints.Tag = 'reflectarraySweepPoints';
bCenterSweep = uibutton(sweepTools,'Text','Center ±10%', ...
    'BackgroundColor',NAVY,'FontColor',WHITE, ...
    'Tag','reflectarrayCenterSweep', ...
    'ButtonPushedFcn',@(~,~)centerSweep());
bCenterSweep.Layout.Row = 1; bCenterSweep.Layout.Column = 7;
bRunSweep = uibutton(sweepTools,'Text','Run sweep', ...
    'BackgroundColor',NAVY,'FontColor',WHITE, ...
    'Tag','reflectarrayRunSweep', ...
    'ButtonPushedFcn',@(~,~)runSweep());
bRunSweep.Layout.Row = 1; bRunSweep.Layout.Column = 9;
sweepInfo = uilabel(sweepTools,'Text','', ...
    'FontColor',MUTED,'FontSize',10, ...
    'Tag','reflectarraySweepInfo');
sweepInfo.Layout.Row = 2; sweepInfo.Layout.Column = [1 9];
sweepPlots = uigridlayout(sweepGrid,[2 2]);
sweepPlots.Layout.Row = 2;
sweepPlots.RowHeight = {'1x','1x'};
sweepPlots.ColumnWidth = {'1x','1x'};
sweepPlots.Padding = [0 0 0 0];
sweepPlots.RowSpacing = 4; sweepPlots.ColumnSpacing = 4;
sweepPlots.BackgroundColor = PANEL;
axSweepDirectivity = uiaxes(sweepPlots, ...
    'Tag','reflectarraySweepDirectivityAxes');
axSweepDirectivity.Layout.Row = 1;
axSweepDirectivity.Layout.Column = [1 2];
styleAxes(axSweepDirectivity);
axSweepTheta = uiaxes(sweepPlots, ...
    'Tag','reflectarraySweepThetaAxes');
axSweepTheta.Layout.Row = 2; axSweepTheta.Layout.Column = 1;
styleAxes(axSweepTheta);
axSweepPhi = uiaxes(sweepPlots, ...
    'Tag','reflectarraySweepPhiAxes');
axSweepPhi.Layout.Row = 2; axSweepPhi.Layout.Column = 2;
styleAxes(axSweepPhi);

status = uilabel(root,'Text','Ready','FontColor',MUTED, ...
    'Tag','reflectarrayStatus', ...
    'FontSize',11,'VerticalAlignment','center');
status.Layout.Row = 3;

createShapeGallery();
fig.AutoResizeChildren = 'off';
fig.SizeChangedFcn = @(~,~)onFigureResize();
onFigureResize();
showMapMode();
showView('Geometry');
refreshAll();
updateWindowTitle();

    function p = section(parent,titleText)
        p = uipanel(parent,'Title',titleText,'BackgroundColor',PANEL, ...
            'ForegroundColor',WHITE,'FontWeight','bold','FontSize',11);
    end

    function gg = sectionGrid(parent,nRows)
        gg = uigridlayout(parent,[nRows 4]);
        gg.ColumnWidth = {70,'1x',70,'1x'};
        gg.RowHeight = repmat({23},1,nRows);
        gg.Padding = [7 3 7 3]; gg.RowSpacing = 3;
        gg.ColumnSpacing = 3; gg.BackgroundColor = PANEL;
    end

    function h = label(parent,str,row,col)
        h = uilabel(parent,'Text',str,'FontColor',WHITE,'FontSize',11);
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = number(parent,value,limits,step,row,col,cb)
        h = uispinner(parent,'Value',value,'Limits',limits,'Step',step, ...
            'ValueChangedFcn',cb,'FontSize',11);
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = dropdown(parent,items,value,row,col,cb)
        h = uidropdown(parent,'Items',items,'Value',value, ...
            'ValueChangedFcn',cb,'FontSize',11);
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function styleAxes(ax)
        ax.Color = PLOT; ax.XColor = WHITE; ax.YColor = WHITE;
        ax.GridColor = MUTED; ax.FontSize = 12;
        grid(ax,'on');
    end

    function showView(name)
        idx = find(strcmp(name,viewNames),1);
        for z = 1:numel(viewNames)
            views(z).Visible = localOnOff(z == idx);
            if z == idx
                viewButtons(z).BackgroundColor = BG;
            else
                viewButtons(z).BackgroundColor = NAVY;
            end
        end
    end

    function showMapMode()
        illuminationOn = strcmp(ddMapMode.Value,'Illumination map');
        phaseGrid.Visible = localOnOff(~illuminationOn);
        illuminationPlots.Visible = localOnOff(illuminationOn);
        if illuminationOn
            mapNote.Text = ['Incident field and effective aperture amplitude ' ...
                '(relative dB)'];
            drawIllumination();
        else
            mapNote.Text = 'Four phase terms at the current synthesis reference';
        end
    end

    function setGridSize()
        newM = round(spM.Value); newN = round(spN.Value);
        spM.Value = newM; spN.Value = newN;
        if newM == S.M && newN == S.N, return; end
        S.M = newM; S.N = newN;
        S.manualEnabled = false; S.manualPhaseDeg = [];
        S.activeMask = [];
        S.xOffsetLambda = zeros(S.M,S.N);
        S.yOffsetLambda = zeros(S.M,S.N);
        if ~isempty(S.selected)
            S.selected = [min(S.selected(1),S.M) min(S.selected(2),S.N)];
            S.selectedCells = S.selected;
        end
        lastDeletedCells = zeros(0,2);
        lastDeletedSelection = zeros(0,2);
        markDirty();
        refreshAll();
    end

    function setSetting(name,value)
        if strcmp(name,'staggerAngleDeg')
            S.rowStaggerLambda = max(-10,min(10, ...
                S.dyLambda*tand(value)));
        else
            S.(name) = value;
        end
        if ismember(name,{'staggerAngleDeg','rowStaggerLambda','dyLambda'})
            S.staggerAngleDeg = atan2d(S.rowStaggerLambda,S.dyLambda);
        end
        markDirty();
        refreshAll();
    end

    function chooseShape(name)
        S.shape = name;
        S.activeMask = [];
        lastDeletedCells = zeros(0,2);
        lastDeletedSelection = zeros(0,2);
        rememberShape(name);
        shapeGallery.Visible = 'off';
        markDirty();
        refreshAll();
    end

    function markDirty()
        dirty = true;
        updateWindowTitle();
    end

    function updateWindowTitle()
        if isempty(projectFile), name = 'Untitled';
        else, [~,name] = fileparts(projectFile); end
        if dirty, suffix = ' *'; else, suffix = ''; end
        fig.Name = ['Reflectarray Designer — ' name suffix];
    end

    function scale = frequencyUnitScale()
        if strcmp(S.freqUnit,'MHz'), scale = 1000;
        else, scale = 1; end
    end

    function syncFrequencyControls()
        scale = frequencyUnitScale();
        for spinner = [spDesignFreq spOperatingFreq]
            spinner.Limits = [0.01 300000];
            spinner.Step = 0.5*scale;
        end
        spDesignFreq.Value = S.designFreqGHz*scale;
        spOperatingFreq.Value = S.opFreqGHz*scale;
        for spinner = [spDesignFreq spOperatingFreq]
            spinner.Limits = [0.01 300]*scale;
        end
        ddFreqUnit.Value = S.freqUnit;
        lblDesignLambda.Text = sprintf( ...
            'Design λ = %.4g mm  ·  dx × dy = %.3g × %.3g mm', ...
            A.designLambdaMM,A.dxMM,A.dyMM);
        for spinner = [spSweepStart spSweepStop]
            spinner.Limits = [0.01 300000];
            spinner.Step = 0.1*scale;
        end
        spSweepStart.Value = sweepStartGHz*scale;
        spSweepStop.Value = sweepStopGHz*scale;
        for spinner = [spSweepStart spSweepStop]
            spinner.Limits = [0.01 300]*scale;
        end
        if ~isempty(sweepResult), drawSweep();
        else, updateSweepInfo(); end
    end

    function setDesignFrequency(value)
        S.designFreqGHz = value/frequencyUnitScale();
        markDirty();
        refreshAll();
    end

    function setOperatingFrequency(value)
        S.opFreqGHz = value/frequencyUnitScale();
        markDirty();
        refreshAll();
    end

    function setFrequencyUnit()
        S.freqUnit = ddFreqUnit.Value;
        markDirty();
        syncFrequencyControls();
    end

    function setSweepLimit(which,value)
        if strcmp(which,'start')
            sweepStartGHz = value/frequencyUnitScale();
        else
            sweepStopGHz = value/frequencyUnitScale();
        end
        invalidateSweep();
    end

    function centerSweep()
        sweepStartGHz = max(0.01,0.9*S.opFreqGHz);
        sweepStopGHz = min(300,1.1*S.opFreqGHz);
        syncFrequencyControls();
        invalidateSweep();
    end

    function updateSweepInfo()
        if S.manualEnabled
            phaseDescription = 'fixed manually assigned phases';
        else
            phaseDescription = sprintf('fixed phases set at %.4g GHz', ...
                A.phaseReferenceGHz);
        end
        sweepInfo.Text = sprintf(['Frequency in %s · %s · ' ...
            'modeled directivity excludes cell dispersion and loss'], ...
            S.freqUnit,phaseDescription);
    end

    function invalidateSweep()
        sweepResult = [];
        for ax = [axSweepDirectivity axSweepTheta axSweepPhi]
            cla(ax);
            grid(ax,'on');
        end
        title(axSweepDirectivity,'Run sweep to evaluate fixed phases', ...
            'Color',WHITE);
        title(axSweepTheta,'Peak θ','Color',WHITE);
        title(axSweepPhi,'Peak φ','Color',WHITE);
        updateSweepInfo();
    end

    function runSweep()
        if sweepStartGHz >= sweepStopGHz
            uialert(fig,'Stop frequency must exceed start frequency.', ...
                'Frequency sweep');
            return;
        end
        if ~any(A.active(:))
            uialert(fig,'Add at least one reflecting cell before sweeping.', ...
                'Frequency sweep');
            return;
        end
        bRunSweep.Enable = 'off';
        sweepInfo.Text = 'Computing fixed-phase frequency sweep…';
        drawnow;
        try
            frequencies = linspace(sweepStartGHz,sweepStopGHz, ...
                round(spSweepPoints.Value));
            sweepResult = reflectarrayFrequencySweep(S,phaseDeg,frequencies);
            drawSweep();
        catch err
            invalidateSweep();
            uialert(fig,err.message,'Frequency sweep');
        end
        bRunSweep.Enable = 'on';
    end

    function drawSweep()
        if isempty(sweepResult), return; end
        scale = frequencyUnitScale();
        f = sweepResult.frequencyGHz*scale;
        for ax = [axSweepDirectivity axSweepTheta axSweepPhi]
            cla(ax); hold(ax,'on');
        end
        plot(axSweepDirectivity,f,sweepResult.peakDbi, ...
            'Color',BLUE,'LineWidth',2.1,'DisplayName','Sampled peak');
        plot(axSweepDirectivity,f,sweepResult.targetDbi, ...
            'Color',[1 0.72 0.25],'LineWidth',1.8, ...
            'DisplayName','Commanded direction');
        ylabel(axSweepDirectivity,'Directivity (dBi)');
        title(axSweepDirectivity,'Fixed-phase directivity', ...
            'Color',WHITE);
        legend(axSweepDirectivity,'Location','best');

        plot(axSweepTheta,f,sweepResult.peakThetaDeg, ...
            'Color',BLUE,'LineWidth',2.1);
        yline(axSweepTheta,S.beamThetaDeg,'--', ...
            'Color',[1 0.72 0.25]);
        ylabel(axSweepTheta,'Peak θ (°)');
        title(axSweepTheta,'Beam pointing · θ','Color',WHITE);

        % Display the equivalent azimuth nearest the commanded azimuth so
        % a peak crossing 0/360° does not appear to jump by a full turn.
        phiDisplay = S.beamPhiDeg + mod( ...
            sweepResult.peakPhiDeg-S.beamPhiDeg+180,360)-180;
        plot(axSweepPhi,f,phiDisplay, ...
            'Color',BLUE,'LineWidth',2.1);
        if S.beamThetaDeg > 0.25
            yline(axSweepPhi,S.beamPhiDeg,'--', ...
                'Color',[1 0.72 0.25]);
        end
        finitePhi = phiDisplay(isfinite(phiDisplay));
        if isempty(finitePhi)
            ylim(axSweepPhi,[0 360]);
        else
            phiExtent = finitePhi;
            if S.beamThetaDeg > 0.25
                phiExtent = [phiExtent S.beamPhiDeg];
            end
            padding = max(3,0.1*(max(phiExtent)-min(phiExtent)));
            lower = min(phiExtent)-padding;
            upper = max(phiExtent)+padding;
            ylim(axSweepPhi,[lower upper]);
        end
        ylabel(axSweepPhi,'Peak φ (°; continuous)');
        title(axSweepPhi,'Beam pointing · φ','Color',WHITE);

        for ax = [axSweepDirectivity axSweepTheta axSweepPhi]
            xlabel(ax,['Frequency (' S.freqUnit ')']);
            xlim(ax,[f(1) f(end)]);
            if S.opFreqGHz >= sweepStartGHz && ...
                    S.opFreqGHz <= sweepStopGHz
                xline(ax,S.opFreqGHz*scale,':','Color',MUTED, ...
                    'DisplayName','Current operating frequency');
            end
            grid(ax,'on'); hold(ax,'off');
        end
        updateSweepInfo();
    end

    function setPhaseMode()
        S.retunePhase = ~cbFixedPhase.Value;
        markDirty();
        refreshAll();
    end

    function setFeed(which,value)
        S.feedMM(which) = value;
        markDirty();
        refreshAll();
    end

    function setIllumination()
        S.illuminationMode = ddIllumination.Value;
        syncIlluminationControls();
        markDirty();
        refreshAll();
    end

    function setIncident(name,value)
        S.(name) = value;
        markDirty();
        refreshAll();
    end

    function syncIlluminationControls()
        ddIllumination.Value = S.illuminationMode;
        spIncidentTheta.Value = S.incidentThetaDeg;
        spIncidentPhi.Value = S.incidentPhiDeg;
        hornOn = strcmp(S.illuminationMode,'Horn');
        hornControls = {lblFeedX,spFx,lblFeedY,spFy,lblFeedZ,spFz,lblFeedQ,spQ};
        planeControls = {lblIncidentTheta,spIncidentTheta,lblIncidentPhi, ...
            spIncidentPhi,planeNote};
        for controlIndex = 1:numel(hornControls)
            hornControls{controlIndex}.Visible = localOnOff(hornOn);
        end
        for controlIndex = 1:numel(planeControls)
            planeControls{controlIndex}.Visible = localOnOff(~hornOn);
        end
        if hornOn
            taperNote.Text = ['Phase correction aligns cells but does not flatten ' ...
                'horn illumination. Hann adds extra amplitude taper.'];
        else
            taperNote.Text = ['Plane-wave illumination is uniform. Hann adds ' ...
                'a mathematical amplitude taper.'];
        end
    end

    function setConvention()
        S.angleConvention = ddConvention.Value;
        markDirty();
        if strcmp(S.angleConvention,'Theta / phi')
            lblTheta.Text = 'Theta (°)'; lblPhi.Text = 'Phi (°)';
            spTheta.Value = S.beamThetaDeg;
            spTheta.Tooltip = 'Theta is measured from +z.';
        else
            lblTheta.Text = 'Elev. (°)'; lblPhi.Text = 'Azim. (°)';
            spTheta.Value = 90-S.beamThetaDeg;
            spTheta.Tooltip = 'Elevation is measured above the xy plane.';
        end
        spPhi.Value = S.beamPhiDeg;
    end

    function setBeam()
        if strcmp(S.angleConvention,'Theta / phi')
            S.beamThetaDeg = spTheta.Value;
        else
            S.beamThetaDeg = 90-spTheta.Value;
        end
        S.beamPhiDeg = spPhi.Value;
        markDirty();
        refreshAll();
    end

    function restoreIdeal()
        if ~S.manualEnabled, return; end
        S.manualEnabled = false;
        S.manualPhaseDeg = [];
        markDirty();
        refreshAll();
    end

    function editSelected(value)
        if isempty(S.selected), return; end
        row = S.selected(1); col = S.selected(2);
        if ~A.active(row,col), return; end
        if ~S.manualEnabled
            S.manualPhaseDeg = A.referencePhaseDeg;
        end
        S.manualEnabled = true;
        S.manualPhaseDeg(row,col) = mod(value,360);
        markDirty();
        refreshAll();
    end

    function onTableEdit(~,event)
        if event.Indices(2) ~= 3, refreshPhaseTable(); return; end
        [cols,rows] = find(A.active.');
        row = rows(event.Indices(1)); col = cols(event.Indices(1));
        if ~isnumeric(event.NewData) || ...
                ~isscalar(event.NewData) || ~isfinite(event.NewData)
            refreshPhaseTable(); return;
        end
        S.selected = [row col];
        S.selectedCells = S.selected;
        editSelected(event.NewData);
    end

    function onTableSelect(~,event)
        if isempty(event.Indices), return; end
        [cols,rows] = find(A.active.');
        chosen = unique(event.Indices(:,1),'stable');
        S.selectedCells = [rows(chosen) cols(chosen)];
        S.selected = S.selectedCells(1,:);
        drawPhaseMap(); updateSelected(); updateGeometrySelection();
        syncEditorControls();
    end

    function onMapClick(sourceAxes)
        point = sourceAxes.CurrentPoint;
        col = round(point(1,1));
        row = round(point(1,2));
        if row >= 1 && row <= S.M && col >= 1 && col <= S.N && ...
                ~A.active(row,col)
            restoreCell(row,col);
        else
            selectCell(row,col);
        end
    end

    function onGeometryCellClick(~,event)
        point = event.IntersectionPoint;
        distance = (A.X-point(1)).^2+(A.Y-point(2)).^2;
        distance(~A.active) = inf;
        [~,index] = min(distance(:));
        if ~isfinite(distance(index)), return; end
        [row,col] = ind2sub(size(A.active),index);
        extend = false;
        try
            modifiers = fig.CurrentModifier;
            extend = any(ismember(modifiers,{'shift','control','command'}));
        catch
        end
        selectCell(row,col,extend);
    end

    function selectCell(row,col,extend)
        if nargin < 3, extend = false; end
        if row < 1 || row > S.M || col < 1 || col > S.N || ...
                ~A.active(row,col), return; end
        S.selected = [row col];
        if extend
            match = find(all(S.selectedCells == [row col],2),1);
            if isempty(match)
                S.selectedCells(end+1,:) = [row col];
            elseif size(S.selectedCells,1) > 1
                S.selectedCells(match,:) = [];
                S.selected = S.selectedCells(1,:);
            end
        else
            S.selectedCells = S.selected;
        end
        drawPhaseMap();
        updateSelected();
        updateGeometrySelection();
        syncEditorControls();
    end

    function refreshAll(throwOnError)
        try
            A = reflectarrayAperture(S);
            if ~S.manualEnabled || ...
                    ~isequal(size(S.manualPhaseDeg),size(A.referencePhaseDeg))
                S.manualEnabled = false;
                S.manualPhaseDeg = A.referencePhaseDeg;
            end
            phaseDeg = S.manualPhaseDeg;
            directivity = reflectarrayDirectivity(A,phaseDeg);
            idealDirectivity = [];
            if any(A.active(:))
                if ~isempty(S.selectedCells)
                    inBounds = S.selectedCells(:,1) <= S.M & ...
                        S.selectedCells(:,2) <= S.N;
                    S.selectedCells = S.selectedCells(inBounds,:);
                    chosen = sub2ind(size(A.active), ...
                        S.selectedCells(:,1),S.selectedCells(:,2));
                    S.selectedCells = S.selectedCells(A.active(chosen),:);
                end
                if isempty(S.selectedCells)
                    S.selected = zeros(0,2);
                elseif isempty(S.selected) || ...
                        ~any(all(S.selectedCells == S.selected,2))
                    S.selected = S.selectedCells(1,:);
                end
            else
                S.selected = zeros(0,2);
                S.selectedCells = zeros(0,2);
            end
            modeText.Text = localModeText();
            btnIdeal.Enable = localOnOff(S.manualEnabled);
            syncFrequencyControls();
            refreshPhaseTable();
            updateSelected();
            patternData = computePattern3D();
            drawGeometry();
            drawPhaseMap();
            refresh2D();
            drawPattern3D();
            invalidateSweep();
            wanted = abs(reflectarrayField(A,phaseDeg, ...
                S.beamThetaDeg,S.beamPhiDeg));
            ideal = abs(reflectarrayField(A,A.idealPhaseDeg, ...
                S.beamThetaDeg,S.beamPhiDeg));
            loss = 20*log10(max(wanted,eps)/max(ideal,eps));
            if isfinite(directivity.factor)
                targetDbi = 10*log10(max(directivity.factor*wanted^2,realmin));
                peakMagnitude = max(patternData.fieldMagnitude(:));
                sampledPeakDbi = 10*log10(max( ...
                    directivity.factor*peakMagnitude^2,realmin));
                status.Text = sprintf(['%d cells  ·  modeled D toward target %.2f dBi' ...
                    '  ·  sampled peak %.2f dBi  ·  target field change vs retuned %.2f dB'], ...
                    nnz(A.active),targetDbi,sampledPeakDbi,loss);
            else
                status.Text = sprintf('%d cells · modeled directivity undefined', ...
                    nnz(A.active));
            end
            if ~any(A.active(:))
                status.Text = '0 active cells · choose a shape or restore a cell to calculate a pattern';
            end
            syncEditorControls();
        catch err
            status.Text = ['Calculation failed: ' err.message];
            if nargin > 0 && throwOnError, rethrow(err); end
            uialert(fig,err.message,'Reflectarray calculation');
        end
    end

    function txt = localModeText()
        if S.manualEnabled
            txt = 'Manual cell phases';
        else
            txt = sprintf('Phases set at %.4g GHz',A.phaseReferenceGHz);
        end
    end

    function refreshPhaseTable()
        [cols,rows] = find(A.active.');
        indices = sub2ind(size(A.active),rows,cols);
        phaseTable.Data = [rows cols round(phaseDeg(indices),2)];
        phaseTable.ColumnName = {'Row','Column','Phase (°)'};
        phaseTable.RowName = {};
        phaseTable.ColumnEditable = [false false true];
        phaseTable.ColumnFormat = {'short','short','bank'};
        phaseTable.ColumnWidth = {65 70 110};
    end

    function updateSelected()
        if ~any(A.active(:))
            selectedText.Text = 'No active cells';
            spCell.Value = 0;
            spCell.Enable = 'off';
            return;
        end
        if isempty(S.selected)
            selectedText.Text = 'No cell selected';
            spCell.Value = 0;
            spCell.Enable = 'off';
            return;
        end
        spCell.Enable = 'on';
        row = S.selected(1); col = S.selected(2);
        if strcmp(A.illuminationMode,'Horn')
            selectedText.Text = sprintf('Cell (%d, %d)   target %.1f°   path %.1f mm', ...
                row,col,A.idealPhaseDeg(row,col),A.R(row,col));
        else
            incidentPhase = mod(-rad2deg(A.k*A.incidentPathMM(row,col)),360);
            selectedText.Text = sprintf('Cell (%d, %d)   target %.1f°   incident %.1f°', ...
                row,col,A.idealPhaseDeg(row,col),incidentPhase);
        end
        if size(S.selectedCells,1) > 1
            selectedText.Text = sprintf('%d cells selected · first: (%d, %d)', ...
                size(S.selectedCells,1),row,col);
        end
        spCell.Value = phaseDeg(row,col);
    end

    function drawGeometry()
        if geometryViewInitialized
            [cameraAz,cameraEl] = view(axGeom);
        end
        cla(axGeom); hold(axGeom,'on');
        axGeom.ZColor = WHITE;
        x = A.X(A.active); x = x(:);
        y = A.Y(A.active); y = y(:);
        nCells = numel(x);
        halfX = 0.34*A.dxMM; halfY = 0.34*A.dyMM;
        cellVertices = zeros(4*nCells,3);
        cellVertices(1:4:end,1:2) = [x-halfX y-halfY];
        cellVertices(2:4:end,1:2) = [x+halfX y-halfY];
        cellVertices(3:4:end,1:2) = [x+halfX y+halfY];
        cellVertices(4:4:end,1:2) = [x-halfX y+halfY];
        cellFaces = reshape(1:4*nCells,4,nCells).';
        patch(axGeom,'Vertices',cellVertices,'Faces',cellFaces, ...
            'FaceColor',BLUE,'EdgeColor',[0.35 0.76 1], ...
            'FaceLighting','none','PickableParts','visible', ...
            'ButtonDownFcn',@onGeometryCellClick,'Tag','reflectarrayCells');
        xEdge = [min(A.X(:))-A.dxMM/2 max(A.X(:))+A.dxMM/2];
        yEdge = [min(A.Y(:))-A.dyMM/2 max(A.Y(:))+A.dyMM/2];
        plot3(axGeom,xEdge([1 2 2 1 1]),yEdge([1 1 2 2 1]), ...
            zeros(1,5),'Color',MUTED,'LineWidth',1.3);
        apertureSize = min(S.M*A.dyMM,S.N*A.dxMM);
        if strcmp(A.illuminationMode,'Horn')
        F = A.feedMM; d = -F/norm(F);
        side = cross(d,[0 1 0]);
        if norm(side) < 1e-8, side = cross(d,[1 0 0]); end
        side = side/norm(side); up = cross(d,side);
        hornLengthScale = 0.65;
        hornWidthScale = 0.43;
        L = hornLengthScale*min(max(0.38*norm(F),0.18*apertureSize), ...
            0.68*norm(F));
        mouthSize = hornWidthScale*min(0.12*apertureSize,0.28*F(3));
        % A rectangular waveguide feeds four planar flare walls. The thin
        % mouth lip and darker inner walls keep the aperture visibly open.
        % F remains the calculation's feed phase centre; this is visual.
        cornerSigns = [-1 -1;1 -1;1 1;-1 1];
        rectangleAt = @(station,halfW,halfH) F(:)+L*station*d(:)+ ...
            halfW*side(:)*cornerSigns(:,1).'+ ...
            halfH*up(:)*cornerSigns(:,2).';
        wallFaces = [1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8];
        mouthHalfW = 1.05*mouthSize;
        mouthHalfH = 0.78*mouthSize;
        neckHalfW = 0.23*mouthHalfW;
        neckHalfH = 0.23*mouthHalfH;
        neckEnd = 0.36;
        backRing = rectangleAt(0,neckHalfW,neckHalfH);
        neckRing = rectangleAt(neckEnd,neckHalfW,neckHalfH);
        mouthRing = rectangleAt(1,mouthHalfW,mouthHalfH);
        neckGold = [0.82 0.53 0.04;1 0.81 0.20; ...
            0.92 0.66 0.07;1 0.76 0.13];
        flareGold = [0.87 0.57 0.04;1 0.84 0.25; ...
            0.96 0.71 0.10;1 0.78 0.16];
        patch(axGeom,'Vertices',[backRing neckRing].','Faces',wallFaces, ...
            'FaceVertexCData',neckGold,'FaceColor','flat', ...
            'EdgeColor','none','Tag','reflectarrayHornSegment');
        patch(axGeom,'Vertices',[neckRing mouthRing].','Faces',wallFaces, ...
            'FaceVertexCData',flareGold,'FaceColor','flat', ...
            'EdgeColor','none','Tag','reflectarrayHornSegment');
        innerMouth = rectangleAt(1,0.87*mouthHalfW,0.87*mouthHalfH);
        patch(axGeom,'Vertices',[mouthRing innerMouth].', ...
            'Faces',wallFaces,'FaceColor',[1 0.87 0.28], ...
            'EdgeColor',[0.62 0.41 0.04],'LineWidth',0.6, ...
            'Tag','reflectarrayHornOpening');
        innerNeck = rectangleAt(neckEnd,0.76*neckHalfW,0.76*neckHalfH);
        patch(axGeom,'Vertices',[innerNeck innerMouth].', ...
            'Faces',wallFaces,'FaceColor',[0.28 0.20 0.04], ...
            'EdgeColor','none','Tag','reflectarrayHornInterior');
        rearPort = rectangleAt(0,0.80*neckHalfW,0.80*neckHalfH);
        patch(axGeom,'Vertices',rearPort.','Faces',1:4, ...
            'FaceColor',[0.10 0.09 0.035], ...
            'EdgeColor',[1 0.88 0.31],'LineWidth',1.2, ...
            'Tag','reflectarrayHornRearPort');
        % Feed-to-corner rays make the change in path length visible.
        corners = [xEdge(1) yEdge(1);xEdge(2) yEdge(1); ...
            xEdge(1) yEdge(2);xEdge(2) yEdge(2)];
        for z = 1:4
            plot3(axGeom,[F(1) corners(z,1)],[F(2) corners(z,2)], ...
                [F(3) 0],':','Color',[0.57 0.48 0.28]);
        end
        else
            drawPlaneWaveGeometry(apertureSize);
        end
        if showGeometryBeam.Value && ~isempty(patternData.x) && ...
                isfinite(directivity.factor)
            beamScale = 0.5*apertureSize;
            surf(axGeom,beamScale*patternData.x, ...
                beamScale*patternData.y,beamScale*patternData.z, ...
                patternData.db,'EdgeColor','none', ...
                'FaceColor','interp','FaceAlpha',1, ...
                'HitTest','off','PickableParts','none', ...
                'Tag','reflectarrayGeometryBeam');
            colormap(axGeom,turbo(256)); clim(axGeom,[-40 0]);
        end
        selectedGeometryHighlight = patch(axGeom, ...
            'Vertices',zeros(4,3),'Faces',1:4,'FaceColor','none', ...
            'EdgeColor',[1 0.85 0.18],'LineWidth',2.5, ...
            'HitTest','off','PickableParts','none', ...
            'Tag','reflectarraySelectedGeometryCell');
        updateGeometrySelection();
        xlabel(axGeom,'x (mm)'); ylabel(axGeom,'y (mm)');
        zlabel(axGeom,'z (mm)');
        if strcmp(A.illuminationMode,'Horn')
            sourceTitle = 'Horn feed';
        else
            sourceTitle = 'Plane-wave illumination';
        end
        if showGeometryBeam.Value
            title(axGeom,sprintf('%s, %d cells, and reflected-field pattern', ...
                sourceTitle,nnz(A.active)),'Color',WHITE);
        else
            title(axGeom,sprintf('%s and %d active reflecting cells', ...
                sourceTitle,nnz(A.active)),'Color',WHITE);
        end
        axis(axGeom,'equal');
        if geometryViewInitialized
            view(axGeom,cameraAz,cameraEl);
        else
            view(axGeom,42,25);
            geometryViewInitialized = true;
        end
        grid(axGeom,'on');
        hold(axGeom,'off');
        axGeom.Interactions = zoomInteraction;
        axGeom.ButtonDownFcn = @onGeometryEmptyClick;
        % Keep the cells selectable while decorative plot objects pass
        % clicks through to the axes for empty-site selection.
        children = axGeom.Children;
        for child = children.'
            if ~isprop(child,'PickableParts'), continue; end
            isSelectableCell = strcmp(child.Tag,'reflectarrayCells');
            if isSelectableCell
                child.PickableParts = 'visible';
            else
                child.PickableParts = 'none';
            end
            child.HitTest = localOnOff(isSelectableCell);
        end
    end

    function drawPlaneWaveGeometry(apertureSize)
        % The red wavefront is perpendicular to the incoming ray.
        % The drawing indicates direction only; cell phase/amplitude come
        % from reflectarrayAperture, not from this graphic.
        theta = A.incidentThetaDeg;
        phi = A.incidentPhiDeg;
        sourceDirection = [sind(theta)*cosd(phi), ...
            sind(theta)*sind(phi),cosd(theta)];
        tangent1 = [-sind(phi),cosd(phi),0];
        tangent2 = cross(sourceDirection,tangent1);
        halfWidth = 0.25*apertureSize;
        corners = [-1 -1;1 -1;1 1;-1 1];
        centre = [0 0 0.35*apertureSize];
        vertices = centre + halfWidth*corners(:,1)*tangent1 + ...
            halfWidth*corners(:,2)*tangent2;
        patch(axGeom,'Vertices',vertices,'Faces',1:4, ...
            'FaceColor',[0.90 0.13 0.16],'FaceAlpha',0.55, ...
            'EdgeColor',[1 0.34 0.31],'LineWidth',1.4, ...
            'FaceLighting','none','HitTest','off', ...
            'PickableParts','none','Tag','reflectarrayPlaneWavefront');
        % One centred arrow shows the incoming propagation direction.
        quiver3(axGeom,centre(1),centre(2),centre(3), ...
            -sourceDirection(1)*0.15*apertureSize, ...
            -sourceDirection(2)*0.15*apertureSize, ...
            -sourceDirection(3)*0.15*apertureSize,0, ...
            'Color',[1 0.34 0.31],'LineWidth',1.8, ...
            'MaxHeadSize',0.55,'HitTest','off', ...
            'PickableParts','none','Tag','reflectarrayPlaneWaveRay');
    end

    function onGeometryEmptyClick(~,event)
        try
            point = event.IntersectionPoint;
        catch
            ray = axGeom.CurrentPoint;
            if abs(ray(2,3)-ray(1,3)) < eps, return; end
            fraction = -ray(1,3)/(ray(2,3)-ray(1,3));
        point = ray(1,:)+fraction*(ray(2,:)-ray(1,:));
        end
        if any(~isfinite(point)) || abs(point(3)) > ...
                max(1e-6,0.01*min(A.dxMM,A.dyMM))
            return;
        end
        distance = (A.X-point(1)).^2+(A.Y-point(2)).^2;
        [~,index] = min(distance(:));
        [row,col] = ind2sub(size(A.active),index);
        if abs(point(1)-A.X(row,col)) > A.dxMM/2 || ...
                abs(point(2)-A.Y(row,col)) > A.dyMM/2
            return;
        end
        if ~A.active(row,col), restoreCell(row,col);
        else, selectCell(row,col); end
    end

    function updateGeometrySelection()
        if ~isgraphics(selectedGeometryHighlight), return; end
        if isempty(S.selectedCells)
            selectedGeometryHighlight.Visible = 'off';
            return;
        end
        selectedGeometryHighlight.Visible = 'on';
        halfX = 0.34*A.dxMM; halfY = 0.34*A.dyMM;
        corners = [-halfX -halfY;halfX -halfY;halfX halfY;-halfX halfY];
        count = size(S.selectedCells,1);
        vertices = zeros(4*count,3);
        for t = 1:count
            row = S.selectedCells(t,1); col = S.selectedCells(t,2);
            vertices(4*t-3:4*t,1:2) = ...
                corners+[A.X(row,col) A.Y(row,col)];
            vertices(4*t-3:4*t,3) = 0.2;
        end
        selectedGeometryHighlight.Vertices = vertices;
        selectedGeometryHighlight.Faces = reshape(1:4*count,4,count).';
    end

    function drawPhaseMap()
        maps = reflectarrayPhaseMaps(A,phaseDeg,S);
        if strcmp(A.illuminationMode,'Horn')
            incidentTitle = 'Spatial delay · synthesis';
        else
            incidentTitle = 'Plane-wave phase · synthesis';
        end
        paintPhaseMap(axSpatial,maps.spatialDeg,A.active, ...
            1:S.N,1:S.M,incidentTitle,'Column','Row',true);
        progressiveDisplay = maps.progressiveDeg;
        progressiveDisplay(abs(progressiveDisplay)<1e-10) = 360;
        paintPhaseMap(axProgressive,progressiveDisplay,A.active, ...
            1:S.N,1:S.M,'Progressive phase · beam','Column','Row',true);
        paintPhaseMap(axPhase,maps.reflectionDeg,A.active, ...
            1:S.N,1:S.M,'Cell reflection phase · editable', ...
            'Column','Row',true);
        paintPhaseMap(axContinuous,maps.continuousDeg, ...
            maps.continuousMask,maps.continuousXMM, ...
            maps.continuousYMM,'Continuous-aperture phase', ...
            'x (mm)','y (mm)',false);
        if any(A.active(:))
            for t = 1:size(S.selectedCells,1)
                row = S.selectedCells(t,1); col = S.selectedCells(t,2);
                rectangle(axPhase,'Position',[col-0.5,row-0.5,1,1], ...
                    'EdgeColor',[0.3 0.85 1],'LineWidth',2.2, ...
                    'HitTest','off');
            end
        end
        hold(axPhase,'off');
        if strcmp(ddMapMode.Value,'Illumination map')
            drawIllumination();
        end
    end

    function drawIllumination()
        if strcmp(A.illuminationMode,'Horn')
            incidentTitle = 'Horn incident field amplitude';
        else
            incidentTitle = 'Plane-wave incident field amplitude';
        end
        paintAmplitudeMap(axIllumination,A.illumination,incidentTitle);
        paintAmplitudeMap(axApertureAmplitude,A.amplitude, ...
            'Effective aperture amplitude after taper');
    end

    function paintAmplitudeMap(ax,values,heading)
        cla(ax); hold(ax,'on');
        if any(A.active(:))
            peak = max(values(A.active));
        else
            peak = 0;
        end
        if peak > 0
            displayDb = max(-40,20*log10(max(values/peak,1e-8)));
        else
            displayDb = -40*ones(size(values));
        end
        im = imagesc(ax,1:S.N,1:S.M,displayDb,[-40 0]);
        im.AlphaData = double(A.active);
        im.ButtonDownFcn = @(~,~)onMapClick(ax);
        ax.ButtonDownFcn = @(~,~)onMapClick(ax);
        ax.Color = PLOT; ax.YDir = 'normal';
        colormap(ax,turbo(256)); clim(ax,[-40 0]);
        axis(ax,'image'); grid(ax,'off');
        xlim(ax,[0.5 S.N+0.5]); ylim(ax,[0.5 S.M+0.5]);
        xlabel(ax,'Column'); ylabel(ax,'Row');
        title(ax,heading,'Color',WHITE,'FontSize',11);
        cb = colorbar(ax); cb.Color = WHITE;
        cb.Ticks = -40:10:0;
        cb.Label.String = 'Relative field amplitude (dB)';
        for t = 1:size(S.selectedCells,1)
            row = S.selectedCells(t,1); col = S.selectedCells(t,2);
            if row <= S.M && col <= S.N && A.active(row,col)
                rectangle(ax,'Position',[col-0.5,row-0.5,1,1], ...
                    'EdgeColor',[0.3 0.85 1],'LineWidth',2.2, ...
                    'HitTest','off');
            end
        end
        hold(ax,'off');
    end

    function paintPhaseMap(ax,values,mask,x,y,heading,xText,yText,selectable)
        cla(ax); hold(ax,'on');
        displayValues = values;
        displayValues(~mask) = 0;
        im = imagesc(ax,x,y,displayValues,[0 360]);
        im.AlphaData = double(mask);
        ax.Color = PLOT; ax.YDir = 'normal';
        colormap(ax,hot(256)); clim(ax,[0 360]);
        axis(ax,'image'); grid(ax,'off');
        if selectable
            xlim(ax,[0.5 S.N+0.5]); ylim(ax,[0.5 S.M+0.5]);
            im.ButtonDownFcn = @(~,~)onMapClick(ax);
            ax.ButtonDownFcn = @(~,~)onMapClick(ax);
        else
            xlim(ax,[x(1) x(end)]); ylim(ax,[y(1) y(end)]);
        end
        xlabel(ax,xText); ylabel(ax,yText);
        title(ax,heading,'Color',WHITE,'FontSize',11);
        cb = colorbar(ax); cb.Color = WHITE; cb.Ticks = [0 180 360];
        hold(ax,'off');
    end

    function data = computePattern3D()
        data = reflectarrayPattern3D(A,phaseDeg, ...
            S.beamThetaDeg,S.beamPhiDeg);
    end

    function refreshCutViews()
        refresh2D();
        drawPattern3D();
    end

    function [angle,theta,phi] = selectedCutDirections()
        if strcmp(ddCut.Value,'Theta cut')
            angle = -90:0.5:90;
            theta = abs(angle);
            phi = spCutPhi.Value + 180*double(angle < 0);
        else
            angle = 0:1:360;
            theta = spCutTheta.Value*ones(size(angle));
            phi = angle;
        end
    end

    function drawPattern3D()
        if patternViewInitialized
            [cameraAz,cameraEl] = view(axPattern);
        end
        cla(axPattern);
        cutOverlayInfo.Text = '';
        if isempty(patternData.x)
            title(axPattern,'No active reflecting cells','Color',WHITE);
            return;
        end
        if ~isfinite(directivity.factor)
            title(axPattern,'No radiated field to normalize','Color',WHITE);
            return;
        end
        surf(axPattern,patternData.x,patternData.y,patternData.z, ...
            patternData.db,'EdgeColor','none', ...
            'FaceColor','interp');
        colormap(axPattern,turbo(256)); clim(axPattern,[-40 0]);
        axPattern.ZColor = WHITE;
        xlabel(axPattern,'x'); ylabel(axPattern,'y'); zlabel(axPattern,'z');
        sampledPeak = max(patternData.fieldMagnitude(:));
        sampledPeakDbi = 10*log10(max( ...
            directivity.factor*sampledPeak^2,realmin));
        if strcmp(ddCut.Value,'Theta cut')
            cutOverlayInfo.Text = sprintf( ...
                'Selected θ cut: φ = %.1f°',spCutPhi.Value);
        else
            cutOverlayInfo.Text = sprintf( ...
                'Selected φ cut: θ = %.1f°',spCutTheta.Value);
        end
        if cbShow3DCut.Value
            [~,cutTheta,cutPhi] = selectedCutDirections();
            cutField = abs(reflectarrayField(A,phaseDeg,cutTheta,cutPhi));
            cutDb = 20*log10(max(cutField,1e-8)) - ...
                20*log10(max(sampledPeak,1e-8));
            cutRadius = max(0.03,1+max(cutDb,-40)/40) + 0.012;
            hold(axPattern,'on');
            if all(cutTheta == 0)
                plot3(axPattern,0,0,cutRadius(1), ...
                    'o','Color',[1 0.92 0.2],'MarkerSize',8, ...
                    'MarkerFaceColor',[1 0.92 0.2], ...
                    'Tag','reflectarray3DCutTrace','HitTest','off');
            else
                plot3(axPattern,cutRadius.*sind(cutTheta).*cosd(cutPhi), ...
                    cutRadius.*sind(cutTheta).*sind(cutPhi), ...
                    cutRadius.*cosd(cutTheta), ...
                    'Color',[1 0.92 0.2],'LineWidth',2.3, ...
                    'Tag','reflectarray3DCutTrace','HitTest','off');
            end
            hold(axPattern,'off');
        end
        title(axPattern,sprintf( ...
            'Relative 3D pattern · sampled peak D %.2f dBi', ...
            sampledPeakDbi),'Color',WHITE);
        cb = colorbar(axPattern); cb.Color = WHITE;
        cb.Label.String = 'Relative field level (dB)';
        axis(axPattern,'equal');
        if patternViewInitialized
            view(axPattern,cameraAz,cameraEl);
        else
            view(axPattern,35,25);
            patternViewInitialized = true;
        end
        grid(axPattern,'on');
    end

    function refresh2D()
        if isempty(A) || ~isgraphics(axCartesian), return; end
        isThetaCut = strcmp(ddCut.Value,'Theta cut');
        lblCutTheta.Visible = localOnOff(~isThetaCut);
        spCutTheta.Visible = localOnOff(~isThetaCut);
        lblCutPhi.Visible = localOnOff(isThetaCut);
        spCutPhi.Visible = localOnOff(isThetaCut);
        spCutTheta.Enable = localOnOff(~cbFollowBeamCut.Value && ~isThetaCut);
        spCutPhi.Enable = localOnOff(~cbFollowBeamCut.Value && isThetaCut);
        if isThetaCut
            cbFollowBeamCut.Tooltip = 'Use the commanded phi for the theta cut.';
        else
            cbFollowBeamCut.Tooltip = ['Use the commanded theta for the phi cut. ' ...
                'At broadside (theta = 0), phi does not change the direction.'];
        end
        if cbFollowBeamCut.Value
            if isThetaCut
                spCutPhi.Value = S.beamPhiDeg;
            else
                spCutTheta.Value = S.beamThetaDeg;
            end
        end
        if ~any(A.active(:))
            cla(axCartesian); cla(axPolar);
            return;
        end
        if isThetaCut
            cutPhi = spCutPhi.Value;
            xlabelText = sprintf('Theta θ (φ = %.1f°)',cutPhi);
            plotTitle = sprintf('Theta cut (φ = %.1f°)',cutPhi);
        else
            cutTheta = spCutTheta.Value;
            xlabelText = sprintf('Phi φ (θ = %.1f°)',cutTheta);
            if cutTheta == 0
                plotTitle = 'Phi cut (θ = 0°; φ has no effect)';
            else
                plotTitle = sprintf('Phi cut (θ = %.1f°; steer φ = %.1f°)', ...
                    cutTheta,S.beamPhiDeg);
            end
        end
        [angle,theta,phi] = selectedCutDirections();
        actual = abs(reflectarrayField(A,phaseDeg,theta,phi));
        ideal = abs(reflectarrayField(A,A.idealPhaseDeg,theta,phi));
        showIdeal = S.manualEnabled || ...
            abs(A.phaseReferenceGHz-S.opFreqGHz) > 1e-12;
        absoluteLevel = cbAbsoluteDbi.Value;
        if absoluteLevel
            if ~isfinite(directivity.factor)
                cla(axCartesian); cla(axPolar);
                axPolar.Visible = 'off'; axCartesian.Visible = 'on';
                title(axCartesian,'Modeled directivity undefined', ...
                    'Color',WHITE);
                return;
            end
            actualDB = 10*log10(max(directivity.factor*actual.^2,realmin));
            sampledPeak = max(patternData.fieldMagnitude(:));
            displayTop = 10*log10(max( ...
                directivity.factor*sampledPeak^2,realmin));
            displayTop = max(displayTop,max(actualDB(:)));
            if showIdeal
                if isempty(idealDirectivity)
                    idealDirectivity = reflectarrayDirectivity(A,A.idealPhaseDeg);
                end
                idealDB = 10*log10(max( ...
                    idealDirectivity.factor*ideal.^2,realmin));
                displayTop = max([displayTop;idealDB(:); ...
                    10*log10(idealDirectivity.factor)]);
            else
                idealDB = [];
            end
            plotFloor = displayTop-40;
            actualDB = max(actualDB,plotFloor);
            idealDB = max(idealDB,plotFloor);
            radialLabels = arrayfun(@(value)sprintf('%.1f',value), ...
                plotFloor+(0:10:40),'UniformOutput',false);
            levelLabel = 'Modeled directivity (dBi)';
        else
            reference = max([actual(:);ideal(:)]);
            if reference < 1e-8
                % A cut through a true null has no useful local peak to
                % normalize against. Show the display floor, not a false 0 dB.
                actualDB = -40*ones(size(actual));
                idealDB = -40*ones(size(ideal));
            else
                actualDB = max(-40,20*log10(max(actual,eps)/reference));
                idealDB = max(-40,20*log10(max(ideal,eps)/reference));
            end
            plotFloor = -40;
            displayTop = 0;
            radialLabels = {'-40','-30','-20','-10','0'};
            levelLabel = 'Relative field level (dB)';
        end
        if strcmp(ddPlotStyle.Value,'Cartesian')
            legend(axPolar,'off');
            cla(axPolar);
            axPolar.Visible = 'off'; axCartesian.Visible = 'on';
            cla(axCartesian); hold(axCartesian,'on');
            if showIdeal
                plot(axCartesian,angle,idealDB,'--','Color',MUTED, ...
                    'LineWidth',1.5, ...
                    'DisplayName','Retuned at operating frequency');
            end
            plot(axCartesian,angle,actualDB,'Color',BLUE, ...
                'LineWidth',2.2,'DisplayName','Assigned phase');
            xlabel(axCartesian,xlabelText);
            ylabel(axCartesian,levelLabel);
            ylim(axCartesian,[plotFloor displayTop+1]);
            xlim(axCartesian,[angle(1) angle(end)]);
            title(axCartesian,plotTitle,'Color',WHITE);
            if showIdeal, legend(axCartesian,'Location','northeast'); end
            grid(axCartesian,'on'); hold(axCartesian,'off');
        else
            legend(axCartesian,'off');
            cla(axCartesian);
            axCartesian.Visible = 'off'; axPolar.Visible = 'on';
            cla(axPolar); hold(axPolar,'on');
            if showIdeal
                polarplot(axPolar,deg2rad(angle),idealDB-plotFloor,'--', ...
                    'Color',MUTED,'LineWidth',1.5, ...
                    'DisplayName','Retuned at operating frequency');
            end
            polarplot(axPolar,deg2rad(angle),actualDB-plotFloor, ...
                'Color',BLUE,'LineWidth',2.2, ...
                'DisplayName','Assigned phase');
            axPolar.RLim = [0 40]; axPolar.RTick = 0:10:40;
            axPolar.RTickLabel = radialLabels;
            if isThetaCut
                axPolar.ThetaZeroLocation = 'top';
                axPolar.ThetaDir = 'clockwise';
                axPolar.ThetaLim = [-90 90];
            else
                axPolar.ThetaZeroLocation = 'right';
                axPolar.ThetaDir = 'counterclockwise';
                axPolar.ThetaLim = [0 360];
            end
            title(axPolar,plotTitle,'Color',WHITE);
            if showIdeal, legend(axPolar,'Location','northeastoutside'); end
            hold(axPolar,'off');
        end
    end

    function buildInlineGeometryControls()
        e = uigridlayout(geometryPanel,[3 1]);
        e.RowHeight = {shapePreviewHeight,159,'1x'};
        e.Padding = [0 0 0 0]; e.RowSpacing = 6;
        e.BackgroundColor = BG;

        shapePanel = section(e,'ARRAY SHAPE');
        shapePanel.Layout.Row = 1;
        editorControls.shapePanel = shapePanel;
        gs = uigridlayout(shapePanel,[1 4]);
        gs.RowHeight = {'1x'};
        gs.ColumnWidth = {'1x','1x','1x',35};
        gs.Padding = [4 4 4 4]; gs.ColumnSpacing = 4;
        gs.BackgroundColor = PANEL;
        for t = 1:3
            shapePreviewButtons(t) = uibutton(gs,'Text','', ...
                'IconAlignment','top','FontSize',11, ...
                'FontColor',WHITE,'Tag',sprintf('reflectarrayShapePreview%d',t), ...
                'ButtonPushedFcn',@(src,~)chooseShape(src.UserData));
            shapePreviewButtons(t).Layout.Column = t;
        end
        arrow = uibutton(gs,'Text','▼','FontSize',16, ...
            'FontColor',WHITE,'BackgroundColor',PANEL, ...
            'Tag','reflectarrayShapeGalleryToggle', ...
            'Tooltip','Show all six array shapes', ...
            'ButtonPushedFcn',@(~,~)toggleShapeGallery());
        arrow.Layout.Column = 4;

        lattice = section(e,'LATTICE'); lattice.Layout.Row = 2;
        gl = uigridlayout(lattice,[3 2]);
        gl.RowHeight = {35,35,35};
        gl.ColumnWidth = {145,'1x'};
        gl.Padding = [7 7 7 7]; gl.RowSpacing = 8;
        gl.BackgroundColor = PANEL;
        label(gl,'Grid angle (°)',1,1);
        editorControls.gridAngle = number(gl,S.gridAngleDeg,[1 179],1,1,2, ...
            @(src,~)setSetting('gridAngleDeg',src.Value));
        editorControls.gridAngle.Tag = 'reflectarrayGridAngle';
        label(gl,'Row stagger (λ)',2,1);
        editorControls.stagger = number(gl,S.rowStaggerLambda,[-10 10],0.05,2,2, ...
            @(src,~)setSetting('rowStaggerLambda',src.Value));
        editorControls.stagger.Tag = 'reflectarrayRowStagger';
        label(gl,'Stagger angle (°)',3,1);
        editorControls.staggerAngle = number(gl,S.staggerAngleDeg,[-89.99 89.99],5,3,2, ...
            @(src,~)setSetting('staggerAngleDeg',src.Value));
        editorControls.staggerAngle.Tag = 'reflectarrayStaggerAngle';
        editorControls.staggerAngle.Tooltip = ...
            'Alternate-row angle atan(row stagger / dy); changes row stagger.';

        cellsPanel = section(e,'CELL EDITING'); cellsPanel.Layout.Row = 3;
        ge = uigridlayout(cellsPanel,[6 4]);
        ge.ColumnWidth = {72,'1x',72,'1x'};
        ge.RowHeight = {29,31,31,31,30,'1x'};
        ge.Padding = [7 3 7 3]; ge.BackgroundColor = PANEL;
        label(ge,'dx (λ)',1,1);
        editorControls.moveX = number(ge,0.1,[-100 100],0.05,1,2,@(~,~)[]);
        editorControls.moveX.Tag = 'reflectarrayMoveX';
        label(ge,'dy (λ)',1,3);
        editorControls.moveY = number(ge,0,[-100 100],0.05,1,4,@(~,~)[]);
        editorControls.moveY.Tag = 'reflectarrayMoveY';
        actions = {'Move selected','Reset to lattice', ...
            'Delete selected','Clear all','Undo last deletion'};
        callbacks = {@moveSelected,@resetSelected,@deleteSelected, ...
            @clearAllCells,@restoreSelected};
        positions = [2 1 2;2 3 4;3 1 2;3 3 4;4 1 4];
        for t = 1:numel(actions)
            callback = callbacks{t};
            b = uibutton(ge,'Text',actions{t}, ...
                'BackgroundColor',NAVY,'FontColor',WHITE, ...
                'ButtonPushedFcn',@(~,~)callback());
            b.Layout.Row = positions(t,1);
            b.Layout.Column = positions(t,2:3);
            b.Tag = ['reflectarrayCellAction' num2str(t)];
            editorControls.cellActions(t) = b;
            if t == 5, editorControls.undoDelete = b; end
        end
        hint = uilabel(ge,'Text', ...
            'Click an empty Phase Map cell to restore it.', ...
            'FontColor',MUTED,'FontSize',10,'WordWrap','on');
        hint.Layout.Row = 5; hint.Layout.Column = [1 4];
    end

    function createShapeGallery()
        % An in-window flyout, like the phased-array ribbon gallery.
        shapeGallery = uipanel(fig,'Title','CHOOSE ARRAY SHAPE', ...
            'Position',[8 340 430 265], ...
            'BackgroundColor',PANEL,'ForegroundColor',WHITE, ...
            'FontWeight','bold','Visible','off', ...
            'Tag','reflectarrayShapeGallery');
        galleryGrid = uigridlayout(shapeGallery,[2 2]);
        galleryGrid.RowHeight = {28,'1x'};
        galleryGrid.ColumnWidth = {'1x',30};
        galleryGrid.Padding = [7 5 7 5];
        galleryGrid.RowSpacing = 4;
        galleryGrid.BackgroundColor = PANEL;
        helper = uilabel(galleryGrid,'Text','Select a surface outline', ...
            'FontColor',MUTED,'FontSize',11);
        helper.Layout.Row = 1; helper.Layout.Column = 1;
        closeButton = uibutton(galleryGrid,'Text','×','FontSize',17, ...
            'BackgroundColor',PANEL,'FontColor',WHITE, ...
            'Tag','reflectarrayShapeGalleryClose', ...
            'ButtonPushedFcn',@(~,~)hideShapeGallery());
        closeButton.Layout.Row = 1;
        closeButton.Layout.Column = 2;
        tiles = uigridlayout(galleryGrid,[2 3]);
        tiles.Layout.Row = 2;
        tiles.Layout.Column = [1 2];
        tiles.RowHeight = {'1x','1x'};
        tiles.ColumnWidth = {'1x','1x','1x'};
        tiles.Padding = [0 0 0 0];
        tiles.ColumnSpacing = 5; tiles.RowSpacing = 5;
        tiles.BackgroundColor = PANEL;
        names = {'Custom','Diamond','Hexagon','Octagon','Circle','Ellipse'};
        for t = 1:6
            name = names{t};
            tile = uibutton(tiles,'Text',name, ...
                'Icon',fullfile(iconRoot,['shape_' lower(name) '.png']), ...
                'IconAlignment','top','FontSize',11, ...
                'FontColor',WHITE,'BackgroundColor',PANEL, ...
                'Tag',['reflectarrayShape' name], ...
                'ButtonPushedFcn',@(src,~)chooseShape(src.Text));
            tile.Layout.Row = ceil(t/3);
            tile.Layout.Column = mod(t-1,3)+1;
        end
    end

    function rememberShape(name)
        shapeRecent(strcmp(shapeRecent,name)) = [];
        shapeRecent{end+1} = name;
        if numel(shapeRecent) > 3
            shapeRecent = shapeRecent(end-2:end);
        end
    end

    function refreshShapePreview()
        for t = 1:3
            name = shapeRecent{t};
            button = shapePreviewButtons(t);
            button.Text = name;
            button.Icon = fullfile(iconRoot,['shape_' lower(name) '.png']);
            button.UserData = name;
            if strcmp(name,S.shape)
                button.BackgroundColor = [0.17 0.43 0.77];
            else
                button.BackgroundColor = PANEL;
            end
        end
        names = {'Custom','Diamond','Hexagon','Octagon','Circle','Ellipse'};
        for t = 1:numel(names)
            tile = findall(shapeGallery,'Tag',['reflectarrayShape' names{t}]);
            if strcmp(names{t},S.shape)
                tile.BackgroundColor = [0.17 0.43 0.77];
            else
                tile.BackgroundColor = PANEL;
            end
        end
    end

    function toggleShapeGallery()
        if strcmp(shapeGallery.Visible,'on')
            hideShapeGallery();
            return;
        end
        positionShapeGallery();
        shapeGallery.Visible = 'on';
    end

    function positionShapeGallery()
        if ~isgraphics(shapeGallery), return; end
        width = 430; height = 265;
        figureWidth = fig.Position(3);
        figureHeight = fig.Position(4);
        shapeBottomFromTop = root.Padding(4) + root.RowHeight{1} + ...
            root.RowSpacing + leftShell.RowHeight{1} + ...
            leftShell.RowSpacing + shapePreviewHeight;
        x = max(8,min(root.Padding(1),figureWidth-width-8));
        y = max(8,min(figureHeight-height-8, ...
            figureHeight-shapeBottomFromTop-height));
        shapeGallery.Position = [x y width height];
    end

    function onFigureResize()
        positionShapeGallery();
    end

    function hideShapeGallery()
        shapeGallery.Visible = 'off';
    end

    function showControlSection(name)
        isDesign = strcmp(name,'Design');
        isGeometry = strcmp(name,'Array Geometry');
        isBeam = strcmp(name,'Beam');
        designPanel.Visible = localOnOff(isDesign);
        geometryPanel.Visible = localOnOff(isGeometry);
        beamPanel.Visible = localOnOff(isBeam);
        hideShapeGallery();
        bDesignSection.BackgroundColor = NAVY;
        bGeometrySection.BackgroundColor = NAVY;
        bBeamSection.BackgroundColor = NAVY;
        if isGeometry
            bGeometrySection.BackgroundColor = BG;
            showView('Geometry');
            syncEditorControls();
        elseif isDesign
            bDesignSection.BackgroundColor = BG;
        else
            bBeamSection.BackgroundColor = BG;
        end
    end

    function syncEditorControls()
        if ~isfield(editorControls,'shapePanel'), return; end
        editorControls.shapePanel.Title = sprintf('ARRAY SHAPE · %d ACTIVE',nnz(A.active));
        refreshShapePreview();
        editorControls.gridAngle.Value = S.gridAngleDeg;
        editorControls.stagger.Value = S.rowStaggerLambda;
        editorControls.staggerAngle.Value = S.staggerAngleDeg;
        editorControls.undoDelete.Enable = ...
            localOnOff(~isempty(lastDeletedCells));
        for actionIndex = 1:3
            editorControls.cellActions(actionIndex).Enable = ...
                localOnOff(~isempty(S.selectedCells));
        end
        editorControls.cellActions(4).Enable = localOnOff(any(A.active(:)));
    end

    function moveSelected()
        if isempty(S.selectedCells), return; end
        if editorControls.moveX.Value == 0 && ...
                editorControls.moveY.Value == 0, return; end
        index = sub2ind(size(A.active),S.selectedCells(:,1),S.selectedCells(:,2));
        S.xOffsetLambda(index) = S.xOffsetLambda(index)+editorControls.moveX.Value;
        S.yOffsetLambda(index) = S.yOffsetLambda(index)+editorControls.moveY.Value;
        markDirty(); refreshAll();
    end

    function resetSelected()
        if isempty(S.selectedCells), return; end
        index = sub2ind(size(A.active),S.selectedCells(:,1),S.selectedCells(:,2));
        if all(S.xOffsetLambda(index) == 0) && ...
                all(S.yOffsetLambda(index) == 0), return; end
        S.xOffsetLambda(index) = 0;
        S.yOffsetLambda(index) = 0;
        markDirty(); refreshAll();
    end

    function deleteSelected()
        if isempty(S.selectedCells), return; end
        lastDeletedCells = S.selectedCells;
        lastDeletedSelection = S.selectedCells;
        if isempty(S.activeMask), S.activeMask = A.active; end
        index = sub2ind(size(A.active),S.selectedCells(:,1),S.selectedCells(:,2));
        S.activeMask(index) = false;
        markDirty(); refreshAll(); syncEditorControls();
    end

    function clearAllCells()
        if ~any(A.active(:)), return; end
        [deletedRows,deletedCols] = find(A.active);
        lastDeletedCells = [deletedRows deletedCols];
        lastDeletedSelection = S.selectedCells;
        S.activeMask = false(S.M,S.N);
        markDirty(); refreshAll(); syncEditorControls();
    end

    function restoreSelected()
        if isempty(lastDeletedCells), return; end
        if isempty(S.activeMask), S.activeMask = A.active; end
        index = sub2ind([S.M S.N], ...
            lastDeletedCells(:,1),lastDeletedCells(:,2));
        S.activeMask(index) = true;
        S.selectedCells = lastDeletedSelection;
        if isempty(S.selectedCells)
            S.selected = zeros(0,2);
        else
            S.selected = S.selectedCells(1,:);
        end
        lastDeletedCells = zeros(0,2);
        lastDeletedSelection = zeros(0,2);
        markDirty(); refreshAll(); syncEditorControls();
    end

    function restoreCell(row,col)
        if row < 1 || row > S.M || col < 1 || col > S.N, return; end
        if A.active(row,col)
            selectCell(row,col);
            return;
        end
        if isempty(S.activeMask), S.activeMask = A.active; end
        S.activeMask(row,col) = true;
        S.selected = [row col];
        S.selectedCells = S.selected;
        markDirty(); refreshAll(); syncEditorControls();
    end

    function ok = saveProject()
        ok = false;
        destination = projectFile;
        if isempty(destination)
            [file,path] = uiputfile('*.mat','Save reflectarray project', ...
                'reflectarray_project.mat');
            if isequal(file,0), return; end
            destination = fullfile(path,file);
        end
        project = S; project.formatVersion = 5;
        try
            save(destination,'project');
            projectFile = destination;
            dirty = false;
            updateWindowTitle();
            status.Text = ['Saved ' destination];
            ok = true;
        catch err
            uialert(fig,err.message,'Save failed');
        end
    end

    function ok = continueAfterChanges()
        ok = true;
        if ~dirty, return; end
        choice = uiconfirm(fig,'Save changes to this reflectarray project?', ...
            'Unsaved changes','Options',{'Save','Discard','Cancel'}, ...
            'DefaultOption',3,'CancelOption',3);
        switch choice
            case 'Save', ok = saveProject();
            case 'Discard', ok = true;
            otherwise, ok = false;
        end
    end

    function newProject()
        if ~continueAfterChanges(), return; end
        S = defaults;
        shapeRecent = {'Diamond','Custom','Circle'};
        projectFile = '';
        lastDeletedCells = zeros(0,2);
        lastDeletedSelection = zeros(0,2);
        geometryViewInitialized = false;
        patternViewInitialized = false;
        showGeometryBeam.Value = false;
        cbShow3DCut.Value = false;
        ddMapMode.Value = 'Phase map';
        showMapMode();
        ddCut.Value = 'Theta cut';
        ddPlotStyle.Value = 'Polar';
        cbAbsoluteDbi.Value = false;
        spCutTheta.Value = S.beamThetaDeg;
        spCutPhi.Value = S.beamPhiDeg;
        cbFollowBeamCut.Value = true;
        spSweepPoints.Value = 21;
        showControlSection('Design');
        showView('Geometry');
        syncControls();
        refreshAll();
        centerSweep();
        dirty = false;
        updateWindowTitle();
        syncEditorControls();
        status.Text = 'New reflectarray project';
    end

    function closeDesigner()
        if ~continueAfterChanges(), return; end
        delete(fig);
    end

    function openProject()
        if ~continueAfterChanges(), return; end
        [file,path] = uigetfile('*.mat','Open reflectarray project');
        if isequal(file,0), return; end
        try
            loaded = load(fullfile(path,file),'project');
            if ~isfield(loaded,'project')
                error('Reflectarray:Project','Not a supported reflectarray project.');
            end
            candidate = reflectarrayUpgradeProject(loaded.project);
            reflectarrayValidateProject(candidate);
            previousS = S;
            previousFile = projectFile;
            previousDirty = dirty;
            previousRecent = shapeRecent;
            S = candidate;
            projectFile = fullfile(path,file);
            try
                syncControls();
                refreshAll(true);
                centerSweep();
            catch renderError
                S = previousS;
                projectFile = previousFile;
                dirty = previousDirty;
                shapeRecent = previousRecent;
                syncControls();
                refreshAll();
                dirty = previousDirty;
                updateWindowTitle();
                rethrow(renderError);
            end
            dirty = false;
            lastDeletedCells = zeros(0,2);
            lastDeletedSelection = zeros(0,2);
            updateWindowTitle();
            syncEditorControls();
            status.Text = ['Opened ' fullfile(path,file)];
        catch err
            uialert(fig,err.message,'Open failed');
        end
    end

    function exportCells(kind)
        if ~any(A.active(:))
            uialert(fig,'There are no active cells to export.','Export');
            return;
        end
        if strcmp(kind,'csv')
            [file,path] = uiputfile('*.csv','Export reflector cells', ...
                'reflectarray_cells.csv');
        else
            [file,path] = uiputfile('*.tsv','Export CST equivalent aperture', ...
                'reflectarray_equivalent_aperture.tsv');
        end
        if isequal(file,0), return; end
        destination = fullfile(path,file);
        try
            count = reflectarrayWriteCells(A,phaseDeg,S,destination,kind);
            status.Text = sprintf('Exported %d cells to %s',count,destination);
        catch err
            uialert(fig,err.message,'Export failed');
        end
    end

    function syncControls()
        spM.Value = S.M; spN.Value = S.N;
        spDx.Value = S.dxLambda; spDy.Value = S.dyLambda;
        rememberShape(S.shape);
        cbFixedPhase.Value = ~S.retunePhase;
        syncFrequencyControls();
        spFx.Value = S.feedMM(1); spFy.Value = S.feedMM(2);
        spFz.Value = S.feedMM(3); spQ.Value = S.feedQ;
        syncIlluminationControls();
        ddConvention.Value = S.angleConvention;
        ddTaper.Value = S.taper;
        setConvention();
    end
end

function value = localOnOff(tf)
if tf, value = 'on'; else, value = 'off'; end
end
