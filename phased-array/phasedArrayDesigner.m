function phasedArrayDesigner
%PHASEDARRAYDESIGNER  Interactive planar phased-array pattern designer.
%
%   Author : Muhammed Said Gökdöl  <saidgokdol@gmail.com>
%
%   Three element-placement modes:
%     1. Uniform grid        - every cell of an MxN lattice is populated
%     2. Sparse (click)      - click cells on the MxN lattice to add/remove
%     3. Sub-position (9-pt) - click a cell, then pick one of 9 offsets
%                              (including the 45 deg diagonals) inside it.
%                              A cell can hold several, one per offset,
%                              so repeated clicks DENSIFY the array
%                              rather than moving one element around.
%
%   Lattice shaping: rectangular (Custom) or trimmed to a Diamond /
%   Hexagon / Octagon / Circle / Ellipse boundary, with an adjustable
%   grid skew angle and per-row stagger for brick lattices.
%
%   Element factor:
%     - Isotropic
%     - cos^q(theta)                q is the POWER-pattern exponent
%     - Short dipole (z-axis)       nulls on the array normal
%     - Dipole (linear pol)         legacy in-plane dipole
%     - Patch (cos^q x lin pol)     built-in analytic linear-polarized patch
%     - Custom (formula)            user-typed E_theta / E_phi expressions
%     - Imported (CST far-field)    a CST Farfield ASCII export used
%                                   directly as the element pattern,
%                                   carrying its own absolute gain and
%                                   radiation efficiency
%   Each element carries its OWN rotation angle, so sequential-rotation
%   (circular polarisation) arrays can be built; rotation drives feed
%   phase only if the "rot angle -> feed phase" box is ticked.
%
%   Amplitude taper: Uniform / Hamming / Hanning / Chebyshev / Taylor /
%   Binomial, or fully manual per-element entry in the table. Chebyshev
%   and Taylor need the Signal Processing Toolbox; the rest do not.
%
%   Readouts: 3D pattern, polar azimuth/elevation cuts, a direction-cosine
%   U cut, and a principal cut (theta or phi), array factor
%   / element factor / total / axial ratio, absolute directivity or gain
%   in dBi, HPBW, FNBW, sidelobe level, a geometry-agnostic grating-lobe
%   search, scan loss vs steering angle, max scan angle, beam squint, and
%   aperture efficiency.
%
%   Import/export: CST far-field ASCII in; element table as CSV and a
%   CST-array .tsv out; full design save/load as .mat.
%
%   All lengths are in wavelengths (lambda = 1). The Beam tab can show
%   directions as theta/phi or azimuth/elevation; elevation = 90 - theta
%   for nonnegative theta. Internal calculations and CST files retain
%   spherical theta from +z and phi from +x in the xy-plane.
%   Elements lie in the z = 0 plane and radiate toward +z, so an imported
%   pattern must have its boresight at theta = 0.

% ------------------------------------------------------------------ state
S.lambda   = 1;
S.k        = 2*pi/S.lambda;
S.M        = 8;          % rows    (y direction)
S.N        = 8;          % columns (x direction)
S.dx       = 0.5;        % column spacing, wavelengths
S.dy       = 0.5;        % row spacing, wavelengths
S.subOff   = 0.25;       % sub-position offset, absolute wavelengths.
                         % Deliberately smaller than the default dx/dy of
                         % 0.5: at 0.5 the nine candidates land exactly on
                         % the NEIGHBOURING cell centres, so the markers
                         % sat on top of the adjacent elements and picking
                         % a corner put two elements in the same place.
                         % A quarter wavelength keeps all nine inside the
                         % cell they belong to. This is an absolute
                         % length, not a fraction of dx, so re-check it if
                         % you move to a very different lattice pitch.
S.theta_s  = 0;          % steering theta, deg
S.phi_s    = 0;          % steering phi, deg
S.angleConvention = 'Theta / phi'; % UI convention; calculations/CST stay spherical
% Cut-plane azimuth for the theta cut, kept separate from the steering
% azimuth. They coincide by default (cutPhiFollow), which is the
% behaviour this app always had -- but they are different questions:
% phi_s is where the BEAM points, cutPhi is which plane you LOOK at.
% Tying them together makes an off-broadside design impossible to
% inspect properly, since you cannot view the H-plane of a beam steered
% in the E-plane without re-steering the beam and changing the pattern
% you were trying to measure. Harmless at broadside (theta_s = 0 makes
% us/vs zero for any phi_s), which is why it went unnoticed.
S.cutPhi       = 0;      % cut-plane azimuth used when not following
S.cutPhiFollow = true;   % cut plane tracks phi_s
% ---- Band sweep (Analysis > Band sweep...) ------------------------------
% Its window, reused while open. Declared here, in the parent scope, so
% the menu item, the window's buttons and the sweep (all siblings) share
% it. The sweep's settings -- band, points, phase-set frequency, phase
% shifters vs true-time delay -- live in that window's own fields: they
% describe an analysis OF the design, not the design, so a saved design
% does not carry them.
bandWin = [];
S.mode     = 'Uniform grid';
S.pending  = [];         % cell awaiting a sub-position pick, [row col]
S.el       = zeros(0,5); % [x  y  amp  phase_deg  rot_deg]
S.elRC     = zeros(0,2); % [row col] each element belongs to, for tapering
S.taper    = 'Uniform';
S.sll      = 30;         % sidelobe level for Chebyshev/Taylor, dB
% ---- result cards and design targets ------------------------------------
% Beam figures measured on the principal cut, kept for the result cards:
% they used to exist only as locals of the cut plot and inside its title.
% measured separates "no sidelobe on this cut" (sll NaN, measured) from
% "nothing measured" (axial-ratio cut, no radiation, out of date); why
% says which, for the cards. Reset by clearCutMetrics.
S.cutMetrics = struct('hpbw',NaN,'fnbw',NaN,'sll',NaN, ...
    'measured',false,'plane','','offBeam',false,'why','');
% Targets the cards are checked against, saved with the design. NaN (and
% false) = no target. minGain applies to whatever the first card shows:
% directivity, realized gain of an imported element, or array gain in
% unit-cell mode. maxSLL is relative to the peak, so it is negative.
S.targets = struct('minGain',NaN,'maxHPBW',NaN,'maxSLL',NaN,'noGrating',false);
% ---- pinned reference trace on the principal cut -------------------------
% cutView describes the cut as last drawn (gain or axial ratio, which
% plane, angle convention, absolute or relative levels); a reference keeps
% the cutView it was pinned under and is drawn only while they agree --
% the same curve over a different x axis or level basis would be a false
% comparison. refCut is [] when nothing is pinned. Not saved with the
% design: it is there to compare designs.
S.cutView = struct('kind','','isPhi',false,'plane',0,'conv','','abs',false);
S.refCut = [];
S.cstOverlay = [];
S.coverageCache = [];
S.coverageWin = gobjects(0);
S.pattern2DBasis = [];
pattern2DWin = gobjects(0);
S.efType   = 'Isotropic';
S.customFormula = 'ct.^2';   % user-editable E_theta(th,ph) formula
S.customFormulaPh = '0';     % user-editable E_phi(th,ph) formula (0 = linear pol)
S.efQ      = 1.5;        % exponent q in cos^q(theta)
S.efBeamAz = 65;        % analytic element azimuth half-power beamwidth (deg)
S.efBeamEl = 65;        % analytic element elevation half-power beamwidth (deg)
S.dualFeedCP = false; % Legacy configuration field; built-in patch is linear.
S.impFF = [];      % Imported (CST far-field) element type: the loaded
                   % far-field, [] until loaded via "Load far-field
                   % pattern..." -- struct of 4 scatteredInterpolant
                   % objects (Re/Im of E_theta, Re/Im of E_phi) built
                   % once at load time, NOT re-parsed on every
                   % elementFactor call. Real/imaginary interpolation
                   % (not magnitude/phase) deliberately, to avoid
                   % phase-wrap artifacts at the +-180 deg boundary.
                   % Persisted with Save/Load config; interpolants are
                   % restored and their power metadata is refreshed.
S.impFFName = '';  % loaded filename, for the button label
% ---- CST path: imported pattern's frequency, macro export scheme --------
S.impFFGHz = NaN;      % frequency the imported pattern was exported at, as
                       % its file name or header states it; NaN = unstated
S.impFFFreqFrom = '';  % where that came from ('file name', ...), for display
S.cstScheme = 'single';% feed scheme of the last CST macro export (session)
% ---- undo / redo: design snapshots ---------------------------------------
% The S fields that make up the DESIGN -- what Undo restores. View settings
% (cut plane, 3D display, angle convention, frequency unit, Auto) are left
% out on purpose, so undoing a taper change never also moves the cut plane
% the user has since switched to. A new saved design field belongs here and
% in syncDesignControls. Parent-scope variables, not S fields, because the
% stacks must survive the S = beforeState rollback in loadConfig.
UNDO_FIELDS = {'M','N','dx','dy','subOff','mode','el','elRC', ...
    'arrayShape','gridAngle','stagger','taper','sll', ...
    'efType','efQ','efBeamAz','efBeamEl','customFormula','customFormulaPh','dualFeedCP', ...
    'impFF','impFFName','impFFGHz','impFFFreqFrom','impNeedsPattern', ...
    'impUnitCell','impTotEffPct','theta_s','phi_s', ...
    'freqGHz','freqOpGHz','squintMode','retunePhase', ...
    'seqPhase','seqBlockM','seqBlockN','portMap','portMapSet'};
UNDO_MAX = 50;       % levels kept; the oldest step is dropped beyond this
undoStack = {};      % {struct('design',snapshot,'what',name)}, newest last
redoStack = {};
undoBase = [];       % the design as of the last checkpoint
undoHold = 0;        % >0 while a restore or a config load is in progress
undoCount = 0;       % steps recorded so far, so a caller can tell if its
                     % action made one (and only then offer "Undo with ⌘Z")
lastUndoUI = [];     % the last Undo/Redo request and where it came from,
                     % to spot one keystroke arriving twice (undoFromUI)
mUndo = []; mRedo = [];  % Edit > Undo / Redo, built with the menu bar
% ---- the open design file: title, unsaved marker, recent files ----------
% Parent-scope variables for the same reason as the undo stacks. "Unsaved"
% means the DESIGN (UNDO_FIELDS) differs from the one last saved or
% opened: view settings travel in the file too, but, as with Undo, a
% changed cut plane is not a change to the design worth a prompt.
APP_NAME = 'Planar Phased Array Designer';
APP_VERSION = '2.1';  % written into every saved design (cfg.appVersion)
docFile = '';         % full path of the open design; '' = not saved yet
docSaved = [];        % that design as saved or opened (a designSnapshot)
docNextStep = [];     % the file a load leaves behind, attached to the load's
                      % undo step: undoing an Open returns to that file, so
                      % Save cannot write the old design over the new file
closeAsking = false;  % the unsaved-changes question is on screen
% Recent files and the last folder live in MATLAB preferences. The group
% can be redirected through the PAD_PREF_GROUP environment variable, so
% the test suites keep their temporary files out of the user's list.
PREF_GROUP = getenv('PAD_PREF_GROUP');
if isempty(PREF_GROUP), PREF_GROUP = 'phasedArrayDesigner'; end
RECENT_MAX = 8;
mRecent = [];         % File > Open Recent, rebuilt by refreshRecentMenu
S.portMapRefresh = {};      % redraw handles, one per OPEN phase map. A
                            % single handle meant a second window
                            % replaced the first (which then went stale
                            % silently) and closing either stopped both.
S.impNeedsPattern = false;  % set when a loaded design named an imported
                            % element the config did not carry: the plots
                            % must say so rather than quietly substituting
% Efficiency below which an imported pattern is treated as a GAIN export
% (loss included) rather than a DIRECTIVITY export. Declared here at the
% top level, not inside a nested function, so loadImportedFF and the
% plotting code genuinely share it -- a variable is only common to two
% nested functions if it also appears in their shared ancestor's own
% code. 0.99 rather than something tighter because integrating a coarse
% CST export grid (5 deg steps, then interpolated) lands ~0.1-0.5% off
% unity even for a perfectly lossless pattern; a genuine gain export
% sits far below this (0.60 on the confirmed real file), and an antenna
% actually at 99% efficiency differs from lossless by 0.04 dB anyway.
EFF_TOL = 0.99;
% Manual coordinates bypass the lattice spinners. This bound is far beyond
% any generated 64x64, 5-lambda lattice while keeping plotting and phase
% arithmetic away from effectively infinite positions.
MAX_POS_LAMBDA = 1e6;
S.fullSphere = false;  % false = upper hemisphere (theta 0-90), the
                   % ground-plane-backed convention this app defaults to;
                   % true = full sphere (theta 0-180). See cbFullSphere.
S.seqBlockM = 2;         % sequential-rotation subarray block size, rows
S.seqBlockN = 2;         % sequential-rotation subarray block size, cols
                          % -- "Apply sequential rotation to all" starts
                          % each tile at its physical BOTTOM-LEFT cell
                          % (0 deg reference) and repeats the SAME
                          % 0/step/2*step/... snake pattern in every
                          % blockM x blockN tile across the array,
                          % rather than incrementing continuously
                          % cell-to-cell. Sets element ORIENTATION only
                          % (S.el(:,5)); whether that also drives feed
                          % phase is a SEPARATE, independent choice --
                          % see S.seqPhase below. Setting blockM=M,
                          % blockN=N makes it one continuously-
                          % incrementing pattern across the whole array.
% ---- physical feed layout of the repeating unit -------------------------
% WHICH physical port carries H, and WHICH SIDE of the element each feed
% sits on, for the repeating unit that tiles the array.
%
% This is a property of the CST model, not of the element positions. The
% tool used to INFER it from column parity and rotation angle, and a
% hand-drawn layout of a real 3-element row proved that inference exactly
% inverted: measured H = 2, 3, 6 where the rule produced H = 1, 4, 5.
% Worse, the sides did not tile at all -- elements 1 and 3 share a block
% column and share H-on-the-even-port, yet element 1 feeds from the LEFT
% and element 3 from the RIGHT. No rule over positions and rotations can
% produce that, because the layout mirrors elements and the tool stores
% no mirror flag. So it is declared, not derived.
%
% hOdd   true  = H is the ODD port (2n-1), false = H is the EVEN port (2n)
% hSide  1 = left,   2 = right
% vSide  1 = bottom, 2 = top
% Internally the map follows S.elRC's bottom-up row numbering. The editor
% presents Block_Row from TOP to BOTTOM and converts automatically.
%
% The default below matches the antenna's declared 2x2 CST topology:
% H port identity is checkerboarded, H feed side is L/R by physical column,
% and V feed side is T/B by physical row. Opposite-side 180 deg corrections
% are derived from the declared physical sides, not guessed from rotation.
%
% The unit is nRow x nCol and is USER-SIZABLE. A fixed 2x2 could not
% express the measured layout at all: elements 1 and 3 of a real row sit
% in the same 2x2 column yet feed from opposite edges, so no entry in a
% 2-periodic table can produce that row. Three elements do not establish
% the true period either, so the period is asked for rather than assumed.
% DEFAULT PHYSICAL 2x2 FEED LAYOUT, viewed from the front:
%   top-left:     H odd,  H side left,  V side top
%   top-right:    H even, H side right, V side top
%   bottom-left:  H even, H side left,  V side bottom
%   bottom-right: H odd,  H side right, V side bottom
% S.elRC row 1 is the physical BOTTOM row, so this structure is stored
% bottom-up internally. The editor displays Block_Row 1 at the TOP.
DEFAULT_PORT_MAP = struct( ...
    'hOdd',  [false true; true false], ...
    'hSide', [1 2; 1 2], ...
    'vSide', [1 1; 2 2]);
S.portMap = DEFAULT_PORT_MAP;
S.efFallbackWarned = '';  % last formula pair already warned about
S.efFallbackWhich = [false false];  % which of E_theta / E_phi failed
S.efFallback = false;     % set when a Custom formula failed on the
                          % REAL grid and was silently substituted
S.portMapSet = false;     % true once the user has edited it, so the map
                          % can distinguish "declared" from "still default"
S.seqPhase = false;       % OFF (default): rotation is purely element
                          % orientation/polarization, never touches feed
                          % phase -- coherently summing rotation-diverse
                          % LINEAR elements (Dipole/Patch) with no phase
                          % compensation cancels at broadside (verified:
                          % four 90-deg-apart orientations in equal
                          % numbers sum to exactly zero field). ON: adds
                          % each element's rotation angle into its feed
                          % phase, which is what actual sequential-
                          % rotation CP arrays do -- physical rotation
                          % and feed-phase rotation together make the
                          % desired circular-polarization sense add
                          % constructively instead of cancelling.
S.largeWarnShown = false; % has the large-array notice already fired?

% ---- CST-style array shape / grid angle (Uniform grid mode only) ----
% Shape trims WHICH cells of the MxN index grid get an element (using
% the orthogonal dx,dy bounding box as the boundary reference, same as
% CST's "Shape" dropdown: Custom/Diamond/Hexagon/Octagon/Circle/Ellipse).
% Grid angle then decides WHERE those surviving cells physically sit:
% row-direction lattice vector = dy*[cosd(gridAngle), sind(gridAngle)],
% col-direction lattice vector = dx*[1,0] -- gridAngle=90 (the CST
% default) reduces exactly to the plain orthogonal grid.
S.arrayShape = 'Custom';   % Custom/Diamond/Hexagon/Octagon/Circle/Ellipse
S.gridAngle  = 90;         % deg, angle between the two lattice vectors
S.stagger    = 0;          % lambda; EVEN rows (2,4,6,...) shift by this
                            % along x, ODD rows (1,3,5,...) stay at x=0
                            % offset -- period-2 "brick wall" lattice,
                            % independent of gridAngle (which instead
                            % skews EVERY row progressively). row1==row3,
                            % row2==row4, etc. in x-alignment.
S.freqGHz = 10;            % Geometry reference: spacing in wavelengths at this frequency.
S.freqOpGHz = 10;          % Steering and radiation are evaluated at this frequency.
S.squintMode = true;       % Physical-frequency scaling is always active: propagation
                           % is evaluated at freqOpGHz while element positions stay
                           % referenced to freqGHz. Not user-switchable; kept as a
                           % field because configs carry it.
S.retunePhase = true;      % TRUE  = the phase shifters are re-set whenever the
                           %         operating frequency moves, so the generated progressive
                           %         steering term remains aimed at the commanded angle.
                           %         This is the "Beam squint" box CLEAR.
                           % FALSE = the phase values are frozen at the DESIGN
                           %         frequency. Propagation still happens at the
                           %         operating frequency, so the beam walks off the
                           %         commanded angle -- classic beam squint:
                           %           sin(theta_beam) = sin(theta_s) * f_design/f_op
                           % The checkbox is the inverse of this field (ticking
                           % "Beam squint" means NOT retuning), and the inversion is
                           % confined to assignSquint() and syncFrequencyMode() so
                           % there is exactly one place each way to get it wrong.
S.cutMode = 'Theta cut (fixed phi_s)';   % or 'Phi cut (fixed theta)'
S.cutFixedTheta = 0;        % theta held fixed for a Phi cut, deg.
                            % 0 by deliberate choice: a neutral starting
                            % point the user dials up from. Note the cut
                            % IS flat at 0 -- the cone collapses to the
                            % +z point, so every phi is the same
                            % direction (measured: 0.000 dB range for
                            % every element factor, including a
                            % sequential-rotation CP array). That is
                            % geometry, not a fault.
S.cutRaw = []; S.radiationValid=false;
S.cutX = []; S.cutDb = [];  % last computed cut curve (Total, dB), for
                            % the click-to-inspect readout on axCut
S.glTheta = NaN; S.glPhi = NaN; S.glRelDb = NaN;  % last grating-lobe search result
S.glInterior = false;    % ...and whether that lobe's PEAK is inside visible
                         % space, which is what separates a real grating lobe
                         % from the horizon cutting through one's skirt
S.DpkTot = NaN;   % last computed TOTAL directivity, for the Metrics readout
S.taperFailed = false;     % taperVec fell back to Uniform this call
% Treat an imported pattern as an infinite-array EMBEDDED element (from a
% unit-cell / periodic-boundary simulation) rather than a free-space one.
% Changes how absolute levels are reported -- see the note in refreshInfo.
% Display unit for the two frequency spinners. S.freqGHz/S.freqOpGHz
% stay in GHz internally no matter what this says -- every physical
% formula in the app (lambdaMM, freqRatio) reads them directly, so a
% unit that leaked into storage would silently rescale the design.
S.freqUnit = 'GHz';
S.impUnitCell = false;
% Total efficiency you supply for an imported unit-cell pattern, percent.
% For such a file the pattern's own integral cannot yield it (see the note
% in loadCstFarfieldASCII), so it is typed in from the solver instead and
% used ONLY to convert the file's realized-gain peak into directivity.
S.impTotEffPct = 100;
S.shapeFellBack = false;   % shapeMaskFn kept the full grid because the
                           % chosen shape would have emptied the array
% Scan-loss window: the one uifigure every Scan loss run draws into, so
% repeated runs replace the curve instead of stacking up windows.
S.scanLossWin = gobjects(0);

% ------------------------------------------------------------------- figure
fig = uifigure('Name',APP_NAME,'Position',[60 10 1480 1020], ...
    'KeyPressFcn',@onKey,'WindowState','maximized');
root = uigridlayout(fig,[3 1]);
% Tall enough for the CONFIGURATION group's FOUR rows plus the tab
% strip above it. The fourth row holds the Beam squint control; keeping
% the older 136 px height can clip that row entirely; 174 leaves some
% platform/theme headroom instead of fitting the four rows edge-to-edge.
% The third row is the status bar, one line of small text.
root.RowHeight = {174,'1x',16};
% Less padding below than above: the status line is the bottom edge, and
% a full 8 px under it would take plot height for nothing.
root.Padding = [8 3 8 8];
root.RowSpacing = 8;

% ============================================ TOP: ribbon with task tabs
% Array / Element / Beam / View are TOP-LEVEL tabs, the way DESIGN and
% ANALYSIS are in MathWorks' own designer apps. Selecting one does two
% things at once: it swaps the ribbon to that task's action groups, and
% it swaps the left panel to that task's parameters. So the whole window
% follows a single choice, instead of the tabs and the ribbon being two
% separate navigations the user has to keep in sync themselves.
%
% Built from ordinary uipanel/uibutton rather than the internal
% matlab.ui.internal.toolstrip package MathWorks' own apps use: that
% package gives the exact native look but is undocumented and free to
% change between releases.
% The tab group holds only TASK-specific groups. Frequency, save/load
% and the exporters sit in a fixed column beside it, outside the tabs
% entirely -- they apply to every task, and the frequency controls in
% particular hold state, so duplicating them per tab would mean four
% spinners disagreeing about one number. Keeping them out here also
% closes the dead band that used to run across the middle of every
% ribbon.
% Declared HERE, in the parent scope, even though commonGroups is what
% actually creates them. MATLAB shares a variable between nested
% functions only if the parent function also uses that name -- a handle
% assigned solely inside commonGroups would be invisible to
% applyFreqUnit, assignFreq and loadConfig, which are siblings, not
% children. The Code Analyzer flags this as "variable might be used
% before it is defined"; it is not a style warning, it is the bug.
spFreq = []; spFreqOp = []; ddFreqUnit = [];
lblFreqOp = []; lblLambdaMM = []; lblSpacingMM = [];
% cbSquint belongs to the same set for the same reason: commonGroups
% creates it, syncFrequencyMode (a sibling, reached from loadConfig)
% writes it. Left out, loadConfig silently stopped restoring the squint
% setting -- and only because syncFrequencyMode uses set() rather than a
% dot assignment did that surface as an error at all instead of writing
% into a throwaway struct.
cbSquint = [];
% Theme bookkeeping, declared here for the same sharing reason. Controls
% left on their theme defaults repaint themselves when the theme
% changes; anything given an explicit colour from pal() does not, so the
% handles carrying one are collected as they are built -- ribbonPanel
% adds its caption to mutedLbls, sep() its rule to sepList, matchTheme()
% each popup window to themedPopups -- and applyTheme walks the lists.
% curTask remembers the open tab so the highlight can be repainted
% without switching tasks.
mutedLbls = gobjects(0);
themedPopups = gobjects(0);
curTask = 'Array';
% Status bar and busy state, shared for the same reason: setStatus,
% beginBusy and the analyses that call them are all siblings. statusMsg
% and statusTone hold the message WITHOUT its tone glyph, so a busy
% section can tell whether the line still shows its own message, and
% applyTheme can repaint the tone colour. busyDepth counts live
% beginBusy tokens, so a busy section inside another (a compute run by
% an analysis) does not hand the pointer back while the outer one runs.
statusBar = []; statusMsg = 'Ready'; statusTone = '';
busyDepth = 0;
% Compute shows the busy pointer only when the wait is noticeable: from
% this many elements, or when the previous compute took half a second
% (an imported or custom element makes a small array slow too). The
% flush that puts the pointer on screen is itself a redraw, and the
% default 8x8 should not pay for it on every spinner click.
BUSY_MIN_ELEMENTS = 400;
lastComputeSec = 0;
defaultDesign = [];
defaultView = [];   % start-up view values, restored by File > New
layoutView = struct('indices',true,'taper',true, ...
    'localAxes',false,'annotation',false);

% The tab strip is drawn from plain buttons rather than a uitabgroup.
% A uitabgroup paints its own framed rectangle around whatever it
% contains, and the always-visible frequency/file/export column has to
% live OUTSIDE it (those controls hold state, so they cannot be
% duplicated per tab) -- which left the tab frame stopping halfway
% across the window with bare background beside it. Driving the ribbon
% body by hand lets it span the full width the way a real toolstrip
% does, with the task groups swapping underneath.
ribbonWrap = uigridlayout(root,[2 1]);
ribbonWrap.RowHeight = {32,'1x'};
ribbonWrap.Padding = [0 0 0 0];
ribbonWrap.RowSpacing = 0;

designTaskNames = {'Array','Element','Beam','View'};
resultTaskNames = {'2D Pattern','Phase Map','Scan Loss','Max Scan', ...
    'Band Sweep','Coverage','CST Compare'};
resultOwners = [4 3 3 3 3 3 4];
taskNames = [designTaskNames resultTaskNames];
resultTaskWidths = {116,108,104,104,114,106,120};
[tabStrip,tabBtns,closeBtns,ribbonHelpIcon] = buildTaskStrip(ribbonWrap, ...
    designTaskNames,resultTaskNames, ...
    @(s,e)selectTask(s.Text),@(s,e)closeResult(s.Tag(6:end)));
ribbonHelpIcon.ImageClickedFcn = @(s,e)showHelp('quick');
lastDesignTask = 1;
resultPanels = cell(1,numel(resultTaskNames));
resultWorkspace = [];

ribbonBody = uipanel(ribbonWrap);
gBody = uigridlayout(ribbonBody,[1 2]);
gBody.ColumnWidth = {'1x',770};   % commonGroups' three widths + gaps
gBody.Padding = [4 4 4 4];
gBody.ColumnSpacing = 6;

% One panel per task, all stacked in the same cell with only the active
% one visible -- the same trick the left-hand sections use.
taskStack = uigridlayout(gBody,[1 1]);
taskStack.Padding = [0 0 0 0];

rbArray = uipanel(taskStack,'BorderType','none');
rbArray.Layout.Row = 1; rbArray.Layout.Column = 1;
gArr = ribbonRow(rbArray,{250,270,130});

gShape = ribbonPanel(gArr,'ARRAY SHAPE');
shapeNames = {'Custom','Diamond','Hexagon','Octagon','Circle','Ellipse'};
shapePreviewOlderModel = 'Diamond';
shapePreviewPreviousModel = 'Circle';
[shapePreviewOlder,shapePreviewPrevious,shapePreviewCurrent,btnShapeGallery] = ...
    buildShapeRibbon(gShape,@ribbonIconFile,@pickShape);
shapeGalleryPopup = gobjects(0);
shapeBtns = gobjects(0);

gLayoutDisplay = ribbonPanel(gArr,'GEOMETRY DISPLAY');
[cbLayoutIndex,cbLayoutTaper,cbLayoutAxes, ...
    cbLayoutAnnotation] = buildGeometryDisplayControls( ...
    gLayoutDisplay,layoutView,@updateLayoutOption);

gElemOps = ribbonPanel(gArr,'ELEMENTS');
gElemOps.ColumnWidth = {'1x','1x'};
uibutton(gElemOps,'Text','Delete','Icon',ribbonIconFile('trash'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnDelete', ...
    'ButtonPushedFcn',@(s,e)deleteSelected(), ...
    'Tooltip','Delete. Apply this setting or action to the current design.');
uibutton(gElemOps,'Text','Clear','Icon',ribbonIconFile('cross'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnClearAll', ...
    'ButtonPushedFcn',@(s,e)clearAll(), ...
    'Tooltip','Clear. Apply this setting or action to the current design.');


% ---- Element ribbon ----------------------------------------------------
rbElem = uipanel(taskStack,'BorderType','none');
rbElem.Layout.Row = 1; rbElem.Layout.Column = 1;
gEl = ribbonRow(rbElem,{150,350});
gImp = ribbonPanel(gEl,'IMPORT');
gImp.ColumnWidth = {'1x'};
bImportRibbon = uibutton(gImp,'Text','Far-field (.txt)', ...
    'Icon',ribbonIconFile('import'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnImportRibbon', ...
    'Tooltip','Load a CST Farfield -> Export -> ASCII .txt export', ...
    'ButtonPushedFcn',@(s,e)loadImportedFF());
gPatterns = ribbonPanel(gEl,'ELEMENT');
gPatterns.ColumnWidth = {'1x','1x',30};
elementPreviewPreviousModel = 'Gaussian';
elementPreviewPrevious = uibutton(gPatterns,'Text','Gaussian', ...
    'Icon',ribbonIconFile('efGaussian'),'IconAlignment','top', ...
    'FontSize',11,'Tag','elementPreviewPrevious', ...
    'Tooltip','Select the previous element pattern', ...
    'ButtonPushedFcn',@(s,e)pickElementPreview(s));
elementPreviewCurrent = uibutton(gPatterns,'Text','Isotropic', ...
    'Icon',ribbonIconFile('efIsotropic'),'IconAlignment','top', ...
    'FontSize',11,'Tag','elementPreviewCurrent', ...
    'Tooltip','Current element pattern; click to open the gallery', ...
    'ButtonPushedFcn',@(s,e)pickElementPreview(s));
uibutton(gPatterns,'Text','▼','FontSize',12,'Tag','btnElementGallery', ...
    'Tooltip','Open the antenna element gallery below this ribbon', ...
    'ButtonPushedFcn',@(s,e)showElementGallery());
elementGalleryPopup = gobjects(0);
elementGalleryBtns = gobjects(0);
elementGalleryModels = {};

% ---- Beam ribbon -------------------------------------------------------
rbBeam = uipanel(taskStack,'BorderType','none');
rbBeam.Layout.Row = 1; rbBeam.Layout.Column = 1;
gBm = ribbonRow(rbBeam,{600});
gAn = ribbonPanel(gBm,'ANALYZE');
gAn.ColumnWidth = {'1x','1x','1x','1x'};
uibutton(gAn,'Text','Scan loss','Icon',ribbonIconFile('wave'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnScanLoss', ...
    'ButtonPushedFcn',@(s,e)showScanLoss(), ...
    'Tooltip','Plot scan loss versus steering angle (°) in the current plane.');
uibutton(gAn,'Text','Max scan','Icon',ribbonIconFile('peak'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnMaxScan', ...
    'ButtonPushedFcn',@(s,e)findMaxScanAngle(), ...
    'Tooltip','Find the largest steering angle (°) before a grating lobe.');
uibutton(gAn,'Text','Band sweep','Icon',ribbonIconFile('band'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnBandSweep', ...
    'ButtonPushedFcn',@(s,e)showBandSweep(), ...
    'Tooltip',['Compare pointing, gain and beamwidth across a frequency ' ...
    'band (GHz), with phase shifters and true-time delay.']);
uibutton(gAn,'Text','Scan coverage','Icon',ribbonIconFile('coverage'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnCoverage', ...
    'ButtonPushedFcn',@(s,e)showCoverage(), ...
    'Tooltip',['Map scan loss (dB) and possible grating lobes over ' ...
    'steering angle (°) and azimuth (°).']);

% ---- View ribbon -------------------------------------------------------
% The cut-plane buttons were here at one point and are not any more: they
% belong beside the cut controls they modify, on the View tab, and having
% them in both places meant two sets of identical buttons in one window.
% What IS here acts on the plots as a whole: the REFERENCE group keeps
% the current principal cut as a dashed trace to compare later designs
% against. Clear starts disabled -- nothing is pinned yet.
rbView = uipanel(taskStack,'BorderType','none');
rbView.Layout.Row = 1; rbView.Layout.Column = 1;
gVw = ribbonRow(rbView,{226,250,140});
gRef = ribbonPanel(gVw,'REFERENCE');
gRef.ColumnWidth = {'1x','1x'};
uibutton(gRef,'Text','Pin as reference','Icon',ribbonIconFile('pin'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnPinRef', ...
    'Tooltip',['Keep the principal cut shown now (Total, dB) as a dashed ' ...
    'reference trace, to compare the next design changes against it'], ...
    'ButtonPushedFcn',@(s,e)pinReference());
btnClearRef = uibutton(gRef,'Text','Clear reference','Icon',ribbonIconFile('cross'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnClearRef','Enable','off', ...
    'Tooltip','Remove the pinned reference trace from the cut', ...
    'ButtonPushedFcn',@(s,e)clearReference());

gOverlay = ribbonPanel(gVw,'CST COMPARISON');
gOverlay.ColumnWidth = {'1x','1x'};
uibutton(gOverlay,'Text','Load CST array','Icon',ribbonIconFile('import'), ...
    'IconAlignment','top','Tag','btnLoadCstOverlay', ...
    'Tooltip','Overlay a whole-array CST realized-gain ASCII export on the cut', ...
    'ButtonPushedFcn',@(s,e)loadCstOverlay());
uibutton(gOverlay,'Text','Clear overlay','Icon',ribbonIconFile('cross'), ...
    'IconAlignment','top','Tag','btnClearCstOverlay', ...
    'Tooltip','Remove the CST array comparison; keep the imported element', ...
    'ButtonPushedFcn',@(s,e)clearCstOverlay());
overlayWin = gobjects(0);
cstComparisonClosed = false;

gPattern2D = ribbonPanel(gVw,'PATTERN');
gPattern2D.ColumnWidth = {'1x'};
uibutton(gPattern2D,'Text','2D pattern','Icon',ribbonIconFile('polar2d'), ...
    'IconAlignment','top','FontSize',11,'Tag','btnPattern2D', ...
    'Tooltip','Open polar azimuth/elevation cuts or a direction-cosine U cut', ...
    'ButtonPushedFcn',@(s,e)showPattern2D(''));

ribbons = [rbArray rbElem rbBeam rbView];

commonRow = uigridlayout(gBody,[1 3]);
commonGroups(commonRow);

% ================================================================ menu bar
% A second route to commands the ribbon already has, plus the place for
% ones with no button. Every item calls the SAME function as its ribbon
% button, never a copy of its logic, so the two cannot drift apart.
% The top-level menus (and the Export submenu) are parent-scope handles
% so later features can append their own items to them. uimenu appends
% at the bottom; 'Position',k puts an item k-th from the top instead
% (the Position of an existing item inserts just above it).
% Opening File also re-reads the recent list, so a design deleted or moved
% since start-up drops out of Open Recent before it can be picked.
mFile = uimenu(fig,'Text','File','Tag','menuFile', ...
    'Tooltip','Open, save and export the design', ...
    'MenuSelectedFcn',@(s,e)refreshRecentMenu());
uimenu(mFile,'Text','New','Tag','menuNew', ...
    'Tooltip','Start a fresh design with the original default settings', ...
    'MenuSelectedFcn',@(s,e)onCloseRequest(true));
uimenu(mFile,'Text','Open…','Tag','menuLoad', ...
    'Tooltip',['Open a design saved as .mat (same as the Open button). ' ...
    'Undo (⌘Z) goes back to the design that was open before.'], ...
    'MenuSelectedFcn',@(s,e)loadConfig());
mRecent = uimenu(mFile,'Text','Open Recent','Tag','menuRecent', ...
    'Tooltip',sprintf(['The last %d designs saved or opened; files that ' ...
    'no longer exist are left out'], RECENT_MAX));
uimenu(mFile,'Text','Save','Tag','menuSave','Separator','on', ...
    'Tooltip',['Save the design to its file (same as the Save button); ' ...
    'a design not saved yet asks for a file name'], ...
    'MenuSelectedFcn',@(s,e)saveConfig());
uimenu(mFile,'Text','Save As…','Tag','menuSaveAs', ...
    'Tooltip',['Save the design to a new .mat file, which then becomes ' ...
    'the open design'], ...
    'MenuSelectedFcn',@(s,e)saveConfigAs());
mExport = uimenu(mFile,'Text','Export','Tag','menuExport','Separator','on', ...
    'Tooltip','Write the element positions and excitations for other tools');
uimenu(mExport,'Text','Element table (CSV)…','Tag','menuExportCSV', ...
    'Tooltip',['Per-element x, y (λ), amplitude, phase and rotation (°) ' ...
    'as CSV (same as the CSV button)'], ...
    'MenuSelectedFcn',@(s,e)exportCSV());
uimenu(mExport,'Text','CST array (.tsv)…','Tag','menuExportTSV', ...
    'Tooltip',['Element positions (m) and excitations as a CST array ' ...
    'task .tsv (same as the CST .tsv button)'], ...
    'MenuSelectedFcn',@(s,e)exportTSV());
uimenu(mExport,'Text','CST macro (.bas)…','Tag','menuExportMacro', ...
    'Tooltip',['CST VBA macro that sets every port''s amplitude and ' ...
    'phase (°) as a Combine Results excitation (same as the CST macro ' ...
    'button and the Phase scheme map''s export)'], ...
    'MenuSelectedFcn',@(s,e)exportCstMacro());

mEdit = uimenu(fig,'Text','Edit','Tag','menuEdit', ...
    'Tooltip','Undo and redo design changes; change the set of elements');
uimenu(mEdit,'Text','Delete selected elements','Tag','menuDelete', ...
    'Separator','on', ...   % below Undo / Redo, added at the top next
    'Tooltip',['Remove the selected elements (same as Delete on the ' ...
    'Array tab, or the Delete key)'], ...
    'MenuSelectedFcn',@(s,e)deleteSelected());
% Undo / Redo head the menu. An Accelerator is one character (Cmd on a
% Mac, Ctrl elsewhere) and cannot carry Shift, so Redo's is Y; onKey adds
% the Mac-standard Shift+Cmd+Z. Both start disabled (nothing to undo yet)
% and name the step they would take back, e.g. "Undo taper".
mUndo = uimenu(mEdit,'Text','Undo','Tag','menuUndo','Position',1, ...
    'Accelerator','Z','Enable','off', ...
    'Tooltip',['Take back the last design change: geometry, element, ' ...
    'taper, steering, frequency or a table edit (up to 50 steps). ' ...
    'View settings are not affected.'], ...
    'MenuSelectedFcn',@(s,e)undoFromUI('undo','menu'));
mRedo = uimenu(mEdit,'Text','Redo','Tag','menuRedo','Position',2, ...
    'Accelerator','Y','Enable','off', ...
    'Tooltip',['Re-apply the design change just undone (also ' ...
    'Shift+Cmd+Z). A new change clears it.'], ...
    'MenuSelectedFcn',@(s,e)undoFromUI('redo','menu'));
uimenu(mEdit,'Text','Clear all elements','Tag','menuClearAll', ...
    'Tooltip','Remove every element (same as Clear on the Array tab)', ...
    'MenuSelectedFcn',@(s,e)clearAll());

mAnalysis = uimenu(fig,'Text','Analysis','Tag','menuAnalysis', ...
    'Tooltip','Pattern calculation and scan analyses');
uimenu(mAnalysis,'Text','Compute pattern','Tag','menuCompute', ...
    'Tooltip',['Recalculate the 3D pattern, the cut and the metrics now ' ...
    '(same as the Compute pattern button)'], ...
    'MenuSelectedFcn',@(s,e)computePattern());
uimenu(mAnalysis,'Text','Scan loss vs steering angle','Tag','menuScanLoss', ...
    'Separator','on', ...
    'Tooltip',['Directivity relative to broadside (dB) while steering ' ...
    '0-80° in the current plane (same as Scan loss on the Beam tab)'], ...
    'MenuSelectedFcn',@(s,e)showScanLoss());
uimenu(mAnalysis,'Text','Max scan angle','Tag','menuMaxScan', ...
    'Tooltip',['Largest steering angle (°) before a grating lobe appears, ' ...
    'in the current plane (same as Max scan on the Beam tab)'], ...
    'MenuSelectedFcn',@(s,e)findMaxScanAngle());
uimenu(mAnalysis,'Text','Phase scheme map','Tag','menuPhaseMap', ...
    'Tooltip',['Feed phase (°) of every port for the current steering ' ...
    '(same as the Phase scheme map button)'], ...
    'MenuSelectedFcn',@(s,e)showPortPhases());
uimenu(mAnalysis,'Text','Band sweep…','Tag','menuBandSweep', ...
    'Separator','on', ...
    'Tooltip',['Beam pointing error (°), gain toward the commanded ' ...
    'direction (dBi) and HPBW (°) across a frequency band (GHz), with ' ...
    'phase shifters and with true-time delay'], ...
    'MenuSelectedFcn',@(s,e)showBandSweep());

% One item per task tab; selectTask ticks the open one.
mView = uimenu(fig,'Text','View','Tag','menuView', ...
    'Tooltip',['Switch the task shown in the ribbon and the left panel, ' ...
    'and pin a reference cut']);
for kMenu = 1:numel(designTaskNames)
    uimenu(mView,'Text',[taskNames{kMenu} ' tab'], ...
        'Tag',['menuTab' taskNames{kMenu}], ...
        'Tooltip',['Show the ' taskNames{kMenu} ' task: its ribbon ' ...
        'groups and its settings on the left (same as the tab)'], ...
        'MenuSelectedFcn',@(s,e)selectTask(taskNames{kMenu}));
end
% The REFERENCE ribbon group's two actions. Tags must not start with
% 'menuTab': selectTask ticks exactly those items.
uimenu(mView,'Text','Pin cut as reference','Tag','menuPinRef','Separator','on', ...
    'Tooltip',['Keep the principal cut shown now as a dashed reference ' ...
    'trace (same as Pin as reference on the View tab)'], ...
    'MenuSelectedFcn',@(s,e)pinReference());
menuClearRef = uimenu(mView,'Text','Clear reference','Tag','menuClearRef', ...
    'Enable','off', ...
    'Tooltip',['Remove the pinned reference trace from the cut (same as ' ...
    'Clear reference on the View tab)'], ...
    'MenuSelectedFcn',@(s,e)clearReference());

uimenu(mView,'Text','Load CST array result…','Tag','menuLoadCstOverlay', ...
    'Tooltip','Load the whole array result, separately from the element pattern', ...
    'MenuSelectedFcn',@(s,e)loadCstOverlay());
uimenu(mView,'Text','Clear CST overlay','Tag','menuClearCstOverlay', ...
    'Tooltip','Remove the CST comparison trace and table', ...
    'MenuSelectedFcn',@(s,e)clearCstOverlay());
menuPattern2D = uimenu(mView,'Text','2D pattern','Tag','menuPattern2D', ...
    'Tooltip','Open an azimuth, elevation, or U-pattern cut in the main workspace');
for patternName = {'Azimuth pattern','Elevation pattern','U pattern'}
    itemName = patternName{1};
    uimenu(menuPattern2D,'Text',itemName, ...
        'Tag',['menuPattern2D' itemName(1)], ...
        'Tooltip',['Open the ' lower(itemName) ' in the main workspace'], ...
        'MenuSelectedFcn',@(s,e)showPattern2D(itemName));
end
% In-app workflow help.

uimenu(mAnalysis,'Text','Scan coverage…','Tag','menuCoverage', ...
    'Tooltip','Map scan loss over azimuth and steering angle with grating-lobe markers', ...
    'MenuSelectedFcn',@(s,e)showCoverage());
mHelp = uimenu(fig,'Text','Help','Tag','menuHelp', ...
    'Tooltip','Guides and information about the app');

body = uigridlayout(root,[1 3]);
body.ColumnWidth = {330,'1x','1.15x'};
body.Padding = [0 0 0 0];
body.ColumnSpacing = 8;
outer = body;

uimenu(mHelp,'Text','Quick start','Tag','menuQuickStart', ...
    'Tooltip','Five steps from an array layout to CST export', ...
    'MenuSelectedFcn',@(s,e)showHelp('quick'));
uimenu(mHelp,'Text','Online guide','Tag','menuUserGuide', ...
    'Tooltip','Open the current phased-array walkthrough on GitHub', ...
    'MenuSelectedFcn',@(s,e)showHelp('guide'));
uimenu(mHelp,'Text','About','Tag','menuAbout', ...
    'Tooltip','App version, author and MATLAB release', ...
    'MenuSelectedFcn',@(s,e)showHelp('about'));
mTemplates = uimenu(mFile,'Text','New from template','Tag','menuTemplates', ...
    'Position',1,'Tooltip','Start a preset design; Undo restores the current design');
uimenu(mTemplates,'Text','Ka-band SATCOM Rx 16×16','Tag','templateKa', ...
    'Tooltip','16×16 square lattice, 0.45λ at 21.4 GHz; operates at 19.45 GHz', ...
    'MenuSelectedFcn',@(s,e)applyTemplate('ka'));
uimenu(mTemplates,'Text','8×8 demo (X-band)','Tag','templateDemo', ...
    'Tooltip','8×8 isotropic array at 10 GHz with half-wavelength spacing', ...
    'MenuSelectedFcn',@(s,e)applyTemplate('demo'));
uimenu(mTemplates,'Text','Blank','Tag','templateBlank', ...
    'Tooltip','Empty 10 GHz sparse grid for placing your own elements', ...
    'MenuSelectedFcn',@(s,e)applyTemplate('blank'));
set(findall(fig,'Tag','menuSave'),'Accelerator','S');
set(findall(fig,'Tag','menuLoad'),'Accelerator','O');
set(findall(fig,'Tag','menuCompute'),'Accelerator','R');

% ============================================================== status bar
% One line under everything: what the app is doing now, how long the
% last calculation took, and why the plots are empty when they are.
statusBar = uilabel(root,'Text','Ready','FontSize',11,'Tag','statusBar', ...
    'Tooltip',['What the app last did or is doing now: calculation ' ...
    'time (s), element count, or why the plots are out of date.']);
statusBar.Layout.Row = 3; statusBar.Layout.Column = 1;

% =================================================== LEFT: control panel
% One panel per ribbon tab, all stacked in the SAME grid cell with only
% the active one visible. A uitabgroup here would put a second row of
% tabs directly under the ribbon's, asking the user to track two
% independent selections that in practice always move together.
%
% Compute and the Metrics readout stay outside the stack: both apply
% whichever task you are on, and a Compute button that vanished when you
% switched tabs would be a trap.
leftG = uigridlayout(outer,[3 1]);
leftG.RowHeight = {'1x', 34, 215};
leftG.Padding = [0 0 0 0];
leftG.RowSpacing = 6;

stack = uigridlayout(leftG,[1 1]);
stack.Padding = [0 0 0 0];

sepList = struct('g',{},'row',{},'h',{});

% ------------------------------------------------------------ tab: Array
% Everything that decides WHERE the elements sit and WHICH WAY they face:
% the lattice, the rotations applied across it, and the per-element
% position edits.
secArray = uipanel(stack,'BorderType','none');
secArray.Layout.Row = 1; secArray.Layout.Column = 1;
c = tabGrid(secArray); r = 1;

lbl(c,r,'Placement mode');
ddMode = uidropdown(c,'Items',{'Uniform grid','Sparse (click cells)', ...
    'Sub-position (9-pt)'},'Value',S.mode,'ValueChangedFcn',@onMode, ...
    'Tag','ddMode','Tooltip','Choose a filled grid, clicked sparse cells, or nine offsets inside each cell.');
ddMode.Layout.Row = r; ddMode.Layout.Column = 2; r = r+1;

lbl(c,r,'Rows M / Cols N');
gMN = uigridlayout(c,[1 2]); gMN.Padding=[0 0 0 0]; gMN.ColumnSpacing=4;
gMN.Layout.Row = r; gMN.Layout.Column = 2;
% RoundFractionalValues, not just Step: Step only sets what the arrows
% add, and a TYPED 8.1 is accepted as-is. That reached
% zeros(S.M*S.N, 5) with a non-integer size and threw "Size inputs must
% be integers" straight out of the callback, leaving the app erroring on
% every later rebuild. Element counts are counts.
spM = uispinner(gMN,'Limits',[1 64],'Value',S.M,'Step',1, ...
    'RoundFractionalValues','on','ValueChangedFcn',@onGeom, ...
    'Tag','spM','Tooltip','Number of array rows, from 1 to 64.');
spN = uispinner(gMN,'Limits',[1 64],'Value',S.N,'Step',1, ...
    'RoundFractionalValues','on','ValueChangedFcn',@onGeom, ...
    'Tag','spN','Tooltip','Number of array columns, from 1 to 64.');
r = r+1;

lbl(c,r,'dx / dy (λ)');
gD = uigridlayout(c,[1 2]); gD.Padding=[0 0 0 0]; gD.ColumnSpacing=4;
gD.Layout.Row = r; gD.Layout.Column = 2;
spDx = uispinner(gD,'Limits',[0.05 5],'Value',S.dx,'Step',0.05,'ValueChangedFcn',@onGeom, ...
    'Tag','spDx','Tooltip','Column spacing in design wavelengths; the mm readout uses the design frequency.');
spDy = uispinner(gD,'Limits',[0.05 5],'Value',S.dy,'Step',0.05,'ValueChangedFcn',@onGeom, ...
    'Tag','spDy','Tooltip','Row spacing in design wavelengths; the mm readout uses the design frequency.');
r = r+1;
lblSpacingMM = uilabel(c,'Text','','Tag','lblSpacingMM','FontSize',10,'WordWrap','on', ...
    'Tooltip','Physical dx and dy in mm at the design frequency');
lblSpacingMM.Layout.Row = [r r+1]; lblSpacingMM.Layout.Column = [1 2];
mutedLbls(end+1) = lblSpacingMM; r = r+2;

sep(c,r); r = r+1;

% CST-style boundary shape: trims which cells of the MxN index grid
% are populated (Uniform grid mode only -- Sparse/Sub-position modes
% keep using the plain orthogonal dx,dy lattice for cell placement).
lbl(c,r,'Array shape');
ddShape = uidropdown(c,'Items',{'Custom','Diamond','Hexagon','Octagon', ...
    'Circle','Ellipse'},'Value',S.arrayShape,'ValueChangedFcn',@(s,e)assignShape(s.Value), ...
    'Tag','ddShape','Tooltip',['Boundary mask applied to the uniform grid. ' ...
    'Ellipse is a horizontal oval at least 4:3 wide; Custom keeps the rectangular array.']);
ddShape.Layout.Row = r; ddShape.Layout.Column = 2; r = r+1;

% CST-style grid angle: angle between the two lattice basis vectors.
% 90 deg (the CST default) is the ordinary orthogonal grid; other
% values skew it into a rhombic/oblique lattice (e.g. 60 deg for a
% triangular/hex-packed lattice).
lbl(c,r,'Grid angle (deg)');
spGridAngle = uispinner(c,'Limits',[1 179],'Value',S.gridAngle,'Step',5, ...
    'ValueChangedFcn',@onGeom, ...
    'Tag','spGridAngle','Tooltip','Angle between lattice basis vectors in degrees; 90 gives a rectangular lattice.');
spGridAngle.Layout.Row = r; spGridAngle.Layout.Column = 2; r = r+1;

% Brick/staggered lattice: even rows shift by this along x, odd rows
% don't -- independent of Grid angle (which skews every row instead of
% just alternating ones). Set Grid angle=90 to use this on its own.
% Row stagger (length) and Stagger angle (deg) are two synced VIEWS of
% the SAME underlying value -- angle = atan2d(stagger,dy), i.e. the
% angle the row-2 diagonal-neighbour direction makes with straight-up.
% Editing either one updates the other; at dy=0.5 (the default), a
% stagger of 0.5 lambda IS exactly 45 deg, since equal x/y offset from
% one row to the next is a 45 deg diagonal by definition.
lbl(c,r,'Row stagger (λ)');
spStagger = uispinner(c,'Limits',[-5 5],'Value',S.stagger,'Step',0.05, ...
    'ValueChangedFcn',@onGeom, ...
    'Tag','spStagger','Tooltip','Offset of alternate rows along x, in design wavelengths.');
spStagger.Layout.Row = r; spStagger.Layout.Column = 2; r = r+1;

lbl(c,r,'Stagger angle (deg)');
% Limits widened to +-89.5: atan2d(stagger,dy) can reach ~89.43 deg at
% the corner of these spinners' own ranges (stagger=5, dy=0.05, its
% Limits' minimum) -- a tighter +-89 would force clampedStaggerAngleDeg
% to actually clamp there, leaving the displayed angle not matching the
% real S.stagger/S.dy ratio. Covering the true reachable range means
% the clamp (kept below as a defensive fallback) should never actually
% need to trigger for any combination the other two spinners allow.
spStaggerAngle = uispinner(c,'Limits',[-89.5 89.5],'Value',atan2d(S.stagger,S.dy),'Step',5, ...
    'ValueChangedFcn',@(s,e)assignStaggerAngle(s.Value), ...
    'Tag','spStaggerAngle','Tooltip','Alternate-row offset angle in degrees; synchronized with row stagger.');
spStaggerAngle.Layout.Row = r; spStaggerAngle.Layout.Column = 2; r = r+1;

sep(c,r); r = r+1;

lbl(c,r,'Rotation (0=+y, +=CW)');
spRot = uispinner(c,'Limits',[-360 360],'Value',45,'Step',15, ...
    'Tag','spRot','Tooltip','Element rotation in degrees; zero points along +y and positive is clockwise.');
spRot.Layout.Row = r; spRot.Layout.Column = 2; r = r+1;

lbl(c,r,'Apply to selected');
gRot = uigridlayout(c,[1 4]); gRot.Padding=[0 0 0 0]; gRot.ColumnSpacing=3;
gRot.Layout.Row = r; gRot.Layout.Column = 2;
uibutton(gRot,'Text','Set','ButtonPushedFcn',@(s,e)setRot(spRot.Value), ...
    'Tag','control17','Tooltip','Set. Apply this setting or action to the current design.');
uibutton(gRot,'Text','-','ButtonPushedFcn',@(s,e)stepRot(-spRot.Value), ...
    'Tag','control18','Tooltip','-. Apply this setting or action to the current design.');
uibutton(gRot,'Text','+','ButtonPushedFcn',@(s,e)stepRot(spRot.Value), ...
    'Tag','control19','Tooltip','+. Apply this setting or action to the current design.');
uibutton(gRot,'Text','0','ButtonPushedFcn',@(s,e)setRot(0), ...
    'Tag','control20','Tooltip','0. Apply this setting or action to the current design.');
r = r+1;

sep(c,r); r = r+1;
lbl(c,r,'Seq. rot. step (deg)');
spSeqStep = uispinner(c,'Limits',[-360 360],'Value',90,'Step',15, ...
    'Tag','spSeqStep','Tooltip','Physical element rotation increment in degrees, applied by the sequential-rotation button.');
spSeqStep.Layout.Row = r; spSeqStep.Layout.Column = 2; r = r+1;

% Block size the rotation pattern repeats over -- default 2x2 is the
% classic CP sequential-rotation subarray (0/90/180/270 tiled
% identically across the array). Set to M x N (the full array size, in
% the Rows M / Cols N controls on the Array tab) to fall back to one
% continuously-incrementing pattern across the whole aperture instead.
lbl(c,r,'Rot. block (rows x cols)');
gSeqBlk = uigridlayout(c,[1 2]); gSeqBlk.Padding=[0 0 0 0]; gSeqBlk.ColumnSpacing=4;
gSeqBlk.Layout.Row = r; gSeqBlk.Layout.Column = 2;
spSeqBlkM = uispinner(gSeqBlk,'Limits',[1 64],'Value',S.seqBlockM,'Step',1, ...
    'RoundFractionalValues','on', ...
    'ValueChangedFcn',@(s,e)assignSeqBlock('m',s.Value), ...
    'Tag','spSeqBlkM','Tooltip','Number of rows in the repeating sequential-rotation block.');
spSeqBlkN = uispinner(gSeqBlk,'Limits',[1 64],'Value',S.seqBlockN,'Step',1, ...
    'RoundFractionalValues','on', ...
    'ValueChangedFcn',@(s,e)assignSeqBlock('n',s.Value), ...
    'Tag','spSeqBlkN','Tooltip','Number of columns in the repeating sequential-rotation block.');
r = r+1;

lbl(c,r,'Sequential rot.');
cbSeq = uicheckbox(c,'Text','extra rot angle -> feed phase','Value',S.seqPhase, ...
    'ValueChangedFcn',@(s,e)assignSeq(s.Value), ...
    'Tag','cbSeq','Tooltip','extra rot angle -> feed phase. Apply this setting or action to the current design.');
cbSeq.Tooltip = ['Advanced extra electrical phase term for the pattern model. ' ...
    'It does not define the physical L/R/B/T sequential-rotation topology and ' ...
    'does not alter dual-feed LP/CP phase-map or CST-export phases.'];
cbSeq.Layout.Row = r; cbSeq.Layout.Column = 2; r = r+1;

bSeqAuto = uibutton(c,'Text','Apply physical sequential rotation to all', ...
    'ButtonPushedFcn',@(s,e)applySeqRot(spSeqStep.Value,spSeqBlkM.Value,spSeqBlkN.Value), ...
    'Tag','bSeqAuto','Tooltip','Apply physical sequential rotation to all. Apply this setting or action to the current design.');
bSeqAuto.Tooltip = ['Starts at the physical bottom-left element of each block. ' ...
    'That element is the 0-degree reference; with a +90-degree 2x2 sequence: ' ...
    'BL=0, BR=90, TR=180, TL=270 degrees.'];
bSeqAuto.Layout.Row = r; bSeqAuto.Layout.Column = [1 2]; r = r+1;

sep(c,r); r = r+1;

lbl(c,r,'Move dx / dy (λ)');
gMv = uigridlayout(c,[1 2]); gMv.Padding=[0 0 0 0]; gMv.ColumnSpacing=4;
gMv.Layout.Row = r; gMv.Layout.Column = 2;
spMvX = uispinner(gMv,'Limits',[-5 5],'Value',0.1,'Step',0.05, ...
    'Tag','spMvX','Tooltip','Displacement of selected elements along x, in design wavelengths.');
spMvY = uispinner(gMv,'Limits',[-5 5],'Value',0,'Step',0.05, ...
    'Tag','spMvY','Tooltip','Displacement of selected elements along y, in design wavelengths.');
r = r+1;

bMove = uibutton(c,'Text','Move selected by (dx,dy)', ...
    'ButtonPushedFcn',@(s,e)moveSelected(spMvX.Value,spMvY.Value), ...
    'Tag','bMove','Tooltip','Move selected elements by the dx/dy offsets entered above, in design wavelengths.');
bMove.Layout.Row = r; bMove.Layout.Column = [1 2]; r = r+1;

bResetLat = uibutton(c,'Text','Reset selected to lattice', ...
    'ButtonPushedFcn',@(s,e)resetSelectedToLattice(), ...
    'Tag','bResetLat','Tooltip','Restore selected elements to their original lattice coordinates.');
bResetLat.Layout.Row = r; bResetLat.Layout.Column = [1 2]; r = r+1;

sep(c,r); r = r+1;

bDel = uibutton(c,'Text','Delete selected element(s)', ...
    'ButtonPushedFcn',@(s,e)deleteSelected(), ...
    'Tag','bDel','Tooltip','Delete selected elements; Undo restores them.');
bDel.Layout.Row = r; bDel.Layout.Column = [1 2]; r = r+1;

bClear = uibutton(c,'Text','Clear all elements','ButtonPushedFcn',@(s,e)clearAll(), ...
    'Tag','bClear','Tooltip','Remove all elements; Undo restores the previous layout.');
bClear.Layout.Row = r; bClear.Layout.Column = [1 2]; r = r+1;

finishTab(c,r-1);

% ---------------------------------------------------------- tab: Element
% What ONE element radiates, and nothing about the array.
secElem = uipanel(stack,'BorderType','none');
secElem.Layout.Row = 1; secElem.Layout.Column = 1;
c = tabGrid(secElem); r = 1;

lbl(c,r,'Element factor');
ddEF = uidropdown(c,'Items',{'Isotropic','cos^q(theta)', ...
    'Cardioid','Gaussian','Sinc','Short dipole (z-axis)','Dipole (linear pol)', ...
    'Crossed dipole (RHCP)','Patch (cos^q x lin pol)', ...
    '3GPP TR 38.901 shape','Custom (formula)', ...
    'Imported (CST far-field)'}, ...
    'Value',S.efType,'ValueChangedFcn',@(s,e)assignEF(s.Value), ...
    'Tag','ddEF','Tooltip',['Radiation model of one element. Short dipole is z-directed ' ...
    'with axial nulls; Dipole (linear pol) is the legacy in-plane model. ' ...
    'CST import uses realized gain and complex Theta/Phi components.']);
ddEF.Layout.Row = r; ddEF.Layout.Column = 2; r = r+1;

lbl(c,r,'q exponent');
spQ = uispinner(c,'Limits',[0 6],'Value',S.efQ,'Step',0.25, ...
    'ValueChangedFcn',@(s,e)assignQ(s.Value), ...
    'Tag','spQ','Tooltip','Power-pattern cosine exponent q; larger values narrow the element pattern.');
spQ.Layout.Row = r; spQ.Layout.Column = 2; r = r+1;

lbl(c,r,'Az HPBW (°)');
spEfBeamAz = uispinner(c,'Limits',[10 180],'Value',S.efBeamAz,'Step',5, ...
    'ValueChangedFcn',@(s,e)assignEfBeamwidth('az',s.Value), ...
    'Tag','spEfBeamAz','Tooltip','Azimuth half-power beamwidth for Gaussian, Sinc and 3GPP-style element patterns.');
spEfBeamAz.Layout.Row = r; spEfBeamAz.Layout.Column = 2; r = r+1;
lbl(c,r,'El HPBW (°)');
spEfBeamEl = uispinner(c,'Limits',[10 180],'Value',S.efBeamEl,'Step',5, ...
    'ValueChangedFcn',@(s,e)assignEfBeamwidth('el',s.Value), ...
    'Tag','spEfBeamEl','Tooltip','Elevation half-power beamwidth for Gaussian, Sinc and 3GPP-style element patterns.');
spEfBeamEl.Layout.Row = r; spEfBeamEl.Layout.Column = 2; r = r+1;

sep(c,r); r = r+1;

lbl(c,r,'Custom E_theta(th,ph)');
efCustom = uieditfield(c,'text','Value',S.customFormula, ...
    'ValueChangedFcn',@(s,e)assignCustomFormula(s.Value), ...
    'Tag','efCustom','Tooltip','Complex E_theta formula; th/ph are in degrees. Use element-wise operators such as .^ and .*.');
efCustom.Layout.Row = r; efCustom.Layout.Column = 2; r = r+1;
lbl(c,r,'Custom E_phi(th,ph)');
efCustomPh = uieditfield(c,'text','Value',S.customFormulaPh, ...
    'ValueChangedFcn',@(s,e)assignCustomFormulaPh(s.Value), ...
    'Tag','efCustomPh','Tooltip','Complex E_phi formula; th/ph are in degrees. Relative phase sets polarization.');
efCustomPh.Layout.Row = r; efCustomPh.Layout.Column = 2; r = r+1;
lblCustomHint = uilabel(c,'Text','vars: th,ph,ct,st (theta,phi also work). NOTE ct = max(cosd(th),0), clamped at 0 -- for back radiation in full-sphere mode use cosd(th) instead. Use dot operators (.* and .^). E_phi=0 gives a linear-pol pattern; set it nonzero for elliptical/CP custom elements. Bad formula -> isotropic', ...
    'FontSize',9,'WordWrap','on');
mutedLbls(end+1) = lblCustomHint;   % coloured by paintChrome at start-up
lblCustomHint.Layout.Row = [r r+2]; lblCustomHint.Layout.Column = [1 2]; r = r+3;

sep(c,r); r = r+1;

% Imported (CST far-field): a single button, a peer to every other
% Element factor option -- load a far-field .txt export (CST Farfield ->
% Export -> ASCII) and it's used directly as the element pattern. This
% app's normal array math (array factor, steering, taper) combines it
% across the array, "Total" shows the real combined result, same as
% Isotropic/cos^q/Dipole/Patch/Custom. No port concept involved.
bImportFF = uibutton(c,'Text','Load far-field pattern (.txt)...', ...
    'ButtonPushedFcn',@(s,e)loadImportedFF(), ...
    'Tag','bImportFF','Tooltip','Load an ELEMENT CST realized-gain ASCII export with complex Theta/Phi components.');
bImportFF.Layout.Row = r; bImportFF.Layout.Column = [1 2]; r = r+1;
% Under the button: WHAT to export from CST before a file is loaded (the
% format was otherwise learned one refusal at a time), and what was read
% once one is -- grid, peak, and the frequency, which the element
% pattern does not scale with. refreshImportInfo writes it.
lblImportInfo = uilabel(c,'Text','','FontSize',10,'WordWrap','on', ...
    'VerticalAlignment','top','Tag','lblImportInfo');
lblImportInfo.Layout.Row = [r r+3]; lblImportInfo.Layout.Column = [1 2];
r = r+4;   % six lines of 10 pt text at the narrowest panel width

% Tick this when the imported pattern came from a UNIT CELL with periodic
% boundaries. Such a file holds the embedded element's pattern SHAPE, but
% its peak is the embedded element GAIN, which the lattice caps at
% 4*pi*A_cell/lambda^2 -- not a free-space gain. The two are set by
% different normalization assumptions. Integrating this exported active
% element pattern as though it were an isolated antenna need not recover
% its efficiency. This mode preserves the supplied realized-gain scale.
cbUnitCell = uicheckbox(c,'Text','Imported unit cell (embedded element)', ...
    'Value',S.impUnitCell,'ValueChangedFcn',@(s,e)assignUnitCell(s.Value), ...
    'Tag','cbUnitCell','Tooltip','Use the embedded-element gain normalization for a periodic unit-cell export.');
cbUnitCell.Layout.Row = r; cbUnitCell.Layout.Column = [1 2]; r = r+1;

% Type in the solver's Total Efficiency for the imported pattern. Used
% only to report directivity as RG - 10log10(eff); it never scales the
% pattern, because the file's peak is ALREADY realized gain and applying
% the efficiency again would count it twice.
lbl(c,r,'Total efficiency (%)');
spTotEff = uispinner(c,'Limits',[0.1 100],'Value',S.impTotEffPct,'Step',1, ...
    'ValueChangedFcn',@(s,e)assignTotEff(s.Value), ...
    'Tag','spTotEff','Tooltip','CST total efficiency in percent; used to derive directivity from imported embedded realized gain.');
spTotEff.Layout.Row = r; spTotEff.Layout.Column = 2; r = r+1;

finishTab(c,r-1);

% ------------------------------------------------------------- tab: Beam
secBeam = uipanel(stack,'BorderType','none');
secBeam.Layout.Row = 1; secBeam.Layout.Column = 1;
c = tabGrid(secBeam); r = 1;

lbl(c,r,'Angle convention');
ddAngleConvention = uidropdown(c,'Items',{'Theta / phi','Azimuth / elevation'}, ...
    'Value',S.angleConvention,'Tag','ddAngleConvention', ...
    'ValueChangedFcn',@(s,e)assignAngleConvention(s.Value), ...
    'Tooltip','Angle convention. Apply this setting or action to the current design.');
ddAngleConvention.Layout.Row = r; ddAngleConvention.Layout.Column = 2; r = r+1;
ddAngleConvention.Tooltip = ['Theta is measured from +z; azimuth is measured in the xy plane. ' ...
    'Elevation = 90 - theta. CST files and field-component formulas retain their required theta/phi definitions.'];

lblSteer = lbl(c,r,'Steer theta / phi');
gS = uigridlayout(c,[1 2]); gS.Padding=[0 0 0 0]; gS.ColumnSpacing=4;
gS.Layout.Row = r; gS.Layout.Column = 2;
spTh = uispinner(gS,'Limits',[-90 90],'Value',S.theta_s,'Step',1,'Tag','spSteerFirst', ...
    'ValueChangedFcn',@(s,e)assignSteer('t',s.Value), ...
    'Tooltip','Steer theta / phi. Edit this value in the units shown.');
spPh = uispinner(gS,'Limits',[-180 360],'Value',S.phi_s,'Step',5,'Tag','spSteerSecond', ...
    'ValueChangedFcn',@(s,e)assignSteer('p',s.Value), ...
    'Tooltip','Steer theta / phi. Edit this value in the units shown.');
r = r+1;
cbCutFollow = uicheckbox(c,'Text','2D cut plane follows steering phi', ...
    'Value',S.cutPhiFollow,'ValueChangedFcn',@(s,e)assignCutPhiFollow(s.Value), ...
    'Tag','cbCutFollow', ...
    'Tooltip',['For elevation and U cuts, keep the observation plane at ' ...
        'the steered phi. Uncheck to inspect another plane.']);
cbCutFollow.Layout.Row = r; cbCutFollow.Layout.Column = [1 2]; r = r+1;

% Phasing scheme map. Lives with the steering controls rather than with
% the element geometry: what it shows is the EXCITATION -- the phase each
% feed is driven with for the current steer angle -- which is a beam
% question. The element table carries one phase per ELEMENT, and that is
% the array phase, not what any individual port sees. Covers a plain
% single-feed array as much as a dual-polarised one, so the scheme is
% picked inside the window rather than tied to the CP checkbox; for a
% 64-element dual-pol build that is 128 numbers the tool could not
% previously state.
bPortMap = uibutton(c,'Text','Phase scheme map','Tag','btnPhaseMap', ...
    'ButtonPushedFcn',@(s,e)showPortPhases(), ...
    'Tooltip','Phase scheme map. Apply this setting or action to the current design.');
bPortMap.Layout.Row = r; bPortMap.Layout.Column = [1 2]; r = r+1;

sep(c,r); r = r+1;
lbl(c,r,'Amplitude taper');
ddTap = uidropdown(c,'Items',{'Uniform','Hamming','Hanning','Chebyshev', ...
    'Taylor','Binomial','Manual (table)'},'Value',S.taper,'ValueChangedFcn',@(s,e)assignTaper(s.Value), ...
    'Tag','ddTap','Tooltip','Amplitude taper. Apply this setting or action to the current design.');
ddTap.Layout.Row = r; ddTap.Layout.Column = 2; r = r+1;
ddTap.Tooltip = ['Amplitude weighting across the aperture. Chebyshev and ' ...
    'Taylor are designed to the Sidelobe level below; the other tapers ' ...
    'have fixed sidelobes.'];

% Design sidelobe level of the Chebyshev and Taylor tapers, in dB below
% the main beam (30 -> sidelobes at -30 dB). Positive, as chebwin and the
% saved configs take it, so a typed "25" means what it says. Disabled,
% not hidden, for the other tapers: it keeps its value for when one of
% the two is picked again. The note beside it states the Taylor n-bar the
% app derives from the level (see taperVec), which used to be invisible.
lbl(c,r,'Sidelobe level (dB)');
gSLL = uigridlayout(c,[1 2]); gSLL.Padding = [0 0 0 0];
gSLL.ColumnSpacing = 8; gSLL.ColumnWidth = {'1x','1x'};
gSLL.Layout.Row = r; gSLL.Layout.Column = 2;
spSLL = uispinner(gSLL,'Limits',[13 80],'Step',1,'Value',S.sll, ...
    'Tag','spSLL','ValueChangedFcn',@(s,e)assignSLL(s.Value), ...
    'Tooltip','Sidelobe level (dB). Edit this value in the units shown.');
spSLL.Tooltip = ['Design sidelobe level for the Chebyshev and Taylor ' ...
    'tapers, in dB below the main beam (30 = sidelobes at -30 dB). ' ...
    'Range 13-80 dB. Lower sidelobes cost taper efficiency and widen ' ...
    'the beam. Only used by those two tapers.'];
lblNbar = uilabel(gSLL,'Text','','FontSize',11,'Tag','lblTaylorNbar');
lblNbar.Tooltip = ['Taylor n̄: how many near-in sidelobes ' ...
    'are held at the design level. Set automatically from the sidelobe ' ...
    'level (the smallest value that can reach it).'];
mutedLbls(end+1) = lblNbar;
r = r+1;
syncSLLControl();

finishTab(c,r-1);

% ------------------------------------------------------------- tab: View
secView = uipanel(stack,'BorderType','none');
secView.Layout.Row = 1; secView.Layout.Column = 1;
c = tabGrid(secView); r = 1;

lbl(c,r,'3D surface shows');
% Canonical list, kept in a variable because refreshPolAvailability()
% removes the Axial Ratio entry for imports that cannot support it and
% has to be able to put it back verbatim -- rebuilding the list by hand
% in two places is how the order drifts.
SHOW_ITEMS = {'Total (EF x AF)','Array factor only', ...
    'Element factor only','Axial Ratio (dB)'};
ddShow = uidropdown(c,'Items',SHOW_ITEMS,'Value','Total (EF x AF)', ...
    'ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','ddShow','Tooltip','Quantity drawn on the 3D surface; this does not change the principal cut.');
ddShow.Layout.Row = r; ddShow.Layout.Column = 2; r = r+1;

lbl(c,r,'Polarization readout');
ddPol = uidropdown(c,'Items',{'Total (any pol)','RHCP component','LHCP component'}, ...
    'Value','Total (any pol)','ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','ddPol','Tooltip','Total field power or a circular-polarization component for pattern readouts.');
ddPol.Layout.Row = r; ddPol.Layout.Column = 2; r = r+1;

lbl(c,r,'3D radius scale');
ddScale = uidropdown(c,'Items',{'Linear power','Linear magnitude','dB (dynamic range)'}, ...
    'Value','dB (dynamic range)','ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','ddScale','Tooltip', ...
    'Linear power shows the radiation-lobe shape; dB expands weak sidelobes for inspection.');
ddScale.Layout.Row = r; ddScale.Layout.Column = 2; r = r+1;

lbl(c,r,'Dynamic range (dB)');
spDR = uispinner(c,'Limits',[10 60],'Value',40,'Step',5, ...
    'ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','spDR','Tooltip','Visible pattern range below the peak, in dB; beam metrics use unfloored data.');
spDR.Layout.Row = r; spDR.Layout.Column = 2; r = r+1;

sep(c,r); r = r+1;

% Label says "level", not "directivity": with an imported CST pattern
% that carries a radiation efficiency, this switches the readout to real
% GAIN, not directivity (see impEffLin). The 3D title and the cut y-axis
% both name which one is actually being shown for the current element.
% Text kept short because a uicheckbox does not wrap -- the old
% parenthetical spelled out "(gain if known, else directivity)" and was
% simply truncated with an ellipsis in a panel this width, so it
% communicated nothing. The 3D title and the cut y-axis already name
% which of the two is being shown.
cbAbs = uicheckbox(c,'Text','Show absolute level in dBi', ...
    'Value',false,'ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','cbAbs','Tooltip','Show absolute dBi when available: directivity for analytic models, realized gain for imported patterns.');
cbAbs.Layout.Row = r; cbAbs.Layout.Column = [1 2]; r = r+1;

% This toggle controls display only. All absolute levels integrate full-sphere
% power; one-sided patch/cosine models explicitly provide zero back radiation.
cbFullSphere = uicheckbox(c,'Text','Show full sphere (theta 0-180)', ...
    'Value',S.fullSphere,'ValueChangedFcn',@(s,e)assignFullSphere(s.Value), ...
    'Tag','cbFullSphere','Tooltip','Display both hemispheres; radiation-power integration always uses the full sphere.');
cbFullSphere.Layout.Row = r; cbFullSphere.Layout.Column = [1 2]; r = r+1;

sep(c,r); r = r+1;

% Bottom 2D plot: sweep theta at fixed phi_s (the usual elevation cut),
% or sweep phi at a fixed theta (an azimuth/conical cut) -- same HPBW/
% SLL/FNBW/grating-lobe analysis code handles either, since it's a
% generic 1D-cut peak search either way; only the observation-angle
% sweep and axis labelling differ (see computePattern).
lbl(c,r,'Cut plane');
ddCutMode = uidropdown(c,'Items',{'Theta cut (fixed phi_s)','Phi cut (fixed theta)'}, ...
    'Value',S.cutMode,'ValueChangedFcn',@(s,e)assignCutMode(s.Value), ...
    'Tag','ddCutMode','Tooltip','Choose a theta/elevation sweep or an azimuth sweep at a fixed theta/elevation.');
ddCutMode.ItemsData = {'Theta cut (fixed phi_s)','Phi cut (fixed theta)'};
ddCutMode.Layout.Row = r; ddCutMode.Layout.Column = 2; r = r+1;

% Which plane the THETA cut is taken in. Disabled while the checkbox
% beside the steering controls keeps it slaved to the steering azimuth,
% so the control always displays the plane actually in use rather than
% a stale value that
% silently disagrees with the plot.
lblCutPhi = lbl(c,r,'Cut phi (deg)');
spCutPhi = uispinner(c,'Limits',[-180 360],'Value',S.phi_s,'Step',5,'Tag','spCutAzimuth', ...
    'Enable','off','ValueChangedFcn',@(s,e)assignCutPhi(s.Value), ...
    'Tooltip','Cut phi (deg). Edit this value in the units shown.');
spCutPhi.Layout.Row = r; spCutPhi.Layout.Column = 2; r = r+1;

% One-click planes. A linear array built along x is read in the phi=0
% plane and one along y in phi=90, and switching between those two is
% the whole reason this control exists -- doing it via "untick, then
% type" every time is friction for the common case. Each button unlinks
% from the steering azimuth and sets the plane in one action. "Follow"
% re-links.
lbl(c,r,'Quick cut plane');
gCutQ = uigridlayout(c,[1 3]); gCutQ.Padding=[0 0 0 0]; gCutQ.ColumnSpacing=4;
gCutQ.Layout.Row = r; gCutQ.Layout.Column = 2;
% Tagged like the ribbon buttons. These captions are exactly the terse,
% rename-prone kind that a text lookup should not depend on; they carried
% tags while a duplicate set lived in the ribbon, and lost them when that
% duplicate was removed.
bCut0 = uibutton(gCutQ,'Text','phi=0',  'Tag','btnCut0', ...
    'ButtonPushedFcn',@(s,e)pinCutPhi(0), ...
    'Tooltip','phi=0. Apply this setting or action to the current design.');
bCut90 = uibutton(gCutQ,'Text','phi=90', 'Tag','btnCut90', ...
    'ButtonPushedFcn',@(s,e)pinCutPhi(90), ...
    'Tooltip','phi=90. Apply this setting or action to the current design.');
uibutton(gCutQ,'Text','Follow', 'Tag','btnCutFollow', ...
    'ButtonPushedFcn',@(s,e)followCutPhi(), ...
    'Tooltip','Follow. Apply this setting or action to the current design.');
r = r+1;

lblCutTheta = lbl(c,r,'Theta for phi-cut (deg)');
spCutTheta = uispinner(c,'Limits',[0 90],'Value',S.cutFixedTheta,'Step',5,'Tag','spCutElevation', ...
    'ValueChangedFcn',@(s,e)assignCutFixedTheta(s.Value), ...
    'Tooltip','Theta for phi-cut (deg). Edit this value in the units shown.');
spCutTheta.Layout.Row = r; spCutTheta.Layout.Column = 2; r = r+1;

% Independent of "3D surface shows" -- swaps the 2D cut plot's usual
% AF/EF/Total gain curves for a single Axial Ratio (dB) curve along the
% same cut. Kept as a separate toggle rather than folded into ddShow
% since AR (0 to ~30 dB) and gain (peak-relative or absolute dBi) don't
% share a sensible y-axis scale -- overlaying them on one axis would
% make one or the other unreadable.
cbARCut = uicheckbox(c,'Text','Cut plot shows Axial Ratio (not gain)', ...
    'Value',false,'ValueChangedFcn',@(s,e)safeCompute(), ...
    'Tag','cbARCut','Tooltip','Show axial ratio in dB on the cut instead of gain; lower values indicate more circular polarization.');
cbARCut.Layout.Row = r; cbARCut.Layout.Column = [1 2]; r = r+1;

finishTab(c,r-1);

sections = [secArray secElem secBeam secView];
selectTask('Array');    % ribbon panel, left section and tab highlight together

% ------------------------------------------------- pinned action + metrics
actG = uigridlayout(leftG,[1 2]);
actG.ColumnWidth = {'1x',96};
actG.Padding = [0 0 0 0]; actG.ColumnSpacing = 6;
uibutton(actG,'Text','Compute pattern','FontWeight','bold', ...
    'Tag','btnCompute','ButtonPushedFcn',@(s,e)computePattern(), ...
    'Tooltip','Compute pattern. Apply this setting or action to the current design.');
cbAuto = uicheckbox(actG,'Text','Auto','Value',true, ...
    'Tooltip','Recompute automatically whenever a control changes', ...
    'Tag','cbAuto');

pInfo = uipanel(leftG,'Title','Details');
gInfo = uigridlayout(pInfo,[1 1]);
gInfo.Padding = [8 6 8 6];
% Scrollable + a fit row: the details grow and shrink with the mode
% (unit-cell adds several rows), and a fixed height would clip the tail
% exactly when there is most to read.
gInfo.Scrollable = 'on';
gInfo.RowHeight = {'fit'};
% A two-column table of labels (quantity | value), filled by setDetails.
% Labels rather than a uitable: values wrap in this narrow column instead
% of being cut off, and a second uitable in the window would be mistaken
% for the per-element table by everything that looks that one up by type.
% detailLbls is the pool of label pairs, grown on demand and reused.
detailsGrid = uigridlayout(gInfo,[1 2],'Tag','detailsTable');
detailsGrid.ColumnWidth = {'fit','1x'};
detailsGrid.RowHeight = {'fit'};
detailsGrid.Padding = [0 0 0 0];
detailsGrid.ColumnSpacing = 10; detailsGrid.RowSpacing = 3;
detailLbls = gobjects(0,2);

% =============================================== MIDDLE: layout + table
mid = uigridlayout(outer,[2 1]);
mid.RowHeight = {'1.5x','1x'};
mid.Padding = [0 0 0 0]; mid.RowSpacing = 8;

pLay = uipanel(mid,'Title','Element layout', ...
    'Tooltip','Click an empty cell to add an element; click an element to select it; press Delete to remove it.');
gLay = uigridlayout(pLay,[1 1]); gLay.Padding=[2 2 2 2];
axLay = uiaxes(gLay);
axLay.ButtonDownFcn = @onGridClick;
axis(axLay,'equal'); grid(axLay,'on'); box(axLay,'on');
xlabel(axLay,'x (λ)'); ylabel(axLay,'y (λ)');

pTab = uipanel(mid,'Title', ...
    'Per-element data  (edit Amp, Offset phase, Rot; shaded = computed)');
gTab = uigridlayout(pTab,[1 1]); gTab.Padding=[4 4 4 4];
% Column ORDER is part of the interface: onTableEdit writes columns 3-5
% straight into S.el(:,3:5), and the tests index tbl.Data by position.
% Only the display changes here.
%
% The two phase columns used to be 'Phase (deg)' and 'Feed phase (deg)',
% which hid the one distinction that matters: column 4 is the offset you
% type, column 6 the total the feed actually gets once steering is added
% (the CSV export once mixed them up). Their names are long, so they wrap
% onto two header lines; weights rather than 'fit' (which sizes to the
% data and cut the headers short) share the width in proportion to what
% each column has to show. At a 1440 px window that is ~70 px for the
% short columns and ~90 px for the total -- -179.99 and both header
% lines fit, where the old default widths clipped -90 to '-9'.
tbl = uitable(gTab, ...
    'ColumnName',{'x (λ)','y (λ)','Amp',['Offset' newline 'phase (°)'], ...
        'Rot (°)',['Total feed' newline 'phase (°)']}, ...
    'ColumnEditable',[false false true true true false], ...
    'ColumnWidth',{'10x','10x','10x','12x','10x','13x'}, ...
    'ColumnFormat',repmat({'shortG'},1,6), ...
    'Tag','tblElements', ...
    'Tooltip',['One row per element, numbered as in the layout plot. ' ...
        'x, y: element centre in wavelengths (λ) at the design frequency. ' ...
        'Amp: relative excitation amplitude (linear). Offset phase: your ' ...
        'manual phase offset (°). Rot: element rotation (°, 0 = +y, ' ...
        'positive clockwise). Total feed phase: offset plus the steering ' ...
        'phase (and the rotation compensation when "extra rot angle -> ' ...
        'feed phase" is ticked), wrapped to -180° ... 180°; this is the ' ...
        'phase the CSV / TSV exports write. Amp, Offset phase and Rot are ' ...
        'editable; the shaded columns are computed.'], ...
    'CellEditCallback',@onTableEdit, ...
    'CellSelectionCallback',@onTableSelect);
elemTblStyles = [];   % column styles owned by styleElementTable
% A uitable with keyboard focus can intercept Delete/Backspace before
% the figure's own KeyPressFcn ever sees them -- attaching the SAME
% handler here too means Delete/Escape work whether focus is on the
% layout axes or the table, instead of silently doing nothing after
% clicking into the table.
tbl.KeyPressFcn = @onKey;
% Row/multi-select only exists on R2023a+; setting it at uitable()
% construction would hard-error on older MATLAB and break the whole
% table, so it's set AFTER construction, wrapped in try/catch, so older
% releases just keep plain cell-selection behavior instead of failing.
% Deliberately NOT also registering SelectionChangedFcn here: it fires
% for the same user action as CellSelectionCallback above, and having
% both meant every click did the refresh work twice, with no guarantee
% which one would "win" if they ever disagreed. One callback, one path.
try
    tbl.SelectionType = 'row';
    tbl.Multiselect = 'on';
catch
end
S.sel = [];   % selected table rows

% =================================================== RIGHT: pattern plots
pPat = uipanel(outer,'Title','Radiation pattern');
% Row 1 is the strip of result cards: the answers (gain, beamwidth,
% sidelobes, grating lobes) read at a glance above the plots
% they come from, each checked against the design targets. The cut gets
% nearly the 3D plot's height: angles and levels are read off it, while
% the 3D surface keeps its aspect ratio and leaves width unused anyway.
gPat = uigridlayout(pPat,[3 1]); gPat.RowHeight = {68,'1.15x','1x'};
gPat.Padding=[4 4 4 4];
% The Targets window, its controls and the verdicts refreshCards last
% reached are shared by sibling nested functions, so declared here.
% TARGET_NUMS names the numeric targets in S.targets and TARGET_LIMITS
% their accepted ranges (dBi, °, dB): the window and sanitizeTargets
% (which screens saved ones) read the same list.
targetsWin = []; tgtCtl = struct(); cardVerdicts = {};
TARGET_NUMS = {'minGain','maxHPBW','maxSLL'};
TARGET_LIMITS = {[-100 100], [0 360], [-200 0]};
cardH = buildResultCards(gPat);
ax3D = uiaxes(gPat);
ax3D.Layout.Row = 2; ax3D.Layout.Column = 1;
axCut = uiaxes(gPat);
axCut.Layout.Row = 3; axCut.Layout.Column = 1;
axCut.ButtonDownFcn = @onCutClick;

% Result tabs reuse the right two columns while the design controls stay
% available at the left.  The panels are created lazily by their actions.
resultWorkspace = uipanel(body,'BorderType','none', ...
    'Tag','resultWorkspace','Visible','off');
resultWorkspace.Layout.Row = 1;
resultWorkspace.Layout.Column = [2 3];
resultHost = uigridlayout(resultWorkspace,[1 1]);
resultHost.Padding = [0 0 0 0];

% Create the dropdown after the main panels so it floats above the layout
% and plots when opened. Its six buttons keep the same tags/callbacks as
% the former always-visible ribbon gallery.
[shapeGalleryPopup,shapeBtns] = buildShapeGalleryPopup( ...
    fig,shapeNames,@ribbonIconFile,@pickShape);
btnShapeGallery.ButtonPushedFcn = ...
    @(s,e)toggleShapeGallery(shapeGalleryPopup,fig,gShape);
shapePreviewCurrent.ButtonPushedFcn = btnShapeGallery.ButtonPushedFcn;

% ------------------------------------------------------------------ start
lblLambdaMM.Text = sprintf('Design λ = %.4g mm', lambdaMM());
% Wired once every control exists but BEFORE the first paint. A window
% that follows the desktop theme is created light and only switches to a
% dark desktop theme once it is on screen -- possibly in the middle of
% the first compute -- so this callback, not the start-up paint alone,
% is what gets a dark desktop the right colours. Guarded because
% releases before figure themes have no such property; there the app
% keeps the colours it started with.
try
    fig.ThemeChangedFcn = @(~,~)applyTheme();
catch
end
rebuildUniform();
paintChrome();          % tab + gallery highlight, dividers, captions
% The import readout compares against the operating frequency and takes
% a theme colour, so it rides along with every redraw (and applyTheme)
% through the same list as the open phase maps.
S.portMapRefresh{end+1} = @refreshImportInfo;
refreshAll();
% The start-up design counts as saved: closing an untouched window asks
% nothing. The close button asks about unsaved changes from here on
% (close(fig,'force') and delete(fig) still close without asking).
docSaved = designSnapshot();
defaultDesign = docSaved;
defaultView = struct('angleConvention',S.angleConvention, ...
    'freqUnit',S.freqUnit,'cutMode',S.cutMode, ...
    'cutFixedTheta',S.cutFixedTheta,'fullSphere',S.fullSphere, ...
    'cutPhi',S.cutPhi,'cutPhiFollow',S.cutPhiFollow, ...
    'targets',S.targets,'show',ddShow.Value,'pol',ddPol.Value, ...
    'scale',ddScale.Value,'dynRange',spDR.Value, ...
    'absLevel',cbAbs.Value,'arCut',cbARCut.Value,'auto',cbAuto.Value, ...
    'layoutView',layoutView);
% JIT warm-up is not evidence that an ordinary edit needs a busy flash.
lastComputeSec = 0;
refreshRecentMenu();
fig.CloseRequestFcn = @(~,~)onCloseRequest();

% ================================================================ helpers
    function h = lbl(parent,row,txt)
        % Returns the handle so callers can show/hide a whole row (label
        % + control) together -- a control hidden without its label
        % leaves an orphaned caption behind.
        h = uilabel(parent,'Text',txt,'FontSize',12);
        h.Layout.Row = row; h.Layout.Column = 1;
    end
    function sep(parent,row)
        % A 1 px rule centred in the (short) separator row, in the
        % theme's divider colour. It used to be the whole row painted
        % light grey, which read as a thick bar -- and a glaring one
        % under the dark theme.
        sg = uigridlayout(parent,[3 1]);
        sg.RowHeight = {'1x',1,'1x'};
        sg.Padding = [0 0 0 0]; sg.RowSpacing = 0;
        sg.Layout.Row = row; sg.Layout.Column = [1 2];
        P = pal();
        h = uilabel(sg,'Text','','BackgroundColor',P.rule);
        h.Layout.Row = 2;
        % Records the GRID as well as the row. There used to be one
        % control grid and a bare list of row numbers was enough; with
        % one grid per tab, row 4 means a different place in each of
        % them, and finishTab has to be able to tell them apart. The
        % rule itself is kept so applyTheme can recolour it.
        sepList(end+1) = struct('g',parent,'row',row,'h',h);
    end
    function [bts, rowG] = ribbonSection(parent, titleTxt, specs)
        %RIBBONSECTION  One labelled group of icon buttons in the ribbon.
        %   specs is an n-by-4 cell: {label, iconName, callback, tag}.
        %   Returns the button handles and the grid holding them, so a
        %   caller can append an extra control (the Auto checkbox) into
        %   the same row.
        %
        %   Each button carries a Tag. The regression suite finds these
        %   by tag rather than by their visible text: ribbon labels are
        %   deliberately terse ("CSV", "Compute") and are the kind of
        %   thing that gets reworded, which would otherwise silently turn
        %   a test into a no-op or a hard failure on a pure UI edit.
        n = size(specs,1);
        rowG = ribbonPanel(parent, titleTxt);
        rowG.ColumnWidth = repmat({'1x'},1,n);

        bts = gobjects(n,1);
        for k = 1:n
            bts(k) = uibutton(rowG,'Text',specs{k,1}, ...
                'Icon',ribbonIconFile(specs{k,2}), ...
                'IconAlignment','top', ...
                'FontSize',11, ...
                'Tag',specs{k,4}, ...
                'ButtonPushedFcn',specs{k,3}, ...
    'Tooltip',[specs{k,1} ' for the current array design.']);
        end
    end

    function g = ribbonRow(tab, widths)
        %RIBBONROW  The horizontal strip of groups inside one ribbon tab.
        %   widths gives the fixed width of each TASK-SPECIFIC group; a
        %   '1x' filler and the two common groups are appended by
        %   commonGroups, so the task groups stay left-packed and the
        %   common ones sit hard right on every tab at the same x.
        g = uigridlayout(tab,[1 numel(widths)+1]);
        % Groups left-pack; the trailing '1x' absorbs the slack so a tab
        % with few groups does not stretch them across the full width.
        g.ColumnWidth = [widths {'1x'}];
        g.Padding = [4 4 4 4];
        g.ColumnSpacing = 6;
    end

    function commonGroups(g)
        %COMMONGROUPS  The always-visible ribbon column, built once.
        %   Frequency, Save/Load and the exporters belong to no single
        %   task: hunting for the right tab to save a session, or to
        %   change the design frequency, is exactly the friction a ribbon
        %   is meant to remove. They live beside the tab group rather
        %   than inside it so they never move and never duplicate.
        g.ColumnWidth = {300,208,236};
        g.Padding = [0 0 0 0];
        g.ColumnSpacing = 6;

gCfg = ribbonPanel(g,'CONFIGURATION');
gCfg.RowHeight  = {22,22,22,22};
gCfg.ColumnWidth = {132,'1x'};
gCfg.RowSpacing = 2;

lblF = uilabel(gCfg,'Text','Design frequency','FontSize',11);
lblF.Layout.Row = 1; lblF.Layout.Column = 1;
% Spinner + unit selector share the row. The unit CONVERTS rather than
% reinterprets: switching GHz -> MHz turns 10 into 10000, it does not
% turn a 10 GHz design into a 10 MHz one. A selector that silently
% redesigned the array would be a trap.
gFreq = uigridlayout(gCfg,[1 2]); gFreq.Padding=[0 0 0 0]; gFreq.ColumnSpacing=4;
gFreq.ColumnWidth = {'1x',62};
gFreq.Layout.Row = 1; gFreq.Layout.Column = 2;
spFreq = uispinner(gFreq,'Limits',[0.001 1e6],'Value',S.freqGHz,'Step',0.5, ...
    'ValueChangedFcn',@(s,e)assignFreq(s.Value), ...
    'Tag','spFreq','Tooltip','Design frequency in the selected GHz or MHz unit; fixes physical spacing and the frozen-phase reference.');
ddFreqUnit = uidropdown(gFreq,'Items',{'GHz','MHz'},'Value',S.freqUnit, ...
    'ValueChangedFcn',@(s,e)assignFreqUnit(s.Value), ...
    'Tag','ddFreqUnit','Tooltip','Display frequencies in GHz or MHz without changing the physical design.');

lblFreqOp = uilabel(gCfg,'Text','Operating frequency','FontSize',11);
lblFreqOp.Layout.Row = 2; lblFreqOp.Layout.Column = 1;
spFreqOp = uispinner(gCfg,'Limits',[0.001 1e6],'Value',S.freqOpGHz,'Step',0.5, ...
    'ValueChangedFcn',@(s,e)assignFreqOp(s.Value), ...
    'Tag','spFreqOp','Tooltip','Operating frequency in the selected GHz or MHz unit; determines wavelength during propagation.');
spFreqOp.Layout.Row = 2; spFreqOp.Layout.Column = 2;
spFreq.Tooltip = 'Reference wavelength for physical element positions and spacing.';
spFreqOp.Tooltip = ['Frequency used for propagation and pattern evaluation. With Beam squint OFF, ' ...
    'steering phases are recalculated at this frequency; with Beam squint ON, ' ...
    'steering phases remain fixed at the design frequency. Imported element data must match the operating frequency.'];
lblLambdaMM = uilabel(gCfg,'Text','','FontSize',10);
lblLambdaMM.Layout.Row = 3; lblLambdaMM.Layout.Column = [1 2];
mutedLbls(end+1) = lblLambdaMM;

% Ticked = the phase values stay as they were computed at the DESIGN
% frequency, so moving the operating frequency walks the beam off the
% commanded angle. Clear = the phases are recomputed at the operating
% frequency and the generated steering term keeps the commanded aim.
%
% Value is ~S.retunePhase: the field says "do retune", the box says "let
% it squint". Spelled out rather than storing a second flag, because two
% fields for one piece of state is how they drift apart.
cbSquint = uicheckbox(gCfg,'Text','Beam squint: keep phases fixed', ...
    'Value',~S.retunePhase,'FontSize',11, ...
    'ValueChangedFcn',@(s,e)assignSquint(s.Value), ...
    'Tag','cbSquint');
cbSquint.Layout.Row = 4; cbSquint.Layout.Column = [1 2];
cbSquint.Tooltip = ['Checked: steering phases stay fixed at the design frequency; ' ...
    'changing the operating frequency can shift the beam.' newline ...
    'Unchecked: phases are recalculated at the operating frequency to ' ...
    'keep the beam aimed at the selected direction.'];

% AFTER lblLambdaMM exists, because applyFreqUnit writes to it. Called
% earlier, that write landed on the placeholder [] this handle is
% pre-declared as -- and MATLAB turns `x.Text = ...` on an empty x into
% a STRUCT rather than erroring, so the assignment was silently thrown
% away when the real label replaced it a line later. Nothing broke only
% because the start-up block sets the text again further down; remove
% that and the readout would go blank with no error anywhere.
applyFreqUnit();          % put both spinners into the startup unit

        gf = ribbonPanel(g,'FILE');
        gf.ColumnWidth = {'1x','1x','1x'};
        uibutton(gf,'Text','New','Icon',ribbonIconFile('new'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnNew', ...
            'Tooltip','Start a fresh design; ask to save unsaved changes first.', ...
            'ButtonPushedFcn',@(s,e)onCloseRequest(true));
        uibutton(gf,'Text','Save','Icon',ribbonIconFile('save'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnSave', ...
            'Tooltip',['Save the design to its file (.mat); the first ' ...
            'save asks for a name. File ▸ Save As… saves under a new name.'], ...
            'ButtonPushedFcn',@(s,e)saveConfig());
        uibutton(gf,'Text','Open','Icon',ribbonIconFile('open'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnLoad', ...
            'Tooltip',['Open a design saved as .mat. File ▸ Open Recent ' ...
            'lists the last ones.'], ...
            'ButtonPushedFcn',@(s,e)loadConfig());

        ge = ribbonPanel(g,'EXPORT');
        ge.ColumnWidth = {'1x','1x','1.3x'};
        uibutton(ge,'Text','CSV','Icon',ribbonIconFile('doc'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnCSV', ...
            'Tooltip',['Element table: x, y (λ), amplitude, phase and ' ...
            'rotation (°) of every element, as CSV'], ...
            'ButtonPushedFcn',@(s,e)exportCSV());
        uibutton(ge,'Text','CST .tsv','Icon',ribbonIconFile('grid'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnTSV', ...
            'Tooltip',['CST array import (.tsv): element positions (m), ' ...
            'amplitudes and phases (°) for a CST array task'], ...
            'ButtonPushedFcn',@(s,e)exportTSV());
        % The macro used to be reachable only from inside the Phase
        % scheme map, three windows deep; it is the last step of the CST
        % round trip, so it sits with the other exports as well.
        uibutton(ge,'Text','CST macro','Icon',ribbonIconFile('macro'), ...
            'IconAlignment','top','FontSize',11,'Tag','btnCstMacro', ...
            'Tooltip',['CST VBA macro (.bas): defines every port''s ' ...
            'amplitude and phase (°) as a Combine Results excitation. ' ...
            'Same export as in the Phase scheme map.'], ...
            'ButtonPushedFcn',@(s,e)exportCstMacro());
    end

    function selectTask(name)
        %SELECTTASK  Switch design settings or an already-open result.
        k = find(strcmp(name, taskNames), 1);
        if isempty(k), k = 1; end
        if ~isempty(elementGalleryPopup) && isgraphics(elementGalleryPopup)
            elementGalleryPopup.Visible = 'off';
        end
        if ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup)
            shapeGalleryPopup.Visible = 'off';
        end
        curTask = taskNames{k};
        isResult = k > numel(designTaskNames);
        if isResult
            % Analysis tabs change only the main workspace. Their source
            % ribbon and settings remain in view, as in the designer UI.
            lastDesignTask = resultOwners(k-numel(designTaskNames));
        else
            lastDesignTask = k;
        end
        for m = 1:numel(ribbons)
            onOff = ternStr(m == lastDesignTask,'on','off');
            ribbons(m).Visible  = onOff;
            sections(m).Visible = ternStr(m == lastDesignTask,'on','off');
        end
        if ~isempty(resultWorkspace) && isgraphics(resultWorkspace)
            resultWorkspace.Visible = ternStr(isResult,'on','off');
            mid.Visible = ternStr(~isResult,'on','off');
            pPat.Visible = ternStr(~isResult,'on','off');
            for j = 1:numel(resultPanels)
                if ~isempty(resultPanels{j}) && isgraphics(resultPanels{j})
                    resultPanels{j}.Visible = ternStr(isResult && j == k-4,'on','off');
                end
            end
        end
        paintResultTabs(tabStrip,tabBtns,closeBtns,k);
        % The View menu ticks the same task. Only its tab items: other
        % features may add View items of their own.
        for hItem = reshape(mView.Children,1,[])
            if startsWith(hItem.Tag,'menuTab')
                hItem.Checked = strcmp(hItem.Tag, ...
                    ['menuTab' taskNames{lastDesignTask}]);
            end
        end
    end

    function pane = resultPane(name,tag)
        idx = find(strcmp(name,resultTaskNames),1);
        pane = resultPanels{idx};
        if isempty(pane) || ~isgraphics(pane)
            pane = makeResultPane(resultHost,tabStrip,tabBtns(idx+4), ...
                closeBtns(idx),numel(designTaskNames),idx, ...
                resultTaskWidths{idx},tag);
            resultPanels{idx} = pane;
        end
        selectTask(name);
    end

    function closeResult(name)
        idx = find(strcmp(name,resultTaskNames),1);
        if isempty(idx) || isempty(resultPanels{idx}), return; end
        pane = resultPanels{idx};
        if idx == 2 && isgraphics(pane)
            mapSlot = pane.UserData;
            if isnumeric(mapSlot) && isscalar(mapSlot) && ...
                    mapSlot <= numel(S.portMapRefresh)
                S.portMapRefresh{mapSlot} = [];
            end
        elseif idx == 7
            cstComparisonClosed = true;
        end
        removeResultPane(pane,tabStrip,tabBtns(idx+4),closeBtns(idx), ...
            numel(designTaskNames),idx);
        resultPanels{idx} = [];
        if strcmp(curTask,name)
            selectTask(designTaskNames{resultOwners(idx)});
        end
    end

    function pickShape(name)
        %PICKSHAPE  Ribbon shape gallery -> the same path as the dropdown.
        %   Routed through ddShape rather than calling assignShape
        %   directly so the dropdown, the gallery highlight and S all
        %   move together; assignShape is also reached from config load
        %   and the dropdown itself, and those must end up in the same
        %   state as a gallery click.
        if ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup)
            shapeGalleryPopup.Visible = 'off';
        end
        ddShape.Value = name;
        assignShape(name);
    end

    function updateLayoutOption(field,value)
        layoutView.(field) = logical(value);
        refreshLayout();
    end

    function syncShapeGallery()
        %SYNCSHAPEGALLERY  Highlight whichever gallery button is active.
        %   uibutton has no "selected" state, so selection is shown by
        %   background colour (see paintSelected).
        if isempty(shapeBtns) || ~isgraphics(shapeBtns(1)), return; end
        activeShape = find(strcmp({shapeBtns.Text}, S.arrayShape), 1);
        paintSelected(shapeBtns, activeShape);
        % Pale dots remain legible on the selected blue tile; the other
        % shapes retain blue dots against the theme's dark button fill.
        for kGallery = 1:numel(shapeBtns)
            if isequal(kGallery,activeShape)
                iconName = ['shp' shapeNames{kGallery} 'Active'];
            else
                iconName = ['shp' shapeNames{kGallery}];
            end
            shapeBtns(kGallery).Icon = ribbonIconFile(iconName);
        end
        if strcmp(shapePreviewPreviousModel,S.arrayShape)
            shapePreviewPreviousModel = 'Circle';
            if strcmp(S.arrayShape,'Circle'), shapePreviewPreviousModel = 'Custom'; end
        end
        if any(strcmp(shapePreviewOlderModel, ...
                {shapePreviewPreviousModel,S.arrayShape}))
            available = shapeNames(~ismember(shapeNames, ...
                {shapePreviewPreviousModel,S.arrayShape}));
            shapePreviewOlderModel = available{1};
        end
        shapePreviewOlder.Text = shapePreviewOlderModel;
        shapePreviewOlder.Icon = ribbonIconFile(['shp' shapePreviewOlderModel]);
        shapePreviewOlder.UserData = shapePreviewOlderModel;
        shapePreviewPrevious.Text = shapePreviewPreviousModel;
        shapePreviewPrevious.Icon = ribbonIconFile(['shp' shapePreviewPreviousModel]);
        shapePreviewPrevious.UserData = shapePreviewPreviousModel;
        shapePreviewCurrent.Text = S.arrayShape;
        shapePreviewCurrent.Icon = ribbonIconFile(['shp' S.arrayShape 'Active']);
        shapePreviewCurrent.UserData = S.arrayShape;
        paintSelected([shapePreviewOlder shapePreviewPrevious shapePreviewCurrent],3);
    end

    function paintSelected(btns, kOn)
        %PAINTSELECTED  Accent the kOn-th button; give the rest back to the theme.
        %   The rest are NOT painted with a remembered "theme colour":
        %   reading BackgroundColor at start-up returns the light value
        %   whatever theme is showing, so the unselected tabs and gallery
        %   buttons came out as white blocks on the dark ribbon, and an
        %   explicit colour would stay put through a theme change anyway.
        %   ColorMode 'auto' hands them back to the theme outright. The
        %   selection still shows as a DIFFERENT BackgroundColor, which
        %   is what the gallery regression test (R10) checks.
        P = pal();
        for kBtn = 1:numel(btns)
            if isequal(kBtn, kOn)
                btns(kBtn).BackgroundColor = P.accent;
                btns(kBtn).FontColor = P.accentText;
            else
                try
                    btns(kBtn).BackgroundColorMode = 'auto';
                    btns(kBtn).FontColorMode = 'auto';
                catch   % releases without ColorMode: nearest fixed match
                    btns(kBtn).BackgroundColor = P.btnBg;
                    btns(kBtn).FontColor = P.btnFg;
                end
            end
        end
    end

    function P = pal()
        P = paletteForFigure(fig);
    end

    function paintChrome()
        %PAINTCHROME  Apply the palette to the window's own furniture.
        %   Tab and gallery highlights, divider rules and caption labels
        %   -- the explicitly coloured parts that are not plots. Run once
        %   at start-up and again by applyTheme.
        P = pal();
        paintResultTabs(tabStrip,tabBtns,closeBtns, ...
            find(strcmp(taskNames, curTask), 1));
        syncShapeGallery();
        syncElementGallery();
        for kSep = 1:numel(sepList)
            if isgraphics(sepList(kSep).h)
                sepList(kSep).h.BackgroundColor = P.rule;
            end
        end
        live = mutedLbls(isgraphics(mutedLbls));
        if ~isempty(live), set(live,'FontColor',P.muted); end
        styleElementTable();
        setStatus(statusMsg, statusTone);   % same text, this theme's colour
    end

    function matchTheme(hWin)
        %MATCHTHEME  Put a popup window on the main window's theme.
        %   A new figure takes the DESKTOP theme, which need not be the
        %   one this window shows (theme(fig,...) can differ), while
        %   pal() answers for the main window -- a popup left alone got
        %   colours meant for the other theme. Only re-themed when the
        %   two differ, so a popup that already follows the desktop in
        %   step with the main window keeps doing so. Remembered so
        %   applyTheme can bring it along on a later theme change.
        if isempty(hWin) || ~isgraphics(hWin), return; end
        try
            want = fig.Theme.BaseColorStyle;
            if ~strcmpi(hWin.Theme.BaseColorStyle, want)
                theme(hWin, want);
            end
        catch
            % No figure themes on this release: nothing to match.
        end
        themedPopups = themedPopups(isgraphics(themedPopups));
        if ~any(themedPopups == hWin), themedPopups(end+1) = hWin; end
    end

    function applyTheme()
        %APPLYTHEME  fig.ThemeChangedFcn: re-apply every explicit colour.
        %   Popups follow the main window, open phase maps redraw with
        %   the new palette, and the pattern plots are RECOLOURED in place
        %   rather than recomputed: a theme change does not touch the
        %   design, so a pattern evaluation (seconds on a large array)
        %   would be pure cost -- and with Auto off it would not run at
        %   all, leaving the old theme's colours on the plots.
        paintChrome();
        live = themedPopups(isgraphics(themedPopups));
        for kWin = 1:numel(live), matchTheme(live(kWin)); end
        refreshPhaseMaps();
        refreshLayout();
        recolorPlots();
        repaintResultTones();
    end

    function recolorPlots()
        %RECOLORPLOTS  Re-apply palette colours to what is already plotted.
        %   Every object drawn with a pal() colour on the cut axes or in
        %   a popup carries a Tag naming its ROLE; this table maps role ->
        %   palette entry for lines/markers and for text. A NEW plotted
        %   object that takes a pal() colour needs a Tag and a row here,
        %   or it keeps the old theme's colour until it is redrawn.
        P = pal();
        roles = { ...   % Tag            line/marker     text
            'cutTotal',    'traceTotal',   ''
            'cutAF',       'traceAF',      ''
            'cutEF',       'traceEF',      ''
            'cutAR',       'traceAR',      ''
            'cutRef',      'refTrace',     ''
            'cutCst',      'overlayTrace', ''
            'coverageGrating','bad', ''
            'coverageRing','gridLines',''
            'coverageContour','ink',''
            'cutARLimit',  'muted',        ''
            'peakLine',    'muted',        ''
            'gratingWarn', 'bad',          ''
            'cutProbe',    'good',         'good'
            'plotCaveat',  '',             'muted'
            'scanCosRef',  'traceAR',      ''
            'scanTotal',   'traceTotal',   ''
            'scanRef3dB',  'muted',        ''
            'scanNow',     'ink',          ''
            'pattern2DTotal','traceTotal', ''
            'pattern2DRHCP','traceTotal', ''
            'pattern2DLHCP','overlayTrace', ''
            'overlayMatlab','traceTotal', ''
            'overlayCst','overlayTrace', ''
            'bandPS',      'traceTotal',   ''
            'bandPeakPS',  'traceTotal',   ''
            'bandTTD',     'traceAR',      ''
            'bandLaw',     'refTrace',     ''
            'bandRefLine', 'muted',        ''};
        % peakLine and gratingWarn are xline markers: their labels take
        % the line colour. ax3D is searched for its caveat subtitle.
        scopes = [{axCut, ax3D, resultWorkspace}, ...
            num2cell(themedPopups(isgraphics(themedPopups)))];
        for kRole = 1:size(roles,1)
            objs = gobjects(0);
            for kScope = 1:numel(scopes)
                objs = [objs; ...
                    findall(scopes{kScope},'Tag',roles{kRole,1})]; %#ok<AGROW>
            end
            for kObj = 1:numel(objs)
                o = objs(kObj);
                if isa(o,'matlab.graphics.primitive.Text')
                    if ~isempty(roles{kRole,3}), o.Color = P.(roles{kRole,3}); end
                elseif ~isempty(roles{kRole,2})
                    lc = P.(roles{kRole,2});
                    if isprop(o,'Color'), o.Color = lc;
                    elseif isprop(o,'LineColor'), o.LineColor = lc; end
                    if isprop(o,'Marker') && ~strcmp(o.Marker,'none')
                        o.MarkerEdgeColor = lc;
                        if ~ischar(o.MarkerFaceColor), o.MarkerFaceColor = lc; end
                    end
                end
            end
        end
    end

    % ------------------------------------------ status bar, busy, progress
    function setStatus(msg, tone)
        %SETSTATUS  Show one line in the status bar.
        %   tone '' (default) is a plain message. 'good' (a finished
        %   save or export), 'warn' and 'bad' take the palette's colour
        %   for it and a leading symbol, so the state does not rest on
        %   colour alone. Nothing clears the line on a timer: a message
        %   stays until the next one replaces it, so it can still be read
        %   after the fact.
        if nargin < 2, tone = ''; end
        statusMsg = char(msg); statusTone = char(tone);
        if isempty(statusBar) || ~isgraphics(statusBar), return; end
        P = pal();
        switch statusTone
            case 'good', txt = ['✓ ' statusMsg]; col = P.good;
            case 'warn', txt = ['⚠ ' statusMsg]; col = P.warn;
            case 'bad',  txt = ['✕ ' statusMsg]; col = P.bad;
            otherwise,   txt = statusMsg;        col = P.muted;
        end
        statusBar.Text = txt; statusBar.FontColor = col;
    end

    function tok = beginBusy(msg)
        %BEGINBUSY  Watch pointer and a status message until tok dies.
        %   Keep tok in a variable of the calling function. It is an
        %   onCleanup object, so the busy state ends when that function
        %   returns OR throws: an error in the wrapped code cannot leave
        %   the window stuck on a watch pointer. endBusy(tok) ends it
        %   early. (An unassigned call therefore ends at once.)
        busyDepth = busyDepth + 1;
        prevMsg = statusMsg; prevTone = statusTone;
        % The token exists before anything below can throw, so even a
        % failure in here gives the depth count and the pointer back.
        tok = onCleanup(@() finishBusy(msg, prevMsg, prevTone));
        setStatus(msg);
        if isgraphics(fig)
            fig.Pointer = 'watch';
            % nocallbacks: puts the pointer and the message on screen,
            % but a click queued meanwhile waits until this work is done
            % instead of starting a second calculation inside it.
            drawnow nocallbacks;
        end
    end

    function endBusy(tok)
        %ENDBUSY  End a busy section before its token goes out of scope.
        if ~isempty(tok) && isvalid(tok), delete(tok); end
    end

    function finishBusy(msg, prevMsg, prevTone)
        % The token's cleanup. The pointer comes back only when the
        % OUTERMOST section ends; the previous message only if the line
        % still shows this section's own one, so a result the section
        % reported ("Pattern computed in 0.4 s") stays on screen.
        busyDepth = max(busyDepth - 1, 0);
        if ~isgraphics(fig), return; end      % window closed meanwhile
        if busyDepth == 0, fig.Pointer = 'arrow'; end
        if strcmp(statusMsg, msg) && isempty(statusTone)
            setStatus(prevMsg, prevTone);
        end
    end

    function [d, closer] = startProgress(titleTxt, msg, cancelable, host)
        %STARTPROGRESS  Progress dialog for a long analysis, or [].
        %   Opens indeterminate; stepProgress turns it into a bar. Keep
        %   the second output in a variable of the caller: like
        %   beginBusy's token it closes the dialog when the caller returns
        %   or throws, so an error cannot leave a modal dialog covering
        %   the window. closeProgress(d) closes it earlier -- before an
        %   alert, which must not open underneath it. uiprogressdlg
        %   refuses a figure that is not visible, so then d is [] and the
        %   analysis runs without one; the other helpers accept [].
        %   host (default: the main window) is the window it covers -- an
        %   analysis run from its own popup shows its progress there.
        d = [];
        if nargin < 3, cancelable = false; end
        if nargin < 4, host = fig; end
        if isgraphics(host) && strcmp(char(host.Visible),'on')
            try
                d = uiprogressdlg(host,'Title',titleTxt,'Message',msg, ...
                    'Indeterminate','on','Cancelable',cancelable);
            catch
                d = [];   % no dialog is better than no analysis
            end
        end
        closer = onCleanup(@() closeProgress(d));
    end

    function stop = stepProgress(d, frac, msg)
        %STEPPROGRESS  Show fraction frac (0-1) done; true if cancelled.
        %   A closed window counts as Cancel. Closing the main window
        %   mid-sweep deletes it and the dialog with it (the native close
        %   button is outside the dialog's modality); read as "not
        %   cancelled", the sweep ran on to the end for a window that no
        %   longer existed, then failed drawing its result or opened an
        %   orphan result window. The same holds for a popup hosting the
        %   dialog, so a deleted dialog is a Cancel too.
        stop = ~isgraphics(fig) || (~isempty(d) && ~isvalid(d));
        if stop || isempty(d), return; end
        try
            d.Indeterminate = 'off';
            d.Value = min(max(frac,0),1);
            if nargin >= 3, d.Message = msg; end
            % limitrate: at most ~20 redraws a second however fast the
            % loop turns, while still taking in a click on Cancel -- or
            % the window close, which is handled right here too.
            drawnow limitrate;
        catch progErr
            % The close can also land inside one of the assignments above
            % (setting Value flushes the dialog), and the next one then
            % throws "Invalid or deleted object". That is the same Cancel;
            % any other error is the analysis's own and goes on up.
            if isgraphics(fig) && isvalid(d), rethrow(progErr); end
        end
        stop = ~isgraphics(fig) || ~isvalid(d) || d.CancelRequested;
    end

    function closeProgress(d)
        %CLOSEPROGRESS  Close a startProgress dialog; [] or closed is fine.
        if ~isempty(d) && isvalid(d), close(d); end
    end

    function g = ribbonPanel(parent, titleTxt)
        %RIBBONPANEL  Empty ribbon group: content area over a grey caption.
        %   Returns the content grid, which the caller reshapes -- a row
        %   of buttons for an action group, or label/field pairs for a
        %   parameter group like CONFIGURATION. Splitting this out of
        %   ribbonSection is what lets the two kinds of group share the
        %   same frame and caption styling.
        pn = uipanel(parent);
        gv = uigridlayout(pn,[2 1]);
        gv.RowHeight = {'1x',13};
        gv.Padding = [3 2 3 2];
        gv.RowSpacing = 2;
        g = uigridlayout(gv,[1 1]);
        g.Padding = [0 0 0 0];
        g.ColumnSpacing = 3;
        P = pal();
        mutedLbls(end+1) = uilabel(gv,'Text',titleTxt, ...
            'HorizontalAlignment','center','FontSize',9,'FontColor',P.muted);
    end

    function p = ribbonIconFile(name)
        % Icons are embedded so this .m file runs without companion PNGs.
        p = embeddedRibbonArtwork(name);
        if ~isempty(p), return; end
        d = fullfile(tempdir,'padRibbonIcons_v26_generated');
        if ~isfolder(d), mkdir(d); end
        p = fullfile(d,[name '.png']);
        if isfile(p), return; end
        [RGB,A] = ribbonGlyph(name);
        if strcmp(name,'efPatch')
            % Match the reference thumbnails' 95 x 75 px white tile.
            rows = round(linspace(15,114,75));
            cols = round(linspace(1,128,95));
            imwrite(RGB(rows,cols,:),p);
        else
            imwrite(RGB,p,'Alpha',A);
        end
    end

    function p = embeddedRibbonArtwork(name)
        % Decode the supplied ribbon artwork on demand into MATLAB's temp
        % folder; uibutton Icon requires a file to preserve transparency.
        switch name
            case 'import'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAACoAAAA0BAMAAAAOMbeDAAAAIGNIUk0AAHomAACAhAAA+gAAAIDoAAB1MAAA6mAAADqYAAAX' ...
                    'cJy6UTwAAAASUExURdb9xNb9xNb9xNb9xNb9xP///2mLDuIAAAAEdFJOUwCAv0BHJ479AAAAAWJLR0QF+G/pxwAAAAd0SU1F' ...
                    'B+oJHg05OgaNF8IAAAAldEVYdGRhdGU6Y3JlYXRlADIwMjYtMDktMzBUMTM6NTc6NTgrMDA6MDA79YCrAAAAJXRFWHRkYXRl' ...
                    'Om1vZGlmeQAyMDI2LTA5LTMwVDEzOjU3OjU4KzAwOjAwSqg4FwAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0wOS0zMFQx' ...
                    'Mzo1Nzo1OCswMDowMB29GcgAAAClSURBVDjLzZPRDYQwDENTjgGK1AHQiQH4yAAFZf+ZaHUcDYnFHwJ/PoEdnED0EgX5DgCz' ...
                    'yALwBHEvItHRT6EyO5wKHR3tCl2xhXeuNKPZwHAMjRnOxjCufp2f+ILm+2kQo3hUoLSol5uy6rbpX0YCBiavNTSBR/caT1km' ...
                    'TxefgMHvHEyWyjNbZpelRrZLTsBgz3ObJ3hRNS97GtCxlzz0G1EX6Rlth0RHmCXQ8ugAAAAASUVORK5CYII=' ...
                    ];
            case 'polar2d'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAIAAAACACAYAAADDPmHLAAAAIGNIUk0AAHomAACAhAAA+gAAAIDoAAB1MAAA6mAAADqYAAAX' ...
                    'cJy6UTwAAAAGYktHRAD/AP8A/6C9p5MAAAAHdElNRQfqCR4OGA/e53hbAAAAJXRFWHRkYXRlOmNyZWF0ZQAyMDI2LTA5LTMw' ...
                    'VDE0OjIxOjI2KzAwOjAwOFmRcwAAACV0RVh0ZGF0ZTptb2RpZnkAMjAyNi0wOS0zMFQxNDoyMToyMSswMDowMIyjF0EAAAAo' ...
                    'dEVYdGRhdGU6dGltZXN0YW1wADIwMjYtMDktMzBUMTQ6MjQ6MTUrMDA6MDBHX94qAAAbEklEQVR42u2debiVVb3HPxyQ8Ygg' ...
                    'IILM4JQoYGpqDuW9OSBF4mxaOKTpTa83GzTvzdnUtK6W1bVMLc0sUXJItHAeSwlBDRWZQRERREDGw/3js97n3Wez9373Ye93' ...
                    '730OfJ9nP2fY+93rXWt9129e62112rhz2YLNF3XVvoEtqC62EGAzxxYCbObYQoDNHG2qfQNVQiugO9AB2ACsApYA66p9Y5XG' ...
                    '5kaAgcARwP7AMKAb0ICT/wbwOPAQML/aN1opbC4E2A44FRgL7AIsA2YAU8L7vYDDgGOBc4HfAbcB71f7xtPG5kCAEcBVuPJn' ...
                    'AhcDjwKzgBXhM/XAIGAk8FXgGuDzwPeAV6vdgTTR0glwKHAzMAD4LfBDYFqOz60CPgD+Dvwe+DZKiz+iRHis2h1JCy3ZC9gP' ...
                    'J78PcClwFrknPxtv4aT/AOgP/BxthhaJlkqAPsCPcOVfgiJ9VROuXx2uj0hwXfjOFoeWSIC2wIXAZ4GbgBuB9ZvwPevDtdeH' ...
                    '77ocaFftzpUbLZEAI4HTgCeBa3E1bypWAz8GJgInAaOq3blyo6URoBdwEbASuILyuHGLgCuRDOdiAKnFoKUR4FhgH+Au4Iky' ...
                    'fu+zGBs4EDih2p0sJ1oSAXphsOc9DOJsKON3rwN+gxLlZIwgtgi0JAJ8DtgDuJc4wldOTAEmAHtjkKhFoKUQYCvgy+jq3Yfx' ...
                    '/XJjHXBHaOOo0GazR0shwC7AQRjJeznFdv4ZXocCu1a70+VASyHAvsD2wF+Bj1Ns5yPgYfQE9q12p8uBlkCA1sDB6Po9U4H2' ...
                    'nsYk0iG0gFxKSyBAL3T9plNcrL9UvIn5gr2AntXufKloCQTYCeiL+v+DCrS3GJgE9AaGVLvzpaKlEKA98Arl9f3zYQMagh2A' ...
                    '3ard+VLREgiwB+r/NyrY5pTQ5qeq3flS0dwJ0BElwIfAnAq2Oy+0uUu4h2aL5k6AbYB+wFyckEphCbAAGAxsW+1BKAXNnQDd' ...
                    'MS4/D1hewXY/RtJtwxYCVBU9gU6YAEoj/JsP61HltAO6VnsQSkFzJ0A3rACaW4W2F6D+36Hag1AKmjsBtsNdPouq0HZUbLJF' ...
                    'AlQR0c6epVVoexmqgs7VHoRS0NwJ0BVYS7zBo5JYgSniZl0i1twJsBVOwidVaPsTtkiAqqMzEmBlFdpei+qnU7UHoRQ0dwK0' ...
                    'xklYU4W2VyAJWlV7EEpBcydAhGpMwjoqG3tIBS2FAJXIAmajDS1g/Jp9B9AQ2+xO9igXWgIBWmM0sNKoo5nr/6gTzRnrURRX' ...
                    'Y9NmPbqhm7LxtGbQ3AmwLPShGgRogxJgWbUHoRQ0dwKsxInoUIW2O6IEaMq5AzWH5k6A5TgJ1QjGtMfx+6jag1AKmjsBPg59' ...
                    'qK9C2/VogC6t9iCUguZOgMXhZzVSsl3Cz2YtAWptZ0svYHfUq6+RXOe3FKNx1diuvT3Fp6K7hn61A6ZiBVNNoJYkwO7AH4AH' ...
                    'cRv2eODfEq5Zgid3VIMAXUPbSxI+dxAwLvTpQdy9PKwK95sTtUKADniA42eBu/F8vuHh968UuG4xpmV7VrgvdUCP0PbiAp8b' ...
                    'i+cV7Bl+3o3nC3yLGtleXisqYHfgcDyQ8Zvo3j2ABzTdgGL24RzXLUN1sR26ZZWqDG6P1cDLyW8DjMaj5j4BzsAziDthDeFI' ...
                    '3F6exkEWTUKtSIChWGI9EQe1AUXlWWhp34QSIRsrcQVui5NSKXRE0i0hdy3CcCRvK+AcJHMDei0TsIpo9wreb17UCgH64ADN' ...
                    'zPr/o7iK+uOBj9nVN5EI7k5lK3M6hzbfZ2MC1KM6648HVD6U9f4MDB/3q+D95kWtEKAbDkqu0q6fof78EjAm670VWJ9fT2Vr' ...
                    '87oS70fIvueT8QiZ+4Ff5rh2NZK9JraW1woBNoRXrgKLlXha5wI8A3Bw1vsL0IjsUcH73YGYAJkYAPwHblO/htw2SUPoa01k' ...
                    'EmuFAMvQIM2nx18GfoEbQbO9ggXh2oEVvN/eaJvMyPr/WLRnbsUzBHKhQ7i2GnsZNkKtEGBJuJdC/vw9aCOcgmf7R5iBBSF9' ...
                    'K3i/Q3AlL8j43zDgdAxg3UL+KqXtkACV3MyaF7VCgHk4YIUMo3fwzP8h+FCHCJEhNggHNm3Uhfv8ONw3KM7HojH7K2B2geuj' ...
                    'CGJNRANrhQCzUV8OpLBuvBOlwFeJj2d5F1fTDlQmKVSPE/0R8fawQXhO4WtosOZDHaqx5VRnP2POG6oFzAMWYnBkmwKfm45S' ...
                    'YCBwfPjfUvQEegFbV+Bet0UCzCc+ku6LaACOo7FayEYXPFRiIbH0qCpqhQCL8fStgSTvtr0LV88YtPyXo1TohsZZ2ugT2pqJ' ...
                    'UqAX8DUU6fclXNsT1cd0knMIFUGtEGA1MBlXV9LBSzMwYjic+Mzet1A096/AvfbGSODb4e+D0QC8F1VAIeyO5Pk7NVJJVCsE' ...
                    'APgHGkd7J3xuPZ4GvgI4EV3At9CIrIQr2D+09SYmdE7CyfwzyRtF9gnXTqJGUEsEeB3F6P4kh3VfBl5ECTACpcJyPLUrzQBL' ...
                    'HfGpZDPQ5z8QeD7cTyF0Dn2bT7KkqBhqiQCzcWJ3A3ZM+OxKTK12xqzbQvQEdiJdQ3BrjER+iGQdg4bdeJIzkUPQyH2ZGjEA' ...
                    'obYIsBbP+t0G6wKS8BRa/yNxY8gsFM9phoS7YsBpevj7yHAPDxdx7eeRLM9Snc2sOZEmAeow6NGX4v3zJ3B1jSK50ndO+PxQ' ...
                    'FMtT0YhMMyI4CAk2DQ2/3TBjmXRGYUesd1iCD7MqBp3QY+hDimotLQL0wyzeM8ALWAp1Hslu2lt45Ote6C8XwjosIWtANTAD' ...
                    'cwlJ6qMU7Ih1fbOAo8P/HiB5d9AuoU+voPGYNHYXYDbxeZR00dNPy440CFAH/BdwNlrq09A4uxFz48eSfy/fCqwK6oqPgEnC' ...
                    'JDwi9iCchDUoDdLC0NDGSuCA0HaS8Qex+J9A/sMs2qJXcx9mP4chqZeEsbyYFCq40iDA9pi7n4R58S+i+LsJxdkd+DDnfImf' ...
                    'TDWQ5A0swke+D8AkywdItjSOb+2IRtxiDOgMwHhE0gnl9fjg6sXh87nQHR9vexsasj/C+MJI4Bh8gPXnQh/LijQI0BXZPgNF' ...
                    '5Sd4uvb5KDYn4cOZf0ruqN9UJME+KDaT8AgGkvZEb2AQ6VQJ90BLfj66nqtR/ydhX3T/niW3+N8Bn098EQbDTgy/T0PPYg4m' ...
                    'wrYlhf0PaRBgFerndjQWWRvQJjgJxdyJmDnLjt6tQb3aEVdAEl7BFRJNSg8ap4vLhX7Ex9KNQKK+knBNK5SA7TGdnV091D+M' ...
                    'wbHoSp6AHkWmTbEV1hCsIYXDsNIgwFLM0A0m90qcA5yJKuEwLPTITgM/ja7WaJKNn6WoBnqin96BdCKCQ5CU9SiKJ5Kc0++D' ...
                    '6u91NOYy0T/0/bAwFmegxMxGl/DZhaSwCykNAixBkT+Y/AbZYuD7aBgehiIw032bhVJgCMU9o28iGpBRnV0a5/gPRynWE1fi' ...
                    'X4u45kjU6Y/QOEvYF+sFDwtj8H3y7y/YGRfBFFLYh5gGARpwQtri49XyYQXwP8j+I9HVySTBeCwVO5HkOMKraJF3D30aQXkN' ...
                    'wU6YyGmNUu2N0GYh1GOkcCnwp4z/90L753B0lf+HwgddjkKp9igpHEaRVhzgOTQCk0T4CrR+70FdeQXxXv8XceXsR3KC6ENU' ...
                    'GxEGU96q2x40ViuRp1II+2Oe4Gni2H8H9IBGIykup/Dk90YPYi7GBMqOtAgwEyd1EMn7+z7CwMeDWPB5ARqPa1EKdMRyq6St' ...
                    'VI8Rx+N7snH1cCkYSBxiXgH8LeHzdbj626Db+wlKxO9hNdOjoZ+LE77nKIw9jKdwmdkmI81Q8L1ouJyOsYFCmI/75d7AgEdU' ...
                    '+fsoRhKPILlO4FWMJIKkKecDnYYSh6bfJFn8D8NVPoX4WYYnoXv3GvCfJJeE9UbD8D3gdlI6kzBNAkzBTZ77UniDZ4TpwHdx' ...
                    'VVyFxt8SJFIPtBMKYRGqnggjKE8MvY7GKuhJkoM/x6EUuj3c18HAleH3C0gOB0ffMTyMYRLhSupcWliPlu5s3CyRFNsHV/wl' ...
                    'OHhXoRv1J6y+OZ3k5/Q9gaoDtJ63oXR0ybj3NSTr4n4YvfsXqsF+wA+RxJeRPxqYiV2Bc1FK3EqKJ5KmnQ6ehiQYgCK+mPP8' ...
                    '7gT+F42/S3G13YV6eEzCtZOIdeUAJFCp6JnxPfNJrub5MhL1L1g0ennoy0/RHkhCW0ycDQrXpLqDuBL1AL/BVXMSxUX2VmMy' ...
                    '5FEstjwDCTAHxWIhe2I+cXKmG0qBUrEr8YOhnqdw6rc3cCrq7bvQeD0ZYwbXUVwdwJdCv59AFZIqKkGARSgC1+OKLiZMuxCN' ...
                    'wXnhZ18MH4/A4Ek+rEMR24Bew4gy3P/eGNZuwEkp5IsfQ3ywxbaozuYCFxLvISiEgRgXWIkucerbxypVETQBgx7DMOpVzF7+' ...
                    'V5Aw3cOgPIcu4xkUToq8TDzYw4tsKx/aE0czF2E1bz70wtX/LrqkF4f7vIziikDboxG8B6rAJ0u476JRKQKsx4jf4yjeivEK' ...
                    'QCPq1+gR7It6dX+MjuXDTEzUgCHhXiXcd09iA3AKZuXy4VCcvPuBL+Dj5W8Ffl9kW6ciuZ8A/o8KnYBeyZrAhcB/Y0z8MnSN' ...
                    'krAK1cdLYYAWYHj4TPJLgRXElnpPSnvC9y7EBHqe/MUcXXHyPkTD7+soia6lOL1/MEqMd8MYVWzn8KYQoC2K5e3Q0GrKMa0v' ...
                    'oFjvgTqumI0cc8Nn22As4A2UBoWkwHNIhI6UdiLXiNC/lRR2/0aFe3oTjbh1qLaKid71Df3rjh5DU0K+bdHW2C6MaZPVXVNK' ...
                    'jPbC/XifQt+2NVrsi0LHX0DWv0NhQ+lOzJBdiL7+mSQ/82cCupPfCW2tx1X2ELm3WE1FVTAU+HToZ1OfKdCa+ByfNzHDmQtd' ...
                    'w72sw5jBThj0KaZYpB4l3IEoLW4v4p4GoGG6F3o5vZGkDegFvYpqZ3IxnSyGAFuheLsYxeE8XJXREa3DUUd/M7z3FBY1PI5i' ...
                    'Pxtr0cjZAwk1HbiawqJyPR66tA+Kyw3EUuB3OT7/PpJxKE5id5q+HTvzIKeXyR/9G4V+/towIU9jejtJh7fBKqkT0La5gfwk' ...
                    '7YI2xigsDesTPvth6GuUmNotfO5klLS3JHWyGAIcij7sYiTC46HBOiRHd0y8HICGzxj0+SchE8ez8UkaC7Fw9H6Mj88imf0L' ...
                    'UUTeg+KujvxSoAGJ+DUsuRpA0wkwBFVUA05qrgmNVn+b8FqEq78Yl+9kXFSvYZAsl94fFMZyFEqydSjd7iQuMfsIF8+GMC6j' ...
                    'UFJegwvyL4VuovWI4z9T6P0OKKL6Ex96GDW4GkX3BxiqnYi++kQ0hPZGv3hkGKjZNC5oWIzq4hCUIK+xMVGyMRulzkEY5++N' ...
                    '0cZc0bJ1SMbuSMaXi5iUTByJNYzzUTznkgCnogqrCxNwPcUFbw7DKN9KrPj9R9b7Q1CiXhfGcC1O+g9QUjwUxjw6pm51mJMl' ...
                    'aDC/iWVmUfxkbb4bSTIC++BETqJxvj0fPkQ35nyUHJfhyriMeG9AZmXro+GzXVBcHZDw/RvQRYqifVvhCszlEWS6g5+maYmh' ...
                    'OiwyBcu5ckX/euPkR1L0SYoT/fujPdMZ/f7M3EBXXL2PZY3doWGciqlDIHznC6gCC9ZFJBGgfXgtpGkPZ2xAa/1SYiJ0wfKn' ...
                    'cWgpR3mB8Siu+qBtkLSxYwEaj0vD3/uRO66Q6Q5GB1EWi22ID6Z8MU/fj6FxkKgY0b8j2jJ9wpjckzHOx2MZ3DVo81wSxu5S' ...
                    '4q1oxWI1EqVjUr+TCLCG+PGom/pYlumhE0ci80fgjp5fhAFsQNF5dRj0n5HsHk7AalpwlZxK7sTPc1iMMZimFYoOCJ9fhSsp' ...
                    'G71RJUbjdxuuzkLoH/q2V+jrjWFsd0Gpdge6rL9FPX45TZ/4CG1xztaSUEmcRICFqJcHESdENhVTMMV5LK6q05Dx56AovwoJ' ...
                    'ciiSo9Aev/VhACNVMByNqmy8gR5Ld5oWD9gjXDMbVUA2jsn4vkmozwuJ/r6hT4ciCa4KfT4H9fnJKK2OR7VSTL1AIfRAd3QO' ...
                    'CVIpiQBLwyDvSHE7dpOwDuv8jsPCiDoMEd+OK+ESzP8fgaKykP6aj5JjefieM9g46reQ2H/fpwn3uW/4OQmjc5kYgLUJdaHt' ...
                    'yNrOh56hL0egyL8s9PV2JE4nLBUbE8ZmLaVj73Cfz5NgMxQTCXwAReHRlO9A5g/CoHwJXcGjQjtHo2v0MK6y6ym83/9BVCeg' ...
                    'mM+WAuuIS7JGED/loxC6EGcRX6BxUKsVGp2R7r8TTwbJh61DH44J93ohqsL7w/8mhL5fT/lKvjuidF0fvr8giiHAi+hLHoUT' ...
                    'Vk5MRv19PorEm3Hb2A1IgpNwxeTbKr4G99FNC3+fzMaVRy+he7QjxdkBA1GSRC5VJoaj7gf1803kD2B1Cvd+UujLjRjnvxXd' ...
                    '6++i8VrM5tKmYBSqmr9SREaxGAKsQL21CgMW5T7lenn4/jHhhs9EAjwWOnEexsrzkeAt4Ce42gdj+VlmgGsGxhu6UVyh6NDw' ...
                    '2bdpbITV4eT3xkn/CZZ95UKncM/nYQXxI0jU0zGQNib8vbTMYzkQ0+2rUb2sSLqg2GTQc2jEfAaZm8bTLl7C1XIFrtYfoG59' ...
                    'CY3HQiT4I3Hs/SvoGkb4kHglF2MHRAWgL9FYf45AsQ2S87d5ro8m/9zwHXPQoh+C1v9XKf+qh7jsfBga08XUHhZNgHUYlRqP' ...
                    'xtb5pEOCRegynhJ+Hxs6NpfCJFiKxtgiGodnIzyLVvqeFN5y3pk4AJQ5Se2ID7h4H1d/rjOBMid/TrhuLNo8pyCp00j1tg5t' ...
                    'nobq5lqKLCRtSjp4SejATAx6nEs6Z/M2INHGoEcwHF2y1hRWBy9g/SG4yyhzQ0pkze9E4RhD//CZBTSu4vkCcUHqLeT2+TPF' ...
                    'fmu0/oehcTsajcU0njMcjcuVqOoupgkHUTe1HmAqPsYlisadR3rPHfoXSpsric/XL0SC9aj3XkRL/lvEVv98jAn0oPAJInuE' ...
                    'z0wlDv92QrukHkmRq1one/JbIZGvQiN3GukgGo+rMdl1Fk3cQ7ApBSFPYwDjvdDwhaT37N4VOLBjiYMjmSTIdhHnExush6Ak' ...
                    'iL4n2jRSKPsVvZcZ/h2NVvUnKFqzff6taTz5hHsdi+RN6+HS7TCWEk3+ORSXr2mETS0JewRTrbNQZ/+E0mrvCqEB/eYTiPfk' ...
                    'RST4ORsHi8aFVxtcudGevufQOh5G7thCPUqANcTh394oSdqhH5/9/J+e4R4yJ/9RDHTdR9OLUIpFXxzzK3AOzsA5aTJKqQl8' ...
                    'GrduP4Wi5/eY6UoLk9GCvhlXY2v0+2+iMflWoRSYj9HL08L/o7DwzuR2ZfvhHoBZxOHfYzCTuACN4MykUN/QzsnhXj4J93I6' ...
                    '6W7m+DTGEs5Gl/JoirT4c6HUotDJKOpuxsH+A4qitM7tfxdX5EXE1vRxbHzUzItYTRxF7gaHa1/BVZvrkW27Y6r6FSTPIOAb' ...
                    '4b1f0dgo7I+2QOQWLkT3+Nvh2jTQBTeV/hmron6K5H6jlC8tR1XwXNRF56G4vgl3xRxIOgccrgltjCW2C47EOEUmCX6NRBiM' ...
                    'LlgDcVg413kDUYzgeTTyzkaJ8AqNHwETHe1yRPh7Oq76n1GeOH426lCy3oVify0usm+zcZ6iyUiqCCoW67Hi5knUuaPDaxuM' ...
                    '1KVhCL2N/v0gnOQdccKiAE70VNFR6No9FQbshDCo44jDuPVYiLEt8RF212Gs4zvEMYFBqPOjyX8EpcQzpIMBWDp3PUYx70b3' ...
                    'ewJlcinLvS9gMro9X0dR+H0MTBRzSuimtncKrsi1WGr1G+JijnEYS+gbBm4eEmenrPvphzmEd5AkZ6Nd8QCx4Tcc8/5HhLZu' ...
                    'w5U/OYV+bY8T/zcsL38/jOlZlNm+KJcEyMQaTME+hEbTQagr/504qvfxJn/7xliBxtBaFO07hjZfxwldgHWJe2AAp2t4/1li' ...
                    'Y+8Q9GrGY4TvCqxZ/CbGAw5GYu2Fxt41SO5yP/VjEOr1a9DgXYbq7nsoXcuuYtIgQIRl4aYnIhH2RoPtMKzUXYKGXDm2QK1D' ...
                    '3T0TdfkQnLS3wj1sjYcybYsFmIcjER8L15+G+YM70KreAxNSd+OKvwWlxlzUvT9Hl7Ic2Cq09w0U9ceH7/4V2lbjSOF0sAhp' ...
                    'EiDCezjQfwsdG0psIwzDyVtE6YcgNmAEbwoO6K5YZPomprM/iyt/GhKkLVY5t0U93wkn+Cso1i9AQ/YXxMe0nYcTUg792x2r' ...
                    'oS/CZNFInOib0KP4E8knkZSMVqeNOzftNhq1h/p2NIq4KCw7GVO/D6JYLlVF7IwFpocTu45tcSXPwMluh5ssNqCBuApFbF/i' ...
                    'nT43oC0wARNgpZZq1ePOqiPDK6qJ/Aeqnz+z6XWAm4RKEyAT3XBVHocrrS+qjX8iGZ5GH3dTn7DZF3XpCWgn/BgNudE46A0Y' ...
                    'yKrDIFZ0vO0E9CS+hSHuqB5/U5/z1x0JeTCGlEdg1nEWSsV7cdt5VZ4illYipxgsRiv7YfStD8bjVfYNv0dHzT+DYdzJaNAV' ...
                    'qyrmol6dgXr7u+F7PkaboA5VwgaM5LVGAkZ59dboEl7bhDbBkqwdcHXvhwbmzkimD9D4vB8N0pmkeP5PMaimBMiFDujTfwYH' ...
                    '7gAs996AtsTraOxNRnE8m+TJaYuu6aUYBVxPTPzoAMeh4efq8N5iTLL8kmRjryO6bbuiobs/ivneoa3ZKFEm4kqfTo08Mg5q' ...
                    'jwCZaI32wp44qHthMKQbDuz7SIhpOJFv42pfRO5SqJEYrcusC4xWX2Y8ZCa6f7n21NVjoGsAxg12wRDyzuH/W+Eqn4p6/VlM' ...
                    'z86jyis9H2qZANmox8kbiittT7Tmu6PuXoWT/w66f2+F36Pa+I/QBrgabY/sMPUGVBEXox3SGUV5/9DOTuinR6eGtkfp8H5o' ...
                    'a3K47lUkUWI9Xi2gmjZAU7EcV9ZU9M87oUgfhBM0DFfjECRIJ5QUy3Hy56M4/gAnJzthtZo4r94Prf9uxKnjiGAzMQT8T2Kp' ...
                    '09StczWD5kSAbKzAwZ9BXCfQGYM90ZbwaOVGzyT+FNoZueoZ22PeYBX64/PQ3ojamBVei0h+RmCzQXMmQC4sC69ZxBVAdbiK' ...
                    't8Z07/aoNgahgRkZau+il7EQpcRHSLKa1N3lQksjQC404GR+RA09sbNWUEtPDt2CKmALATZzbCHAZo4tBNjMsYUAmzn+Hzm9' ...
                    'ehBaFgCcAAAAAElFTkSuQmCC' ...
                    ];
            case 'ef3GPP'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF4AAABLCAIAAAB2q6cRAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAABu9JREFU' ...
                    'eNrtm09sFFUcx79vt91td7cbdkspUhKCBcOfUA3FSE/IgRjTRLwAgZiYeFUTLxiC8VCPopwIYtTERGJMipGI4WIIeAALRWwX' ...
                    'KdDSSmkXoX9m7bSd6ezuzPPwtsN0trPtvJndnZr9pofX2Zk3v/nk+/u9N29myMaNG1HRYiK9vb0tLS3lDsNzSiQShFJa7jA8' ...
                    'Kl+5A/CuKmgsVUFjqQoaS1XQWKqCxlJVZTy3OnVRFR9SOq6JDwFQUQLwK1HZrxfpa40k1UimRuSXmuqk93e9XeLwyjCvSY+c' ...
                    'prRPHR1f9FcdzXX68iTiAKTZ9bMz6wGEA0Igteqbt178H6JRHnWoo30FdkiB3iSaCQ2A8ae7WSM9vEaVgqtDT777YEuxoy1R' ...
                    'rUmPnJauHSrMxaTNZFBvhyOjrBHYMOYPKRPS2jePDf98bXplo2FQsiO/OekkFB4NBETWrmqYAqCQ6JefKx1nu1YqGlXrtAUl' ...
                    'BqK34xDqIRjp5CIOKdWrRdQR7QX/jZ82fX3p9spDk1U7VO2c3aOaDXSMqg6IRuMwOvQ53/nvG1cSGkr7smqHRvsAkCh/P8Zy' ...
                    'A4NxGB1/SKHrfFpt1cHjg7a7LheazDwXAKTO3rExaplT1QFRr8cAAhvGANB1vlnZAf5SosmqHQtO0GQTzcKEyjeOnlYAghvG' ...
                    'AGgN1cWoxy6jUbVO3S/cdJqtizEWplWuJAO3hPWeRqNqnYvWXbtonqc+Ix2TcUxplSs6roNxEY0VFz46RsUh5KeViU5mJvLV' ...
                    'jb88i6bQOO3UOHhg2sFYdFhapYS7XkSjap1Ln8luPaYL6vFu0p1PR29XNUz1iM3eRLP01M7XZI9ODMRUj/OLzqrYs5LvC4vL' ...
                    '77xEaJZjGT46+WlVgI6Urv/2+imvobFxN+Brsjc/tkXn8uQOD6Gh1MY6A5Pf5krL89QXN/xbgM60Ej92ebkWLjoajd7hOMq/' ...
                    '1d7+rdS/JJ2Gxq5wZLT/3y2dCXfoOF3lS2cP8R1Ip6HaHG3/IKpxXiwgPkCb9ZVAJrZaerj53oGWAw7ROHINRzbpInU83jGN' ...
                    'Wa+Q7vzZYENj16WnLjwOcOSawjPg5XZyD9TOsDtEtEGYYx7ApgG6YF7TSFKftH24stEA0JLQkjb2T4EOEU3I2z6ATQB0Rg7p' ...
                    'eAINBx0sfPxgkoA4q0HxYPidnR/zheQoJ53UGpPYfEdL2kiuGMg+6l80v+IQ4mwpQwGlfYRs4wnJrWtzLlIH/xaeG9F91L9r' ...
                    '4dzHqL6+T/niceQazT3X6OKwD4AYSCv1A0iBpghNGXz0QJzZXno0RRKzD522DYgxMt2yc8tRQvm4cniZYoD8Wx09kwBAoiG+' ...
                    'A73omgUX5sBBDuV1NEwMEMDDiCDMd1JHaAjZhiJU4kJnnDcRkJsHLYmJRBvKgMZHtqtwZ8pnFxAMixt0OgeIITPB8pUFjUdE' ...
                    '6iyekcoqhHTi6u3Wdp5uHY1QhGwr6iDFKVlFUqaDM/SxfP6HZGv7F3zdOHWNlryHmjRq/aj1l5+IkKZzqr7hTs8c8QW4+3OK' ...
                    'pqr+SGb0DFIAQGIBAIjzR8OHA4CRCNPjv9M3eydjMc7hCa68y5cefsMUWY5RMawkqwCscBi53Lqt3B9MHj/D/36JG6thoR10' ...
                    'rse4habSAJiVAJCaeUbLhyWrzxqyWhjEolwmU+Lm5jWOrsuVN0CV+1xjgB5ETY7X8q9/UWVk+s9wJtGvTKbECWHaiWXg1qKE' ...
                    'v36fk8PpnMr+nHQiCepIfzrRr0iyMiFMO7QMXHxvOP3oMJVLe5NjgCJNan3DaVHSJFkZeTwBwKFl4OJSVlX9kdJDych0Kpl9' ...
                    'MpLtujtn5OLcMnD3bfPM+Eea0FMaIpKgZiQqStroeFaUNABGLgeO/u4tNChmWmVkCoARAWCEUgwu7qMBkOren5meq64lAKpr' ...
                    'fazBjSMjaxmZMhxMJigA2HjE2s5LTBHRaLMXJq6fMl4MgOpQDlABUswXAEzHMhwATEQASLIymRIlOQ0gGgm+95mbKyTF+rJl' ...
                    '7Ep7ZibjsBNR0sRZTZQ0ExETFMbl3RPn+B6qlBoNgIELrwsTSi70sA9ANLT4gKhfOQNh3JIvExQUwS9FRwPgzx/3Dw2lnPdj' ...
                    'BQWu1t2SogHQ+8vBxK3hUG2QjwUAeU6RZMVEBEA0Emzbu4d7Oab8aJg6T7QNDI4BCNUGdEy1NWZe8pyiQ8lnYdTm5rUHjl4t' ...
                    'aswl/cBQB8StaCS4e+/OXe1nSxBtGT5LPXfyVaoqthhFI0EARU0fT6DRNdB9sudKJ4CnYyIoFWfmh7PIs0Rr27sHQCmJeAKN' ...
                    'x+Whl0i8pgoaS1XQWKqCxlIVNJaqoLFUBY2lKmgs5UskEuWOwYtKJBL/ARk5xEjp9Et7AAAAAElFTkSuQmCC' ...
                    ];
            case 'efCardioid'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF8AAABLCAIAAACZacwvAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAADINJREFU' ...
                    'eNrtnGt0FOUZx593Zmd2Z7K72VtISAIYuUXhJNgYwXIRWw/aqlXQg1ZrT6FHqkcFLXjt5Whtz6GWWqEX2lo4nlpa7yCIWsqh' ...
                    'CAqoVITKPQmE7EJi9pLsZWbn3g/vZjKZ3azZ7CabevL/NJndvPu8v/0/z8x7mUU1NTVbtmypq6uDUfXVkSNHUE1NTUtLS7Ej' ...
                    'GaEiih3AiNYonWwapZNNo3SyaZRONo3SyaZROtlkKXYA4D/XdHDvGxdaPkp2nbFRQCIVn5dVhBDiBVWzVpSNHf+1hQ9Xj580' ...
                    'zLEV7W7w/d1vfbzjD15rGCGtdpIdAFgbwi/FoyoA+Dtl4/tD4RgAdCadV9308Jz5N3w56Zz3n9n8lx8i/mz9pc7pDbWUFWi7' ...
                    'Dwga5G4AUJMBTY4DgMRpEq8BwKlTgpEUxws8L57yJyY33nbjoqWV1TVfEjqdHYEX1z9+zeXC1IZGkuWB4jUyAgBIcQMAqB6k' ...
                    'eAAA5LAaPw4AarI9RSSkcmHFH5SNbgqFY8fOdl1Uv/DOJSvLyqv+v+m88uK6Sd6P62ZP0WznNTLc39sIcSKSJmFGmhLW+I4s' ...
                    'jELhWFtH9yk/d8tdTyy+a/lQhE263e4VK1YMHZdAW/O6n3xz6bKqiqmSZm0Fgs/yZo2MIACkeoBgEMEAoSJAONcoFlEMYVXA' ...
                    'WUJ0disAwDJWikRORjvx2f5d/959aV2js9RT2OCH1jsfvLe9kt1VXcdl8Uu6eh2k8poUAAA11owZAUC3X5Z4TTcRxwtt54Ph' ...
                    'uHY6ID/y5PrZV11fwPiH8H5n62sbptXsqrrMnxMaAFDpZo1qAgAgGERVAQBiyvVXS6stFIOqfRYnSwAAy1jHVfo8djS5yvLT' ...
                    'VXdsfW1DAbswVJn16qbfLZjXZKvwD+7fNTJCYPsgCjQeIQJRdk2M4FdtToILq2UuMsqpgqRRlAUB0lTRamW2b9/Cljim1V0x' ...
                    'cum8+erz113VQo1py6cRpHqQxgAAQhSoMUTQeg0CAIohhJhqpZFegxAgCxIJ0vKvHe+43L7aaQ35d6TwmbV319Y5decsZefy' ...
                    'bEejmntiZIBgoG9+USyiGORkiUsn0PiM1+MAgCofU+Gzr3tm1d5dW0ccndYzJ/0nXnBeFCxss73h2ir0Y9ZDAoCTJXABAgCf' ...
                    '20mRylivw0aTv3zyntYzJ0cWnbWrH7r1O9W5luGMMjaCyNSlOt0+AFBdlhotej0OlrE6GLW63MHzibWrHxpBdDZtXPOzxxs1' ...
                    'tqOAbaYLWez6cbp9vG4HRSoOlnE5qMOffLBp45oRQSfQ1tx0YINtjDSkaExKtw+Wg1HLvS4bTW5c/3SgrXmQrReQzt82rln5' ...
                    'xGKVHnwo5sjEiRnPG5NLl+4dlrGyjJUiFavV7nZQOLAi02lpOgrBvTZftFBoBi6cXEZAWBSplLA2G03ueOvvLU1Hi0ln2+sb' ...
                    'ln6/AV93CybVMGhS+S98u7Okt/QAQAktWOkSbJ9trw/yBroAdEQh2XToTdeYGjwjUSilJjQGLGNy4QPSQttoAgDe3bZJFJLF' ...
                    'obN75+ZvzakgrPaCXMhTYfUtOprWr3coFmU+Tyo9gEhRSO7eubk4dPbteXtKrb2waZUao+syZpaUGGAjFKmQJIWTa9+et4tD' ...
                    'J3B6f1llaQHRmK9W8oAsaarKWBYyNc74z4FdRaBz4ujBqeNYgrbn2Y4upHhMxtGUPnT0gSiWxGnpjeilBwBw6eG4+ImjB4eb' ...
                    'zsljn04Z7yApBLnX0cwBJRv7/J1mHDMdPkUnyqn6SY4XjO+x0SQOdbjpnG055nYUbFHMjCbdOHxuwxTSQhtDHW46gbYWj4/S' ...
                    '5DioPJIm5tMUIU40u29gFQfL6B2TbFYChzrcdDo7AnhxTtN4pHgGnVxEstF8nZLD6cbRlyh0cWHFfKYnrSSFNIU63HTc1hjg' ...
                    'WqDyg7MPUjxEstGMVeVNaDKKC/X6JZpIHfO8qJ9U5N7jrkjOs0750hHF1BelyXFNCeOu5vDx4sTMaCTz95zROEaZMstkHABI' ...
                    'xHMeBuZbUDtCiWhClTjNakuAxQ4qj8BDJq7VqKbs4/XeZRmTMqEBw9KoUXpa+YO9y6TBSFSnIyu93tG0DNf+oaUjarbUZ8tx' ...
                    'BOWaEkZEFQAgaRKheoAIa2QEjzB0gyBpYr/lKa3WpBrPdKkyppUuvBkho9iSnG/K8qVT6vJGOZULKxQb1/gOxJRrUiC1CKV4' ...
                    'QPGgAU+HaVIg41i8v5wy1mPThg0ASIhWAFAUCQCSgopDzbV3+dYdX1klxwuhkAK6+VVeE5oGMufQq/7/pV80BuOkpxUnpm5z' ...
                    'jFXZV1Y53HQqq2sURYtyKr6jV2OpWqNJAZDDX8BI5UEOa0JTxkKTHU1G45jSShBSd9VJUcGh5tq7fDNrwsW1+956p2oMy4WV' ...
                    'Utaiyan8AnybqwAQDEIMAPQO4lU+NSORlV1/aCRO64MmzTjQk1ZYkbikhzrcdCbX1r+8oTvKVYVCCushKRapyXbCOPur8hrw' ...
                    'AADKQNvEiE3jKV3dgV4cpg0r+EBPK1HoM9cxubY+197lm1nT6mZ2xZU4J/g7Zf0rVZPteorlKo3vMO64MMl0nTIWY5Nx9LTq' ...
                    'ikkAgBAxrW5mrsHk6x2Kohtmzve3H7ezVlybS6stAKDJcSVymLBVZFxCyABFjoOUyH6/h/em6H8ea+2tuG3nU/fBJuPoadUw' ...
                    'cz5F0ZCjCjD7deXc6w6dDHK84O+UJV4zfr1qsl2JHNb4jv7G1poc1+S4GmtWY805ofEHZf3mOBSO4bEVJ9Im4xiDHETXCjD5' ...
                    'MO/rN//2V4+cu9DNMlZ/UK4GAADW28s91e1kO/RdyewvfUySOM1YayCt3KQXY904OK1wkMWh4/GWX73glkP7t44fW+rvBABI' ...
                    'B5QrEV2mi7cJDccLOpounsUH6ca5esEtHu+AEtykwqxn3bBoSSQmnTwbBgB/p4zvnoOnpYzTmgOUxGnB05IJzbFW0YjGWG7w' ...
                    'wEoQ4tg4SVHVjXPDoiWDi6EwdGY0zJ05e8H+/37efC6U6kNQBoDugMyF1IwDoixQuJDa7ZdN2RTl1GOtol5rTGhMOQUAkR40' ...
                    'M2cvmNEwd3D9Ktik5x1LVn74wY6Pjna47LTX48DfcLXPgr98LqywHpJiUH/LT9hlXFgxll5d6ZuV9YQyouG41HJjJC7h+2Mc' ...
                    '2KA7VTA60+tn3bx42ZZX/rz30465MwAD8nfK1WWWap8F+g4a8caJXjR8vwkY5VScql+IRhDieFRlLMY3L142vX7WoDtVyB25' ...
                    'oijc/e0r/eeaayc46yf78FY1LJ3RABXl1CinRhOqaaUhFInpE6Mm12A0SVG9EEotClePn/j8P/bTtHXgn2tSIZ+xoWnr8kef' ...
                    'feS+m060RkkCVUSiPrcTM9J9BH03Iw0ECpbRMpJCJkSrPvtnRKOXGwBY/uiz+aCBodjN/cZL63//68cAwONkp05wjS+3MQxt' ...
                    'XH7DSgfU36JCFssosiiIiXTXAMB9K1cvuv3ePPtS+OezFt1+b6iz/aW/PheOcodOCpLkLi0Bn9sJPdtCs7PoDwqkWUaRRb0M' ...
                    'm9Dc/t0H80czJHQA4O4HnhIEfvPLf0qKysETwQmVXk5MsrQYjEQxpoxuwiB4XuSSAvRdz0znolsG+pZhAFh42w/ufuCpgnRk' ...
                    'qJ7tu3/VMzRtffnFdQDQej7UHaPKvS6Pg7oQTFCkAgPb6GOCkoLYU2Wgp9DoF28AuO2u5cuWP12oXgzhk4/Llj/t8Zavf+5H' ...
                    'ANAVk7pinS4HVe51Wa0OilRoUoaeXTZYFKlgELismKCY/AJplgGAex/8xa133l/ALgz581kf7du5dvVD7Rd6t767HJTPVUKS' ...
                    'lIWkjevc6cIsTFCSosqLiolLxdjxKx77zRVfvaawwQ/H02tcIvbHtT/evvkF40mXgwIAt50yAiJJCnoWEoxEMBQAMOUR1vUL' ...
                    'v3fPip+zJQ4otIbvycdPPn5v08Y1nx7cYzpvo0m8C4ChyZ4zhI6DFxW83pIOBQBmXD7vzqWrvtJ41RDFPNxPzR54/5/b3th4' ...
                    'YO+7ebYza+51Ny5aOmvOtUMabXGeuA60Ne/euXn/nneOf5bbfqxLpl9+5bxvzL9mYdW4vHbDjGg6uroinZ8d/vDU8UNnW46f' ...
                    'bzvT+XkgHuvWX7U7SsvGVFWOq7no4kumXHLZ9PqZLnfZcIY34n67SZYlnosDAMPaLRaquMEU/5cezAFZKIfTXewoUhr9lZBs' ...
                    'GqWTTaN0smmUTjYRAHDkyJFihzESlfo1xmKHMXL1P86PGApnoCp7AAAAAElFTkSuQmCC' ...
                    ];
            case 'efCosine'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAGAAAABLCAIAAAAAkRWaAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAABjhJREFU' ...
                    'eNrtm09s01Ycx3/PSZMmJFlawugfxuoUqKYOJjTEgcuEhDQxjRyZBAekSWMTAm3SJMSJcdgBDa6gXXbYAQ5oGxJCTBob0jiw' ...
                    'IYVDA2Ui0DRAG0ja2Gn+2HFi++3gNHUcN3ZS+TlU+Z6cZ/vlq09+v997fn4BmqanpqZwT02ampqiaZoCgF27dkFPTVKwUHbb' ...
                    '6Hb1ABmoB8hAPUAG6gEyUA+QgZx2G2iQkDyNq2nMV3BRrjciH+UIvecaOWuLpa4AJAqXMBfD/hxGIC9oz+KiLBanxeRnzrFJ' ...
                    '8phsTjGZe1hlv5Qdd7E/BwDUKFCjq14sJqeFmROEHdoJSBQuia7vFTQmJaWzhBnZBkgsnpUdd3UMtQwiAJDS2erCT8R82gNI' ...
                    'LJ6V+590fvurP4lZtQGQIZ3WEQQA0VJ1+tEXZNySBiRzD83EjiEjhsuRMUwcUOXO2jthABgJv4qfImCYKCBp8YqU/8uUrdUj' ...
                    'iAWsHDD8IgHPb96jRgLVJtnPuCqBryMHSFq8ImavIo9jLZ2wgBnVRwLjfZdGEC7ot9fDR1E8/a/VTsgBErNX19jDAyQxjS0z' ...
                    'mRGrbRMCJC1eaet6eV7bkkCyhk50NvJH/BOrnXdjiuEC4HxDCwt4ZnnwUhSdjUSTERG31XEnIgSorfzShA8LONpYelK5iWgy' ...
                    'AgDP2QmrnZNeD8K8hAZWP1sA6T9to6YwA0B0NqIcoLhktWESESRzD81cthodTelJ5SZSuQkAgAJGBctzjEQEYS628oHX/83l' ...
                    'eZ3C3Fx6QB0+r2T3AG+1+a4o0rp0QC+5NOFTepG02psdKdYURLrTQs2kWZE6fMqZp/v3J9YDIK3KDYCaB/UaIKRTX9ThIyw8' ...
                    'pVxD6xFQo3STC1RP7Vo6y+EzNBk6/e3HVtsjkmJ8Q4phpmLmLqblWZFjfvnV8mk02BZBrClGGqWWp4XS6yzllDvooQPZn2K6' ...
                    'as4vu2QPIHWWIX8nPWzcsqZ1JfOyL4I6yjJFzg2DxGySAER5djY31oPI8AVGXXvoG/Xj7JzlT2E1e2S+Rl9s24wU4WGKcoVM' ...
                    'PuKtUXYCMjneq7Vn7AYAYD/qf3t75BgJk3aPYmwF2o8g8CM8TGF5+NR39602SKQGeXeudgozlbYY7aFv1IJohMJ+lHg0bnWi' ...
                    'kQCEvK02qmOmAryk2dQxAKgFo9qNOxx4mPr6B+GNB2QonOKhreFsbJnRCFXq67PUm80pVlczo/GWQTQSrO2A6Bssv/GATImX' ...
                    'cIqnBiWTcRTZfUFh9MLidXtCgHTnilo1MgpjA2+R3RdIOCfwHQDgCB01eSVO8ZS3ojBqkWWKCDAiFUEmylBdmKkgxFOb5TCm' ...
                    'WjMaCT45YjEjcjXIVJbVxUtI5KhNwvgG1JrR1oHO9zqask0EDkA7WVYXkquUu7QtVB13ocFVljc8krVPreTerLaVZWohubrd' ...
                    'VwWAbLV/powYzbapkrVLa0SHeefGIx3fK5d8QTb4QaF/x8sce2+em553Zgp+h+z0DVvrmSQgR+hoZ7uE5JIPl3y8sJRhco+f' ...
                    'Vyosn3paAMgEfO6TF3+z1DPpiWJ7pRpALvmkzJCaTpbNLzK1N43vbtlkuWHCgPq2nk+/NFU11GjYpZcZJje3IAJAnQ4AHPrm' ...
                    'b6sN2/B3qHgC/rmfPvjRuCtYQX0ra2a46sIVFwDgkk9p4YUlACiX83ML4tyiCABZduUl7OhQkIBbGwAdOH7r2vl912/H905u' ...
                    'HAx4PO631GerUhmgxkVpUdNRh8/hz8+tT0AAcPjMvR/P7Lx+Z3Z0c2DvJB/w6mR6npPzJVlB00xn29hmz9ZDBKwimqYTCcv3' ...
                    'SOjq53Mfzr/OAUBo0E+PBpXGgJfKc3KeW3lxyvFCls1z/Eoyjg4Fj517QMBhOBy2ExAAXDu/71kyrRx7PS6vx+3pd9fP8mWB' ...
                    '4wU1GgAI+NwnLz4mY89+QBpGhiJJBwDC4bD9f+o9fOYeANy8fDAWi7e4zOtxvT+x5cBXtwnbsz+C1Lp5+SCXX+LLQqHI54sC' ...
                    'AAR8bqfTsS38zoHjt8j76YoIUuvTE7/bbUGrrlmT7lb1ABmoB8hAPUAG6gEyUA+QgXqADNQDZCAKAGKx2Jr7WYdSsCCapu12' ...
                    '0tX6H56RvLRKAMPyAAAAAElFTkSuQmCC' ...
                    ];
            case 'efCrossed'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF8AAABLCAIAAACZacwvAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAACgRJREFU' ...
                    'eNrtm1ts29Ydxj9eJVGiHFGWY3vKxbXlXGsnzdxLNrQrhg5LkyUthgRoMmBth+xhSIsB3UPWZntYh3QBhgIttoesxdw9pMiG' ...
                    'Ag0WNENaZDCaZPMSNGnk3Bo5V7u52ZJiWZQsURT3cGyKupiSJXlyMn0PxjkkD3n48//7n8NDimpra0Nd04gFcODAga6urlr3' ...
                    'ZM7J7/dTbW1tV65cqXVP5qjoWndgTqtOx0x1Omaq0zFTnY6Z6nTMxNa6AzOTlgq+ceD6sDBOqtGEZHcMC/bhZjr8jOVch3Xh' ...
                    'vIZV1tYf/t/R+ejswf0nOuXb84GHLIvu0kKCbOe4CAAlrQUnZE09u/LeyVOH342rdo/v2a6ndlR40fvAWdpEYNPOoT//fq18' ...
                    'e37+Xo6PAAhCAhBSGEFili21OZiJ/sPv7921Bunog0xnyxtXN2xzpQK2okeGCKAULUhM51K+2ydRaXXXT7oOvvejB5POltev' ...
                    'yBfFGTUZjLMABIlpW8R1+6TGeZbDn/bt3dXzoNFZ93zf2KlUpu6g9GJatuYfr5uLVAWJ8XrYFR3zW92Wk/7ht3aseHDorHu+' ...
                    'b/zCjRk1CWmuyUKKBsDZKMFNez1soyS6RO7yjVAZETQX6azfdqR0NErSmbOFmAuAIDHNC9hun+RbJLlE7qR/+ItPfnF/0/n4' ...
                    '8NDYF18X2BHV9GI6ZsnfT5yVI+Kvbp/kEjmXyPV+sB+ach/TeefXZ4oeoxaiA33YUhhirklABn9NJNS3Xll1v9L51Z8+TY2F' ...
                    '9GpKDuplalwzHqmHj6I488+jmwuAIDGCm+72SXoC0lKj9yWdY39JlHikPmwZ805Aa5/ueK+HBeASOQC/+/lT9x+dM1dGlWDY' ...
                    '7Ijs8CFKJgvETq65ssMnFB6fW3RGI+fSsQHzY/b0ns/ZkoqFjFWjuZTRDBQ9fIKQQlO52WguTKVnAC6RC48rB99/qZRuz9ZT' ...
                    'aDo2EIxeGxk5GIyPkBlaAB0BrV2T3ZHwwnTU0SFGVkrztn63U29y659jOSdRDXknX6mRBtYzBiAmexv48zogCSFMhY/EpjOA' ...
                    '3HS3T/q3P+ESlbNfnvxBTeikYwPq6L50fMAJ3FXYkMKR7T4MhuAK2pGM+tIW+nSg8/S4tv+vt9pWjb/zs04A6VSxQI5qxS6O' ...
                    'gNbuowZJeTDOPiomjXu9HlawWVxiEkgXPRWq66x0bEC5sVMZ2pmOT5qow5bqsGXmFz7qMgC7Y5gWElxjBCKluK2DR6SNP731' ...
                    '9+NhdaxApEzcDejlnGFLN1cy6TTmZt1c+dkHQLdPAkClUyhBVaOjju4zctFlBCQh5KMuC/Zhu2OY9YwRQFoLrQa5vW+r0tqV' ...
                    'Ra9C3cz6n6dGGkghJnv1jcaRKzf7TM195HhJc8Lq0FFH96WCH063t8OW+Uf5MKiXWc8YIyS0VlproQEkbjocnQty2iZGAln1' ...
                    'aPHwCUIKoIOUc8KHs9EAnAJts5R041WgY45mClCWv0j4EEAACCDWLilhq33Jwpy2JubC9OFTcPDibBSAtm/MK/HWKqVTChpk' ...
                    '+8uHQZKAAEwmIAJIpKxNvolhOR+QUTnm0pWTfXR/hRQma+rspr0e1tvcMOt0SkSjAzJWjeHDCAkAxF9Wj0/+6kbT6swiaY65' ...
                    'qFtZdJRRZ8HwMfprMM6FsgdE0VFghajKdEpHMwWoQPhgyl8kQ7N2ibG7g5eyUBrNhbzwMWafHH8ZAE2GD0k9K5d4ZpeOOrpv' ...
                    'pk2mCx9aSEyGTyuttdBWj0+Vg/bOjL/MwweG7CNHvTn+0h/cCSCSerwrS3rUKp/OTANnClAmfIzb+UV3SUFrpVm7ZPH45EtZ' ...
                    'C2BFw0d/ar8XXm4E1K/1EECDcU6PIEvz92aRThmBky8yeOlVkp4BaC20tcnH2N3islZ9b374TDf3IYCMFuvXeojFSAISPCJt' ...
                    'Kz63Kp9OeYGDQnMfYi7o2cfgr9iNrAWN3PDJ9pcasxgByVFvwRw0GGfF9udK7G05dIo+bRcDlDv30at6+ABg7RKABSsyPUyM' ...
                    'BIr6qyigkMIw8zbMIh0t5q+MTu4zjkn4jDmyLJAYCaTkzLJGvr+KAhKcqyjWNYt0KowdABKnkgIxF3kXTsQIGTexdonKW0WY' ...
                    'KJaA8gGN3HmcMHrMMvTMil2l97M2a4PGZRcfdZnjIzw/CSgrfESKYpty2qpyMD8BmQMijMbCy19cvnVG/SwrduKVxk6+ubgp' ...
                    'OvrcBwAcVFpl8psXSEAlAFrXfIm2LZ11OtWVGyEAxtxsDJ/pWpUIKHl9MvSeXnDhx2tKWi2dE3SMiz6koOdm2pB6yMNX2YDU' ...
                    'mCV5velbzbde+eaWMjpZzsopbXu4cnMZ5UYo500mIyT0V3pMg7vgsiEBlIqFrB4fGf4JIEQ1rYWGSAGwN9/52wtp2lbSzLjA' ...
                    'nVbxJmckY+qRqDCmMRffEus7vlZcNu2ahioH5Wv9WcP8uEZfUjGurfZd+ujlJbRtWdmdrH3ewVTqAcAbcjMpCKJMsa5/fPyd' ...
                    'o2eeXLh5IdsgFTyDfK1fdxm9KDV/ZfjQb6ndmx+tsGPlOItp3JYe2lk5lA6bMhjnYEg9HB/R394Rc3U2xEiVsiz+8M3FeBPv' ...
                    '/uF0/51I+KwWv5n56o1tkjg2ITbdeO314Z4FG4HmyrtXJh1aeLgq186XcVpI9EhHR86WV3esftXsHI9UsT9lOou2zQogfdYD' ...
                    'gPWM8Q55U5djNi5U6m2W18w/oB3tu1fhtSUuM/q688Z1AN7G4EzPOSforFm/52Lg3tG+e7GQWj4dw/MEGbayeiYk/ri1xr+p' ...
                    'K3/MeuLpDSdO3/7PicjNq8myT5IvPfW88FDxz5xmW+XTWbN+j9NhPXXhzmefj3z5r3gspFYSRxk6fATAE6lPNneX/51xtVTR' ...
                    'VwbtXd8/1ncIofEzAQyPOr2NrNfDCm4aUy+tZ6rFqZEn4+cunfvA0b661mQAoNLfhR7q3X786DE7nwTQKIlul9PbyALwelhO' ...
                    'oMgLAPKShJSNUuLaiQgPQEyprMI3KeMDl6MXr4WdTZ6Xf/l5rclUgw6A3t3rLnw1ZOeTHKMKNl6wWdwuJwCCyWmnnYKZfyOx' ...
                    'NIDhkdTVr++NhsZbF7e8uPNYrbFUjw4B5D9/S+ATJIiMjHQRRuQvIUJ0OxgPhiOxeBLAnEJTNToADvVu/+zISQACnwCgew2A' ...
                    'zWoRbFnf0Mbik49RhIuiMnKS7/l2z+btvbUGMjt0iN5+7bGhO1MPkHwCAM+oHFN4LFNUBoCc5K1W/jd7T9UaxezTyWdExDG5' ...
                    'S6WKygKw8ulNz214fOPuWnP4H9IhOtS7/eJZf1SOhaOZ0d3GpwEwDP3E2p5nX3qv1rdfOzoPgObE6tecVZ2Omep0zFSnY6Y6' ...
                    'HTPV6ZipTsdMdTpmqtMxU52Omep0zFSnYyYagN9f0VeSD6r8fv9/Ac+zhOFmhbFqAAAAAElFTkSuQmCC' ...
                    ];
            case 'efCustom'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF4AAABLCAIAAAB2q6cRAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAABxlJREFU' ...
                    'eNrtnEtsG1UUhv/xeMaOHduJH2mahIJT3ISWtqi09KEAolCBAJUKpApYIDawoYuKh1hVLCohXmJBWVRiU1GeolJ5lCJ1gaCk' ...
                    'qIWqUItUadKEpE2cNLHjZFxPYs+LxSTOeDzj1q3vOAn+F9HMzczo+NN/z733zLWp8+fPr1u3DlXlKxqNUoqiVDqMBSpbpQNY' ...
                    'uKqiMVUVjamqaExlr3QAuBw74EyecTeJABiPAkDhIFIUx3tmbI7w8o8rFVhlRqi/YydD/NG6wIjKorgSV73Lmz/5X6Dp79nT' ...
                    'HBynvCXcUhE6lqK5OvBWnbe7JCg5CSnaXf+FZaHCSjQjwy8FlnFm/1U4AFBSJlF6AGBi2lLvWITGkIvCQY7Nc7kRDTXuWtX6' ...
                    '3NJBo+OiErlxHFolbMxtWz6zBg3xeY2Wi8JB6obUfZNcAARkoaf/y6WApqf/S5WLwkH845ag5FQz9uNSQLOC/jbnlHLJMuMQ' ...
                    'RDMy/BLlKY9TtOqn5L8Sp0hzAdGFgm+Qk8r6wCSl9EFJQnGKo4SxAORc09+zp7xQzlLyWchJzI6nf8dOkkZDyjX8NUdZiABQ' ...
                    'nUIaRKFIoWmdGVL9D6Dwg9WDAuA3v33C6K6c/vht66Gj3qPvLUI0/5x7DZTcZ/7ZkrPIbka/f9d26tuQ87ZpsmDIuaaPQBe4' ...
                    '0h049V3ble4gWSRzIoMmc7m8zyuE4lxxJ0ksAAk0ca6rXJYZHm39/Zu2gXOMrt0RipDmQgRNbPh47ngC9Qn4J5TZhOsRRCHr' ...
                    'BbCMHQTQyA5obxzhIopsAxCbXBVLtg8cUzL9PYXPd4QizoYIyI9Z5UczIDMAerGyV8nz/FRydTZrUMWSeYc47pP4vME+29lt' ...
                    'yIV2B5wNEQC1TcLiQwPgjLIpkT80m3HJDjbooBTnUnvHZtJECKIZljwJuLUtfLqlkIvMOzKDDQawbozLV/ubFh+aaboRiGm5' ...
                    'pK+16K4Rx31CvIQS8Vx+mdXaZyYB4kM48fdQQoFfinCR/o0XWsZ9xxa7e757elZn33me+MhNBM3FqZj2VNeVivtFGo7rWnxr' ...
                    'HreAgqHIlrIKLVO8H4nxhPbUcP5iQZaxBI2gt0zx6yUuD43dHbCGgqGq2wFMVX40bb55wzNMXu1T5ksr4syM9xY2vnPEihIf' ...
                    'ETRFZA9NlXS9lE6I6QldY+dPtDXRkkXDsLdaMU8PnJ4Zy/OOMkk9uy92s88rQcRdw2ro2FwZ2pUpcjHtNci7mfFeHZ3UBdYC' ...
                    'OuVHE/E2ak9d7iHtafE+5VjfZtieGe+d6jquBZS6wC4+NG3evHkHw3I64zhuHzO7lw4HDY2TA6SlQzofW5GGdcaxuTJM0DQH' ...
                    'mRmnkA7pfFx+NKt8y3UtDMvV1V/QtthDU2beocNBR+uqG6SzyNAA2Oq8el06Nlem5q4rhvZhO9pd92+7Lh1lklp8aO5pfCBC' ...
                    'XSqkE1p22l2rz8qGgOhw0PPCziL2yYz3Fk55yitSRYkI+kBBVwOFJu9oizi5YUu3+GQ72tmO9mxntxhP6JZXFojUrqyv/9zj' ...
                    'FEfVsnkhIAB8ukUHSJW6BC1coEv/xqXhuJYR7fP/euZJcmhIueZuibkE+JH0I1nEPi73kJD1CoJXyHrVyo7qIHtoSuYdctqZ' ...
                    'w0SHg3Q4qE5msp3dAKQEWR8Rq/L5N2zs6brodaUYewR9EaqvFysBA0YMyzEsl6smq6QAoHaupXnOQRkAaGCHlt/XH2JHGgKt' ...
                    'RNEQ3OZ44cQuX5aPuRyxmrwFdy9WTij+RLHdAHoFMAEgQl3ya16UU96OnWv2Lko0XT2H63oOq8eFgHKYdC3q+zw/NaHiAOA3' ...
                    '2Tfw1NYj5LiQRQPg5KePOiShJWQvDqhUeQSRcq996N79RNGQXShsfuI9jpeHxkX1tInPbExwTdOZpunMzT3QI4htXLqN4z0z' ...
                    'ZCc1sGBLdefxF6929QNoCdpz9lGVYugUYweQss8eGLIA4BEljyB6hPm9gc1Pnlj0aAAcPfiYnOLV45agHYCOUUnieDlR29Lx' ...
                    '+KGlgAbAgddXO1g2UD8/kSuVkdoxOV6ua/Y9/PwPFsRs3TdbPti7NjvDB/0eLSBVXpfN6zLIehwv5/7mrtzxyi/WBGzp96E+' ...
                    'f3vz4OU4gKDfA6CQURHx0xkAjdvu2/HIgSWIBsC7r27iJq8BcLNZAK4a1lXjMMTEz41iiSQHgJ/Otm/f8fTug5aFWoEvGP76' ...
                    '8/7Tx74fHRNcbAYAS0sMbbwtXZBm63gU7Xhw9+4Ht++zMs6K/XDC4Y92xQaGR8f0m6sYWhSkvNzc2MC8+eGf1kdY+d+UeP+N' ...
                    'LbIkpdOiIEgz2dlk7GRlAOu3bHj25UOVCqzyaBasqtsBTFVFY6oqGlNV0ZiqisZUVTSmqqIxVRWNqapoTFVFYypbNBqtdAwL' ...
                    'UdFolAqHw5UOY4HqPyUQHCzmh73IAAAAAElFTkSuQmCC' ...
                    ];
            case 'efDipole'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF8AAABLCAIAAACZacwvAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAAB4pJREFU' ...
                    'eNrtnG9sG2cdx7/PnRvHTuwmaZylTf8oTdomrA0RGZBqjUQBsaGAKiFlL9Aktjfbi03AC4oEEtoqELzYCyRgIAR7wYtNIAXR' ...
                    'F+tWNBhstKxV96fxumzt0mxpHLdxnHPsOHe2757n4cU5rhPbF9/54nPRfV5EF989zz3PR7/nued57rFJb28vXCrgAXD27Nmh' ...
                    'oSGnS9JwhMNh0tvbOzs763RJGhTB6QI0NK4dI1w7Rrh2jHDtGOHaMcK1Y4RrxwjXjhGuHSNcO0a4dozwOHHTLE2+zvgijVzm' ...
                    'UAEVoDylFk6TIAEEEvQBzUKgW9x5hJC9BAcI6fm/tcP5dS35mhZ5m6eULa5McYDyVBpIM8Q1XBN6AEDoCQrkMwI5RoQRgvY6' ...
                    'lHnbVzA4v6MlJ7XIZZ7K1Z6b0IO8KTIsCGMCOXGv2uF8hrJXtMhFtmBzzgVHBLsE4WFRGAfEe8YO5zHKJhl/I//vKngKAIo1' ...
                    'KRntdjy9JMnSSiYmyZmslslp3iZPc5PY2e7v2Nkc6vB3d7YEWpoq3UUcBAnojtpF4VuC8LV7wA5lL1P2EkDLnl2clz++kJh6' ...
                    'b0mRtWpyE0Xhs0dCffvaeu4LGAgCQMigR3iUkP4GtcOxTOkfGX+37NnITPq/r0ZvXU9Zy7xrl/+LQ3v69rVtKH0Q4sBGX8Kj' ...
                    'ovDNhrPD+BSlz3Mky569dP72hZdt6H6GB7q+9IX9G3QUhY+OQMY84tO22LFnNMj4mxr9eSU1b70atUUNgKsfxf55aW7DrUsy' ...
                    'Zvw/Kj0DpBvCDmP/1ujzlc6uLGUvnovaokbn/RtLt25v0Tw5n1bpL2oXVKsdxq9o7HcGF9y6YbGjMcqzyA4JlL+G8xmN/tJJ' ...
                    'O5xHDaLGcRi/prHf15JDTXYoewHYYlqw/3DQ9mrv353PszAsrCiIvc7YvxywQ9krjF/b8rK2kPfB8T02qjl2OFRsZ0s09iKw' ...
                    'Vmc7Wcr+WuWlx7++58Q37JleDw90fWX0gH4sDlaZaJWyv1m7ncU5OmXnTT0RRh/evbc/YNdokAQh9FTsj8uV9pwonAKqTrCO' ...
                    'xdGgSr/H+R0LlVycl2emElcvLClpizOJLfuasgjZL3tanqyHHc6nVXrGgppilLQW/TQdm5eX72RiEVlZ0zJrmtfn8bWIoR5/' ...
                    'R3dzV4+/O9Tij96dhZoNmQ1E0bT/L2YTWWlZjE/VqAaAr9XTd7St72ibwTV8FSwNoAYpOgrlisKxTLDLVDorvTLjH9VupxpI' ...
                    'AOIAxIHa1ABQKAAa/4PZdFbscD5nIZWD8EQOQDYZMZvQtB2O5JYjwMZCyq/Yri2YtmO+3+FJ00kcRKF64MgSBTOd2kLLUs0n' ...
                    'cY71wJGXWXK1qjFEMRbs+J2ucdVIOZ6h0AMH6Oz2mc3AtB1C6vEiyQaknN6mAMjLDMC7U6aH6RZip9nsqMEBitVIFMD7N9PS' ...
                    'iukXalae6IQ09vb4jWr0wLm1KB8cPFofO9XOjuuNQrGglKqZnsvFpdWR8V/Vw45Ahp3WUA4px6OK3g0Xq4ksaZ8srPQf7ACa' ...
                    'zGZpLXb2EnLEaRklahJ3u5ViNZG4FpdW++638kUQi+s7Ijmp8etOKynjpVTNciIFYGT8BQt5W1wbFISTJNPmvJeb6S3VxKXV' ...
                    'Q/27rd3B+v4dkjnA5Rg6TDdme7wkNj+eVYXLElVlXlAjK9m4tBps9U784IK1+1hfdfd0/piQPVio74y0XLyoCk8uaMmIpso8' ...
                    'JbPpuVwkrgGYj8YBHD/5oPU61lJUT+czuRuPkwUFHU3wbcsOmoIUrC9EbKLQlLAeMvrxfHQJwKG+0Mi46WUde+wQ0uXZ+x0t' ...
                    '8idEFQCkvcnmhlZZiqpwVWEFLymZRZa0lMwKamQl13+wfeL0pVruX+u+QbHlEdYxzaQr+WokcqRZzDuyFk0K1ZfyykrBxv6l' ...
                    '1IusZJcTKVnJBVu9Ez88X2PtbNhVuSP0rIqfMCm/bYdnaD6UmsW8o4KmTb4UevcvAIUWxnKVpACo5AWA/oQCEGz1PvXc3wk6' ...
                    'nbcDYEfopyqeYdLbxR/mq5qhSFjMVlU4gOJIqeSlEDLrav5BiA0vYG3bkbsjdIb6X9IiL9aSia5DVZiq8E1GykrRKYQMgGDA' ...
                    '+/Rzly282NteOwDElm+TffuUuV/z7NavSVWFFXSgqL2UGgFQVkpxvOgc6uuaOP2WjTWyeTe34B9rGRxbufrEWuRTazkUdBSO' ...
                    'Syn1Egx4j3/1kZGHnrW3Otu1X5nJV5Mf/ExeLLMcV1zn1Bor/qSSDt0IgE1SdA71dU+cvmh7FbB93wQQ/MPtn5/8ePKxlPTJ' ...
                    'tQ9jfp8XgK/ZC0A/NkZ3AUDJZGUlW2oEQLDVOzp27IFTpt//Vk89fsuAyZcmf/Pdmdnl0lN+34bRY1kLZbwEvKMnhh449edt' ...
                    'LXad7BR459xTNz+4EoulU+ms2bTBVi8IGT1xdOTUb+u2sO3Y72C8c+77sx++x5kSi+UfcAVlwdb1pkfIfV07ew8fFJrv/9xD' ...
                    'P6p/Id1fCTHC/W6fEa4dI1w7Rrh2jHDtGOHaMcK1Y4RrxwjXjhGuHSMEAOFw2OliNCLhcPh/PPH09YJSKswAAAAASUVORK5C' ...
                    'YII=' ...
                    ];
            case 'efGaussian'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF8AAABLCAIAAACZacwvAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAABqBJREFU' ...
                    'eNrt3F9sU1UcB/Dvuf239q7d2rKBbFO6rUJJ+BPHP02MKBiDZEF9EyKJLzyY6BOGhyWK8d+rERESI4koPplAHAZjiCERIkQ0' ...
                    '23QTt1FhNPxZtw663Xa7t73Hh9t1t/9Ou9vbtZL7fSBdb2/3O5/8zrnn3iUQn8935syZ9evXw0hmBgYGiM/nCwaD1a6kRsNV' ...
                    'u4CajqHDiqHDiqHDiqHDiqHDiqHDiqHDiqHDiqHDiqHDiqHDiqHDirnaBeicgx9fSshi+kchCdBEbNIhUxmgc6I9Pi6s3kQ/' ...
                    '6Xn24dR558h5k7v+RrA+cksU7sYT41M0QQAkhckSv+Hq39h19cK509sfBp33vrn8z7XG6VgjfcDR+xywMX3IZoWtFQkhkhAm' ...
                    'zQ7PXHikxO+MjQmlfKwWdX4Y6h8cG740/CQFpJl6YA0cQJKiGWgGuSOTaZoxBt5j5j0A6pr9s+MjpRjZlrlKqaSGng1+9fuX' ...
                    'P9/rTFIaE71yzAZAFuqUQ9JEzmCmaS5TOkWNLl7bX0pJ1e+dn/rfvTLn/DfWBtTLlhgBeGtIsrpE0cU55lJVNj1IhBsymJyE' ...
                    'Ok0UwDTlhpNZ31nX7AdQCKhr3/ISa6ta75wfPfZrRAglGhifiQmtkugSxYzGSYQbcluJDCdz+2jmxpXcpdoZeLSU9VhJFXpn' ...
                    'aPDtoZnpP5I+AF5EJuEp9EkHHwIPSXRJkkuYaU1V3PQAOXONPm7CbZnckdVv1jX5hUwdZ6CtdJql1hkcfGM0Og7AAmwlE+pD' ...
                    'I+gEMEI7cs+yWKMWaxSAGkiZa2ojupIDoAYy8x4T71W3z7nTJW1zllon7VIofowC8JPREXRGqDu3oRx8yMGHYkKr2kiO2ZIx' ...
                    'mxooq33U2fF652LLrvi6I976PHzrgvJ6ilAAU6CRYmcxWkkNBEC82awGIpnzK33xauxaefbUzsUWX8HekYVe8e738r2oG0R5' ...
                    'x01J+ugUaJDIhZiUVgLJA+TgQ8icZcmbzemjedvHGWg7e2pxc0pJpe5CE8IRafokzFFSYNvlBumipuepqQOk0Jf4Meon13Pf' ...
                    'd/AhqzWaGoBjzjR/4U8BPZIxKGegbbHLTWV1pPAB2XaRLOO4FpjWwBQAKbw1bafcJsoVum75MbqN/JYXKP1auYrljWV5o2Ya' ...
                    '/XVkoVcKH6DujHKJs4iR0keFmsiDSC6QxRpltI8S6iQ/fhcoZzh66shCb0I4kUVTulE75RhAuVNM3T55aT46dL/MEemmI4V7' ...
                    'EsIJeKzsj7GNGEB+jHqRsYhb5nsHAKfunRkKYMtzf20MdKC86KMjhXsoHSpKk2XEtSwSKKd9FiYXP6t+f/Oe/sP7nyl/XDro' ...
                    'yEIvjf1ZOs3C726BKd+yUAjIg+wNQLp90nfzAFq6bh/ev6P8cemjI4WOo8Wu7VziLAiU9yrmRZGNpH/V0BdvrtOFRgcdcexV' ...
                    '4l5015QIlGfw+bY/C6c8dvPTvU/pRVOujhTuQd2shjlVCpAbhLFRzK5kwtW1of/ovg060pSlQ+mQHOkrnyYNlHeRLiWJcMOu' ...
                    'rb98sFuHZVg3ncTEqTLnVHYpLdlAeSeXOsrd1s6O/rd2dOtOA813oVK4R4706auTNx0g15H/4XFMaOWtkydfbCakIjQod1W2' ...
                    'm3SuptjkUj/3cZHZr3dvJ2RthWigWUeO9FVCBzlA6oceACLUrbzYahk7uvPlyrko0TKzpHBP6lU8WQkgRpTeeaVu4IUnjizB' ...
                    'r6v+X2wyEk+SuQRgy3twBJ1eRF5a6dmyailooE2HxueftEZEzbvkXBdERDqbpFLBZvSYLQc3H1oaF+066dDZJImIZW155lEK' ...
                    'HQ+S1GNQf/2KPeveX0qacnUA0CkxtWwuyqgYihLl+XynjQusfo3jK3XZrqCOAgQAUyJxW1OLdNZSHU+m/o0nAbBRaGLh3KBJ' ...
                    '9jtWrF332dK7aNch9nYa78vPNFV2RfPrDuf1PL36WLVclGjRMTm2pfY7eodKJiqZOBex+vZVZSplRYsOx3cT+7c0HtVwbpFv' ...
                    'dvFW395acFGicd0xe/dKoeM61kHsLrO3hlxSw9R2Gsd3W1pRPlBtoiyUV+bf0ZWbdW0oAGrWRUm5V3RL04dogiz0JmOXaTxY' ...
                    'aDEidhcAYm83ObbVuIieOko4vls9ZlnoVR+q9hi1pyJ3of9rkYyBVLuAmo6hw4qhw4qhw4qhw4qhw4qhw4qhw4qhw4qhw4qh' ...
                    'w4qhw4qhw4qhw4qhw4qhwwoHYGBgoNpl1GJS/xtjtcuo3fwHslJ7wB/UxYcAAAAASUVORK5CYII=' ...
                    ];
            case 'efIsotropic'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAGAAAABLCAIAAAAAkRWaAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAABjtJREFU' ...
                    'eNrtm3tMFEcYwL97gICkemCPhwVdXiK15zW11CqSaqzQVhRtbUupbX1E1ASa1Ni0miitRYtV0xBrNFHTWh+UqE1NI7bGaqUa' ...
                    'FQWPahCRhygcKMchHOxyz/5xdHs+uN3bmdmlZH/hj7uE/b6ZX3ZmZ76bBYqiDAaDS+YxDAYDRVFKANDpdCDzGG4tSqmbMdiR' ...
                    'BXEgC+JAFsSBLIgDWRAHamnTN7XdM9TeOnKpotNiae+l7XaH0+VygUupUKpVKn+1Wj86PEqrnZ38QkL0M5K0UEFRVH19vchZ' ...
                    'D54pKy2vbGxtZdR+fBtqt1NPh87S65akzxStnTExMWILyt2591pDgwVhaCsdjrQJ4z99e35wUOCQEpS7c++Vmps2/2FYojmt' ...
                    '1td1z25Y9P5QEFRRU/vZDwfNDif2yCqnc9ErKTlvpJETRHySXla0q7K5hVBwh1K5++z5tra2dYs/IJSC4GPe4XTO/2ozOTss' ...
                    'JRcuZ63b0Hy/nURwUoJMD7qmr/3yTncPQTEelFfXZK8vuFbfiD0yEUGmB11zNm6lnS7SXtwEajQBmpC79+8vLtj8d10D3uD4' ...
                    'BTkcznlfb7OK4+ZhOrq6Vn7z7e3WNowx8QtasGmLaPcOizqwf01kNJk+KdrhcmFrAGZBy4p2ijbveOIXEMB+NtTWrd21B1dk' ...
                    'nIIOnSmrbDaK6cWTAE0I+7nk1Jmf//wLS1icgnacOCWyFC8U7j9koWn0ONgEbf7pKINv5AuAnYbctHc+KCo5ih4Wm6Aj5VfE' ...
                    'l+Kdvb+WNrXdQwyCR1Bh8RGnSuLSkuc8zbKv9DfEsHgEHas0iG+ED8UnTzN9SGsyDIIOnD5rBYXUKp4MY7X+UnYOJQIGQfvP' ...
                    'npfagzdOXCxHuRyDILNFgpUhf8quVvXQjODLUQUdv1zhUA6Kn0ZszIAWLl6vFhwWtW+37hIv96BztfaW4GtRBV2pa5S6+/3Y' ...
                    'B14337jdJDgsqqA75k6JhDyKnRlQUIOxVXBYVEFWm00qI/xpNXUIvhZVUJ/UnXdjYxgvQ4zu62OsApeLg+IBhA5j5rhHaKHr' ...
                    '6SEiyM5Z2RBaaUAVpHTYpRDyEF5WQCzD/PkeAni0g4iNU6tU4ht5BM7x5adWBz1pr88HVEH+vI9nEML79OwmLEQjOD6qoNEj' ...
                    'nxJfiiectw8AjAkPExwfVVDoiBEiG/GENpu5p2eAcdFRglOgClqRNkNkKSw2huFz+wDAxLhYwVlQBSWOiVbYpXmQ8bQDAC8m' ...
                    'JQrOgmEdNHZUCHoQX+k2tvAZXADwfEK8VjNScCIMgmIjI8XS0g/PqcdN2kuTUHJhELTpwywRpLDQZjP/wQUAGSlTUNLh2Wro' ...
                    'IyOISmHx1U5m6lSURRDgErRl6UKiXgDAxjDm+jqf7ADAwvRZiHnxCBoxfLguMpyQGgCgzWZLS7OvV82dNnVivPAHvBtsu/nd' ...
                    'ecuDmF4++0Zf1fg6rFhyF8xDbwA2QQqF4vOsBZaW5m5jCxZNrBphdlZnvzs2AsNNjbMelD45OTttpp2mWU3CTCGqAYBUvS4n' ...
                    'czaWTmE+cfDF0o+qG5sqam7aadpCN8O/p1Lcp5v8Bqg5sDtywUY80Wo0BTlLcPUI/5GMrbnL38svMLab3F/dPXfL+i+rx1ke' ...
                    '/ks+nmzLWxExKhRXNPwl16gw7XerPh4ZHOzlf+w0zf7hzb59Vd7kCUkYAxKpSeviYvasWa3VIK3QBLB9VV765GS8MUkV7SfG' ...
                    'xx7IX/NcbAxxKwAAoNVofly/BrsdIPqrBhUZcbhg/ZvTUwmKAQCAVL3u8Mb8l7GOLBay5+ZUKlXhymXJ4xML9x/q6OomkWJ1' ...
                    '9js5mRnkuiDSC3VdPT1FJUe/P456YtCTOdOm5L41jyK5Txb7lcxGY+u+0t+LT/5hRStCzp02deFrr+rj40g3WIJ3VgGgh2GO' ...
                    'lZ0/ceHSuaprPl2oT4hLS56UkTIlPFSkGqY0gli6e3svXq++Wlt343ZTY4vRaOro8zgrolapwkI0Y8LDxkVH6eJik5MSESs7' ...
                    '/z9Bj9PLMIzV6nLBMH+/4EDiLzVzIsY7qz4RFBAg+DdiQgyR0x3kkAVxIAviQBbEgSyIA1kQB7IgDpQAUFVVJXUzBiNuLQqK' ...
                    'oqRuyaDmH1sZc2d3xmShAAAAAElFTkSuQmCC' ...
                    ];
            case 'efSinc'
                encoded = [ ...
                    'iVBORw0KGgoAAAANSUhEUgAAAF4AAABLCAIAAAB2q6cRAAABQGlDQ1BpY2MAACiRfZGxS8NAFMa/loJU6iB2rJixg0pRySIO' ...
                    'bcWiOIRYweqUpkksJPFIUkpn/xXB2U0ROjsoCIKjODuJiGv8LrWkKvrg5f3uu3fvO3JAtmQI4eYqgOdHgd6oKQetQ2XqBXnM' ...
                    'o4gFqIYZiqqm7YIxrt/j4xEZWR+W5Kzf+/9GvmOFJusrs2OKIAIyTbLWj4TkM3Ix4KXIl5KdEd9Kbo/4Oelp6nXyO7ncnmBn' ...
                    'gj23Z375yhsXLH9/j3WLWcImuggh4MLAAAo0rPLb4MpDDxG5z44Ix6QQOndqpCYCdvhUbFhk/OGxlnjUcUKHAfu6cDhJTq1S' ...
                    'ka4WeZuTTCxjkbyCClOV7/HzP6fa6ROwMYzj+DrVdobAhQpMX6VaeR2YLQA3d8IIjETKMbO2DbydAzMtYO6eZ47GD/MJR8ZZ' ...
                    '1F33I0wAAAAgY0hSTQAAeiYAAICEAAD6AAAAgOgAAHUwAADqYAAAOpgAABdwnLpRPAAAADhlWElmTU0AKgAAAAgAAYdpAAQA' ...
                    'AAABAAAAGgAAAAAAAqACAAQAAAABAAACGKADAAQAAAABAAAEFAAAAACaSJMDAAAABmJLR0QA/wD/AP+gvaeTAAAAB3RJTUUH' ...
                    '6goDFTQiV6n2gwAAAYp6VFh0UmF3IHByb2ZpbGUgdHlwZSBpY2MAADiNlVNbbsQgDPznFD2C8TM5ToBE6v0vUGNgla7SSmsp' ...
                    'QYzNeGxD+q41fXUThgTdMFfNBgZagSkgbXoaGwqyMSLIJrscCGDn4e4dIDf/2L8LgDRpVjIy4CwgwBWmve//s8uzdkV5AY2w' ...
                    'vZR9aOnD+F1ZxUhpapmwcvLCwNA49lnG6kWZeYdg4cfEc+9adGtsy7bwpHZ31Loc+uvACTci9gEMRfmaCXBL+ofDEzzjTavJ' ...
                    'Ku1lefWoqugpIjOAdfpdITuJ+rjNfdb74Zj6gMwxjDb4+CmNjbcjnOp3w9zTe9Nxs0nSybaJOaHgOMM4ydiJDCfJOYLU1agM' ...
                    'PC7mYynvlczSaE2ll9LJrlFGlx+JfFUbqpZazjMmyk0jezj2oapnR1+xDpz8R2vV11O52ePNfgrcqMX8qZXzCivlPS5K46NE' ...
                    'YMk0pnWU/YkQNyrRgAxx3+pWmz+nbKj5/kR4xatqEF542JAqexBQGxGHlLfmQ/oBDZzlIubF5lIAAAAldEVYdGRhdGU6Y3Jl' ...
                    'YXRlADIwMjYtMTAtMDNUMjA6NDE6MjcrMDA6MDALNVYkAAAAJXRFWHRkYXRlOm1vZGlmeQAyMDI2LTEwLTAzVDIwOjI1OjEy' ...
                    'KzAwOjAwR4Gs5QAAACh0RVh0ZGF0ZTp0aW1lc3RhbXAAMjAyNi0xMC0wM1QyMTo1MjozNCswMDowMCcQaAIAAAASdEVYdGV4' ...
                    'aWY6RXhpZk9mZnNldAAyNlMbomUAAAAYdEVYdGV4aWY6UGl4ZWxYRGltZW5zaW9uADUzNoN1HkcAAAAZdEVYdGV4aWY6UGl4' ...
                    'ZWxZRGltZW5zaW9uADEwNDQdZCpxAAAAKHRFWHRpY2M6Y29weXJpZ2h0AENvcHlyaWdodCBBcHBsZSBJbmMuLCAyMDI249l7' ...
                    'hQAAADN0RVh0aWNjOmRlc2NyaXB0aW9uAERpc3BsYXkgUDMgR2FtdXQgd2l0aCBzUkdCIFRyYW5zZmVyJzf6egAAC8RJREFU' ...
                    'eNrtmnlwG9Udx7+7K62OtSRLsnwfsRMbx0kMDXFiCCHJlCskIZ22FCgDHYahTOkMM0xnKNCWhjC0FP4ozZQpncBADyiEa8gB' ...
                    '5cgUnDSDQ+LExkl8x44c27Gt09KutKvd1z8kS7KsdezYOdzoO/rj7Xu/93tvP/OO374nqrm5uba2FhlNVEtLC0UIudTduExF' ...
                    'X+oOXL7KoFFVBo2qMmhUlUGjqgwaVWXQqCqDRlUZNKrKoFFVBo2qMmhUlUGjqgwaVf2foPn2uPTtcWlufc4zNP/q+e9eZ9NZ' ...
                    'wZeS/+Y7/BNP++aWznxC45eEL4dO7HIeefrozoahk8lFy5ZowyNdf3yh/QpFY9Ya8gwWgS8YOVu/vVF4trE5XuR1BULDHR2N' ...
                    'Tc9t674S0QDYVLJcUbQAFEV7dID+w+Fj0fwHH7BFEx3NfVcompU5i67OyY+mZS93oMn21CedAPQGzdKVCwF0H+8X+MictMVs' ...
                    '3br1Ur/vOfTmO/yuj4WWVglAUSGzuiC/2d0/OqYJ9+YpvG5oiOsLDd24MLu7T3PicA+AM8P69Tc5Zt/uPLhR2PiDUQCit1/y' ...
                    'nkEkePudtY//uvqxz5rbDpfFLEJk80rxkU156+r+HRkbthcVfLTv5tm3Ow8m1IM/4QAIA62RoCsSDu3656Gb6/cs7bPnl3gB' ...
                    'UH0KfVze+zrT2sPfuqUKgOvMYDg0B3Pqoo4aWQnz4ZGw6JXkYE/AMyqGDIzGzJrtensul2/R29UqfrhL+NO2LyOBkeTMoqpy' ...
                    '8y2FA/+xRB/pbLL7Ncfa5btl3nPLj+qf3lY1P9A4R77qHtjt43sB+GE6Sq4OgptsFhHsUsjhMIolJsYcyd+8tMSRFTM7csjz' ...
                    '/NbmwZ7TKVVyapdKcml4pDM03MmVl6yuK/ps59d55aXvf7JuHqA5efqtjjMfRNODyD9Krp5sI0kmn2cxIbEJLg1ZI54sAIw+' ...
                    'pB0S7Tmh9Stx77qq/Q2ev25v623tTK5rsOeJgiTzbgCLblnS9flJEOVA2/3zAE2b8+32/vei6UNkxShSJ04kYvS4Evfu0tns' ...
                    'iNsU61+vTLnGe8iAsUo6qremSuPtYpu+bCWKHK+ltRRIvkFQ9OIblp3c3/zki7du3Jx3uaMBMOhu7BrY5R5r78TCTrIopVQS' ...
                    'zV5PTfxRaC+GQgGgRhWqT5lg6R/inU0AKJrR2q2WCiMTEoaaY2uQJsuhiEFaxypiZNXaihe3XzMP0ETlDXQPe499OjLYKuhS' ...
                    'igJjCwQ+FsuF+3IVXgeA8hGqS04287d9QWQxpS6t1evL7Dqr6PlmBICGs7F2hzbk2ttw+7xBkyxZkZ3+nsHAmWFh1B0e84jB' ...
                    'k94CX6A0IuuVECuedhCZBkB1y5Q30cNAz0FZ8Kq+DM1oHTnGUq00HCI8+Wz/xnmJRl3KzpbdvChEFFEm4sG272C8g15nrtw7' ...
                    'YbhFfO6QsyvxrIHewQEUAEWirfnhd147/9jvkqEZHnz3lDB23DskiMGzJBsgEUKcwdpgoHgCJ14X7sud0OMOmRpL9Dk03Bke' ...
                    '6Uxxbsh1yCLR2vWlFs2Ot268fNEMD75LiN/laXYFRwDilhOhaiOpc8GWYu/z1Iiiefp0Ar2NXI1em6P1HfQDAGG1FoPgPFW9' ...
                    'orrtcJuGs/3uz+uuvy7rskAj+xoI6Zf6G4Aw8fM9lNIN1VYm00nZsNLQGSN0h2wwjkl5g7797ojPHS9hjHbGqCeiKPlHVt20' ...
                    '7OiBDjEUpjS6F1/bVL+Kw0w0Z2gIGSDkePjE28QfmFw6BR03bF+TupRMPlicMrPEvlyZjy00Js7Jt+QSLxPobZSDrmQugKwr' ...
                    'z+WPdwC4du2S+htyXn7uKwDZuY49DRtm9Eaz/LyUFXI4ouyQIo9I8mMR5VWmOkAXpbGrIPRCUGld2OCupFKP5oxcf0qOxuGz' ...
                    '2Jy6LD/VIQcbCoiXAaB3VCZzkXmXaZE5MhKOvRvD3HNfmaO0CIDIFs305FhzfkiU4B6Z30f0QzAoKUV0ESgz5JOpVSoI7aFk' ...
                    'dzpvlehyw5oyrVjWH19x2Ego7Db5xnIBoIBQY7FgR8PZGM4uB11RLmxOYYinFd4fLTVnswDe+3j92++HACxbor2waOTgG3Lw' ...
                    'cxjCsDJqNpQJzOL0dNyUkrZKJdXtIhPQGLl+UazJs7RC0jtPXh+fTTBRxETFl2GN0QYCmXcxRrM4OqDPXsBmlxFZplmuZqkF' ...
                    'AKOh773LONPXnBEaRRp5RnEfpqwsbCzATG1NmUCZQfwTMq2gbIAb05KW9duNvWd9SwFoHD45aRkmBXR84CiSIPMuAJRWC4Bm' ...
                    'KEqry6pYXVTI3PVj63kQiWtaa40cfCPcvllxH6YKDbCx03TNVKfJrCDpW7TBbU+CpkfI56lx8QtivTSGtTlJmE2xZSs03Cl5' ...
                    '+wHY65ZGfC6Nxa5hKOFMS35ZYOuvzJidzjlqQuLph4jgBkAVGmBgpuEzofMbOBaM+WCyWE8kBzg0F8LohLeNB3umkrJAbxjA' ...
                    '4i0FHiGfoqqDnFiYP7OuTtZUo4ZgVDx9/3lzAaC2W6U1ju5TxRjwIXYiYbGeSLgyhtNyYTizpb5Sn1vJLaj/y1PLZ4ljmmhk' ...
                    '6fTPiRAEQFnZ8+CC8RVn+qqmOmqpb+upb+I5XFZiF2fG6YgH2sa52K/7xdKxEywAa10i3g2HznPnnRYaaeQJIowHb9NeX6Yj' ...
                    'q0qAY4BQgVOYGOkkBzjRgSMeaAv3dERzHv3lNa37YmfDa5PQSO70TcwBGgKX4o4NZso6Ky5p55QtnaUBQjxdicT3NMsmlqtk' ...
                    'LlseWPupLBEvBYDOJj+7wwzAXCMqVYy8cGYhzAzQKMG9s3c9S00OkZUubZzLinU1VRvszlPZxEQRE7V6Q9KJl2kOhow6Gn4u' ...
                    '/3JwTk0n/Ii0mEi3lltQz3D2FetqXnplxd8aAgBIFaNfIT3xw9gJYUmFN5p49dDABUFDGUsTD4I8TV/pXZmmhEKBZRbudN83' ...
                    'ucg+cYuXW7IAaDjbTx9f/9IrK3acaBL8+mhRXa0nblaXb9OVDRsWO4vzZntLl34lp1AQT5OQjO7AeBA8E0WZukXAkBYKr8t9' ...
                    '3bve560u1Z+cwo0ompWzsaaf32aJfgrtOyVoc/zSqNlgDj154+Jk++hq/c2Q+7ayUsxC6dEw3B2y4e9ESKyLxCMmZrCBSbOX' ...
                    'R0EIcjRBQqpjLcLCjYKPxVoI0GgDAAQ5zdDqJAsB8MHiKvvA7793c+saOf59+GzjMVE0axw+jcO34SoJyJ8NgpmhAUDbN8j9' ...
                    'HyTnEM/4Ub4H0xeRxiGyyrBRed+3qENMXLbQdITL6h8JFPOy2chMCJxdsK22Zt29cg2ryQKwbEli7h8diPksMCkP1UwI824r' ...
                    'K32z4xNJNOebDJidVNFouAeJrUlx986yAUrHUjpm2Gg65nXBR+kQTjEwcv2KrN3lerjK0ORg+1mNcmN5fZF1+RZjSVqHu51H' ...
                    'gNj8uqYwzdg0cv3gkG+abWQ8VdSodWyX8Oj50aH0OlBWjf1Omruta+Cj433/iOabMDbZOMt8SlHYs2L2QvuK+5fVTe15j7Mp' ...
                    '22rmg8VLcriUIRPVVZbCdt9styec6/OS0TpelvAM4VtIKDQtIgYrY7+JQgXNJc7x+85+EU9nw1eO3lNYMKEWyHdzdHdX331O' ...
                    '/7udRwBoWb+FPfGbVQ+ltak057f7Bq4yF15QNACgdfwWgBz8UOEPAT4ijAIKQI3/jJShjDbWUnDQ3Jq0HmRlwn3jYqpdIbQT' ...
                    'xQpoDSIVBs09lRsLs4rO2RMAnf6haGJTiep82VxyLYAqS8F0HE6hi3HZcqj9hUH3ocn5xTlrFhXeYeHKp+/q4YM7olyi739B' ...
                    'NQdfqOdUVdH3x/j+QGgAAEOzNlN1nnV5of06A2ufqatNJcs7/UMXgQsu5u2lFAkoRNZpLRenudnrMrzzvlw0D/7meKmUQaOq' ...
                    'DBpVZdCoKoNGVRk0qsqgUVUGjaoyaFSVQaOqDBpVZdCoim5pabnUfbgc1dLSQpWXz+Ak6YrS/wAppxWm/vqQGQAAAABJRU5E' ...
                    'rkJggg==' ...
                    ];
            otherwise
                p = '';
                return;
        end
        d = fullfile(tempdir,'padRibbonIcons_v1_embedded');
        if ~isfolder(d), mkdir(d); end
        p = fullfile(d,[name '.png']);
        if isfile(p), return; end
        bytes = matlab.net.base64decode(encoded);
        fid = fopen(p,'wb');
        assert(fid>0,'phasedArrayDesigner:iconWrite', ...
            'Could not create a temporary icon file.');
        written = fwrite(fid,bytes,'uint8');
        fclose(fid);
        assert(written==numel(bytes),'phasedArrayDesigner:iconWrite', ...
            'Could not write a complete temporary icon file.');
    end
    function g = tabGrid(parent)
        % Every tab gets the same two-column shape (label | control) so
        % controls line up when you switch tabs. Rows are allocated
        % generously and trimmed by finishTab -- the count cannot be
        % known until the tab's controls have been placed.
        g = uigridlayout(parent,[60 2]);
        g.Scrollable = 'on';
        g.RowHeight = repmat({26},1,60);
        g.ColumnWidth = {150,'1x'};
        g.Padding = [8 8 8 8];
        g.RowSpacing = 5;
    end
    function finishTab(g,nUsed)
        % Trim to the rows actually used, thin out the separators, and
        % add one '1x' filler so the controls sit at the TOP of the tab
        % instead of being stretched apart to fill its height.
        rh = repmat({26},1,nUsed);
        for k = 1:numel(sepList)
            if isequal(sepList(k).g, g), rh{sepList(k).row} = 8; end
        end
        g.RowHeight = [rh {'1x'}];
    end

    function assignStaggerAngle(v)
        % User edited "Stagger angle (deg)" directly: convert to the
        % equivalent length (stagger = dy*tand(angle), the inverse of
        % angle = atan2d(stagger,dy)) and write it into the Row stagger
        % spinner, then run the SAME onGeom path that spinner already
        % uses -- onGeom re-derives S.stagger from spStagger.Value and,
        % at its end, resyncs this angle spinner back to the same value
        % (closed loop, so this round-trips exactly).
        %
        % CLAMPED to spStagger's own Limits before assignment: the two
        % spinners' ranges were set independently ([-5,5] for length,
        % [-89.5,89.5] for angle) without checking they stay compatible --
        % at the default dy=0.5, an angle as mild as ~85 deg already
        % implies stagger=0.5*tand(85)=5.7, outside the length spinner's
        % range. Without clamping, this either produces an inconsistent
        % pair of fields or an invalid spinner assignment depending on
        % MATLAB version. Clamping first, then letting onGeom's own
        % end-of-function resync set the angle field back from the
        % CLAMPED value, keeps both fields honest about what's actually
        % applied, rather than the angle field silently promising a
        % stagger that never actually got used.
        newStagger = S.dy*tand(v);
        newStagger = min(max(newStagger, spStagger.Limits(1)), spStagger.Limits(2));
        spStagger.Value = newStagger;
        onGeom([],[]);
    end
    function assignShape(v)
        % Shape changes the SET of populated cells (not just spacing),
        % so -- unlike dx/dy/gridAngle in onGeom -- there is no sensible
        % "preserve manual perturbation" path here: which perturbation
        % would a newly-included cell inherit? Always a full rebuild.
        % Still one click: the rebuild discards per-element edits, so it
        % is an undo step and the status line says so.
        stepsWas = undoCount;
        if ~strcmp(v,S.arrayShape)
            shapePreviewOlderModel = shapePreviewPreviousModel;
            shapePreviewPreviousModel = S.arrayShape;
        end
        S.arrayShape = v;
        S.shapeFellBack = false;
        if strcmp(S.mode,'Uniform grid'), rebuildUniform(); end
        refreshAll();
        % Reported only from here (the one place the user actively PICKS
        % a shape), not from rebuildUniform -- which also runs on every
        % M/N/spacing change and would otherwise re-alert repeatedly for
        % a shape the user already knows is inactive.
        if S.shapeFellBack && ~strcmp(v,'Custom')
            S.arrayShape = 'Custom'; ddShape.Value = 'Custom';
            uialert(fig, sprintf(['A %s boundary has no lattice point inside it ' ...
                'on a %d x %d grid, so applying it would have deleted every ' ...
                'element.\n\nThe full rectangular grid has been kept and the shape ' ...
                'set back to Custom. Add more rows (M >= 3) and re-pick the shape ' ...
                'if you want it trimmed.'], v, S.M, S.N), 'Shape not applicable', ...
                'Icon','warning');
        end
        syncShapeGallery();
        announceUndoable(stepsWas, sprintf('%s shape: %d %s', S.arrayShape, ...
            size(S.el,1), plural(size(S.el,1),'element')));
    end
    function assignSteer(w,v)
        if useAzEl()
            % Derive BOTH spherical coordinates from the visible pair.
            % A negative stored theta displays at phi+180; reusing its
            % old phi while editing elevation (or retaining its sign
            % while editing azimuth) would steer to the opposite side.
            if w=='t', az = v; el = spPh.Value;
            else, az = spTh.Value; el = v; end
            S.theta_s = 90-el;
            S.phi_s = az;
            [shownAz,shownEl] = azEl(S.theta_s,S.phi_s);
            spTh.Value = shownAz; spPh.Value = shownEl;
        else
            if w=='t', S.theta_s = v; else, S.phi_s = v; end
        end
        % Keep the linked cut-plane spinner showing the plane actually in
        % use. Display only -- cutPhiVal() reads S.phi_s directly while
        % linked, so the cut is correct either way; this just stops the
        % control from displaying a stale number next to the plot.
        if S.cutPhiFollow, spCutPhi.Value = cutPhiVal(); end
        % refreshTable() needed here specifically: the Feed phase column
        % depends on theta_s/phi_s (via effectivePhaseDeg), and this was
        % the actual bug behind "steering angle doesn't change the phase
        % column" -- this callback only ever called maybeCompute()
        % (the pattern plots), never refreshTable(), so the table was
        % never told to redraw regardless of what it would have shown.
        refreshTable();
        maybeCompute();
    end
    function tf = useAzEl()
        tf = strcmp(S.angleConvention,'Azimuth / elevation');
    end
    function [az,el] = azEl(th,ph)
        % A negative signed theta is the opposite azimuth, not negative
        % elevation. Keep this identity when switching display modes.
        az = mod(ph + 180*(th<0),360);
        el = 90-abs(th);
    end
    function s = directionText(th,ph)
        if useAzEl()
            [az,el] = azEl(th,ph);
            s = sprintf('azimuth %.1f deg, elevation %.1f deg',az,el);
        else
            s = sprintf('theta %.1f deg, phi %.1f deg',th,ph);
        end
    end
    function s = scanPlaneText()
        if useAzEl()
            s = sprintf('azimuth held at %.0f deg',scanPhiVal());
        elseif S.theta_s < 0
            s = sprintf('phi_s=%.0f deg (negative-theta side)',S.phi_s);
        else
            s = sprintf('phi_s=%.0f deg',S.phi_s);
        end
    end
    function s = scanSteerText(th)
        if useAzEl(), s = scanDirectionText(th);
        else, s = sprintf('theta_s = %.0f deg',scanThetaSign()*th); end
    end
    function s = scanDirectionText(th)
        if ~useAzEl() && S.theta_s < 0
            s = directionText(-th,S.phi_s);
        else
            s = directionText(th,scanPhiVal());
        end
    end
    function sg = scanThetaSign()
        sg = 1;
        if S.theta_s < 0, sg = -1; end
    end
    function p = scanPhiVal()
        % A convention switch changes labels, never the scan plane.
        % Positive polar-theta samples on a negative signed-theta command
        % run at phi+180 in both UI conventions.
        [p,~] = azEl(S.theta_s,S.phi_s);
    end
    function s = cutDirectionText(x)
        % A point of the cut as a direction, in the "°" form of the other
        % result readouts (angleText) -- the cut plot's labels only.
        if strcmp(S.cutMode,'Phi cut (fixed theta)')
            s = angleText(S.cutFixedTheta,mod(x,360));
        else
            s = angleText(abs(x),mod(cutPhiVal()+180*(x<0),360));
        end
    end
    function s = cutPointText(x)
        % A point of the cut the way its own x axis reads it: theta or phi
        % in the theta/phi convention; elevation (plus the azimuth on the
        % az+180 half) or azimuth in az/el. Short enough for plot labels.
        if strcmp(S.cutMode,'Phi cut (fixed theta)')
            if useAzEl(), s = sprintf('az %.1f°',mod(x,360));
            else, s = sprintf('φ %.1f°',mod(x,360)); end
        elseif useAzEl()
            s = sprintf('el %.1f°',90-abs(x));
            if x < 0, s = sprintf('%s, az %.0f°',s,mod(cutPhiVal()+180,360)); end
        else
            s = sprintf('θ %.1f°',x);
        end
    end
    function labelCutAngleAxis(ax,isPhiCut,thFix)
        % Every view starts from AUTOMATIC ticks and labels. The az/el
        % elevation view below writes XTickLabel by hand, which silently
        % switches the axes to manual labelling, and nothing else ever
        % switched it back. So after one visit to that view, every other
        % cut kept its frozen elevation list: the theta cut read
        % 10,30,...,90,...,10 over ticks -80..80 (a beam at theta = 30
        % sat under "70"), the azimuth and phi cuts laid it over 0..350,
        % and the full-sphere theta cut read -60,-10,40,90,40,-10,-60.
        % The title always named the right quantity, which is what hid it.
        % Measured in the real app (MATLAB R2025b): 6 of 10 views wrong
        % before this line, across both the gain and axial-ratio cuts.
        ax.XTickMode = 'auto'; ax.XTickLabelMode = 'auto';
        if isPhiCut
            if useAzEl()
                xlabel(ax,sprintf('Azimuth (°)   [elevation fixed at %.0f°]',90-thFix));
            else
                xlabel(ax,sprintf('φ (°)   [θ fixed at %.0f°]',thFix));
            end
        elseif useAzEl()
            % The existing full great-circle cut includes both azimuth
            % half-planes. Its x positions stay physical and its tick
            % labels report each point's elevation; the caption names
            % which half-plane each side belongs to.
            %
            % The ticks are PINNED as well as the labels. Labelling alone
            % left XTickMode on auto, so any zoom or pan re-ticked the axis
            % while the label list stayed put -- zooming to [20 40] laid 9
            % stale labels over 11 new ticks, exactly when someone zooms in
            % to read a beamwidth. Pinned, each label stays welded to its
            % own tick: a deep zoom shows fewer ticks, never wrong ones.
            ticks = ax.XTick;
            ax.XTick = ticks;
            ax.XTickLabel = arrayfun(@(v)sprintf('%.0f',90-abs(v)), ...
                ticks,'UniformOutput',false);
            xlabel(ax,sprintf(['Elevation (°) along cut; left az %.0f°, ' ...
                'right az %.0f°'],mod(cutPhiVal()+180,360),mod(cutPhiVal(),360)));
        else
            xlabel(ax,sprintf('θ (°)   [negative side at φ = %g°]', ...
                mod(cutPhiVal()+180,360)));
        end
    end
    function syncAngleControls()
        if useAzEl()
            [az,el] = azEl(S.theta_s,S.phi_s);
            spTh.Limits = [-180 360];
            spPh.Value = 0; spPh.Limits = [0 90];
            spTh.Value = az; spPh.Value = el;
            lblSteer.Text = 'Steer az / el';
            spTh.Tooltip = 'Steering azimuth in degrees, measured from +x toward +y.';
            spPh.Tooltip = 'Steering elevation in degrees, measured up from the xy plane.';
            ddCutMode.Items = {'Elevation cut (az and az+180)', ...
                'Azimuth cut (fixed elevation)'};
            lblCutPhi.Text = 'Cut azimuth (deg)';
            lblCutTheta.Text = 'Az-cut elevation (deg)';
            bCut0.Text = 'az=0'; bCut90.Text = 'az=90';
            cbCutFollow.Text = '2D cut plane follows steering azimuth';
            cbCutFollow.Tooltip = ['For elevation and U cuts, keep the ' ...
                'observation plane at the steered azimuth. Uncheck to ' ...
                'inspect another plane.'];
            cbFullSphere.Text = 'Show full sphere (elevation -90 to 90)';
            cutLimits = [0 90];
            if S.fullSphere, cutLimits = [-90 90]; end
            spCutTheta.Value = min(max(spCutTheta.Value,cutLimits(1)),cutLimits(2));
            spCutTheta.Limits = cutLimits;
            spCutTheta.Value = 90-S.cutFixedTheta;
        else
            spTh.Value = 0; spTh.Limits = [-90 90];
            spPh.Limits = [-180 360];
            spTh.Value = S.theta_s; spPh.Value = S.phi_s;
            lblSteer.Text = 'Steer theta / phi';
            spTh.Tooltip = 'Signed steering theta in degrees from +z. Negative theta points at phi + 180 degrees.';
            spPh.Tooltip = 'Steering phi in degrees, measured from +x toward +y.';
            ddCutMode.Items = {'Theta cut (fixed phi_s)','Phi cut (fixed theta)'};
            lblCutPhi.Text = 'Cut phi (deg)';
            lblCutTheta.Text = 'Theta for phi-cut (deg)';
            bCut0.Text = 'phi=0'; bCut90.Text = 'phi=90';
            cbCutFollow.Text = '2D cut plane follows steering phi';
            cbCutFollow.Tooltip = ['For theta and U cuts, keep the ' ...
                'observation plane at the steered phi. Uncheck to ' ...
                'inspect another plane.'];
            cbFullSphere.Text = 'Show full sphere (theta 0-180)';
            cutLimits = [0 90];
            if S.fullSphere, cutLimits = [0 180]; end
            spCutTheta.Value = min(max(spCutTheta.Value,cutLimits(1)),cutLimits(2));
            spCutTheta.Limits = cutLimits;
            spCutTheta.Value = S.cutFixedTheta;
        end
        ddAngleConvention.Value = S.angleConvention;
        ddCutMode.Value = S.cutMode;
        spCutPhi.Value = cutPhiVal();
        if isgraphics(pattern2DWin), configurePattern2D(); end
    end
    function assignAngleConvention(v)
        S.angleConvention = v;
        syncAngleControls();
        safeCompute();
        refreshInfo();
        refreshPhaseMaps();
    end
    function ready = patternReady(announce)
        ready = ~S.impNeedsPattern && ...
            ~(strcmp(S.efType,'Imported (CST far-field)') && isempty(S.impFF));
        if ready, return; end
        if isgraphics(overlayWin)
            set(findall(overlayWin,'Tag','overlayNote'),'Text', ...
                'Imported pattern is unavailable; recompute after loading it.');
        end
        S.cutX = []; S.cutDb = []; S.cutRaw=[]; S.radiationValid=false; S.DpkTot = NaN;
        S.glTheta = NaN; S.glPhi = NaN; S.glRelDb = NaN;
        clearCutMetrics();
        blankPlots('IMPORTED PATTERN REQUIRED - results unavailable', ...
            'Import the missing pattern or explicitly select a built-in element');
        setStatus('Imported pattern required: no results until it is loaded.','warn');
        if announce
            uialert(fig, ['This design needs its imported far-field pattern. ' ...
                'Import it, or explicitly choose a built-in element before calculating or exporting phases.'], ...
                'Imported pattern required');
        end
    end

    function assignEF(v)
        if ~strcmp(S.efType,v), elementPreviewPreviousModel = S.efType; end
        S.efType = v;
        syncElementGallery();
        % Selecting the imported element with nothing imported lands in
        % elementFactor's fallback and quietly draws an ISOTROPIC pattern.
        % The config-load path already flags that state; this one did not,
        % so the same substitution was announced on one route and silent
        % on the other.
        S.impNeedsPattern = strcmp(v,'Imported (CST far-field)') && isempty(S.impFF);
        refreshPolAvailability();   % availability follows the element type
        refreshLayout(); didCompute = maybeCompute();
        % The table shows effectivePhaseDeg(), which folds in
        % seqRotSign() -- and THAT changes with the element. Without
        % this the panel kept the old sign's numbers while the
        % pattern beside it used the new one.
        refreshTable();
        % patternReady() clears the plots when an import is missing, but
        % does not update the Metrics panel. Without this, selecting the
        % unavailable imported element leaves the previous gain and lobe
        % numbers visible beside empty plots.
        if ~didCompute, refreshInfo(); end
    end
    function assignEfBeamwidth(which,v)
        if strcmp(which,'az'), S.efBeamAz = v;
        else, S.efBeamEl = v; end
        maybeCompute();
    end

    function showElementGallery()
        % A tiled dropdown on the main figure, anchored below the Element
        % ribbon. It retains the exact reference thumbnails but does not
        % make the user switch to a separate gallery window.
        if ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup)
            shapeGalleryPopup.Visible = 'off';
        end
        if ~isempty(elementGalleryPopup) && isgraphics(elementGalleryPopup)
            if strcmp(elementGalleryPopup.Visible,'on')
                elementGalleryPopup.Visible = 'off';
            else
                positionElementGalleryPopup(elementGalleryPopup,fig,gPatterns);
                elementGalleryPopup.Visible = 'on';
            end
            return;
        end
        [elementGalleryPopup,elementGalleryBtns,elementGalleryModels] = ...
            buildElementGalleryPopup(fig,@ribbonIconFile,@pickElementPattern);
        syncElementGallery();
        positionElementGalleryPopup(elementGalleryPopup,fig,gPatterns);
        elementGalleryPopup.Visible = 'on';
    end

    function pickElementPreview(button)
        v = button.UserData;
        if strcmp(v,S.efType), showElementGallery(); return; end
        ddEF.Value = v;
        assignEF(v);
    end

    function pickElementPattern(button)
        kGallery = find(elementGalleryBtns == button,1);
        if isempty(kGallery), return; end
        v = elementGalleryModels{kGallery};
        elementGalleryPopup.Visible = 'off';
        ddEF.Value = v;
        assignEF(v);
    end

    function syncElementGallery()
        if isgraphics(elementPreviewPrevious) && isgraphics(elementPreviewCurrent)
            if strcmp(elementPreviewPreviousModel,S.efType)
                elementPreviewPreviousModel = 'Gaussian';
                if strcmp(S.efType,'Gaussian'), elementPreviewPreviousModel = 'Isotropic'; end
            end
            setElementPreviewButton(elementPreviewPrevious,elementPreviewPreviousModel,@ribbonIconFile);
            setElementPreviewButton(elementPreviewCurrent,S.efType,@ribbonIconFile);
            paintSelected([elementPreviewPrevious elementPreviewCurrent],2);
        end
        if ~isempty(elementGalleryPopup) && isgraphics(elementGalleryPopup) && ...
                ~isempty(elementGalleryBtns) && ...
                all(isgraphics(elementGalleryBtns))
            paintSelected(elementGalleryBtns,find(strcmp(elementGalleryModels,S.efType),1));
        end
    end

    function assignCustomFormula(v)
        % refreshTable() for the same reason as assignEF: a custom formula
        % IS the element, so it can move seqRotSign() and with it every
        % Feed phase the table prints.
        S.customFormula = v; warnIfBadFormula(v,'E_theta');
        refreshTable(); maybeCompute();
    end
    function assignCustomFormulaPh(v)
        S.customFormulaPh = v; warnIfBadFormula(v,'E_phi');
        refreshTable(); maybeCompute();
    end
    function warnIfBadFormula(expr, whichOne)
        %WARNIFBADFORMULA  Say so when a typed formula cannot be used.
        %   elementFactor falls back to isotropic (E_theta) or zero
        %   (E_phi) when eval fails, which is the right thing to compute
        %   but the wrong thing to do silently: a mistyped formula gave a
        %   completely different antenna with nothing on screen to say
        %   so, and every number in the window then described that other
        %   antenna. The taper path already warns in exactly this
        %   situation (see assignTaper's toolbox check).
        %
        %   Checked HERE, once, where the user just typed it -- not
        %   inside elementFactor, which runs on every recompute and would
        %   pop the same dialog over and over. Evaluated on a tiny
        %   stand-in grid, since only whether it RUNS and what SHAPE it
        %   returns is in question, not its values.
        % The stand-in grid INCLUDES theta = 90 and 180, where ct is
        % clamped to zero. Without those the check would miss exactly the
        % formulas that matter -- 1./ct and ct./ct are perfectly finite
        % at 0/30/45/60 and only blow up at grazing and behind.
        TH = [0 45; 90 180]; PH = [0 90; 180 270];
        th = TH; ph = PH; st = sind(TH); ct = max(cosd(TH),0); %#ok<NASGU>
        theta = TH; phi = PH; %#ok<NASGU>
        bad = '';
        try
            V = eval(expr);
            V = formulaFieldValue(V,size(TH)); %#ok<NASGU>
        catch err
            bad = err.message;
        end
        if isempty(bad), return; end

        uialert(fig, sprintf(['The %s formula "%s" cannot be used:\n\n%s\n\n' ...
            'Calculation is refused until the formula is corrected.\n\nAvailable variables: th, ph (degrees), ' ...
            'ct = max(cos(th),0), st = sin(th). Use element-wise operators ' ...
            '(.* ./ .^) and avoid transposes.'], ...
            whichOne, expr, bad), 'Formula not usable', 'Icon','warning');
    end
    function assignQ(v),     S.efQ = v;    maybeCompute(); end
    function assignFullSphere(v)
        S.fullSphere = v;
        % The phi-cut's fixed-theta spinner may reach into the back
        % hemisphere only while the back hemisphere is actually being
        % shown. Order matters: MATLAB errors if a spinner's current
        % Value falls outside the Limits being assigned, so when
        % NARROWING the range the value has to be brought inside first
        % (and S.cutFixedTheta kept in step, since setting .Value
        % programmatically does not fire its ValueChangedFcn).
        if ~v, S.cutFixedTheta = min(S.cutFixedTheta,90); end
        syncAngleControls();
        safeCompute();
    end

    function s = polTag()
        % Short marker naming which polarization the plotted curve is,
        % for titles. Empty for the default total-power readout so the
        % common case stays uncluttered.
        %
        % This matters because HPBW / FNBW / SLL are all measured off
        % whatever polCombine returned, so they change with the
        % dropdown. On a linear array the circular components are just
        % the total scaled by 1/sqrt(2), so the shape -- and every
        % metric -- is unchanged. On a CP array they are genuinely
        % different patterns: a 2x2 sequential-rotation design measured
        % HPBW 49.0 / FNBW 180 on Total but HPBW 31.0 / FNBW 91.3 /
        % SLL 0.0 dB on LHCP, because LHCP there is the SUPPRESSED
        % sense and those numbers describe cross-pol leakage, not a
        % beam. Without this tag a screenshot of that gives no clue.
        switch ddPol.Value
            case 'RHCP component', s = 'RHCP';
            case 'LHCP component', s = 'LHCP';
            otherwise,             s = '';
        end
    end

    function tf = magOnlyMixedRot()
        %MAGONLYMIXEDROT  The one combination that makes an array result
        %   not merely approximate but WRONG.
        %
        %   A magnitude-only ("Abs") CST export carries no split between
        %   the two field components, so this app has to put the whole
        %   magnitude into E_theta with E_phi = 0 -- an invented linear
        %   polarisation. That is harmless while every element shares one
        %   orientation, because the same invented component then
        %   multiplies the array factor everywhere and divides out of
        %   every ratio.
        %
        %   It stops being harmless the moment elements carry DIFFERENT
        %   rotations. Rotating an element rotates its polarisation, and
        %   the coherent sum then combines components that do not exist
        %   as written. Measured on the analytic Patch exported both
        %   ways: with the Theta/Phi split the 8x8 sequentially-rotated
        %   array reads 23.29 dBi, matching the analytic element exactly;
        %   the SAME element exported Abs-only reads 14.70 dBi, and in
        %   unit-cell mode it collapses to the -120 dB build-up floor.
        %   Nothing about the plot looks broken -- it just quietly loses
        %   8 to 138 dB.
        %
        %   The import already warns about this once, at load. That is
        %   not enough: the rotation is usually applied later, long after
        %   the dialog was dismissed, so the warning is gone by the time
        %   the state it describes actually exists. This predicate is
        %   what lets the plots say so at the moment it is true.
        tf = strcmp(S.efType,'Imported (CST far-field)') && ~isempty(S.impFF) ...
            && isfield(S.impFF,'noComponents') && S.impFF.noComponents ...
            && ~isempty(S.el) && numel(unique(S.el(:,5))) > 1;
    end

    function tf = impBasisOK()
        % False only when the CURRENT element factor is an imported file
        % whose two field components are not the theta/phi spherical
        % pair. Everything analytic in this app is defined directly in
        % theta/phi, so those are always fine.
        tf = true;
        if strcmp(S.efType,'Imported (CST far-field)') && ~isempty(S.impFF) ...
                && isfield(S.impFF,'compThetaPhi') && ~S.impFF.compThetaPhi ...
                && ~isempty(S.impFF.compNames{1})
            tf = false;
        end
    end

    function assignTotEff(v), S.impTotEffPct = v; refreshAll(); end
    function assignUnitCell(v)
        S.impUnitCell = v;
        refreshAll();
    end

    function e = impEffLin()
        % Linear efficiency factor that converts a computed DIRECTIVITY
        % into real GAIN. Only an imported CST pattern carries this (in
        % its absolute levels -- see ffS.eff); every analytic element
        % factor in this app is a pure normalized pattern shape with no
        % loss model attached, so 1.0 (gain == directivity) there.
        % Deliberately NOT applied to the Array-factor-only display: the
        % AF is the isotropic-element array pattern, a mathematical
        % construct with no physical antenna to be lossy.
        e = 1;
        if ~strcmp(S.efType,'Imported (CST far-field)') || isempty(S.impFF)
            return;
        end
        % A unit-cell file's peak IS the embedded element GAIN -- realized
        % gain, with mismatch and loss already applied (its peak sits
        % 10log10(eff) below the lossless 4*pi*A_cell/lambda^2 ceiling).
        % Applying any efficiency on top would count it twice, so neither
        % the file-derived value nor the typed override is used here.
        if S.impUnitCell
            return;
        end
        e = NaN;
        if isfield(S.impFF,'eff') && isfinite(S.impFF.eff) && S.impFF.eff > 0
            e = S.impFF.eff;
        end
    end
    function refreshPolAvailability()
        %REFRESHPOLAVAILABILITY  Gate every polarization readout at once.
        %   Axial ratio and the RHCP/LHCP split are derived from the
        %   theta/phi components. A magnitude-only import has no such
        %   split (the whole magnitude was put into E_theta), and a file
        %   stored in another orthogonal basis has one that means
        %   something else -- in both cases these readouts produce
        %   confident numbers from a polarisation the file never carried.
        %
        %   ONE function, called from the element-type change, the
        %   import, and the config load. The first attempt gated these
        %   only inside the import path, so switching from a
        %   magnitude-only file to a built-in Patch left them disabled
        %   with nothing to turn them back on -- a control that can never
        %   be re-enabled is a worse failure than the one being guarded
        %   against. The 3D "Axial Ratio (dB)" surface option is included
        %   here too; it was reachable for exactly the imports the other
        %   two were being hidden for.
        ok = true;
        if strcmp(S.efType,'Imported (CST far-field)') && ~isempty(S.impFF)
            ok = S.impFF.compThetaPhi && ~S.impFF.noComponents;
        end
        if ok
            ddPol.Enable = 'on'; cbARCut.Enable = 'on'; ddShow.Enable = 'on';
            % Put the AR entry back. Items first, Value second: the
            % order only matters in the other branch, but keeping both
            % assignments in the same order in both places is what stops
            % a later edit from reintroducing the throw.
            if ~ismember('Axial Ratio (dB)', ddShow.Items)
                ddShow.Items = SHOW_ITEMS;
            end
        else
            ddPol.Value = 'Total (any pol)'; ddPol.Enable = 'off';
            cbARCut.Value = false;           cbARCut.Enable = 'off';
            % The surface dropdown stays usable -- only its AR entry is
            % unsupported -- so the entry is REMOVED rather than the whole
            % control being locked.
            %
            % Moving the selection off the entry was not enough. The item
            % stayed in the list, so the user could simply pick it again
            % straight after the import and get a full Axial Ratio surface
            % computed from an invented polarisation split -- and for a
            % file whose component columns are unnamed, impBasisOK() is
            % true, so it rendered with no warning at all. Taking the item
            % out is what actually makes the state unreachable.
            %
            % It also closes the config-load route: cfgItem() validates a
            % saved view against ddShow.Items, so with the entry gone a
            % config carrying vShow = 'Axial Ratio (dB)' can no longer
            % restore it either.
            %
            % Value must be moved BEFORE Items shrinks -- a dropdown
            % throws if its Value is not in its Items.
            if strcmp(ddShow.Value,'Axial Ratio (dB)')
                ddShow.Value = 'Total (EF x AF)';
            end
            ddShow.Items = SHOW_ITEMS(~strcmp(SHOW_ITEMS,'Axial Ratio (dB)'));
        end
    end


    function [ph, phWithRot, rotDeg, swapHV, hSide, vSide, feedOffset] = portPhaseFor(scheme)
        %PORTPHASEFOR  Phase driven into each feed, for one feed scheme.
        %
        %   Returns one column per PORT: [nEl x 1] for a single-feed
        %   element, [nEl x 2] (H then V) for a dual-fed one.
        %
        %   ph         the SCHEME phase, always WITHOUT the rotation
        %              term -- see below.
        %   phWithRot  the same plus each element's rotation feed phase.
        %   rotDeg     that rotation contribution on its own [nEl x 1].
        %   swapHV     [nEl x 1] logical, true where the DECLARED feed layout
        %              assigns H to the even physical port rather than the odd one.
        %
        %   Every scheme starts from the same array phase -- the
        %   progressive phase that steers the beam, plus any manual
        %   per-element offset. The scheme only adds to it.
        %
        %   ROTATION IS DELIBERATELY EXCLUDED from ph, whatever the
        %   global extra-rotation checkbox says. For the dual-feed map the
        %   declared L/R/B/T topology and the LP/CP scheme already define
        %   the physical port reference corrections. Adding S.el(:,5)
        %   again double-counts that sequential-rotation convention.
        %   rotDeg/phWithRot are retained only as diagnostic alternatives
        %   for the general pattern/single-feed model.
        %
        %   SINGLE FEED. One port, nothing to add. This is the plain
        %   steered array.
        %
        %   DUAL FEED, LP. Alternate elements are laid out mirrored so
        %   their probe-induced cross-polar contributions cancel in the
        %   far field. Mirroring also reverses the co-polar field, so the
        %   feed puts it back with 180 deg. H mirrors across COLUMNS and V
        %   across ROWS, which is why the two ports do not share a
        %   correction.
        %
        %   DUAL FEED, CP. The LP scheme with a uniform 90 deg added to V
        %   at every element. Uniform, so it changes the polarisation
        %   without moving the beam.
        %
        %   Wrapped to [0,360): that is the range a phase shifter is set
        %   in, and negative feed phases are awkward to transcribe.
        if isempty(S.el)
            ph = zeros(0,1); phWithRot = zeros(0,1);
            rotDeg = zeros(0,1); swapHV = false(0,1);
            hSide = zeros(0,1); vSide = zeros(0,1); feedOffset = zeros(0,1);
            return;
        end
        % Array phase WITHOUT rotation, regardless of S.seqPhase -- the
        % same expression effectivePhaseDeg() evaluates with its rotation
        % term switched off. Written out rather than calling that
        % function so this map does not change under a checkbox that is
        % about the pattern, not about the feed scheme.
        x = S.el(:,1); y = S.el(:,2); ph0 = S.el(:,4);
        us = sind(S.theta_s)*cosd(S.phi_s);
        vs = sind(S.theta_s)*sind(S.phi_s);
        base = rad2deg(-S.k*phaseFreqRatio()*(x*us + y*vs) + deg2rad(ph0));
        % The feed phase the app ACTUALLY applies, sign included -- not the
        % raw rotation angle. seqRotSign() picks -1 for an LHCP-dominant
        % element, and this window is where the numbers get transcribed
        % into a solver, so printing +90 where the pattern was computed
        % with -90 would hand over a design that does not match the plot
        % beside it.
        %
        % Note this is the optional rotation-derived FEED-PHASE term only.
        % Which physical port carries H is NOT inferred from S.el(:,5); it is
        % read from the declared S.portMap below, because mirrored feed
        % topology and physical rotation are separate pieces of information.
        rotDeg = mod(seqRotSign()*S.el(:,5),360);
        % WHICH PHYSICAL PORT CARRIES H, and which side each feed is on:
        % read straight out of the DECLARED repeating unit (S.portMap),
        % tiled over the array. Nothing here is inferred from positions
        % or rotation angles -- see S.portMap's note for the measured
        % layout that ruled that out.
        % Tile the declared unit at ITS OWN period, which need not be 2.
        [nbR,nbC] = size(S.portMap.hOdd);
        bi = sub2ind([nbR nbC], mod(S.elRC(:,1)-1,nbR)+1, mod(S.elRC(:,2)-1,nbC)+1);
        swapHV = ~S.portMap.hOdd(bi);   % H on the EVEN port = swapped
        swapHV = swapHV(:);
        hSide  = S.portMap.hSide(bi); hSide = hSide(:);
        vSide  = S.portMap.vSide(bi); vSide = vSide(:);
        if strcmp(scheme,'single')
            feedOffset = ph0;
            ph = mod(base,360);
            phWithRot = mod(base + rotDeg,360);
            swapHV = false(size(swapHV));   % one port: nothing to exchange
            return;
        end
        % The 180 deg mirror corrections are DERIVED from the declared
        % feed sides, not from column/row index.
        %
        % They exist because a feed on the opposite edge presents a
        % reversed co-polar field, so the element needs 180 deg to put it
        % back. That is a statement about WHERE THE FEED IS -- which is
        % now declared -- and it used to be hardcoded as 180*colPar /
        % 180*rowPar, applied whatever the layout said. The two disagreed
        % as soon as the layout was declared: a unit with every H feed on
        % the same edge is not mirrored at all, yet still got 0/180/0...
        % across columns, which cancels at broadside for elements that
        % are in truth identically fed. The default unit declares the
        % alternating sides, so the historical scheme is reproduced
        % exactly -- it is now a consequence of the declaration instead
        % of an assumption underneath it.
        %
        % Fixed physical-side convention used by this antenna/CST map:
        %   H-left = 0 deg, H-right = +180 deg
        %   V-bottom = 0 deg, V-top = +180 deg
        % Do NOT reference these corrections to whichever cell happens to
        % appear first in the repeating unit. Doing that independently for H
        % and V can introduce a GLOBAL 180-deg shift on only one polarization
        % when the first cell is edited, which flips the CP H/V relationship.
        hFlip = 180*(hSide == 2);   % right relative to left
        vFlip = 180*(vSide == 2);   % top relative to bottom
        h = base + hFlip;
        v = base + vFlip;
        if strcmp(scheme,'cp'), v = v + 90; end
        feedOffset = [ph0+hFlip, ph0+vFlip+90*strcmp(scheme,'cp')];
        ph = mod([h v],360);
        phWithRot = mod([h v] + rotDeg,360);
    end

    function showPortPhases()
        if ~patternReady(true), return; end
        if isempty(S.el), uialert(fig,'Place at least one element first.','No elements'); return; end
        if ~isempty(resultPanels{2}) && isgraphics(resultPanels{2})
            selectTask('Phase Map');
            refreshPhaseMaps();
            return;
        end
        mapSchemes={'Single feed','Dual feed LP','Dual feed CP'};
        mapFig=resultPane('Phase Map','phaseMapView');
        mapOuter=uigridlayout(mapFig,[2 1]); mapOuter.RowHeight={36,'1x'};
        mapTop=uigridlayout(mapOuter,[1 5]); mapTop.Layout.Row=1;
        mapTop.ColumnWidth={80,190,150,110,'1x'}; mapTop.Padding=[4 0 4 0];
        uilabel(mapTop,'Text','Feed scheme');
        mapScheme=uidropdown(mapTop,'Items',mapSchemes,'Value',mapSchemes{1}, ...
    'Tag','mapScheme','Tooltip','Physical feed scheme: single feed, dual linear polarization, or dual circular polarization.');
        mapScheme.Tooltip='Port excitation scheme; does not change the imported element pattern.';
        mapEdit=uibutton(mapTop,'Text','Edit feed layout...', ...
    'Tag','mapEdit','Tooltip','Edit feed layout. Apply this setting or action to the current design.');
        mapCSV=uibutton(mapTop,'Text','Export CSV', ...
    'Tag','mapCSV','Tooltip','Export CSV. Apply this setting or action to the current design.');
        mapCST=uibutton(mapTop,'Text','Export CST macro', ...
    'Tag','mapCST','Tooltip','Export CST macro. Apply this setting or action to the current design.');
        mapBody=uigridlayout(mapOuter,[1 2]); mapBody.Layout.Row=2;
        mapBody.ColumnWidth={'1.8x','1x'};
        mapAxes=uiaxes(mapBody); mapAxes.Layout.Row=1; mapAxes.Layout.Column=1;
        mapTable=uitable(mapBody, ...
    'Tag','mapTable','Tooltip','Review the displayed values; editable cells update the current settings.'); mapTable.Layout.Row=1; mapTable.Layout.Column=2;
        mapTable.ColumnWidth={75,'auto','auto'};
        mapTable.RowName={};
        mapTable.Tooltip=['Final scheme phases. Dual-feed LP/CP includes the declared opposite-side 180-degree correction, ' ...
            'the LP/CP port relation, steering/manual phase and negative-weight phase shifts, and deliberately ignores the extra rotation term. ' ...
            'Single-feed follows the extra rotation-phase checkbox, matching the main pattern and CST single-feed export.'];
        fullMapTable=table();
        mapScheme.ValueChangedFcn=@(~,~)redrawMap();
        mapEdit.ButtonPushedFcn=@(~,~)editPortMap();
        mapCSV.ButtonPushedFcn=@(~,~)exportMapCSV();
        mapCST.ButtonPushedFcn=@(~,~)exportCstExcitation(mapKey(),fig,true);
        S.portMapRefresh{end+1}=@redrawMap;
        mapFig.UserData = numel(S.portMapRefresh);
        redrawMap();
        function key=mapKey()
            if strcmp(mapScheme.Value,mapSchemes{1})
                key='single';
            elseif strcmp(mapScheme.Value,mapSchemes{2})
                key='lp';
            else
                key='cp';
            end
        end
        function exportMapCSV()
            redrawMap();
            exportPortCSV(fullMapTable,fig);
        end
        function redrawMap()
            if ~isgraphics(mapFig), return; end
            P=pal();
            fullMapTable=table(); mapTable.Data=[]; cla(mapAxes);
            colorbar(mapAxes,'off');   % back below, once there are tiles
            if isempty(S.el), title(mapAxes,'No elements placed'); return; end
            if ~patternReady(false), title(mapAxes,'Import an element pattern first'); return; end
            key=mapKey();
            [basePh,rotPh,~,swapped,hSides,vSides]=portPhaseFor(key);
            % Dual-feed LP/CP is defined by the declared physical L/R/B/T
            % topology plus the LP/CP scheme itself, so adding the element
            % rotation angle again would double-count sequential rotation.
            %
            % Single-feed is different: there is no dual-feed side/CP scheme,
            % so the optional "extra rot angle -> feed phase" setting must be
            % reflected here exactly as it is in the main pattern and the
            % single-feed CST export.
            useMapRotation = strcmp(key,'single') && S.seqPhase;
            finalPh=appliedPortPhases(basePh,rotPh,S.el(:,3),useMapRotation);
            if any(~isfinite(finalPh(:)))
                title(mapAxes,'Phase calculation unavailable'); return;
            end
            rowIds=S.elRC(:,1); colIds=S.elRC(:,2); count=size(S.el,1);
            nRows=max(rowIds); nCols=max(colIds);
            [~,~,cellIds]=unique([rowIds colIds],'rows');
            cellCount=accumarray(cellIds,1); firstInCell=accumarray(cellIds,(1:count)',[],@min);
            axis(mapAxes,'equal'); hold(mapAxes,'on');
            mapAxes.XTick=[]; mapAxes.YTick=[];
            % Transparent, so the tiles sit on the window's own background
            % in either theme; the tile edges use that same colour.
            mapAxes.Color='none'; mapAxes.XColor='none'; mapAxes.YColor='none';
            xlim(mapAxes,[0.4 nCols+0.6]); ylim(mapAxes,[0.4 nRows+0.6]);
            % The tiles are coloured from the same cyclic map the colour
            % bar shows, so the bar reads them directly. Dual-feed tiles
            % are coloured by the H port.
            colormap(mapAxes,P.phaseMap); clim(mapAxes,[0 360]);
            mapBar=colorbar(mapAxes,'Ticks',0:90:360);
            if strcmp(key,'single'), mapBar.Label.String='Phase (°)';
            else, mapBar.Label.String='H-port phase (°)'; end
            nMapCol=size(P.phaseMap,1);
            % Map texts give angles in degrees with the ° sign;
            % directionText's 'deg' wording is shared with other readouts.
            degText=@(th,ph) strrep(directionText(th,ph),' deg','°');
            if S.retunePhase
                mapFreqTxt=sprintf(['Phases retuned at %.4g GHz (operating ' ...
                    '%.4g GHz): the steering phase is recalculated for the ' ...
                    'commanded angle.'],phaseReferenceGHz(),S.freqOpGHz);
            elseif size(S.el,1) <= 1
                % A single element has no progressive array phase, hence no
                % array-factor beam squint. Its element pattern may still be
                % frequency dependent, but that is a different effect.
                mapFreqTxt=sprintf(['Phases fixed at %.4g GHz, operating ' ...
                    '%.4g GHz. A single element has no array-factor beam ' ...
                    'squint.'],phaseReferenceGHz(),S.freqOpGHz);
            else
                sMap=squintSinTheta(S.theta_s);
                if abs(sMap)<=1
                    mapFreqTxt=sprintf(['Beam squint: phases fixed at %.4g ' ...
                        'GHz, operating %.4g GHz, so the steering term ' ...
                        'points to %s.'],phaseReferenceGHz(),S.freqOpGHz, ...
                        degText(asind(sMap),S.phi_s));
                else
                    mapFreqTxt=sprintf(['Beam squint: phases fixed at %.4g ' ...
                        'GHz, operating %.4g GHz; the steering term ' ...
                        'predicts a beam outside visible space.'], ...
                        phaseReferenceGHz(),S.freqOpGHz);
                end
            end
            [unitRows,unitCols]=size(S.portMap.hOdd);
            if S.portMapSet, layoutStatus='Declared'; else, layoutStatus='Default: check your CST port layout'; end
            if S.portMapSet, layoutTxt='declared';
            else, layoutTxt='default, check it against your CST ports'; end
            if strcmp(key,'single')
                sideLegend=sprintf(['Single-feed map: the H/V feed-side ' ...
                    'correction is not used. Feed layout: %s.'],layoutTxt);
            else
                sideLegend=sprintf(['H (L/R) and V (B/T) name the physical ' ...
                    'feed edge. Opposite-side feeds carry an ASSUMED 180° ' ...
                    'correction, and the dual-feed scheme does not add the ' ...
                    'rotation angle again: verify against the CST port ' ...
                    'references. Feed layout (%d × %d unit): %s.'], ...
                    unitRows,unitCols,layoutTxt);
            end
            title(mapAxes,sprintf('%s · %d × %d · steered to %s', ...
                mapScheme.Value,S.M,S.N,degText(S.theta_s,S.phi_s)));
            mapScheme.Tooltip=[mapFreqTxt newline sideLegend];
            mapEdit.Tooltip=sprintf('%s %dx%d layout. H side is Left/Right; V side is Bottom/Top. Opposite sides use an assumed 180-degree reference correction; verify this against the actual CST port reference phases.',layoutStatus,unitRows,unitCols);
            showLabels=max(nRows,nCols)<=12 && count<=144;
            for elementId=1:count
                groupId=cellIds(elementId);
                if elementId~=firstInCell(groupId), continue; end
                rr=rowIds(elementId); cc=colIds(elementId);
                if cellCount(groupId)>1
                    % Several elements share the cell, so no single phase
                    % describes it: a neutral tile rather than a hue.
                    shade=P.cellNeutral;
                else
                    colIdx=1+floor(mod(finalPh(elementId,1),360)/360*nMapCol);
                    shade=P.phaseMap(min(colIdx,nMapCol),:);
                end
                ink=inkOn(P,shade);
                rectangle(mapAxes,'Position',[cc-.46 rr-.46 .92 .92], ...
                    'FaceColor',shade,'EdgeColor',P.panelBg,'LineWidth',1);
                if ~showLabels, continue; end
                if cellCount(groupId)>1
                    text(mapAxes,cc,rr,sprintf('%d elements',cellCount(groupId)), ...
                        'HorizontalAlignment','center','FontSize',8,'Color',ink,'Interpreter','none');
                elseif strcmp(key,'single')
                    text(mapAxes,cc,rr,sprintf('P%d  %g°',elementId,round(finalPh(elementId),1)), ...
                        'HorizontalAlignment','center','FontSize',9,'Color',ink,'Interpreter','none');
                else
                    hp=2*elementId-1+double(swapped(elementId)); vp=2*elementId-double(swapped(elementId));
                    hSideTxt = sideName(hSides(elementId),'h');
                    vSideTxt = sideName(vSides(elementId),'v');
                    % Always show the PHYSICAL CST port side directly on the map.
                    % Example: H17(L) -> 0 deg, V18(T) -> 180 deg.
                    % H uses L/R; V uses B/T. The displayed phase already contains
                    % the declared opposite-side 180-degree correction.
                    text(mapAxes,cc,rr+.14, ...
                        sprintf('H%d(%s)  ->  %g°',hp,hSideTxt,round(finalPh(elementId,1),1)), ...
                        'HorizontalAlignment','center','FontSize',8,'Color',ink,'Interpreter','none');
                    text(mapAxes,cc,rr-.14, ...
                        sprintf('V%d(%s)  ->  %g°',vp,vSideTxt,round(finalPh(elementId,2),1)), ...
                        'HorizontalAlignment','center','FontSize',8,'Color',ink,'Interpreter','none');
                end
            end
            hold(mapAxes,'off');
            elementIds=(1:count)'; amps=abs(calculationWeights());
            if strcmp(key,'single')
                ports=elementIds; feeds=repmat("Single",count,1); sides=repmat("",count,1);
                idx=elementIds;
            else
                ports=[2*elementIds-1+double(swapped);2*elementIds-double(swapped)];
                feeds=[repmat("H",count,1);repmat("V",count,1)];
                sides=[string(arrayfun(@(v)sideName(v,'h'),hSides,'UniformOutput',false)); ...
                    string(arrayfun(@(v)sideName(v,'v'),vSides,'UniformOutput',false))];
                idx=[elementIds;elementIds];
            end
            phaseRefCol=repmat(phaseReferenceGHz(),numel(idx),1);
            opFreqCol=repmat(S.freqOpGHz,numel(idx),1);
            squintCol=repmat(~S.retunePhase,numel(idx),1);
            fullMapTable=table(ports,idx,feeds,finalPh(:),rowIds(idx),colIds(idx), ...
                amps(idx),phaseRefCol,opFreqCol,squintCol,basePh(:),rotPh(:), ...
                repmat(S.seqPhase,numel(idx),1),repmat(useMapRotation,numel(idx),1),sides,180*(S.el(idx,3)<0), ...
                'VariableNames',{'Port','Element','Feed','Phase_deg','Lattice_Row','Lattice_Col', ...
                'Amplitude','Phase_reference_GHz','Operating_frequency_GHz','Beam_squint', ...
                'Scheme_phase_deg','With_rotation_deg','Extra_rotation_checkbox','Rotation_applied_to_scheme', ...
                'Feed_side','Weight_phase_deg'});
            fullMapTable=sortrows(fullMapTable,'Port');
            % One row per physical element; finalPh columns are H then V,
            % independently of the declared CST port numbering.
            if strcmp(key,'single')
                mapTable.Data=[elementIds,finalPh(:,1)];
                mapTable.ColumnName={'Element','Phase (°)'};
                mapTable.ColumnWidth={75,'auto'};
            else
                % Show hardware port identity and physical feed edge together
                % with the final phase. A phase value by itself is ambiguous
                % for mirrored elements because opposite-side feeds need a
                % 180-degree reference correction.
                hpVec = 2*elementIds-1+double(swapped);
                vpVec = 2*elementIds-double(swapped);
                hSideTxt = string(arrayfun(@(v)sideName(v,'h'),hSides,'UniformOutput',false));
                vSideTxt = string(arrayfun(@(v)sideName(v,'v'),vSides,'UniformOutput',false));
                mapTable.Data = table(elementIds,hpVec,hSideTxt,finalPh(:,1), ...
                    vpVec,vSideTxt,finalPh(:,2), ...
                    'VariableNames',{'Element','H_Port','H_Side','H_Phase_deg', ...
                    'V_Port','V_Side','V_Phase_deg'});
                mapTable.ColumnName={'Element','H port','H side','H phase (°)', ...
                    'V port','V side','V phase (°)'};
                mapTable.ColumnWidth={60,55,55,82,55,55,82};
            end
        end
    end

    function phaseValues=appliedPortPhases(basePh,rotPh,weights,useRotation)
        if useRotation, phaseValues=rotPh; else, phaseValues=basePh; end
        phaseValues=mod(phaseValues+180*(weights<0),360);
    end

    function nm = sideName(code,which)
        %SIDENAME  Human name for a feed-side code.
        %   1/2 mean left/right for an H feed and bottom/top for a V feed,
        %   so the same code reads differently depending on which feed it
        %   describes -- hence the second argument.
        if which == 'h'
            if code == 2, nm = 'R'; else, nm = 'L'; end
        else
            if code == 2, nm = 'T'; else, nm = 'B'; end
        end
    end

    function editPortMap()
        %EDITPORTMAP  Declare the repeating unit's physical feed layout.
        %
        %   One displayed row per cell of the user-sized repeating unit.
        %   Everything here is a fact about the CST model that the tool
        %   cannot see, so it is typed in rather than guessed -- an earlier
        %   version inferred it from column parity and rotation and came
        %   out exactly inverted against a real layout.
        fe = uifigure('Name','Feed layout of the repeating unit', ...
            'Position',[200 200 640 330]);
        matchTheme(fe);
        gl = uigridlayout(fe,[4 1]); gl.RowHeight = {96,32,'1x',36};
        uilabel(gl,'WordWrap','on','Text', ...
            ['One row per cell of the unit that repeats across the array. ' ...
             'In this editor Block_Row is numbered from TOP to BOTTOM, so ' ...
             'Block_Row 1 / Block_Col 1 is the top-left cell. The app converts ' ...
             'this automatically to its internal bottom-up lattice row numbering.' newline newline ...
             'H port: which of an element''s two ports (odd = 2n-1, even = 2n) carries H. ' ...
             'H side / V side: which physical edge that feed sits on. These are read from your CST ' ...
             'model and are independent of the optional extra electrical rotation-phase term.']);
        % Period of the repeating unit. Asked for, never inferred: three
        % elements of a real row cannot establish it, and a wrong period
        % silently forces distant elements to share a declaration.
        szBar = uigridlayout(gl,[1 5]);
        szBar.ColumnWidth = {150,70,150,70,'1x'}; szBar.Padding=[0 0 0 0];
        uilabel(szBar,'Text','Repeat every (rows)');
        spBR = uispinner(szBar,'Limits',[1 8],'Value',size(S.portMap.hOdd,1), ...
            'Step',1,'RoundFractionalValues','on', ...
    'Tag','spBR','Tooltip','Number of rows in the repeating physical feed-mapping block.');
        uilabel(szBar,'Text','and (columns)');
        spBC = uispinner(szBar,'Limits',[1 8],'Value',size(S.portMap.hOdd,2), ...
            'Step',1,'RoundFractionalValues','on', ...
    'Tag','spBC','Tooltip','Number of columns in the repeating physical feed-mapping block.');
        uilabel(szBar,'Text','');
        t = uitable(gl, ...
    'Tag','t','Tooltip','Review the displayed values; editable cells update the current settings.');
        fillTable(S.portMap);
        spBR.ValueChangedFcn = @(s,e)resizeUnit();
        spBC.ValueChangedFcn = @(s,e)resizeUnit();
        bar = uigridlayout(gl,[1 3]); bar.ColumnWidth = {'1x',110,110};
        uilabel(bar,'Text','');
        bCan = uibutton(bar,'Text','Cancel', ...
    'Tag','bCan','Tooltip','Cancel. Apply this setting or action to the current design.');
        bOk  = uibutton(bar,'Text','Apply', ...
    'Tag','bOk','Tooltip','Apply. Apply this setting or action to the current design.');
        bCan.ButtonPushedFcn = @(s,e)close(fe);
        bOk.ButtonPushedFcn  = @(s,e)applyIt();
        function fillTable(pm)
            [nR,nC] = size(pm.hOdd);
            % Display rows TOP-to-BOTTOM to match the physical/CST drawing.
            % Internally S.elRC and portMap rows are BOTTOM-to-TOP.
            [displayRow,cc] = ndgrid(1:nR,1:nC); nn = nR*nC;
            internalRow = nR - displayRow + 1;
            hp = strings(nn,1); hs = strings(nn,1); vs = strings(nn,1);
            for k = 1:nn
                ir = internalRow(k); ic = cc(k);
                hp(k) = ternStr(pm.hOdd(ir,ic),"odd (2n-1)","even (2n)");
                hs(k) = ternStr(pm.hSide(ir,ic)==2,"right","left");
                vs(k) = ternStr(pm.vSide(ir,ic)==2,"top","bottom");
            end
            t.Data = table(displayRow(:), cc(:), categorical(hp,["odd (2n-1)","even (2n)"]), ...
                categorical(hs,["left","right"]), categorical(vs,["bottom","top"]), ...
                'VariableNames',{'Block_Row','Block_Col','H_port','H_side','V_side'});
            t.ColumnEditable = [false false true true true];
        end
        function resizeUnit()
            % Grow/shrink while KEEPING the declarations in the same
            % TOP-DOWN Block_Row positions shown to the user. readTable()
            % returns bottom-up internal rows, so tiling that matrix directly
            % makes a 2->3 row resize visibly swap the old top and bottom
            % rows. Tile in DISPLAY coordinates, then convert each destination
            % row back to the internal bottom-up convention.
            nR = spBR.Value; nC = spBC.Value;
            cur = readTable();
            pm.hOdd  = false(nR,nC); pm.hSide = ones(nR,nC); pm.vSide = ones(nR,nC);
            [oR,oC] = size(cur.hOdd);
            for displayR = 1:nR
                srcDisplayR = mod(displayR-1,oR)+1;
                srcInternalR = oR-srcDisplayR+1;
                dstInternalR = nR-displayR+1;
                for c0 = 1:nC
                    sc = mod(c0-1,oC)+1;
                    pm.hOdd(dstInternalR,c0)  = cur.hOdd(srcInternalR,sc);
                    pm.hSide(dstInternalR,c0) = cur.hSide(srcInternalR,sc);
                    pm.vSide(dstInternalR,c0) = cur.vSide(srcInternalR,sc);
                end
            end
            fillTable(pm);
        end
        function pm = readTable()
            D = t.Data;
            nR = max(D.Block_Row); nC = max(D.Block_Col);
            pm.hOdd = false(nR,nC); pm.hSide = ones(nR,nC); pm.vSide = ones(nR,nC);
            for q = 1:height(D)
                displayRow = D.Block_Row(q); c0 = D.Block_Col(q);
                r0 = nR - displayRow + 1; % top-down display -> bottom-up internal
                pm.hOdd(r0,c0)  = (string(D.H_port(q)) == "odd (2n-1)");
                pm.hSide(r0,c0) = 1 + (string(D.H_side(q)) == "right");
                pm.vSide(r0,c0) = 1 + (string(D.V_side(q)) == "top");
            end
        end
        function applyIt()
            S.portMap = readTable();
            S.portMapSet = true;
            close(fe);
            refreshPhaseMaps();
            undoCheckpoint();   % no compute follows, so no checkpoint either
        end
    end

    function outStr = ternStr(c,a,b)
        if c, outStr = a; else, outStr = b; end
    end

    function exportPortCSV(T,parentFig)
        if ~isgraphics(parentFig), return; end
        if ~patternReady(false) || isempty(T) || ~any(S.el(:,3)~=0)
            uialert(parentFig,'No valid active phase data to export.','Export unavailable'); return;
        end
        figure(parentFig);
        [f2,p2]=uiputfile('phase_scheme.csv','Export phase scheme');
        if ~isgraphics(parentFig), return; end
        figure(parentFig);
        if isequal(f2,0), return; end
        try
            writetable(T,fullfile(p2,f2));
        catch err
            uialert(parentFig,err.message,'Export failed');
        end
    end

    function [cPorts,cAmp,cPh,cN,cUseRot,cOffsets] = cstPortWeights(scheme)
        [cPh,cRot,cRotDeg,cSwap,~,~,cOffsets] = portPhaseFor(scheme);
        % For dual-feed LP/CP, the declared physical feed-side topology and
        % scheme phases are already the complete CST port excitation
        % convention. Adding S.el(:,5) again double-counts the sequential
        % rotation (the exact error visible as 0/270 instead of
        % H=0/180, V=90/270 at broadside).
        %
        % Preserve the optional rotation-derived phase only for single-feed
        % export, where no dual-feed side/CP scheme exists.
        cUseRot = S.seqPhase && strcmp(scheme,'single');
        if cUseRot, cPh = cRot; cOffsets=cOffsets+cRotDeg; end
        cAmp = S.el(:,3);
        if any(~isfinite(cPh(:))) || any(~isfinite(cAmp)) || ~any(cAmp ~= 0)
            error('PAD:Excitation','Export needs finite phases and at least one nonzero weight.');
        end
        % Use nonnegative, relative wave amplitudes. Absorb signed weights
        % into phase; common normalization preserves every excitation ratio.
        cPh = appliedPortPhases(cPh,cRot,cAmp,cUseRot);
        cOffsets=mod(cOffsets+180*(cAmp<0),360);
        cAmp = abs(cAmp)/max(abs(cAmp));
        cN = numel(cAmp); cPorts = (1:cN)';
        if ~strcmp(scheme,'single')
            cH = 2*cPorts-1+double(cSwap);
            cV = 2*cPorts-double(cSwap);
            cPorts = [cH;cV];
            cAmp = [cAmp;cAmp]; % Equal H/V wave amplitudes, not watts.
        end
        cPh = cPh(:);
        [cPorts,cOrder] = sort(cPorts);
        cAmp = cAmp(cOrder); cPh = cPh(cOrder);
        cOffsets=cOffsets(:); cOffsets=cOffsets(cOrder);
    end

    function cText = buildCstMacro(cPorts,cMode,cAmp,cPh,cRef,cLabel,scheme,cUseRot,cTheta,cPhi,cFreq,cGeometryFreq,cUnit,cParam)
        if nargin < 14, cParam = []; end
        cCount=numel(cPorts);
        [macroAz,macroEl] = azEl(cTheta,cPhi);
            cLines = { ...
                ''' Generated by phasedArrayDesigner; CST beam excitation export.'; ...
                ''' Paste this complete file into a new CST VBA macro and run Main.'; ...
                ''' Existing geometry and solver settings are not modified.'; ...
                ''' Waveguide ports only. Verify port numbering and mode first.'; ...
                ''' Relative wave amplitudes, not power in watts.'; ...
                sprintf(''' Scheme %s; rotation compensation %d; steering theta %.12g phi %.12g deg.', ...
                    scheme,cUseRot,cTheta,cPhi); ...
                sprintf(''' Equivalent azimuth %.12g elevation %.12g deg; CST parameters remain Theta/Phi.', ...
                    macroAz,macroEl); ...
                sprintf(''' Phase reference %.12g GHz; geometry reference %.12g GHz; CST unit %s.',cFreq,cGeometryFreq,cUnit); ...
                ''' Dual feeds use equal amplitudes; CP adds +90 deg to V.'; ...
                ''' Opposite-side 180 deg correction is an uncalibrated preset.'; ...
                ''' The chosen H/V scheme must agree with the excitation behind an imported combined pattern.'; ...
                ''' Negative weights are represented by positive amplitude and +180 deg.'; ...
                'Option Explicit'; ...
                'Sub Main()'; ...
                '    On Error GoTo ExportFailed'; ...
                '    Dim existingLabels As Variant'; ...
                '    Dim existingLabel As Variant'; ...
                '    CombineResults.GetDefinedLabels existingLabels'; ...
                '    If IsArray(existingLabels) Then'; ...
                '        For Each existingLabel In existingLabels'; ...
                ['            If CStr(existingLabel) = "' cLabel '" Then']; ...
                '                MsgBox "This combination name already exists. Export again with a different name, such as a _v2 suffix."'; ...
                '                Exit Sub'; ...
                '            End If'; ...
                '        Next existingLabel'; ...
                '    End If'; ...
                '    With CombineResults'; ...
                '        .Reset'; ...
                '        .SetMonitorType "frequency"'; ...
                '        .SetOffsetType "phase"'; ...
                sprintf('        .SetReferenceFrequency %.15g',cRef); ...
                '        .EnableAutomaticLabeling False'; ...
                ['        .SetLabel "' cLabel '"']; ...
                '        .SetNone'};
            if ~isempty(cParam)
                % Quoted expressions must reach CST, not be evaluated by VBA.
                cThetaName=[cParam.prefix '_Theta_deg'];
                cPhiName=[cParam.prefix '_Phi_deg'];
                cXName=[cParam.prefix '_PhaseX_deg_per_mm'];
                cYName=[cParam.prefix '_PhaseY_deg_per_mm'];
                if isfield(cParam,'freezePhaseFrequency') && cParam.freezePhaseFrequency
                    % Beam-squint export: theta/phi can remain parametric, but the
                    % frequency inside the phase gradient MUST stay at the design
                    % reference. Making it a CST parameter here would silently turn
                    % a fixed-phase squint case back into frequency-retuned steering.
                    cInit={sprintf('    If MsgBox("Create/update %s steering-angle parameters? Phase frequency stays frozen for beam squint.", vbOKCancel) <> vbOK Then Exit Sub',cParam.prefix); ...
                        sprintf('    StoreParameter "%s", "%.15g"',cThetaName,cTheta); ...
                        sprintf('    StoreParameter "%s", "%.15g"',cPhiName,cPhi); ...
                        sprintf('    StoreParameter "%s", "-360*%.15g/299.792458*sin(%s*pi/180)*cos(%s*pi/180)"',cXName,cFreq,cThetaName,cPhiName); ...
                        sprintf('    StoreParameter "%s", "-360*%.15g/299.792458*sin(%s*pi/180)*sin(%s*pi/180)"',cYName,cFreq,cThetaName,cPhiName)};
                    cLines=[{sprintf(''' Parametric steering angles; BEAM SQUINT is ON, so phase frequency is frozen at %.12g GHz.',cFreq); ...
                        ''' Edit only *_Theta_deg and *_Phi_deg to change the commanded design-frequency steering direction.'; ...
                        ''' Change the solver/monitor operating frequency in CST without changing these phase expressions to observe squint.'; ...
                        ''' Positions and feed offsets are fixed at export; angles are degrees converted to radians in sin/cos.'; ...
                        ''' Check CST retains expressions after running this macro and updating a parameter.'}; cLines];
                else
                    cFreqName=[cParam.prefix '_Freq_GHz'];
                    cInit={sprintf('    If MsgBox("Create/update %s beam parameters? Combinations using this same prefix will share the updated settings.", vbOKCancel) <> vbOK Then Exit Sub',cParam.prefix); ...
                        sprintf('    StoreParameter "%s", "%.15g"',cFreqName,cFreq); ...
                        sprintf('    StoreParameter "%s", "%.15g"',cThetaName,cTheta); ...
                        sprintf('    StoreParameter "%s", "%.15g"',cPhiName,cPhi); ...
                        sprintf('    StoreParameter "%s", "-360*%s/299.792458*sin(%s*pi/180)*cos(%s*pi/180)"',cXName,cFreqName,cThetaName,cPhiName); ...
                        sprintf('    StoreParameter "%s", "-360*%s/299.792458*sin(%s*pi/180)*sin(%s*pi/180)"',cYName,cFreqName,cThetaName,cPhiName)};
                    cRefLine=find(startsWith(cLines,'        .SetReferenceFrequency '),1);
                    cLines{cRefLine}=sprintf('        .SetReferenceFrequency "%s*%.15g"',cFreqName,cRef/cFreq);
                    cLines=[{''' Parametric steering: edit the exported *_Freq_GHz, *_Theta_deg and *_Phi_deg parameters.'; ...
                        ''' Frequency changes recalculate the steering phases, matching Beam squint OFF in the MATLAB tool.'; ...
                        ''' Positions and feed offsets are fixed at export; angles are degrees converted to radians in sin/cos.'; ...
                        ''' Feed/rotation compensation is frozen at export, not re-estimated from a frequency-dependent pattern.'; ...
                        ''' Check CST retains expressions after running this macro and updating a parameter.'}; cLines];
                end
                cBefore=find(strcmp(cLines,'    With CombineResults'),1);
                cLines=[cLines(1:cBefore-1); cInit; cLines(cBefore:end)];
            end
            for cI = 1:cCount
                if isempty(cParam)
                    cLines{end+1,1} = sprintf( ...
                        '        .SetExcitationValues "port", "%d", %d, %.15g, %.15g', ...
                        cPorts(cI),cMode,cAmp(cI),cPh(cI)); %#ok<AGROW>
                else
                    % Preallocated: at most three terms (x, y, offset).
                    cTerms=cell(1,3); nTerm=0;
                    if cParam.xyMM(cI,1)~=0
                        nTerm=nTerm+1;
                        cTerms{nTerm}=sprintf('(%.15g)*%s',cParam.xyMM(cI,1),cXName);
                    end
                    if cParam.xyMM(cI,2)~=0
                        nTerm=nTerm+1;
                        cTerms{nTerm}=sprintf('(%.15g)*%s',cParam.xyMM(cI,2),cYName);
                    end
                    if cParam.offsetDeg(cI)~=0
                        nTerm=nTerm+1;
                        cTerms{nTerm}=sprintf('(%.15g)',cParam.offsetDeg(cI));
                    end
                    cTerms=cTerms(1:nTerm);
                    if isempty(cTerms), cExpr='0'; else, cExpr=strjoin(cTerms,' + '); end
                    cLines{end+1,1}=sprintf( ...
                        '        .SetExcitationValues "port", "%d", %d, %.15g, "%s"', ...
                        cPorts(cI),cMode,cAmp(cI),cExpr); %#ok<AGROW>
                end
            end
            cLines = [cLines; { ...
                '        .AddToExcitationList'; ...
                '    End With'; ...
                ['    MsgBox "Created ' cLabel '. Select it in Combine Excitation and check the port table. No solver was started."']; ...
                '    Exit Sub'; ...
                'ExportFailed:'; ...
                '    MsgBox "Could not create the combination: " & Err.Description'; ...
                'End Sub'}];
            cText = strjoin(cLines,sprintf('\r\n'));
    end

    function cParam = cstParametricData(cPorts,cPh,scheme,cPrefix)
        % Match physical positions to sorted CST port IDs, including H/V swaps.
        if strcmp(scheme,'single'), cElem=cPorts; else, cElem=ceil(cPorts/2); end
        cParam.xyMM=S.el(cElem,1:2)*lambdaMM();
        if any(~isfinite(cParam.xyMM(:)))
            error('PAD:CSTGeometry','Physical coordinates exceed numeric range. Check element positions and design frequency.');
        end
        % Read offsets directly, avoiding cancellation against wrapped steering phases.
        [cCheckPorts,~,cCheckPh,~,~,cParam.offsetDeg]=cstPortWeights(scheme);
        if ~isequal(cPorts,cCheckPorts) || ~isequal(cPh,cCheckPh)
            error('PAD:CSTSnapshot','Excitation changed while preparing the export. Export again.');
        end
        cParam.prefix=cPrefix;
        cParam.freezePhaseFrequency=~S.retunePhase;
    end

    function closeIfOpen(fid)
        try
            if ~isempty(fopen(fid)), fclose(fid); end
        catch
        end
    end

    function tf = validPortMode(v)
        tf=isnumeric(v)&&isscalar(v)&&isreal(v)&&isfinite(v)&&v>=1&&v==fix(v)&&v<=2147483647;
    end

    function exportCstExcitation(scheme,cParent,fromMap)
        % CST 2025 CombineResults API, from the supplied local Help pages.
        % A definition only: never starts the solver or combines old results.
        % All export messages belong to the window that launched the export.
        %
        % The ONE implementation behind both routes: the Phase scheme
        % map's "Export CST macro" (its scheme, its window) and EXPORT >
        % CST macro in the main window (exportCstMacro). Same design, same
        % answers in the form -> the same file, byte for byte.
        %
        % The answers are collected in one form window (openCstMacroDialog)
        % instead of the old chain uiconfirm -> inputdlg -> file dialog.
        % Both of those waited for the user INSIDE this call, which a
        % batch session refuses outright, and inputdlg took the frequency
        % unit as free text that could only be checked after the fact.
        if ~isgraphics(cParent), return; end
        figure(cParent);
        if ~patternReady(false)
            uialert(cParent,'Load the required imported pattern before exporting.', ...
                'Imported pattern required'); return;
        end
        if isempty(S.el)
            uialert(cParent,'Place at least one element first.','No elements'); return;
        end
        % Refused before the form opens: nothing the user types there can
        % repair an excitation that cannot be exported.
        try
            [cPorts,~,cPh]=cstPortWeights(scheme);
            cstParametricData(cPorts,cPh,scheme,'Beam');
        catch exportErr
            uialert(cParent,exportErr.message,'Invalid excitation'); return;
        end
        openCstMacroDialog(scheme,cParent,fromMap);
    end

    function refreshImportLabel()
        %REFRESHIMPORTLABEL  Make the import buttons report S.impFF.
        %   Called from the import itself AND from config load. It used to
        %   live inline in loadImportedFF, so a design reopened from a
        %   config restored its pattern correctly but went on naming
        %   whatever file had been imported by hand in this session --
        %   the data was right and the label lied about which data it was.
        % The empty-state caption is the one the button is built with:
        % a design loaded without a pattern used to rename it to a second
        % wording of the same action.
        if isempty(S.impFF)
            bImportFF.Text = 'Load far-field pattern (.txt)...';
        else
            % Only the name: grid, peak and frequency are spelled out in
            % the readout underneath, where they do not get clipped.
            bImportFF.Text = sprintf('Loaded: %s', S.impFFName);
        end
        refreshImportInfo();
    end

    function loadImportedFF()
        [f,p] = uigetfile({'*.txt;*.dat;*.csv','CST farfield ASCII export'}, ...
            'Load far-field pattern (CST ASCII export)');
        figure(fig);   % native file dialogs can leave the uifigure behind
                       % other windows on macOS, especially with real
                       % processing time between the dialog closing and
                       % the next screen update -- force it back to front.
        if isequal(f,0), return; end
        try
            ffS = loadCstFarfieldASCII(fullfile(p,f));
        catch err
            uialert(fig, err.message, 'Far-field import failed');
            setStatus(['Far-field import refused: ' f],'bad');
            return;
        end
        S.impFF = ffS; S.impFFName = f;
        [S.impFFGHz,S.impFFFreqFrom] = cstFileFrequency(f,fullfile(p,f));
        S.impNeedsPattern = false;   % supplied now
        % Reports the file's own peak gain and where it points. Peak far
        % from theta=0 means the CST model is oriented wrong for this
        % app (elements face +z here) -- re-orient in CST and re-export.
        % The efficiency tag is built inside refreshImportLabel now, so
        % the import path and the config-restore path cannot word it
        % differently.
        refreshImportLabel();
        % Every concern about the loaded file is COLLECTED and shown in
        % ONE dialog. uialert is non-blocking -- it returns immediately
        % rather than waiting for the user -- so firing two or three of
        % them back to back means the later calls land while the first
        % is still up, and they overlay or replace each other. The user
        % then sees one warning and never learns about the rest. This is
        % reachable with a real file: a Ludwig-3 export that ALSO points
        % away from broadside raises two at once.
        notes = {};
        % Checked before (and instead of) the basis note below: this file
        % carries no polarisation data at all, so telling the user to
        % "re-export in the Theta/Phi basis" would be the wrong fix --
        % there is no basis to change, the decomposition is simply absent.
        if isfield(ffS,'dupWarn') && ~isempty(ffS.dupWarn)
            notes(end+1,:) = {'Duplicate directions disagree', ffS.dupWarn};
        end
        if ffS.noComponents
            notes(end+1,:) = {'Magnitude-only export (Abs) loaded', ...
                ['This is a CST "Abs" export: it carries the total magnitude ' ...
                'only, with no split into polarisation components (both component ' ...
                'columns repeat the total and every phase entry is zero).' newline newline ...
                'The pattern has been loaded from the total column, so the beam ' ...
                'shape, gain levels, directivity, HPBW and sidelobe level are ' ...
                'correct and will match CST.' newline newline ...
                'ARRAY results are correct ONLY while every element shares the ' ...
                'same orientation. The file carries no phase and no polarisation ' ...
                'split, so this app has had to put the whole magnitude into ' ...
                'E_theta with E_phi zero -- an invented polarisation. That is ' ...
                'harmless when the elements are identical and identically ' ...
                'oriented, because the same invented component then multiplies ' ...
                'the array factor everywhere. It is NOT harmless once elements ' ...
                'carry different rotations: the coherent sum then combines ' ...
                'components that do not exist as written, so sequential rotation ' ...
                'and any mixed-rotation array will be wrong.' newline newline ...
                'Also unavailable for the same reason: Axial Ratio and the ' ...
                'RHCP/LHCP dropdown, which are disabled for this pattern. ' ...
                'Re-export with a component pair (Theta/Phi) selected if you ' ...
                'need polarisation or a rotated array.']}; 
        end
        if ~ffS.noComponents && ~ffS.compThetaPhi && ~isempty(ffS.compNames{1})
            notes(end+1,:) = {'Component basis is not Theta/Phi', ...
                sprintf(['This export stores its field as %s / %s components ' ...
                'rather than the Theta / Phi spherical pair.\n\nTotal power is the same ' ...
                'in any orthogonal basis, so the GAIN pattern, directivity and beam ' ...
                'shape are all correct. But the polarization-sensitive readouts -- ' ...
                'Axial Ratio, and the RHCP/LHCP polarization dropdown -- assume ' ...
                'Theta/Phi and are NOT valid for this file.\n\nIf you need those, ' ...
                're-export from CST with the far-field component basis set to ' ...
                'Theta/Phi (spherical) instead.'], ...
                ffS.compNames{1}, ffS.compNames{2})}; 
        end
        if ~ffS.phiFull
            notes(end+1,:) = {'Partial azimuth coverage', ...
                ['This export does not span a full 360 deg in phi, so angles ' ...
                'outside the exported azimuth range are clamped to the nearest ' ...
                'exported value rather than interpolated, and no efficiency (gain ' ...
                'vs directivity) can be determined from it.' newline newline ...
                'Re-export from CST over the full phi range for correct array results.']}; 
        end
        if ffS.peak(2) > 30
            notes(end+1,:) = {'Pattern not pointing at broadside', ...
                sprintf(['This pattern peaks at theta = %g deg (phi = %g deg), ' ...
                'not near theta = 0.\n\nThis app places every element in the z = 0 plane ' ...
                'facing +z, so the element boresight must be theta = 0 (CST''s +z axis). ' ...
                'Re-orient the model in CST so it radiates along +z, then re-export.\n\n' ...
                'The pattern has been loaded exactly as exported -- nothing was rotated.'], ...
                ffS.peak(2), ffS.peak(3))};
        end
        if importFreqMismatch()
            notes(end+1,:) = {'Pattern frequency differs', importFreqWarning(true)};
        end
        if size(notes,1) == 1
            uialert(fig, notes{1,2}, notes{1,1}, 'Icon','warning');
        elseif size(notes,1) > 1
            parts = cell(1,size(notes,1));
            for i = 1:size(notes,1)
                parts{i} = sprintf('%d) %s\n\n%s', i, notes{i,1}, notes{i,2});
            end
            uialert(fig, strjoin(parts, sprintf('\n\n%s\n\n', repmat('-',1,40))), ...
                sprintf('Imported pattern: %d things to check', size(notes,1)), ...
                'Icon','warning');
        end
        % Selecting the imported pattern is the whole point of loading
        % it. Leaving Element factor alone meant the file was parsed,
        % reported on, and then ignored until you also changed a
        % dropdown -- which is indistinguishable from a failed import,
        % and cost real time before it was understood.
        S.efType = 'Imported (CST far-field)';
        ddEF.Value = S.efType;
        refreshPolAvailability();
        % AFTER the element type and the availability gate, not before.
        % The table prints effectivePhaseDeg(), which folds in
        % seqRotSign() -- and seqRotSign() reads S.efType. Refreshing
        % while efType still named the OLD element computed the
        % compensation for that element: importing an LHCP pattern over
        % a built-in RHCP patch left +rot phases in the table while the
        % pattern was computed with -rot.
        refreshTable();
        refreshLayout();
        maybeCompute();
        % After the compute, which writes its own line: what was read.
        if importFreqMismatch()
            importStatus(sprintf('Imported %s: %s',f,importFreqWarning(false)),'warn');
        else
            importStatus(sprintf('Imported %s: %s',f,importReadText()),'good');
        end
    end
    function assignSeqBlock(w,v)
        v=round(v);
        if w=='m'
            S.seqBlockM=v; spSeqBlkM.Value=v;
        else
            S.seqBlockN=v; spSeqBlkN.Value=v;
        end
        % A saved design setting that recomputes nothing, so it never
        % reaches maybeCompute's checkpoint: record the step here.
        undoCheckpoint();
    end
    function assignSeq(v), S.seqPhase = v; refreshTable(); maybeCompute(); end
    function assignTaper(v)
        S.taper = v;
        applyTaper(); syncSLLControl(); refreshAll();
    end

    function w = calculationWeights()
        w=S.el(:,3);
        if isempty(w), return; end
        wm=max(abs(w));
        if wm>0, w=w/wm; end
    end
    function invalidateRadiation(reason)
        if isgraphics(overlayWin)
            set(findall(overlayWin,'Tag','overlayNote'),'Text', ...
                'Pattern is out of date; recompute to update the comparison.');
        end
        S.radiationValid=false; S.cutX=[]; S.cutDb=[]; S.cutRaw=[];
        S.pattern2DBasis=[];
        S.DpkTot=NaN; S.glTheta=NaN; S.glPhi=NaN; S.glRelDb=NaN;
        clearCutMetrics();
        blankPlots(reason, 'No current radiation result');
        if isgraphics(pattern2DWin), drawPattern2D(); end
        showDetailsMessage(reason);
        setStatus(reason,'warn');
    end
    function refreshPhaseMaps()
        for rm=1:numel(S.portMapRefresh)
            if isempty(S.portMapRefresh{rm}), continue; end
            try
                S.portMapRefresh{rm}();
            catch mapErr
                warning('PAD:MapRefresh','Phase map could not refresh: %s',mapErr.message);
            end
        end
    end
    function computePattern()
        nElem = size(S.el,1);
        if nElem >= BUSY_MIN_ELEMENTS || lastComputeSec >= 0.5
            busyTok = beginBusy(sprintf( ...
                'Computing the pattern of %d elements…', nElem)); %#ok<NASGU>
        end
        tCompute = tic;
        S.radiationValid=false; S.cutRaw=[]; clearCutMetrics();
        S.pattern2DBasis=[];
        try
            if isempty(S.el)||~any(S.el(:,3)~=0)
                invalidateRadiation('No radiation: no active elements.'); return;
            end
            computePatternCore();
            S.radiationValid=~isempty(S.cutX) && ~S.efFallback && ~S.impNeedsPattern;
            if isgraphics(pattern2DWin), drawPattern2D(); end
            redrawCstOverlay();
            % An empty cut means the core refused (patternReady: the
            % imported pattern is missing) and has already put the reason
            % in the status bar; "Pattern computed" would overwrite it.
            if isempty(S.cutX), return; end
            lastComputeSec = toc(tCompute);
            if S.efFallback
                setStatus(sprintf(['Pattern computed in %.2f s with a ' ...
                    'SUBSTITUTED element: the Custom formula failed on ' ...
                    'part of the grid.'], lastComputeSec), 'warn');
            else
                setStatus(sprintf('Pattern computed in %.2f s · %d elements', ...
                    lastComputeSec, nElem));
            end
        catch patternErr
            invalidateRadiation(['Calculation unavailable: ' patternErr.message]);
            setStatus(['Pattern calculation failed: ' patternErr.message],'bad');
            uialert(fig,patternErr.message,'Pattern calculation refused');
        end
    end
    function showScanLoss()
        busyTok = beginBusy('Computing scan loss…'); %#ok<NASGU>
        try
            showScanLossCore();
        catch scanErr
            setStatus(['Scan loss failed: ' scanErr.message],'bad');
            % Not on a window closed mid-sweep: an alert on a deleted
            % figure is a second error, and this one would escape.
            if isgraphics(fig)
                uialert(fig,scanErr.message,'No valid scan-loss curve');
            end
        end
    end

    function showPattern2D(mode)
        % Show the current vector-field cut in the main result workspace.
        % These controls change observation directions, never steering.
        % The ribbon button passes an empty mode: Elevation is the first
        % view, but an already-open 2D pane keeps the user's chosen cut.
        if isempty(mode)
            if ~isempty(pattern2DWin) && isgraphics(pattern2DWin)
                currentMode = findall(pattern2DWin,'Tag','ddPattern2DMode');
                if ~isempty(currentMode)
                    mode = currentMode.Value;
                end
            end
            if isempty(mode), mode = 'Elevation pattern'; end
        end
        if ~S.radiationValid || isempty(S.pattern2DBasis)
            computePattern();
        end
        if ~S.radiationValid || isempty(S.pattern2DBasis), return; end
        newWindow = false;
        if isempty(pattern2DWin) || ~isgraphics(pattern2DWin)
            newWindow = true;
            pattern2DWin = resultPane('2D Pattern','pattern2DWin');
            layout2D = uigridlayout(pattern2DWin,[3 1]);
            layout2D.RowHeight = {76,'1x',52};
            layout2D.Padding = [12 8 12 8];
            controls2D = uigridlayout(layout2D,[2 7]);
            controls2D.Layout.Row = 1;
            controls2D.ColumnWidth = {85,175,105,100,62,175,'1x'};
            controls2D.RowHeight = {32,32};
            controls2D.Padding = [0 0 0 0];
            c2 = uilabel(controls2D,'Text','2D view');
            c2.Layout.Row = 1; c2.Layout.Column = 1;
            c2 = uidropdown(controls2D, ...
                'Items',{'Azimuth pattern','Elevation pattern','U pattern'}, ...
                'Value',mode,'Tag','ddPattern2DMode', ...
                'ValueChangedFcn',@(s,e)configurePattern2D());
            c2.Layout.Row = 1; c2.Layout.Column = 2;
            c2 = uilabel(controls2D,'Text','Fixed angle','Tag','lblPattern2DFixed');
            c2.Layout.Row = 1; c2.Layout.Column = 3;
            c2 = uispinner(controls2D,'Limits',[0 360],'Value',0,'Step',5, ...
                'Tag','spPattern2DFixed', ...
                'ValueChangedFcn',@(s,e)drawPattern2D());
            c2.Layout.Row = 1; c2.Layout.Column = 4;
            c2 = uilabel(controls2D,'Text','Pattern');
            c2.Layout.Row = 1; c2.Layout.Column = 5;
            c2 = uidropdown(controls2D, ...
                'Items',{'Total (EF x AF)','Array factor only','Element factor only'}, ...
                'Value','Total (EF x AF)','Tag','ddPattern2DSource', ...
                'Tooltip','Show the total field, isotropic-element array factor, or first element pattern.', ...
                'ValueChangedFcn',@(s,e)drawPattern2D());
            c2.Layout.Row = 1; c2.Layout.Column = 6;
            c2 = uicheckbox(controls2D,'Text','Absolute dBi','Value',true, ...
                'Tag','cbPattern2DAbs', ...
                'Tooltip','For dB and U views, use each selected pattern’s full-sphere gain/directivity basis.', ...
                'ValueChangedFcn',@(s,e)drawPattern2D());
            c2.Layout.Row = 1; c2.Layout.Column = 7;
            c2 = uilabel(controls2D,'Text','Polarization');
            c2.Layout.Row = 2; c2.Layout.Column = 1;
            initialPol = ddPol.Value;
            c2 = uidropdown(controls2D, ...
                'Items',{'Total (any pol)','RHCP + LHCP', ...
                         'RHCP component','LHCP component'}, ...
                'Value',initialPol,'Tag','ddPattern2DPol', ...
                'Tooltip','Plot total power, one circular-polarization component, or both components.', ...
                'ValueChangedFcn',@(s,e)drawPattern2D());
            c2.Layout.Row = 2; c2.Layout.Column = 2;
            c2 = uilabel(controls2D,'Text','Polar radius');
            c2.Layout.Row = 2; c2.Layout.Column = 3;
            c2 = uidropdown(controls2D, ...
                'Items',{'Linear power','dB (dynamic range)'}, ...
                'Value','dB (dynamic range)','Tag','ddPattern2DScale', ...
                'Tooltip','Linear power reveals lobe shape; dB expands weak radiation and enables absolute dBi ticks.', ...
                'ValueChangedFcn',@(s,e)drawPattern2D());
            c2.Layout.Row = 2; c2.Layout.Column = 4;
            plotHost = uigridlayout(layout2D,[1 1]);
            plotHost.Layout.Row = 2;
            plotHost.Padding = [2 2 2 2];
            plotHost.Tag = 'pattern2DPlotHost';
            note2D = uilabel(layout2D,'Text','','WordWrap','on', ...
                'Tag','lblPattern2DNote');
            note2D.Layout.Row = 3;
            mutedLbls(end+1) = note2D;
        else
            selectTask('2D Pattern');
        end
        dd = findall(pattern2DWin,'Tag','ddPattern2DMode');
        if ~strcmp(dd.Value,mode) || newWindow
            dd.Value = mode;
            configurePattern2D();
        else
            drawPattern2D();
        end
    end

    function configurePattern2D()
        if isempty(pattern2DWin) || ~isgraphics(pattern2DWin), return; end
        dd = findall(pattern2DWin,'Tag','ddPattern2DMode');
        sp = findall(pattern2DWin,'Tag','spPattern2DFixed');
        lbl2D = findall(pattern2DWin,'Tag','lblPattern2DFixed');
        if strcmp(dd.Value,'Azimuth pattern')
            sp.Value = 0;
            if useAzEl()
                lbl2D.Text = 'Elevation (°)';
                if S.fullSphere, sp.Limits = [-90 90];
                else, sp.Limits = [0 90]; end
            else
                lbl2D.Text = 'Theta (°)';
                if S.fullSphere, sp.Limits = [0 180];
                else, sp.Limits = [0 90]; end
                sp.Value = 90;
            end
        else
            if useAzEl(), lbl2D.Text = 'Azimuth (°)';
            else, lbl2D.Text = 'Phi (°)'; end
            sp.Value = 0;
            sp.Limits = [0 360];
            sp.Value = mod(cutPhiVal(),360);
        end
        drawPattern2D();
    end

    function syncPattern2DFixedAngle()
        if isempty(pattern2DWin) || ~isgraphics(pattern2DWin), return; end
        mode = findall(pattern2DWin,'Tag','ddPattern2DMode');
        fixedAngle = findall(pattern2DWin,'Tag','spPattern2DFixed');
        if isempty(mode) || isempty(fixedAngle), return; end
        follows = S.cutPhiFollow && ...
            ~strcmp(mode.Value,'Azimuth pattern');
        if follows, fixedAngle.Value = mod(cutPhiVal(),360); end
        fixedAngle.Enable = ternStr(~follows,'on','off');
    end

    function drawPattern2D()
        if isempty(pattern2DWin) || ~isgraphics(pattern2DWin), return; end
        syncPattern2DFixedAngle();
        if ~S.radiationValid || isempty(S.pattern2DBasis)
            host = findall(pattern2DWin,'Tag','pattern2DPlotHost');
            delete(host.Children);
            note = findall(pattern2DWin,'Tag','lblPattern2DNote');
            note.Text = 'Pattern is out of date. Click 2D pattern to recompute.';
            return;
        end
        dd = findall(pattern2DWin,'Tag','ddPattern2DMode');
        sp = findall(pattern2DWin,'Tag','spPattern2DFixed');
        cb = findall(pattern2DWin,'Tag','cbPattern2DAbs');
        pol2D = findall(pattern2DWin,'Tag','ddPattern2DPol');
        radius2D = findall(pattern2DWin,'Tag','ddPattern2DScale');
        source2D = findall(pattern2DWin,'Tag','ddPattern2DSource');
        sourceChoice = source2D.Value;
        host = findall(pattern2DWin,'Tag','pattern2DPlotHost');
        delete(host.Children);
        isU = strcmp(dd.Value,'U pattern');
        isDb = isU || strcmp(radius2D.Value,'dB (dynamic range)');
        if isU, radius2D.Enable = 'off';
        else, radius2D.Enable = 'on'; end
        if isDb, cb.Visible = 'on';
        else, cb.Visible = 'off'; end
        if strcmp(sourceChoice,'Array factor only')
            pol2D.Enable = 'off';
        else
            pol2D.Enable = 'on';
        end
        fixed = sp.Value;
        switch dd.Value
            case 'Azimuth pattern'
                angle = 0:0.5:360;
                if useAzEl(), thetaFixed = 90-fixed;
                else, thetaFixed = fixed; end
                th = thetaFixed*ones(size(angle));
                ph = angle;
                if useAzEl()
                    head = sprintf('Azimuth cut (elevation = %.1f°)',fixed);
                else
                    head = sprintf('Phi cut (theta = %.1f°)',fixed);
                end
            case 'Elevation pattern'
                if S.fullSphere, angle = -180:0.5:180;
                else, angle = -90:0.25:90; end
                th = abs(angle);
                ph = fixed + 180*(angle<0);
                if useAzEl()
                    head = sprintf('Elevation cut (azimuth = %.1f°)',fixed);
                else
                    head = sprintf('Theta cut (phi = %.1f°)',fixed);
                end
            otherwise
                angle = -1:0.002:1;
                th = asind(abs(angle));
                ph = fixed + 180*(angle<0);
                if useAzEl()
                    head = sprintf('U cut (azimuth = %.1f°)',fixed);
                else
                    head = sprintf('U cut (phi = %.1f°)',fixed);
                end
        end
        bothCP = ~strcmp(sourceChoice,'Array factor only') && ...
            strcmp(pol2D.Value,'RHCP + LHCP');
        if bothCP
            choices = {'RHCP component','LHCP component'};
        elseif strcmp(sourceChoice,'Array factor only')
            choices = {'Total (any pol)'};
        else
            choices = {pol2D.Value};
        end
        levels = cell(size(choices));
        try
            for j2 = 1:numel(choices)
                levels{j2} = pattern2DLevels(th,mod(ph,360), ...
                    cb.Value && isDb,choices{j2},sourceChoice);
            end
        catch cutErr
            note = findall(pattern2DWin,'Tag','lblPattern2DNote');
            note.Text = ['2D pattern unavailable: ' cutErr.message];
            setStatus(['2D pattern unavailable: ' cutErr.message],'warn');
            return;
        end
        DR = spDR.Value;
        top = 5*ceil(max(cellfun(@max,levels))/5);
        bottom = top-DR;
        P = pal();
        if isU
            ax = uiaxes(host,'Tag','pattern2DUAxes');
            hold(ax,'on');
            for j2 = 1:numel(levels)
                [lineColor,lineTag] = pattern2DTraceStyle(choices{j2},P);
                plot(ax,angle,levels{j2},'Color',lineColor,'LineWidth',2, ...
                    'Tag',lineTag);
            end
            grid(ax,'on'); xlim(ax,[-1 1]); ylim(ax,[bottom top+1]);
            xlabel(ax,'u = signed sin(θ) along the cut plane');
            ylabel(ax,pattern2DLevelLabel(cb.Value,sourceChoice));
            title(ax,head);
        else
            ax = polaraxes(host);
            ax.Tag = 'pattern2DPolarAxes';
            hold(ax,'on');
            for j2 = 1:numel(levels)
                [lineColor,lineTag] = pattern2DTraceStyle(choices{j2},P);
                if isDb
                    rho = max(levels{j2}-bottom,0);
                else
                    % Keep one FULL-SPHERE reference for every cut and
                    % both CP components. Renormalizing each cut or each
                    % polarization separately would disguise their loss.
                    switch sourceChoice
                        case 'Array factor only'
                            peakDb = S.pattern2DBasis.peakDbAF;
                        case 'Element factor only'
                            peakDb = S.pattern2DBasis.peakDbEF;
                        otherwise
                            peakDb = S.pattern2DBasis.peakDb;
                    end
                    rho = min(10.^((levels{j2}-peakDb)/10),1);
                end
                polarplot(ax,deg2rad(angle),rho, ...
                    'Color',lineColor,'LineWidth',2,'Tag',lineTag);
            end
            if isDb
                rlim(ax,[0 DR]);
                ticks = linspace(0,DR,5);
                rticks(ax,ticks);
                ax.RTickLabel = arrayfun(@(v)sprintf('%g',v+bottom), ...
                    ticks,'UniformOutput',false);
            else
                rlim(ax,[0 1]);
                rticks(ax,0:0.25:1);
            end
            if strcmp(dd.Value,'Elevation pattern')
                % In the signed theta cut, 0 is array broadside (+z).
                % Put broadside at the top, with positive cut angles to
                % its right and negative ones to its left. Leave the
                % azimuth plot at the conventional +x/right origin.
                ax.ThetaZeroLocation = 'top';
                ax.ThetaDir = 'clockwise';
            else
                ax.ThetaZeroLocation = 'right';
                ax.ThetaDir = 'counterclockwise';
            end
            ax.ThetaTick = 0:30:330;
            labels = mod(ax.ThetaTick+180,360)-180;
            labels(ax.ThetaTick==180) = 180;
            ax.ThetaTickLabel = arrayfun(@num2str,labels,'UniformOutput',false);
            title(ax,head);
        end
        if bothCP
            legend(ax,{'RHCP','LHCP'},'Location','best');
        end
        ax.Layout.Row = 1; ax.Layout.Column = 1;
        note = findall(pattern2DWin,'Tag','lblPattern2DNote');
        if isU
            axisNote = 'U is signed direction cosine; -1 and +1 are the two horizon directions.';
            levelNote = pattern2DLevelLabel(cb.Value,sourceChoice);
        elseif isDb
            axisNote = 'Polar radius is the pattern level in dB; inward means lower radiation.';
            levelNote = pattern2DLevelLabel(cb.Value,sourceChoice);
        else
            axisNote = 'Polar radius is normalized linear power (0–1), referenced to this factor’s full-sphere peak.';
            levelNote = '';
        end
        if strcmp(sourceChoice,'Array factor only')
            polNote = 'Array factor assumes isotropic elements; polarization does not apply.';
        elseif strcmp(sourceChoice,'Element factor only')
            polNote = sprintf('First element pattern (rotation %.1f°); %s.', ...
                S.el(1,5),pol2D.Value);
        elseif bothCP
            polNote = 'Blue: RHCP; orange: LHCP.';
        elseif strcmp(pol2D.Value,'Total (any pol)')
            polNote = 'Total polarization.';
        else
            polNote = [pol2D.Value '.'];
        end
        note.Text = strtrim(sprintf('%s %s %s',levelNote,axisNote,polNote));
        note.FontColor = P.muted;
    end

    function [colour,tag] = pattern2DTraceStyle(choice,P)
        switch choice
            case 'RHCP component'
                colour = P.traceTotal; tag = 'pattern2DRHCP';
            case 'LHCP component'
                colour = P.overlayTrace; tag = 'pattern2DLHCP';
            otherwise
                colour = P.traceTotal; tag = 'pattern2DTotal';
        end
    end

    function label = pattern2DLevelLabel(isAbsolute,sourceChoice)
        if ~isAbsolute
            label = 'Relative pattern level (dB).';
        elseif strcmp(sourceChoice,'Array factor only')
            label = 'Array factor directivity (dBi).';
        elseif S.impUnitCell && strcmp(S.efType,'Imported (CST far-field)')
            label = 'Array realized gain (dBi).';
        elseif strcmp(S.efType,'Imported (CST far-field)')
            label = 'Realized gain (dBi).';
        else
            label = 'Directivity (dBi).';
        end
    end

    function level = pattern2DLevels(th,ph,isAbsolute,polChoice,sourceChoice)
        % Evaluate the SAME steered vector sum used by computePatternCore,
        % at fine 1D samples. Its full-sphere normalization was cached by
        % that computation, so changing the observation cut cannot
        % silently renormalize a weak or off-beam plane to 0 dB.
        b = S.pattern2DBasis;
        fallbackBefore = S.efFallback;
        fallbackWhichBefore = S.efFallbackWhich;
        w = calculationWeights();
        us = sind(S.theta_s)*cosd(S.phi_s);
        vs = sind(S.theta_s)*sind(S.phi_s);
        u = sind(th).*cosd(ph);
        v = sind(th).*sind(ph);
        fr = freqRatio();
        rot = S.el(:,5);
        switch sourceChoice
            case 'Array factor only'
                arrayField = zeros(size(th));
                for k2 = 1:size(S.el,1)
                    xn = S.el(k2,1); yn = S.el(k2,2);
                    feed = -S.k*phaseFreqRatio()*(xn*us+yn*vs) + ...
                        deg2rad(S.el(k2,4));
                    if S.seqPhase, feed = feed+b.seqSign*deg2rad(rot(k2)); end
                    arrayField = arrayField + w(k2)*exp(1j*feed)* ...
                        exp(1j*S.k*fr*(xn*u+yn*v));
                end
                field = abs(arrayField);
            case 'Element factor only'
                [ethTotal,ephTotal] = elementFactor(th,ph,rot(1));
            otherwise
                [uniqueRot,~,rotIndex] = unique(rot);
                ethByRot = cell(numel(uniqueRot),1);
                ephByRot = cell(numel(uniqueRot),1);
                for k2 = 1:numel(uniqueRot)
                    [ethByRot{k2},ephByRot{k2}] = ...
                        elementFactor(th,ph,uniqueRot(k2));
                end
                ethTotal = zeros(size(th)); ephTotal = zeros(size(th));
                for k2 = 1:size(S.el,1)
                    xn = S.el(k2,1); yn = S.el(k2,2);
                    feed = -S.k*phaseFreqRatio()*(xn*us+yn*vs) + ...
                        deg2rad(S.el(k2,4));
                    if S.seqPhase, feed = feed+b.seqSign*deg2rad(rot(k2)); end
                    wave = w(k2)*exp(1j*feed)*exp(1j*S.k*fr*(xn*u+yn*v));
                    ethTotal = ethTotal + wave.*ethByRot{rotIndex(k2)};
                    ephTotal = ephTotal + wave.*ephByRot{rotIndex(k2)};
                end
        end
        if S.efFallback && ~fallbackBefore
            S.efFallback = fallbackBefore;
            S.efFallbackWhich = fallbackWhichBefore;
            error('PAD:2DCutFormula','The custom element formula failed at a 2D cut angle.'); %#ok<CTPCT>
        end
        if ~strcmp(sourceChoice,'Array factor only')
            switch polChoice
                case 'RHCP component'
                    field = abs(ethTotal+1j*ephTotal)/sqrt(2);
                case 'LHCP component'
                    field = abs(ethTotal-1j*ephTotal)/sqrt(2);
                otherwise
                    field = hypot(abs(ethTotal),abs(ephTotal));
            end
        end
        switch sourceChoice
            case 'Array factor only'
                reference = b.ref; floorValue = b.floorAF;
                Prad = b.PradAF; efficiency = 1;
            case 'Element factor only'
                reference = 1; floorValue = b.floorEF;
                Prad = b.PradEF; efficiency = impEffLin();
            otherwise
                reference = b.ref; floorValue = b.floor;
                Prad = b.Prad; efficiency = impEffLin();
        end
        level = 20*log10(max(abs(field)/reference,floorValue));
        if isAbsolute
            offset = safeOff(efficiency,Prad);
            uc = unitCellFigures();
            if ~isempty(uc)
                switch sourceChoice
                    case 'Element factor only'
                        offset = uc.rg;
                    case 'Total (EF x AF)'
                        pIn = sum(w.^2);
                        if pIn > 0
                            offset = uc.rg + 10*log10(b.ref^2/pIn);
                        else
                            offset = uc.rg;
                        end
                end
            end
            if isfinite(offset), level = level+offset; end
        end
    end

    % ------------------------------------------------ analyses: shared
    function reportAnalysis(msg)
        %REPORTANALYSIS  Status line for an analysis that ended (done or
        %   cancelled). With Auto off, an edit clears the main plots and
        %   the status line says why ("Settings changed: click Compute
        %   pattern."). An analysis message written over that line hid
        %   the only explanation of the empty plots, so while they are out
        %   of date the message carries the warning with it.
        if S.radiationValid
            setStatus(msg);
        else
            setStatus([msg ' · main plots out of date: click Compute ' ...
                'pattern.'],'warn');
        end
    end

    % ------------------------------------------------------- band sweep
    % Beam pointing, gain toward the commanded direction and beamwidth
    % across a frequency band, for the two ways of steering: phase
    % shifters set at one frequency (the beam squints away from it,
    % sin(theta) = sin(theta_s)*f_set/f) and true-time delay (no squint).
    % It answers the wideband receive-array question the single
    % Operating-frequency view cannot: how much of the band a scanned
    % beam can serve before its pointing loss is too large.
    %
    % Built on the app's own squint machinery -- freqRatio(f),
    % squintSinTheta and the exact radiated-power integrals, each given
    % the swept frequency explicitly -- rather than on a copy of it. The
    % main window's Operating frequency and Beam squint settings are
    % neither read nor changed.
    function showCoverage()
        if isempty(S.coverageWin) || ~isgraphics(S.coverageWin)
            S.coverageWin = resultPane('Coverage','coverageWin');
            cg = uigridlayout(S.coverageWin,[3 1]); cg.RowHeight={65,'1x',65};
            controls = uigridlayout(cg,[2 4]);
            controls.ColumnWidth={'1x','1x','1x',125};
            uilabel(controls,'Text','Maximum θ (°)');
            uilabel(controls,'Text','θ step (°)');
            uilabel(controls,'Text','Azimuth step (°)');
            uilabel(controls,'Text','');
            uispinner(controls,'Limits',[1 88],'Value',60,'Tag','coverageMax', ...
                'Tooltip','Maximum steering angle from broadside, in degrees');
            uispinner(controls,'Limits',[1 45],'Value',10,'Tag','coverageThetaStep', ...
                'Tooltip','Steering-angle sample interval, in degrees');
            uispinner(controls,'Limits',[5 90],'Value',30,'Tag','coveragePhiStep', ...
                'Tooltip','Azimuth sample interval around the full circle, in degrees');
            uibutton(controls,'Text','Run coverage','Tag','coverageRun', ...
                'Tooltip','Compute scan loss and grating-lobe flags; unchanged designs use the cache', ...
                'ButtonPushedFcn',@(s,e)runCoverage());
            uiaxes(cg,'Tag','coverageAxes');
            uilabel(cg,'Text','Choose the scan cone and sampling steps, then Run coverage.', ...
                'Tag','coverageNote','WordWrap','on');
            S.coverageWin.UserData = struct('result',[],'cacheHit',false);
        else
            selectTask('Coverage');
        end
    end

    function runCoverage()
        win = S.coverageWin;
        if ~isgraphics(win) || ~patternReady(true), return; end
        tm = findall(win,'Tag','coverageMax'); dt = findall(win,'Tag','coverageThetaStep');
        dp = findall(win,'Tag','coveragePhiStep');
        settings = [tm.Value dt.Value dp.Value];
        design = designSnapshot();
        if ~isempty(S.coverageCache) && isequaln(S.coverageCache.design,design) && ...
                isequal(S.coverageCache.settings,settings)
            win.UserData = struct('result',S.coverageCache.result,'cacheHit',true);
            drawCoverage(); reportAnalysis('Scan coverage: reused cached result'); return;
        end
        th = unique([0:dt.Value:tm.Value tm.Value]); ph = 0:dp.Value:360;
        ph = ph(ph<360); [tt,pp] = meshgrid(th,ph);
        samples = [0 0; tt(:) pp(:)];
        b = findall(win,'Tag','coverageRun'); b.Enable = 'off';
        coverageDone = onCleanup(@()coverageFinished(b));
        tok = beginBusy('Computing scan coverage…'); %#ok<NASGU>
        try
            R = showScanLossCore(samples);
            if isempty(R) || ~isgraphics(win), return; end
            if ~isequaln(design,designSnapshot())
                reportAnalysis('Scan coverage cancelled: design changed during the sweep'); return;
            end
            R.theta = th; R.phi = ph;
            R.lossDb = reshape(R.lossDb(2:end),size(tt));
            R.grating = reshape(R.grating(2:end),size(tt));
            R.samples = [tt(:) pp(:)];
            R.freqGHz = S.freqOpGHz;
            S.coverageCache = struct('design',design,'settings',settings,'result',R);
            win.UserData = struct('result',R,'cacheHit',false);
            drawCoverage(); reportAnalysis(sprintf('Scan coverage: %d directions',numel(tt)));
        catch err
            if isgraphics(fig)
                setStatus(['Scan coverage failed: ' err.message],'bad');
                uialert(fig,err.message,'Scan coverage');
            end
        end
    end

    function coverageFinished(button)
        if isgraphics(button), button.Enable = 'on'; end
    end

    function drawCoverage()
        win = S.coverageWin; if ~isgraphics(win), return; end
        R = win.UserData.result; if isempty(R), return; end
        ax = findall(win,'Tag','coverageAxes'); cla(ax); P=pal();
        [rr,ang] = meshgrid(R.theta,[R.phi 360]);
        xx=rr.*cosd(ang); yy=rr.*sind(ang);
        zz=[R.lossDb;R.lossDb(1,:)];
        surf(ax,xx,yy,zeros(size(xx)),zz,'EdgeColor','none','Tag','coverageHeat');
        view(ax,2); axis(ax,'equal'); hold(ax,'on');
        if max(zz(:))-min(zz(:))>1e-8
            % Rounded levels are readable on a coarse engineering map.
            step=max(.5,ceil((max(zz(:))-min(zz(:)))/6*2)/2);
            levels=ceil(min(zz(:))/step)*step:step:floor(max(zz(:))/step)*step;
            contour(ax,xx,yy,zz,levels,'LineColor',P.ink, ...
                'ShowText','on','LabelSpacing',300,'Tag','coverageContour');
        end
        a=linspace(0,360,361);
        for radius=linspace(R.theta(end)/3,R.theta(end),3)
            plot(ax,radius*cosd(a),radius*sind(a),':','Color',P.gridLines,'Tag','coverageRing');
        end
        [rt,pt]=meshgrid(R.theta,R.phi);
        plot(ax,rt(R.grating).*cosd(pt(R.grating)),rt(R.grating).*sind(pt(R.grating)), ...
            'x','Color',P.bad,'LineWidth',1.5,'Tag','coverageGrating','MarkerSize',7);
        hold(ax,'off'); colormap(ax,parula);
        cbCoverage = colorbar(ax); cbCoverage.Label.String='Change from broadside (dB)';
        title(ax,sprintf('Scan coverage at %g GHz',R.freqGHz));
        xlabel(ax,'θ cos(azimuth) (°)'); ylabel(ax,'θ sin(azimuth) (°)');
        lim=R.theta(end)*1.08; xlim(ax,[-lim lim]); ylim(ax,[-lim lim]);
        txt = ['Radius = θ from broadside; azimuth runs counterclockwise from +x. ' ...
            'Contours show dB relative to broadside; crosses mark possible grating lobes ' ...
            '(including grazing). Rectangular complete lattices use exact visible-space ' ...
            'replicas; other layouts use the main pattern’s numerical lobe classifier.'];
        if R.unitCell, txt=[txt ' Levels preserve the embedded gain scale.'];
        else, txt=[txt ' Levels compare intended-direction directivity.']; end
        set(findall(win,'Tag','coverageNote'),'Text',txt);
    end

    function flag = coverageGrating(th,ph,af,thetaGrid,phiGrid)
        % Exact reciprocal-lattice visibility avoids a sampled-grid error
        % at grazing. With skew, holes or manual phases, use the same
        % numerical classifier as the main pattern instead of pretending
        % that nominal dx/dy describe the actual array.
        regular = strcmp(S.mode,'Uniform grid') && S.M>1 && S.N>1 && ...
            S.gridAngle==90 && S.stagger==0 && size(S.el,1)==S.M*S.N && ...
            all(S.el(:,3)>0) && ~S.seqPhase && max(S.el(:,4))-min(S.el(:,4))<1e-9;
        if regular
            [gx,gy]=latticePosAt(S.elRC(:,1),S.elRC(:,2),S.dx,S.dy,90,0);
            regular=max(abs(S.el(:,1)-gx))<1e-9 && max(abs(S.el(:,2)-gy))<1e-9;
        end
        if regular
            fr=freqRatio(); u=phaseFreqRatio()/fr*sind(th)*cosd(ph);
            v=phaseFreqRatio()/fr*sind(th)*sind(ph);
            % Only the nearest lattice orders can minimize distance to
            % the visible disk. Bound the work even at extreme f/f0.
            mx=unique([round(-u*S.dx*fr)+(-1:1) -1:1]);
            my=unique([round(-v*S.dy*fr)+(-1:1) -1:1]);
            [ix,iy]=meshgrid(mx,my);
            visible=(u+ix/(S.dx*fr)).^2+(v+iy/(S.dy*fr)).^2<=1+1e-10;
            flag=any(visible & (ix~=0 | iy~=0),'all');
        else
            keep=thetaGrid<=90; db=20*log10(max(abs(af(:,keep)),realmin));
            [~,~,gd,~,~,id]=findGratingLobe(db,thetaGrid(keep),phiGrid);
            flag=gd>-6 || id>-6;
        end
    end

    function showBandSweep()
        %SHOWBANDSWEEP  Analysis > Band sweep...: open or raise it, and run.
        %   Runs on every call, as Scan loss does, so the window always
        %   shows the current design; a window already open keeps the
        %   band and settings typed into it.
        if ~isgraphics(fig), return; end
        if isempty(bandWin) || ~isgraphics(bandWin)
            try
                buildBandWindow();
            catch buildErr
                % Drawing a new window takes in events, so it can be
                % closed (or the app with it) before it is finished; the
                % next line of the build then fails on the deleted window.
                % Nothing is left to run then. Anything else is a bug.
                if isgraphics(bandWin) && isgraphics(fig)
                    rethrow(buildErr);
                end
                return;
            end
        else
            selectTask('Band Sweep');
        end
        runBandSweep();
    end

    function buildBandWindow()
        %BUILDBANDWINDOW  Settings rows and three full-height plots.
        kUnit = freqScale(); unitTxt = S.freqUnit;
        win = resultPane('Band Sweep','bandSweepWin');
        bandWin = win;   % at once: showBandSweep checks it if the build fails
        % The unit is fixed per window: the fields are read with the
        % scale they were written in, even if the app's unit changes.
        win.UserData = struct('kUnit',kUnit,'unit',unitTxt, ...
            'table',[],'result',[]);
        gW = uigridlayout(win,[3 1]);
        gW.RowHeight = {'fit','fit','1x'};
        gW.Padding = [8 5 8 5]; gW.RowSpacing = 4;

        gA = uigridlayout(gW,[1 8]);
        gA.ColumnWidth = {'fit',96,'fit',96,'fit',64,'fit',96};
        gA.Padding = [0 0 0 0]; gA.ColumnSpacing = 6;
        uilabel(gA,'Text',sprintf('From (%s)',unitTxt));
        uispinner(gA,'Tag','bandFmin','Limits',[0.001 1e6]*kUnit, ...
            'Step',0.1*kUnit,'ValueDisplayFormat','%.5g', ...
            'Tooltip',sprintf('Lowest frequency of the sweep (%s).',unitTxt));
        uilabel(gA,'Text',sprintf('To (%s)',unitTxt));
        uispinner(gA,'Tag','bandFmax','Limits',[0.001 1e6]*kUnit, ...
            'Step',0.1*kUnit,'ValueDisplayFormat','%.5g', ...
            'Tooltip',sprintf('Highest frequency of the sweep (%s).',unitTxt));
        uilabel(gA,'Text','Points');
        uispinner(gA,'Tag','bandPoints','Limits',[2 401],'Step',1, ...
            'RoundFractionalValues','on', ...
            'Tooltip',['Number of frequencies, evenly spaced from the ' ...
            'lowest to the highest (2-401).']);
        uilabel(gA,'Text',sprintf('Phases set at (%s)',unitTxt));
        bwPhase = uispinner(gA,'Tag','bandFphase', ...
            'Limits',[0.001 1e6]*kUnit,'Step',0.1*kUnit, ...
            'ValueDisplayFormat','%.5g', ...
            'Tooltip',sprintf(['Frequency (%s) the phase shifters are ' ...
            'set at. The beam points at the commanded direction there ' ...
            'and squints elsewhere: sin θ = sin θs · f_set / f. ' ...
            'Default: the design frequency.'],unitTxt));

        gB = uigridlayout(gW,[1 6]);
        gB.ColumnWidth = {'fit','fit','1x',100,84,110};
        gB.Padding = [0 0 0 0]; gB.ColumnSpacing = 8;
        uicheckbox(gB,'Text','Phase shifters','Tag','bandUsePS', ...
            'Value',true, ...
            'Tooltip',['Steer with phase shifters set at the phase-set ' ...
            'frequency: fixed phases, so the beam squints with frequency.'], ...
            'ValueChangedFcn',@(s,e)set(bwPhase,'Enable', ...
            ternStr(s.Value,'on','off')));
        uicheckbox(gB,'Text','True-time delay','Tag','bandUseTTD', ...
            'Value',true, ...
            'Tooltip',['Steer with true-time delays: they hold the beam on ' ...
            'the commanded direction at every frequency (no squint).']);
        uilabel(gB,'Text','');
        uibutton(gB,'Text','Run sweep','Tag','bandRun', ...
            'Tooltip','Sweep the current design with these settings', ...
            'ButtonPushedFcn',@(s,e)runBandSweep());
        uibutton(gB,'Text','Defaults','Tag','bandDefaults', ...
            'Tooltip',['Band, points and phase-set frequency for the ' ...
            'current design frequency'], ...
            'ButtonPushedFcn',@(s,e)bandDefaults(win));
        uibutton(gB,'Text','Export CSV…','Tag','bandExport', ...
            'Enable','off', ...
            'Tooltip',['Save the swept values as CSV: frequency (GHz), ' ...
            'angles and HPBW (°), levels (dBi), pointing loss (dB)'], ...
            'ButtonPushedFcn',@(s,e)exportBandCSV());

        % Tiled layout aligns the actual plot boxes, even when their Y-axis
        % tick labels have different widths. It also shares the available
        % height more efficiently than three separate grid cells.
        plotPane = uipanel(gW,'BorderType','none');
        plots = tiledlayout(plotPane,3,1, ...
            'TileSpacing','compact','Padding','compact');
        ax = uiaxes(plots,'Tag','bandAxPointing'); ax.Layout.Tile = 1;
        ax = uiaxes(plots,'Tag','bandAxGain'); ax.Layout.Tile = 2;
        ax = uiaxes(plots,'Tag','bandAxHPBW'); ax.Layout.Tile = 3;
        bandDefaults(win);
        drawBandResults(win,[]);
    end

    function h = bandCtl(win, tag)
        h = findall(win,'Tag',tag);
    end

    function bandDefaults(win)
        %BANDDEFAULTS  Band fields for the current design frequency.
        %   The Ka-band SATCOM receive band, 17.7-21.2 GHz in 0.1 GHz
        %   steps, when the design frequency is at or near it -- a receive
        %   array's geometry reference often sits at or just above the
        %   top channel (21.4 GHz in the Ka Rx design this app was built
        %   for). Otherwise the design frequency +-10 % in 1 % steps.
        %   Phases are set at the design frequency, the reference the main
        %   window's Beam squint mode uses too.
        KA_RX_GHZ = [17.7 21.2];
        kUnit = win.UserData.kUnit;
        f0 = S.freqGHz;
        if f0 >= 0.95*KA_RX_GHZ(1) && f0 <= 1.05*KA_RX_GHZ(2)
            fLo = KA_RX_GHZ(1); fHi = KA_RX_GHZ(2); nFreq = 36;
        else
            fLo = 0.9*f0; fHi = 1.1*f0; nFreq = 21;
        end
        set(bandCtl(win,'bandFmin'),'Value',fLo*kUnit);
        set(bandCtl(win,'bandFmax'),'Value',fHi*kUnit);
        set(bandCtl(win,'bandPoints'),'Value',nFreq);
        set(bandCtl(win,'bandFphase'),'Value',f0*kUnit,'Enable','on');
        set(bandCtl(win,'bandUsePS'),'Value',true);
        set(bandCtl(win,'bandUseTTD'),'Value',true);
    end

    function runBandSweep()
        %RUNBANDSWEEP  Run the sweep: busy state, errors, one at a time.
        if isempty(bandWin) || ~isgraphics(bandWin) || ~isgraphics(fig)
            return;
        end
        win = bandWin;
        bwRun = bandCtl(win,'bandRun');
        if strcmp(bwRun.Enable,'off'), return; end   % already running
        % This window is not covered by the main window's progress dialog,
        % so its Run button can fire while Scan loss or Max scan is
        % mid-sweep (their progress steps take in clicks). Running inside
        % them would reset the element-fallback flag they are collecting.
        if busyDepth > 0
            uialert(fig,['Another calculation is running in the main ' ...
                'window. Run the band sweep when it has finished.'], ...
                'Band sweep');
            return;
        end
        bwRun.Enable = 'off';
        fig.Pointer = 'watch';
        runDone = onCleanup(@()bandRunEnded(win));
        busyTok = beginBusy('Running the band sweep…'); %#ok<NASGU>
        try
            runBandSweepCore(win);
        catch bandErr
            setStatus(['Band sweep failed: ' bandErr.message],'bad');
            if isgraphics(fig)
                uialert(fig,bandErr.message,'Band sweep failed');
            end
        end
    end

    function bandRunEnded(win)
        if ~isgraphics(win), return; end
        if isgraphics(fig), fig.Pointer = 'arrow'; end
        set(bandCtl(win,'bandRun'),'Enable','on');
    end

    function runBandSweepCore(win)
        if ~patternReady(true), return; end
        % beginBusy and patternReady both redraw, which takes in events:
        % either window may have been closed since the caller checked.
        if ~isgraphics(win) || ~isgraphics(fig), return; end
        fu = win.UserData; kUnit = fu.kUnit;
        fLo = bandCtl(win,'bandFmin').Value/kUnit;
        fHi = bandCtl(win,'bandFmax').Value/kUnit;
        nFreq = round(bandCtl(win,'bandPoints').Value);
        fPh = bandCtl(win,'bandFphase').Value/kUnit;
        usePS = logical(bandCtl(win,'bandUsePS').Value);
        useTTD = logical(bandCtl(win,'bandUseTTD').Value);
        % Invalid settings need a direct message now that the result pane
        % gives its full height to the plots.
        if ~(fHi > fLo)
            uialert(fig,'The highest frequency must be above the lowest.', ...
                'Invalid band');
            setStatus('Band sweep not run: the band is empty.','warn');
            return;
        end
        if ~usePS && ~useTTD
            uialert(fig,'Select Phase shifters, True-time delay, or both.', ...
                'No steering mode selected');
            setStatus('Band sweep not run: no steering chosen.','warn');
            return;
        end
        if isempty(S.el)
            uialert(fig,'Place at least one element first.','No elements');
            return;
        end
        w = calculationWeights();
        if ~(sum(abs(w)) > 0)
            uialert(fig,['Every element amplitude is zero - nothing to ' ...
                'sweep. Set at least one amplitude to a nonzero value.'], ...
                'Zero amplitude');
            return;
        end
        % Cleared ahead of every element evaluation below, and checked at
        % the end -- see showScanLossCore for why the order matters.
        S.efFallback = false; S.efFallbackWhich = [false false];
        tBand = tic;
        % The progress dialog blocks edits to this design during the sweep.
        selectTask('Band Sweep');
        [prog, progDone] = startProgress('Band sweep', ...
            'Evaluating the element pattern…',true,fig); %#ok<ASGLU>

        fGHz = linspace(fLo,fHi,nFreq).';
        isPS = [true false]; isPS = isPS([usePS useTTD]);
        nM = numel(isPS);
        thS = S.theta_s; phS = S.phi_s;
        us = sind(thS)*cosd(phS); vs = sind(thS)*sind(phS);
        sgn = scanThetaSign();     % + = away from broadside, either side
        % The cut through the commanded direction, sampled exactly as the
        % main window's theta cut: signed theta, the negative side on
        % phi+180. At a frequency and phase setting the main window can
        % show, the two curves are the same numbers.
        cutX = -90:0.1:90; nCut = numel(cutX);
        uc = (sind(cutX)*cosd(phS)).'; vc = (sind(cutX)*sind(phS)).';
        thEF = abs(cutX); phEF = phS + 180*(cutX < 0);
        x = S.el(:,1); y = S.el(:,2); rot = S.el(:,5);
        % Element patterns by rotation, once for the whole sweep: they do
        % not depend on frequency (an imported pattern is single-frequency
        % data and is held fixed -- the note under the plots says so).
        [uRot,~,rIdx] = unique(rot); nUR = numel(uRot);
        EthC = zeros(nCut,nUR); EphC = zeros(nCut,nUR);
        eth0 = zeros(1,nUR); eph0 = zeros(1,nUR);
        for q = 1:nUR
            [a,b] = elementFactor(thEF,phEF,uRot(q));
            EthC(:,q) = a(:); EphC(:,q) = b(:);
            [eth0(q),eph0(q)] = elementFactor(abs(thS), ...
                phS+180*(thS < 0),uRot(q));
        end
        ph0 = deg2rad(S.el(:,4));
        if S.seqPhase, ph0 = ph0 + seqRotSign()*deg2rad(rot); end
        steerGeo = -S.k*(x*us + y*vs);   % steering phase per unit ratio

        % Level basis, as the main window reports it: directivity for an
        % analytic element, realized gain for an imported one (its
        % efficiency), and the unit-cell relation G = G_emb*|sum|^2/sum|w|^2
        % for an embedded-element import, whose pattern integral is not
        % its radiated power (see unitCellFigures).
        ucFig = unitCellFigures();
        effLin = 1; basis = 'Directivity';
        if ~isempty(ucFig)
            basis = 'Realized gain';
            ucOff = ucFig.rg - 10*log10(sum(w.^2));
        elseif strcmp(S.efType,'Imported (CST far-field)')
            effLin = impEffLin();
            if isfinite(effLin) && effLin > 0
                basis = 'Realized gain';
            else
                effLin = 1;
            end
        end

        thPk = nan(nFreq,nM); lvlCmd = thPk; lvlPk = thPk; hp = thPk;
        edge = false(nFreq,nM); law = nan(nFreq,1);
        for iF = 1:nFreq
            % The window itself is checked too: without a dialog (it can
            % fail to open) stepProgress cannot see this window close.
            if ~isgraphics(win) || stepProgress(prog,(iF-1)/nFreq, ...
                    sprintf('Frequency %d of %d (%s)',iF,nFreq, ...
                    bandFreqText(fu,fGHz(iF))))
                reportAnalysis('Band sweep cancelled: no result');
                return;
            end
            fr = freqRatio(fGHz(iF)); kOp = S.k*fr;
            % Phase ratio per mode, the role phaseFreqRatio() plays in the
            % main window: phase shifters keep the phases computed at the
            % phase-set frequency; a true-time delay is the retuned case
            % at every frequency.
            phRatio = zeros(1,nM);
            phRatio(isPS) = freqRatio(fPh); phRatio(~isPS) = fr;
            W = w .* exp(1j*(steerGeo*phRatio + ph0));   % N x modes
            % Both modes share one propagation matrix per block of
            % elements; blocks bound its memory on very large arrays.
            Ect = zeros(nCut,nM); Ecp = Ect; Et0 = zeros(1,nM); Ep0 = Et0;
            for q = 1:nUR
                idx = find(rIdx == q);
                A = zeros(nCut,nM);
                for b0 = 1:512:numel(idx)
                    blk = idx(b0:min(b0+511,numel(idx)));
                    A = A + exp(1j*kOp*(uc*x(blk).' + vc*y(blk).')) ...
                        * W(blk,:);
                end
                Ect = Ect + A.*EthC(:,q); Ecp = Ecp + A.*EphC(:,q);
                a0 = exp(1j*kOp*(x(idx)*us + y(idx)*vs)).' * W(idx,:);
                Et0 = Et0 + a0*eth0(q); Ep0 = Ep0 + a0*eph0(q);
            end
            kernels = {};   % pair kernels depend on frequency, not mode
            for m = 1:nM
                if ~isempty(ucFig)
                    offDb = ucOff;
                else
                    switch S.efType
                        case 'Isotropic'
                            Prad = exactIsotropicPrad(W(:,m),fr);
                        case 'cos^q(theta)'
                            [Prad,kernels] = exactCosQPrad(W(:,m),kernels,fr);
                        case 'Short dipole (z-axis)'
                            [Prad,kernels] = exactZDipolePrad(W(:,m),kernels,fr);
                        case {'Dipole (linear pol)','Patch (cos^q x lin pol)'}
                            [Prad,kernels] = exactVectorPrad(W(:,m),rot, ...
                                kernels,fr);
                        otherwise
                            Prad = quadraturePrad(W(:,m),fr);
                    end
                    offDb = NaN;
                    if Prad > 0, offDb = 10*log10(4*pi*effLin/Prad); end
                end
                % Total power, both polarisations, as Scan loss uses: the
                % pointing and level of the beam, not of one readout.
                Pcut = fieldPower(Ect(:,m),Ecp(:,m));
                if ~any(Pcut > 0), continue; end   % nothing radiates
                cutDb = 10*log10(Pcut);
                % Start from the steering term's own direction and climb
                % to the lobe there. The global maximum is not the beam
                % once a grating lobe enters at the top of the band.
                sBeam = squintSinTheta(thS,phRatio(m),fr);
                if abs(sBeam) <= 1
                    thPred = asind(sBeam);
                    if isPS(m), law(iF) = thPred; end
                else
                    thPred = sign(sBeam)*90;
                end
                [~,i0] = min(abs(cutX - thPred));
                ipk = bandClimb(cutDb,i0);
                [thPk(iF,m),pkDb] = bandRefinePeak(cutX,cutDb,ipk);
                edge(iF,m) = abs(sBeam) > 1 || ipk == 1 || ipk == nCut;
                [okR,ra] = crossOut(cutDb,cutX,ipk,+1,nCut,false,cutDb(ipk)-3);
                [okL,la] = crossOut(cutDb,cutX,ipk,-1,nCut,false,cutDb(ipk)-3);
                if okR && okL, hp(iF,m) = ra - la; end
                lvlCmd(iF,m) = offDb + 10*log10(fieldPower(Et0(m),Ep0(m)));
                lvlPk(iF,m) = offDb + pkDb;
            end
        end
        closeProgress(prog);
        % A window closed during the last step leaves nothing to draw in.
        if ~isgraphics(fig) || ~isgraphics(win), return; end
        if S.efFallback
            drawBandResults(win,[]);
            setStatus(['Band sweep refused: the Custom formula failed ' ...
                'on its angles.'],'warn');
            uialert(fig,['No valid band sweep.' newline newline ...
                'The Custom element formula could not be evaluated at ' ...
                'every angle this sweep needs, so an element was ' ...
                'substituted for it. Every curve would describe that ' ...
                'substitute rather than your design, so none is drawn.'], ...
                'Band sweep','Icon','warning');
            return;
        end
        R = struct('f',fGHz,'fPhase',fPh,'isPS',isPS,'thS',thS,'phS',phS, ...
            'thPk',thPk,'err',sgn*(thPk - thS),'law',law, ...
            'errLaw',sgn*(law - thS),'lvlCmd',lvlCmd,'lvlPk',lvlPk, ...
            'loss',lvlPk - lvlCmd,'hp',hp,'edge',edge,'basis',basis, ...
            'nEl',size(S.el,1),'efType',S.efType, ...
            'unitCell',~isempty(ucFig), ...
            'when',char(datetime('now','Format','HH:mm')));
        drawBandResults(win,R);
        reportAnalysis(sprintf('Band sweep: %d frequencies, %s to %s, in %.1f s', ...
            nFreq,bandFreqText(fu,fLo),bandFreqText(fu,fHi),toc(tBand)));
    end

    function iPk = bandClimb(db,iPk)
        %BANDCLIMB  Walk uphill from sample iPk to the top of its lobe.
        n = numel(db);
        while true
            best = iPk;
            if iPk > 1 && db(iPk-1) > db(best), best = iPk-1; end
            if iPk < n && db(iPk+1) > db(best), best = iPk+1; end
            if best == iPk, return; end
            iPk = best;
        end
    end

    function [xPk,dbPk] = bandRefinePeak(xs,db,iPk)
        %BANDREFINEPEAK  Parabola through the peak sample and its neighbours.
        %   The cut is sampled every 0.1 degree; squint moves the beam by
        %   fractions of that between neighbouring frequencies, so the raw
        %   sample index would draw the pointing curve as a staircase.
        xPk = xs(iPk); dbPk = db(iPk);
        if iPk <= 1 || iPk >= numel(db), return; end
        a = db(iPk-1); b = db(iPk); c3 = db(iPk+1);
        den = a - 2*b + c3;
        if ~all(isfinite([a b c3])) || den >= 0, return; end
        d = 0.5*(a - c3)/den;             % offset in samples, |d| <= 0.5
        xPk = xs(iPk) + d*(xs(iPk+1) - xs(iPk));
        dbPk = b - 0.25*(a - c3)*d;
    end

    function s = bandFreqText(fu,fGHz)
        % fu: the window's unit (its UserData), copied by the sweep so a
        % window closed mid-sweep leaves nothing here to read from.
        s = sprintf('%.4g %s',fGHz*fu.kUnit,fu.unit);
    end

    function drawBandResults(win,R)
        %DRAWBANDRESULTS  Plots and export table; R [] clears.
        P = pal(); fu = win.UserData;
        axE = bandCtl(win,'bandAxPointing');
        axG = bandCtl(win,'bandAxGain');
        axH = bandCtl(win,'bandAxHPBW');
        for ax = [axE axG axH]
            % cla leaves graphics with HandleVisibility='off' in place,
            % including the old phase-set marker. Clear those too before
            % drawing a sweep with a different phase-set frequency.
            delete(allchild(ax));
            cla(ax); legend(ax,'off'); grid(ax,'on');
        end
        title(axE,'Beam offset from target (+ = farther from broadside)', ...
            'FontSize',11);
        ylabel(axE,'Offset (°)');
        title(axG,'Toward the commanded direction','','FontSize',11);
        ylabel(axG,'dBi');
        title(axH,'Half-power beamwidth in the scan plane','FontSize',11);
        ylabel(axH,'HPBW (°)');
        xlabel(axH,sprintf('Frequency (%s)',fu.unit));
        % Only the bottom chart needs frequency labels. Suppressing the
        % repeated labels gives each plot more vertical drawing space.
        axE.XTickLabel = [];
        axG.XTickLabel = [];
        bwExport = bandCtl(win,'bandExport');
        if isempty(R)
            win.UserData.table = []; win.UserData.result = [];
            bwExport.Enable = 'off';
            return;
        end
        fx = R.f*fu.kUnit;
        iPS = find(R.isPS); iTT = find(~R.isPS);
        for ax = [axE axG axH], hold(ax,'on'); end
        yline(axE,0,':','Color',P.muted,'Tag','bandRefLine', ...
            'HandleVisibility','off');
        hLeg = gobjects(0);
        if ~isempty(iPS)
            hLeg(end+1) = plot(axE,fx,R.errLaw,'--','Color',P.refTrace, ...
                'LineWidth',1.2,'Tag','bandLaw', ...
                'DisplayName','sin-law prediction');
            psName = sprintf('Phase shifters, set at %s', ...
                bandFreqText(fu,R.fPhase));
            hLeg(end+1) = bandLine(axE,fx,R.err(:,iPS),'bandPS',psName);
            plot(axG,fx,R.lvlPk(:,iPS),'--','Color',P.traceTotal, ...
                'LineWidth',1.1,'Tag','bandPeakPS', ...
                'DisplayName','Beam peak, phase shifters (gap = pointing loss)');
            bandLine(axG,fx,R.lvlCmd(:,iPS),'bandPS',psName);
            bandLine(axH,fx,R.hp(:,iPS),'bandPS',psName);
            if R.fPhase >= R.f(1) && R.fPhase <= R.f(end)
                xline(axE,R.fPhase*fu.kUnit,':','phases set', ...
                    'Color',P.muted,'Tag','bandRefLine', ...
                    'LabelVerticalAlignment','bottom', ...
                    'HandleVisibility','off');
            end
        end
        if ~isempty(iTT)
            hLeg(end+1) = bandLine(axE,fx,R.err(:,iTT),'bandTTD', ...
                'True-time delay');
            bandLine(axG,fx,R.lvlCmd(:,iTT),'bandTTD','True-time delay');
            bandLine(axH,fx,R.hp(:,iTT),'bandTTD','True-time delay');
        end
        for ax = [axE axG axH]
            hold(ax,'off'); xlim(ax,[fx(1) fx(end)]);
        end
        % cla does not reliably restore automatic Y limits after a previous
        % sweep. A different design can otherwise make HPBW disappear.
        ylim(axG,'auto');
        ylim(axH,'auto');
        % A flat zero curve (broadside, or true-time delay alone) would
        % otherwise be autoscaled to a 1e-14 degree axis of noise.
        errAll = [R.err(:); R.errLaw(:)];
        yMax = max([0.5; 1.15*abs(errAll(isfinite(errAll)))]);
        ylim(axE,[-yMax yMax]);
        legend(axE,hLeg,'Location','best','AutoUpdate','off');
        title(axG,sprintf('%s toward the commanded direction',R.basis), ...
            '','FontSize',11);
        win.UserData.result = R;
        win.UserData.table = bandTable(R);
        bwExport.Enable = 'on';
    end

    function h = bandLine(ax,fx,yv,tag,name)
        % One swept mode's curve; the Tag is its colour role.
        P = pal();
        if strcmp(tag,'bandPS'), col = P.traceTotal; mk = 'o';
        else, col = P.traceAR; mk = 's'; end
        h = plot(ax,fx,yv,'-','Marker',mk,'MarkerSize',3, ...
            'MarkerFaceColor',col,'Color',col,'LineWidth',1.8, ...
            'Tag',tag,'DisplayName',name);
    end

    function T = bandTable(R)
        %BANDTABLE  One row per frequency and steering mode, for export.
        modeNames = {'True-time delay','Phase shifters'};
        nF = numel(R.f); T = table();
        for m = 1:numel(R.isPS)
            fSet = nan(nF,1);
            if R.isPS(m), fSet(:) = R.fPhase; end
            Tm = table(R.f,repmat(modeNames(R.isPS(m)+1),nF,1),fSet, ...
                repmat(R.thS,nF,1),R.thPk(:,m),R.err(:,m),R.lvlCmd(:,m), ...
                R.lvlPk(:,m),R.loss(:,m),R.hp(:,m), ...
                'VariableNames',{'Frequency_GHz','Steering', ...
                'PhaseSet_GHz','Commanded_theta_deg','Beam_theta_deg', ...
                'Pointing_error_deg','Toward_commanded_dBi', ...
                'Beam_peak_dBi','Pointing_loss_dB','HPBW_deg'});
            T = [T; Tm]; %#ok<AGROW>
        end
    end

    function exportBandCSV()
        %EXPORTBANDCSV  The last sweep's table as CSV (Export CSV... button).
        if isempty(bandWin) || ~isgraphics(bandWin), return; end
        T = bandWin.UserData.table;
        if isempty(T), return; end
        [fn,pn] = uiputfile('*.csv','Save the band sweep as CSV', ...
            'band_sweep.csv');
        if isgraphics(bandWin), selectTask('Band Sweep'); end
        if isequal(fn,0) || ~isgraphics(bandWin), return; end
        try
            writetable(T,fullfile(pn,fn));
            setStatus(sprintf('Saved the band sweep (%d rows) to %s', ...
                height(T),fn),'good');
        catch expErr
            setStatus(['Band sweep export failed: ' expErr.message],'bad');
            uialert(fig,expErr.message,'Export failed');
        end
    end

    function findMaxScanAngle()
        showMaxScanResult('Searching the steering directions…');
        busyTok = beginBusy('Searching for the max scan angle…'); %#ok<NASGU>
        try
            findMaxScanAngleCore();
        catch scanErr
            setStatus(['Max scan failed: ' scanErr.message],'bad');
            if isgraphics(fig), showMaxScanResult(scanErr.message); end
        end
    end

    function showMaxScanResult(message)
        pane = resultPanels{4};
        if isempty(pane) || ~isgraphics(pane)
            pane = resultPane('Max Scan','maxScanView');
            mg = uigridlayout(pane,[2 1]);
            mg.RowHeight = {42,'1x'};
            bar = uigridlayout(mg,[1 2]);
            bar.ColumnWidth = {'1x',120};
            uilabel(bar,'Text','Maximum scan angle','FontWeight','bold', ...
                'FontSize',17);
            uibutton(bar,'Text','Run again','Tag','maxScanRun', ...
                'ButtonPushedFcn',@(~,~)findMaxScanAngle());
            uitextarea(mg,'Tag','maxScanResult','Editable','off');
        else
            selectTask('Max Scan');
        end
        set(findall(pane,'Tag','maxScanResult'),'Value',splitlines(message));
    end

    function safeCompute()
        refreshControlVisibility();
        % recompute for display-only changes, without alerting on an
        % empty / zero-amplitude array
        if ~isempty(S.el) && sum(abs(S.el(:,3))) > 0, computePattern(); end
    end

    function assignCutMode(v), S.cutMode = v; safeCompute(); end
    function assignCutFixedTheta(v)
        if useAzEl(), S.cutFixedTheta = 90-v;
        else, S.cutFixedTheta = v; end
        safeCompute();
    end
    function assignCutPhi(v)
        S.cutPhi = v;
        if useAzEl(), spCutPhi.Value = mod(v,360); end
        safeCompute();
    end
    function assignCutPhiFollow(v)
        %ASSIGNCUTPHIFOLLOW  Link/unlink the cut plane from the steer phi.
        %   Unticking must LEAVE THE CUT WHERE IT IS. S.cutPhi still held
        %   whatever azimuth was last set by hand, which could be from a
        %   different design entirely, so releasing the link made the cut
        %   jump to that old value instead of staying on the plane the
        %   user was looking at.
        if ~v && S.cutPhiFollow
            S.cutPhi = cutPhiVal();      % freeze at the displayed plane
            spCutPhi.Value = S.cutPhi;
        end
        S.cutPhiFollow = v;
        if S.cutPhiFollow
            spCutPhi.Enable = 'off'; spCutPhi.Value = cutPhiVal();
        else
            spCutPhi.Enable = 'on';
        end
        syncPattern2DFixedAngle();
        safeCompute();
    end
    function pinCutPhi(v)
        % Unlink and set, in one action -- the two-step the buttons exist
        % to replace. Goes through the same state the checkbox and spinner
        % write, so nothing here is a second code path.
        S.cutPhiFollow = false; cbCutFollow.Value = false;
        spCutPhi.Enable = 'on';
        S.cutPhi = v; spCutPhi.Value = v;
        syncPattern2DFixedAngle();
        safeCompute();
    end
    function followCutPhi()
        S.cutPhiFollow = true; cbCutFollow.Value = true;
        spCutPhi.Enable = 'off'; spCutPhi.Value = cutPhiVal();
        syncPattern2DFixedAngle();
        safeCompute();
    end

    function p = cutPhiVal()
        % The plane the theta cut is taken in.
        if S.cutPhiFollow
            p = S.phi_s;
            if useAzEl() && S.theta_s<0, p = mod(p+180,360); end
        else
            p = S.cutPhi;
        end
        if useAzEl(), p = mod(p,360); end
    end

    function did = maybeCompute()
        refreshControlVisibility();
        refreshSpacingMM();
        if ~isempty(S.coverageCache) && ~isequaln(S.coverageCache.design,designSnapshot())
            S.coverageCache = [];
            if isgraphics(S.coverageWin)
                set(findall(S.coverageWin,'Tag','coverageNote'),'Text', ...
                    'Design changed. This map is the previous result; Run coverage to update it.');
            end
        end
        % Keep the plots in step with the controls, unless the user has
        % turned Auto off (large arrays can make every edit sluggish).
        % patternReady(false) has a deliberate side effect when an imported
        % pattern is missing: it clears the plots and prints the SPECIFIC
        % blocking reason. Preserve that message instead of immediately
        % overwriting it with the generic 'Settings changed' text.
        % Every design change arrives here (refreshAll ends here too), so
        % this is the one place a changed design becomes an undo step --
        % first, so a compute that throws cannot lose the step.
        undoCheckpoint();
        hasActive = ~isempty(S.el) && any(S.el(:,3)~=0);
        ready = patternReady(false);
        did = cbAuto.Value && hasActive && ready;
        if did
            computePattern();
        elseif ~hasActive
            invalidateRadiation('No radiation: no active elements.');
            refreshPhaseMaps();
        elseif ~ready
            % patternReady already invalidated the radiation result and
            % explained what is missing. Only keep open phase maps in sync.
            refreshPhaseMaps();
        else
            invalidateRadiation('Settings changed: click Compute pattern.');
            refreshPhaseMaps();
        end
    end

    function onMode(src,~)
        prevMode = S.mode;      % captured BEFORE the overwrite -- see below
        S.mode = src.Value;
        S.pending = []; S.sel = [];
        syncTableSelection();
        if strcmp(S.mode,'Uniform grid')
            rebuildUniform();
        elseif strcmp(prevMode,'Uniform grid') && ...
                ~isempty(S.el) && (S.gridAngle ~= 90 || S.stagger ~= 0)
            % Gated on the PREVIOUS mode, not just the new one. The test
            % used to read only "the new mode is not Uniform", but S.mode
            % had already been overwritten one line up, so the condition
            % was true for EVERY transition into a non-Uniform mode --
            % Sparse -> Sub-position and Sub-position -> Sparse included.
            % Those elements are already on the orthogonal lattice, so
            % the block computed their offsets against the skewed lattice
            % a second time and shifted the whole array again: a measured
            % 0.25 lambda walk per switch, silent, and cumulative if the
            % user toggled back and forth. The conversion is only ever
            % meaningful once, on the way OUT of Uniform grid.
            %
            % Leaving Uniform grid mode with a skewed/staggered lattice
            % active: Sparse/Sub-position always place elements on the
            % plain orthogonal dx,dy grid (see onGridClick/onGeom), and
            % refreshLayout's background dots already switch to that
            % plain grid outside Uniform mode -- so without this, any
            % elements carried over from a skewed Uniform-grid session
            % would sit off the (now-orthogonal) reference dots.
            % Re-reference onto the orthogonal lattice while PRESERVING
            % each element's manual perturbation offset (moveSelected), the same "compute offset from old
            % lattice, reapply to new lattice" approach onGeom already
            % uses for a dx/dy/gridAngle/stagger change -- rather than
            % snapping every element exactly onto the lattice and
            % silently discarding any perturbation study in progress.
            [oldGX,oldGY] = latticePosAt(S.elRC(:,1),S.elRC(:,2),S.dx,S.dy,S.gridAngle,S.stagger);
            offx = S.el(:,1) - oldGX;
            offy = S.el(:,2) - oldGY;
            [newGX,newGY] = latticePosAt(S.elRC(:,1),S.elRC(:,2),S.dx,S.dy,90,0);
            S.el(:,1) = newGX + offx;
            S.el(:,2) = newGY + offy;
        end
        refreshAll();
    end

    function ang = clampedStaggerAngleDeg(stag,ddy)
        % Defensive fallback -- spStaggerAngle's Limits ([-89.5 89.5])
        % are already widened to cover every angle atan2d(stagger,dy)
        % can reach given spStagger/spDy's own Limits, so this clamp
        % should never actually trigger in practice. Kept anyway:
        % assigning an out-of-Limits value to a uispinner throws an
        % error and breaks the UI, and reading Limits from the spinner
        % itself (rather than hardcoding) keeps this correct even if
        % those Limits ever change.
        ang = atan2d(stag,ddy);
        ang = min(max(ang, spStaggerAngle.Limits(1)), spStaggerAngle.Limits(2));
    end

    function [gx,gy] = latticePosAt(row,col,ddx,ddy,ang,stag)
        % Skewed + staggered lattice position for lattice cell (row,col),
        % given explicit dx,dy,gridAngle,stagger -- NOT read from S, so
        % this can be evaluated at the OLD spacing/angle/stagger before S
        % is overwritten, and again at the NEW values, to compute/
        % preserve a manual perturbation offset across a geometry change.
        % gridAngle=90, stag=0 reduces exactly to the plain orthogonal
        % grid (col-1)*dx, (row-1)*dy used throughout the rest of the file.
        % Stagger: EVEN rows (2,4,6,...) get +stag along x; ODD rows
        % (1,3,5,...) get +0 -- independent of the gridAngle skew, which
        % applies to every row via the ddy*cosd/sind term below.
        rowIsEven = mod(row-1,2)==1;   % row=2,4,6,... -> true
        gx = (col-1).*ddx + (row-1).*ddy.*cosd(ang) + rowIsEven.*stag;
        gy = (row-1).*ddy.*sind(ang);
    end

    function onGeom(~,~)
        oldDx = S.dx; oldDy = S.dy; oldAngle = S.gridAngle; oldStag = S.stagger;
        oldM = S.M; oldN = S.N;
        % A full rebuild (new M or N) writes amplitude 1, phase 0, rot 0 into
        % every cell. Recorded here so the status line can say that hand
        % edits went with it, and that Undo brings them back.
        stepsWas = undoCount; lostEdits = false;
        hadEdits = ~isempty(S.el) && (any(S.el(:,4:5) ~= 0,'all') || ...
            strcmp(S.taper,'Manual (table)'));
        S.M = spM.Value; S.N = spN.Value;
        S.dx = spDx.Value; S.dy = spDy.Value;
        S.gridAngle = spGridAngle.Value;
        S.stagger = spStagger.Value;
        S.pending = [];
        if strcmp(S.mode,'Uniform grid')
            % If M,N are unchanged, most geometry edits only move the
            % existing cells. Circle/Hexagon/Octagon use the smaller
            % half-extent; Ellipse caps its height at 0.75*halfW. Thus a
            % dx/dy ASPECT RATIO change can change WHICH lattice cells
            % survive even when M and N do not change. Diamond alone is
            % normalized independently by halfW/halfH.
            countSame = (S.M==oldM) && (S.N==oldN);
            oldAspect = oldDx/oldDy;
            newAspect = S.dx/S.dy;
            aspectChanged = abs(oldAspect-newAspect) > ...
                1e-12*max([1 abs(oldAspect) abs(newAspect)]);
            maskDependsOnAspect = ismember(S.arrayShape, ...
                {'Circle','Hexagon','Octagon','Ellipse'});
            needRemask = countSame && ~isempty(S.el) && aspectChanged && maskDependsOnAspect;

            if needRemask
                % Snapshot the per-cell state and each manual xy perturbation,
                % rebuild the correct mask, then restore the state on cells
                % that still exist. Newly admitted cells keep the normal
                % rebuild defaults/taper; removed cells are intentionally gone.
                oldEl = S.el; oldRC = S.elRC;
                [oldGX,oldGY] = latticePosAt(oldRC(:,1),oldRC(:,2), ...
                    oldDx,oldDy,oldAngle,oldStag);
                oldOffX = oldEl(:,1)-oldGX;
                oldOffY = oldEl(:,2)-oldGY;

                rebuildUniform();

                [survives,oldIdx] = ismember(S.elRC,oldRC,'rows');
                if any(survives)
                    srcIdx = oldIdx(survives);
                    S.el(survives,3:5) = oldEl(srcIdx,3:5);
                    S.el(survives,1) = S.el(survives,1) + oldOffX(srcIdx);
                    S.el(survives,2) = S.el(survives,2) + oldOffY(srcIdx);
                end
            elseif countSame && ~isempty(S.el)
                % Same populated-cell set: preserve manual perturbations
                % instead of snapping elements back to lattice centres.
                [oldGX,oldGY] = latticePosAt(S.elRC(:,1),S.elRC(:,2),oldDx,oldDy,oldAngle,oldStag);
                offx = S.el(:,1) - oldGX;
                offy = S.el(:,2) - oldGY;
                [newGX,newGY] = latticePosAt(S.elRC(:,1),S.elRC(:,2),S.dx,S.dy,S.gridAngle,S.stagger);
                S.el(:,1) = newGX + offx;
                S.el(:,2) = newGY + offy;
            else
                rebuildUniform();
                lostEdits = hadEdits;
            end
        else
            % Sparse / Sub-position modes stay on the plain orthogonal
            % dx,dy lattice for cell placement (grid angle and shape
            % masking apply to Uniform grid mode only -- see the state
            % comment above S.arrayShape for the reasoning).
            if ~isempty(S.el)
                offx = S.el(:,1) - (S.elRC(:,2)-1)*oldDx;
                offy = S.el(:,2) - (S.elRC(:,1)-1)*oldDy;
                S.el(:,1) = (S.elRC(:,2)-1)*S.dx + offx;
                S.el(:,2) = (S.elRC(:,1)-1)*S.dy + offy;
            end
            keep = S.elRC(:,1)<=S.M & S.elRC(:,2)<=S.N;
            S.el = S.el(keep,:); S.elRC = S.elRC(keep,:);
            S.sel = [];
            % taper weights (Hamming/Hanning/Chebyshev/Taylor/Binomial)
            % depend on each element's (row,col) relative to the CURRENT
            % M,N -- rebuildUniform's path already recomputes these via
            % its own applyTaper() call, but this Sparse/Sub-position
            % path didn't, leaving surviving elements with stale weights
            % from whatever M,N was active when they were placed. No-op
            % for Uniform taper (all weights 1) and for Manual (table),
            % which applyTaper() itself skips.
            applyTaper();
            checkLargeArrayWarning('');   % M/N shrinking here can drop the count back down
        end
        % Keep the Stagger angle spinner in sync with the underlying
        % S.stagger/S.dy -- runs regardless of which control triggered
        % onGeom (dx/dy/gridAngle/stagger all funnel through here), so
        % editing dy alone (which changes the angle for a FIXED stagger
        % length) also updates the displayed angle correctly.
        spStaggerAngle.Value = clampedStaggerAngleDeg(S.stagger, S.dy);
        refreshAll();
        if lostEdits
            announceUndoable(stepsWas, sprintf(['%d × %d grid rebuilt: the ' ...
                'per-element amplitude, phase and rotation edits were reset'], ...
                S.M, S.N));
        end
    end

    function keep = insideRegularPolygon(xr,yr,nSides,apothem)
        % True where (xr,yr) lies inside a regular n-gon centred at the
        % origin with the given apothem (inradius), via the standard
        % half-plane intersection test: for a convex polygon, a point is
        % interior iff its projection onto every edge's outward normal
        % does not exceed the apothem. Edge normals are placed at
        % multiples of 360/nSides starting from 0 deg (Octagon, n=8,
        % this way comes out as an ordinary "stop sign" octagon).
        keep = true(size(xr));
        for kk = 0:nSides-1
            ang = kk*360/nSides;
            proj = xr*cosd(ang) + yr*sind(ang);
            keep = keep & (proj <= apothem*(1+1e-9));
        end
    end

    function keep = shapeMaskFn(xr,yr,halfW,halfH,shape)
        % (xr,yr): position relative to the array's own centre, in the
        % ORTHOGONAL (unskewed) dx,dy bounding box -- shape is decided
        % before gridAngle skews the surviving cells' physical position,
        % so the boundary shape itself is never distorted by the skew.
        % (halfW,halfH): half-extents of that orthogonal bounding box.
        %
        % Degenerate line array (M=1 or N=1) guard: Circle/Hexagon/
        % Octagon all use min(halfW,halfH) as their radius/apothem --
        % for a 1-row or 1-column array that minimum is exactly 0,
        % which used to reject every element except the exact centre
        % (confirmed: an 8-element line array with Circle selected kept
        % 0 of 8 elements). A true 1D line has no meaningful 2D circular
        % or polygonal boundary to apply, so these shapes fall back to
        % keeping the full line rather than deleting it. Ellipse also
        % uses halfW to cap its height, so it needs the same guard when
        % the grid has only one column. Diamond's per-axis normalization
        % already keeps a full degenerate line.
        if (halfW==0 || halfH==0) && ...
                ismember(shape,{'Circle','Hexagon','Octagon','Ellipse'})
            keep = true(size(xr));
            return;
        end
        switch shape
            case 'Diamond'
                hw = max(halfW,eps); hh = max(halfH,eps);
                keep = (abs(xr)/hw + abs(yr)/hh) <= 1+1e-9;
            case 'Circle'
                % true round shape: uses the SMALLER half-extent, so it
                % always fits fully inside a non-square MxN bounding box
                R = min(halfW,halfH);
                keep = (xr.^2+yr.^2) <= R^2*(1+1e-9);
            case 'Ellipse'
                % A fitted ellipse became identical to Circle on the
                % default square 8x8 lattice: both kept the same 32
                % cells. Widen the horizontal radius slightly so the
                % outer columns survive near the middle (an even grid
                % has no row at y=0), and cap the vertical radius to
                % keep a clearly horizontal oval within the grid.
                hw = max(1.15*halfW,eps);
                hh = max(min(halfH,0.75*halfW),eps);
                keep = ((xr/hw).^2+(yr/hh).^2) <= 1+1e-9;
            case 'Hexagon'
                keep = insideRegularPolygon(xr,yr,6,min(halfW,halfH));
            case 'Octagon'
                keep = insideRegularPolygon(xr,yr,8,min(halfW,halfH));
            otherwise   % 'Custom' -- full rectangular grid, no trimming
                keep = true(size(xr));
        end
        % Outcome-based backstop for the same failure the halfW/halfH==0
        % guard above catches: if the shape would keep NOTHING, keep the
        % full grid instead of silently deleting the user's whole array.
        %
        % That guard only fires when a half-extent is exactly zero (M=1
        % or N=1). It misses the NEARLY degenerate case, which fails just
        % as completely: at M=2 with an even N, every non-Custom shape
        % kept 0 of 16 elements. The two occupied rows sit exactly on the
        % boundary (|yr| == halfH, so the shape's whole budget is spent
        % on y) while an even N puts no column at xr == 0, so no lattice
        % point satisfies the test. Measured across M=1..6, N=8:
        %   M=1  all shapes keep 8   (handled by the guard above)
        %   M=2  Diamond/Hexagon/Octagon/Circle/Ellipse ALL keep ZERO
        %   M=3+ every shape keeps at least some
        % Testing the RESULT rather than enumerating shapes and sizes
        % also covers odd/even N, non-square lattices and any shape added
        % later. assignShape reports it so the dropdown doesn't end up
        % claiming a shape the array does not actually have.
        if ~any(keep(:))
            keep = true(size(xr));
            S.shapeFellBack = true;
        end
    end

    function checkLargeArrayWarning(prefix)
        % One-time, non-blocking notice: computePattern recomputes over
        % the full angle grid on every change and gets noticeably slow
        % well before ~1500 elements. Re-arms whenever the count drops
        % back to/below 1500, so deleting elements down and later adding
        % (or loading) past the threshold again still warns -- shared by
        % every path that can change size(S.el,1): rebuildUniform,
        % addOrSelect, loadConfig, clearAll, deleteSelected, and the
        % Sparse/Sub-position M/N-shrink filtering in onGeom.
        msg = largeArrayNotice(prefix);
        if ~isempty(msg)
            uialert(fig, msg, 'Large array', 'Icon','info');
        end
    end

    function msg = largeArrayNotice(prefix)
        % The decision and the message, separated from the act of
        % showing it, so a caller that ALREADY has something to say can
        % merge the two into a single dialog instead of firing two.
        % uialert does not block, so back-to-back calls overlay or
        % replace each other and the user sees only one -- loadConfig
        % could hit exactly that by loading a version-skewed config
        % (missing-fields warning) that also holds >1500 elements.
        nEl = size(S.el,1);
        msg = '';
        if nEl > 1500 && ~S.largeWarnShown
            S.largeWarnShown = true;
            msg = sprintf(['%s%d elements. computePattern recomputes ' ...
                'over the full angle grid on every change and can get ' ...
                'slow at this size -- consider turning off Auto-recompute.'], ...
                prefix, nEl);
        elseif nEl <= 1500
            S.largeWarnShown = false;
        end
    end

    function rebuildUniform()
        S.sel = []; S.pending = [];
        % This flag describes THIS rebuild only. Without resetting it, one
        % earlier no-cell fallback could remain latched during a later valid
        % rebuild and make the UI/state report the wrong shape.
        S.shapeFellBack = false;
        Ntot = S.M*S.N;
        elAll = zeros(Ntot,5); rcAll = zeros(Ntot,2);
        xoAll = zeros(Ntot,1); yoAll = zeros(Ntot,1);   % orthogonal coords, for shape test only
        % The skew/stagger belong to Uniform grid mode ONLY. Sparse and
        % Sub-position are documented to sit on the plain orthogonal
        % dx,dy lattice (see onGridClick's note below, and the background
        % dots refreshLayout draws outside Uniform mode), and onMode
        % converts the array onto it when leaving Uniform.
        %
        % This function did not know that. It always placed cells with
        % S.gridAngle/S.stagger, which is right for its original caller
        % but wrong for the one loadConfig added: a geometry-only config
        % carrying a skew, loaded while in Sparse mode, rebuilt a fully
        % skewed array underneath an orthogonal reference grid, with
        % every element off its own dot.
        if strcmp(S.mode,'Uniform grid')
            useAng = S.gridAngle; useStag = S.stagger;
        else
            useAng = 90; useStag = 0;
        end
        n = 0;
        for row = 1:S.M
            for col = 1:S.N
                n = n+1;
                xo = (col-1)*S.dx; yo = (row-1)*S.dy;
                [xg,yg] = latticePosAt(row,col,S.dx,S.dy,useAng,useStag);
                elAll(n,:)  = [xg, yg, 1, 0, 0];
                rcAll(n,:)  = [row col];
                xoAll(n) = xo; yoAll(n) = yo;
            end
        end
        halfW = (S.N-1)*S.dx/2; halfH = (S.M-1)*S.dy/2;
        if strcmp(S.mode,'Uniform grid')
            keepMask = shapeMaskFn(xoAll-halfW, yoAll-halfH, halfW, halfH, S.arrayShape);
        else
            % Sparse/Sub-position geometry uses the complete orthogonal grid.
            keepMask = true(Ntot,1);
        end
        S.el   = elAll(keepMask,:);
        S.elRC = rcAll(keepMask,:);
        if S.shapeFellBack && ~strcmp(S.arrayShape,'Custom')
            % shapeMaskFn kept the full bounding grid because the requested
            % boundary contained no lattice point. Keep the state/UI honest:
            % the geometry on screen is rectangular, therefore its shape is
            % now Custom rather than silently continuing to claim Circle/
            % Hexagon/etc. assignShape() still sees shapeFellBack=true and
            % provides its dedicated explanation when this came from a click.
            S.arrayShape = 'Custom';
            if exist('ddShape','var') && isgraphics(ddShape)
                ddShape.Value = 'Custom';
            end
            % The ribbon gallery is the THIRD place this state is shown,
            % and it was missed: syncShapeGallery ran only from startup,
            % assignShape and loadConfigCore, never on the onGeom path.
            % So choosing Circle at 8x8 and then shrinking to 2x2 (where
            % the circle encloses no lattice point) left the dropdown
            % reading Custom while the gallery still highlighted Circle.
            % syncShapeGallery guards its own handles, so this is safe at
            % startup too.
            syncShapeGallery();
        end
        applyTaper();
        % Uses the ACTUAL post-shape-mask element count (checked inside
        % checkLargeArrayWarning via size(S.el,1)), not M*N -- a Circle/
        % Hexagon/etc mask can drop a large fraction of the bounding MxN
        % grid, and warning off the un-trimmed M*N would fire (or stay
        % armed) well before the real element count justifies it.
        checkLargeArrayWarning('');
    end

    function clearAll()
        % One click, no confirmation (a dialog in front of every Clear is
        % the kind users learn to click through); it is an undo step
        % instead, and the status line says how to get the elements back.
        nWas = size(S.el,1); stepsWas = undoCount;
        S.el = zeros(0,5); S.elRC = zeros(0,2); S.pending = []; S.sel = [];
        checkLargeArrayWarning('');   % re-arms the flag (count is now 0)
        refreshAll();
        % Kept with the compute's warning (why the plots are empty):
        % "⚠ Cleared 64 elements — No radiation: … Undo with ⌘Z".
        announceUndoable(stepsWas, sprintf('Cleared %d %s', nWas, ...
            plural(nWas,'element')));
    end

    % --------------------------------------------------- cut-plot probe
    function cutSteer = signedSteerOnCut(thetaSteer,phiSteer,phiCut)
        % A signed elevation cut contains directions at phiCut and phiCut+180.
        if thetaSteer==0
            cutSteer=0;
        elseif abs(sind(phiSteer-phiCut))<=1e-10
            cutSteer=thetaSteer*sign(cosd(phiSteer-phiCut));
        else
            cutSteer=NaN; % steer direction does not lie on this cut
        end
    end

    function onCutClick(~,evt)
        % Click anywhere on the cut plot (theta or phi axis) to read off
        % the exact Total-curve dB level at that point -- nearest-sample
        % lookup against S.cutX/S.cutDb, which computePattern refreshes
        % every time it runs. Click again elsewhere to move the probe;
        % A successful recompute clears it with the rest of the cut plot.
        % Inactive or blocked designs clear the plot via
        % invalidateRadiation() or patternReady(), so a marker never
        % survives beside results from a different design state.
        if isempty(S.cutX), return; end
        xClick = evt.IntersectionPoint(1);
        [~,idx] = min(abs(S.cutX - xClick));
        xAt = S.cutX(idx); dbAt = S.cutDb(idx);
        delete(findall(axCut,'Tag','cutProbe'));
        P = pal();
        hold(axCut,'on');
        plot(axCut, xAt, dbAt, 'o', 'MarkerSize',8, 'Color',P.good, ...
            'MarkerEdgeColor',P.good, 'MarkerFaceColor',P.good, ...
            'Tag','cutProbe','HitTest','off');
        probeAngle=cutPointText(xAt);
        probeLabel=sprintf('%s, %.1f dB',probeAngle,dbAt);
        if ~isempty(S.cutRaw) && numel(S.cutRaw)==numel(S.cutX)
            rawAt=S.cutRaw(idx);
            if isinf(rawAt) && rawAt<0
                probeLabel=sprintf('%s: zero field (marker at display floor)',probeAngle);
            elseif isfinite(rawAt) && rawAt<dbAt-1e-6
                probeLabel=sprintf('%s, %.1f dB (below display floor)',probeAngle,rawAt);
            end
        end
        % Above the marker, or below it near the top of the axes: the
        % band above the beam holds the peak label, and a label above the
        % peak used to run into the title.
        yLim = axCut.YLim;
        if dbAt <= yLim(2) - 0.3*diff(yLim)
            yLbl = dbAt + 1.5; vAl = 'bottom';
        else
            yLbl = dbAt - 1.5; vAl = 'top';
        end
        % Kept clear of the side edges as well.
        xLim = axCut.XLim; hAl = 'center';
        if xAt < xLim(1) + 0.15*diff(xLim), hAl = 'left';
        elseif xAt > xLim(2) - 0.15*diff(xLim), hAl = 'right'; end
        text(axCut, xAt, yLbl, probeLabel, ...
            'Color',P.good, 'FontWeight','bold', 'FontSize',9, ...
            'HorizontalAlignment',hAl, 'VerticalAlignment',vAl, ...
            'Tag','cutProbe','HitTest','off');
        hold(axCut,'off');
    end

    % ------------------------------------------------- layout interaction
    function onGridClick(~,evt)
        p = evt.IntersectionPoint(1:2);
        isSub = strcmp(S.mode,'Sub-position (9-pt)');

        % ---- sub-position mode, second click: pick one of the 9 offsets ----
        if isSub && ~isempty(S.pending)
            row = S.pending(1); col = S.pending(2);
            cx = (col-1)*S.dx;  cy = (row-1)*S.dy;
            ox = S.subOff*[-1 0 1];
            oy = S.subOff*[-1 0 1];
            best = inf; bx = cx; by = cy;
            for a = 1:3
                for b = 1:3
                    px = cx+ox(a); py = cy+oy(b);
                    d = hypot(p(1)-px, p(2)-py);
                    if d < best, best = d; bx = px; by = py; end
                end
            end
            % a click nowhere near any candidate cancels instead of placing
            cancelTol = max(0.7*S.subOff, 0.15*min(S.dx,S.dy));
            if best > cancelTol
                S.pending = [];
                refreshLayout(); refreshPick();
                return;
            end
            addSubPosition(bx,by,row,col);
            S.pending = [];
            refreshAll();
            return;
        end

        % ---- sub-position mode, first click ----
        % Land ON an element and you SELECT it; land anywhere else in a
        % cell and you arm that cell's 9 candidates.
        %
        % This used to arm unconditionally, ahead of the select branch,
        % so an existing element could never be picked by clicking it --
        % the click armed a cell and cleared the selection instead, which
        % made a just-placed element impossible to select or delete. The
        % old comment justified that by saying a click would otherwise
        % "just select and the 9 candidates would never appear"; that is
        % only true for a click landing exactly on an element, which the
        % tight tolerance below is what distinguishes.
        %
        % Tolerance is tied to subOff, not to the lattice pitch as in the
        % other modes: candidates sit only subOff apart, so the 0.3*dx
        % used elsewhere would swallow neighbouring candidates and make
        % "which one did I click" ambiguous.
        if isSub
            row = []; col = [];
            if ~isempty(S.el)
                selTol = min(0.35*S.subOff, 0.3*min(S.dx,S.dy));
                [dmin,nHit] = min(hypot(S.el(:,1)-p(1), S.el(:,2)-p(2)));
                if dmin < selTol
                    % BOTH, not either. Selecting and arming are not in
                    % conflict, and making one click choose between them
                    % is what went wrong twice: arming only left a
                    % just-placed element impossible to select or
                    % delete, and selecting only meant clicking a cell
                    % centre never opened its candidates.
                    S.sel = toggleSel(S.sel, nHit, isAdditiveClick());
                    syncTableSelection();
                    % Arm the element's OWN cell, from its stored
                    % identity. Deriving it from the click position
                    % would round an offset element into a NEIGHBOURING
                    % cell -- at the default pitch an element at
                    % (0.25,0.25) rounds to cell (2,2), not its own.
                    row = S.elRC(nHit,1); col = S.elRC(nHit,2);
                end
            end
            if isempty(row)
                col = round(p(1)/S.dx)+1;  row = round(p(2)/S.dy)+1;
            end
            if col<1||col>S.N||row<1||row>S.M, return; end
            S.pending = [row col];
            refreshLayout(); refreshPick();
            return;
        end

        % ---- other modes: clicking an existing element selects it ----
        % Ctrl (Windows/Linux) or Cmd (Mac) held during the click ADDS to
        % the current selection instead of replacing it.
        if ~isempty(S.el)
            [dmin,nHit] = min(hypot(S.el(:,1)-p(1), S.el(:,2)-p(2)));
            if dmin < 0.3*min(S.dx,S.dy)
                S.sel = toggleSel(S.sel, nHit, isAdditiveClick());
                syncTableSelection();
                refreshLayout(); refreshInfo();
                return;
            end
        end

        % ---- otherwise: add at the nearest lattice cell ----
        % Uniform grid mode with a skewed (Grid angle != 90) or
        % staggered lattice can't invert click position -> (row,col) by
        % simple division (that only works for a plain orthogonal grid)
        % -- so for that case, find the nearest MISSING cell by its
        % actual computed position instead (brute-force over the full
        % MxN grid, excluding already-occupied cells). This is the same
        % class of bug as resetSelectedToLattice's fix above: any place
        % that derives element position from (row,col) needs to go
        % through latticePosAt once gridAngle/stagger can be nonzero.
        % Sparse mode intentionally stays on the plain orthogonal
        % lattice (documented scoping decision), so keeps the simple
        % division-based lookup.
        if strcmp(S.mode,'Uniform grid')
            [GR,GC] = ndgrid(1:S.M,1:S.N);
            [GX,GY] = latticePosAt(GR,GC,S.dx,S.dy,S.gridAngle,S.stagger);
            occupied = false(S.M,S.N);
            if ~isempty(S.elRC)
                occupied(sub2ind([S.M,S.N],S.elRC(:,1),S.elRC(:,2))) = true;
            end
            distAll = hypot(GX-p(1), GY-p(2));
            distAll(occupied) = Inf;
            [minD,idxMin] = min(distAll(:));
            if isinf(minD) || minD > 0.6*min(S.dx,S.dy), return; end
            [row,col] = ind2sub([S.M,S.N],idxMin);
            addOrSelect(GX(idxMin),GY(idxMin),row,col);
        else
            col = round(p(1)/S.dx)+1;  row = round(p(2)/S.dy)+1;
            if col<1||col>S.N||row<1||row>S.M, return; end
            addOrSelect((col-1)*S.dx,(row-1)*S.dy,row,col);
        end
        refreshAll();
    end

    function addSubPosition(x,y,row,col)
        %ADDSUBPOSITION  Place an element at a chosen 9-point offset.
        %   Deliberately NOT addOrSelect. That one refuses to add when
        %   the (row,col) cell identity is already claimed, which made
        %   this whole mode inert: coming from a uniform grid every cell
        %   is occupied, so the nine candidates opened, the click landed,
        %   and nothing happened at all.
        %
        %   Here a cell may hold SEVERAL elements, one per offset, so the
        %   array gets denser as you click rather than just rearranged.
        %   Elements sharing a cell also share its (row,col), and that is
        %   fine for everything that consumes it: applyTaper assigns
        %   wy(row)*wx(col) per element, so co-cell elements get the same
        %   aperture weight -- correct, they sit at nearly the same
        %   place -- and applySeqRot derives its rotation from
        %   mod(row-1,blockM)/mod(col-1,blockN) per element, with no sort
        %   key to collide. The dx/dy rescale path preserves each
        %   element's own offset from its cell, so they stay distinct
        %   when the lattice pitch changes.
        %
        %   The one thing worth blocking is stacking duplicates at the
        %   SAME offset: clicking one candidate twice should select what
        %   is there, not pile a second element on the same point.
        tol = min(0.01*min(S.dx,S.dy), 0.2*S.subOff);
        hit = [];
        if ~isempty(S.el)
            hit = find(hypot(S.el(:,1)-x, S.el(:,2)-y) < tol, 1);
        end
        if isempty(hit)
            S.el(end+1,:)   = [x y 1 0 0];
            S.elRC(end+1,:) = [row col];
            applyTaper();
            S.sel = toggleSel(S.sel, size(S.el,1), isAdditiveClick());
            checkLargeArrayWarning('');
        else
            S.sel = toggleSel(S.sel, hit, isAdditiveClick());
        end
        syncTableSelection();
    end

    function addOrSelect(x,y,row,col)
        if isempty(S.el)
            hit = [];
        else
            % capped by subOff too: with a small sub-offset, adjacent
            % 9-point candidates can sit closer together than 1% of the
            % lattice spacing, which would otherwise misidentify them as
            % duplicates of each other and silently block placement
            tol = min(0.01*min(S.dx,S.dy), 0.2*S.subOff);
            hit = find(hypot(S.el(:,1)-x, S.el(:,2)-y) < tol, 1);
            if isempty(hit)
                % No element physically near this click point -- but the
                % (row,col) CELL IDENTITY might already be claimed by an
                % element that was since moved/randomized elsewhere
                % (moving an element updates its x,y but deliberately
                % preserves its elRC identity, e.g. for taper/sequential-
                % rotation assignment). Without this check, clicking the
                % now-vacant original spot would add a SECOND element
                % sharing the same (row,col) identity as the displaced
                % one -- corrupting tapering (both get the same weight)
                % and sequential-rotation ordering (duplicate sort key).
                hit = find(S.elRC(:,1)==row & S.elRC(:,2)==col, 1);
            end
        end
        additive = isAdditiveClick();
        if isempty(hit)
            S.el(end+1,:)  = [x y 1 0 0];
            S.elRC(end+1,:) = [row col];
            applyTaper();
            S.sel = toggleSel(S.sel, size(S.el,1), additive);  % select the new element
            % Sparse/Sub-position build up their element count one click
            % at a time instead of one rebuildUniform() call, so it needs
            % its own check here too.
            checkLargeArrayWarning('');
        else
            S.sel = toggleSel(S.sel, hit, additive);           % existing element
        end
        syncTableSelection();
    end

    function syncTableSelection()
        % Wait, rather than fail loudly, when the table has not caught up
        % with S.el yet: on an add the caller sets S.sel to the new
        % element's index before refreshTable has grown tbl.Data.
        % refreshTable re-applies the selection immediately after it sets
        % Data, so nothing is lost by skipping here. The catch below is
        % kept for genuinely unexpected failures, which this no longer
        % drowns out.
        if ~isempty(S.sel) && max(S.sel) > size(tbl.Data,1)
            refreshPick(); return;
        end
        try
            tbl.Selection = S.sel(:).';   % must be 1-by-N under row-selection mode
            if ~isempty(S.sel), scroll(tbl,'row',S.sel(1)); end
        catch err
            fprintf(2,'[table sync] %s\n', err.message);
        end
        refreshPick();
    end

    function refreshPick()
        if ~isempty(S.pending)
            pLay.Tooltip = sprintf(['Cell (row %d, col %d) armed - click one of ' ...
                'the 9 orange candidates to place an element there (Esc cancels).'], ...
                S.pending(1), S.pending(2));
            return;
        end
        if strcmp(S.mode,'Sub-position (9-pt)') && isempty(S.sel)
            pLay.Tooltip = 'Sub-position mode: click a lattice cell to open its 9 offsets.';
            return;
        end
        if isempty(S.sel) || isempty(S.el)
            pLay.Tooltip = 'Click an empty cell to add an element; click an element to select it; press Delete to remove it.';
            return;
        end
        n = S.sel(1);
        if n > size(S.el,1)
            pLay.Tooltip = 'Click an empty cell to add an element; click an element to select it; press Delete to remove it.';
            return;
        end
        extra = '';
        if numel(S.sel) > 1
            extra = sprintf('   |   (+%d more selected, Ctrl/Cmd-click to add)', numel(S.sel)-1);
        end
        pLay.Tooltip = sprintf(['Element #%d  (table row %d)   |   ' ...
            'x = %.3f λ,  y = %.3f λ   |   Amp = %.3f   |   ' ...
            'Phase = %.1f°   |   Rot = %.1f°%s'], ...
            n, n, S.el(n,1), S.el(n,2), S.el(n,3), S.el(n,4), S.el(n,5), extra);
    end

    function onKey(keySrc,evt)
        switch lower(evt.Key)
            case {'delete','backspace'}
                deleteSelected();
            case 'escape'
                if ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup) && ...
                        strcmp(shapeGalleryPopup.Visible,'on')
                    shapeGalleryPopup.Visible = 'off';
                    return;
                end
                if ~isempty(elementGalleryPopup) && isgraphics(elementGalleryPopup) && ...
                        strcmp(elementGalleryPopup.Visible,'on')
                    elementGalleryPopup.Visible = 'off';
                    return;
                end
                S.sel = []; S.pending = [];
                refreshAll();
            case {'z','y'}
                % Cmd+Z undo, Shift+Cmd+Z and Cmd+Y redo (Ctrl elsewhere).
                % Shift+Cmd+Z exists only here, since a menu Accelerator
                % cannot hold Shift; the other two are the Edit menu's
                % Accelerators as well, and undoFromUI lets one keystroke
                % that arrives both ways act once.
                keyMods = evt.Modifier;
                if ~any(ismember(keyMods,{'command','control'})), return; end
                if strcmpi(evt.Key,'z') && ~any(strcmp(keyMods,'shift'))
                    keyAct = 'undo';
                else
                    keyAct = 'redo';
                end
                keyFrom = 'window';
                if isequal(keySrc,tbl), keyFrom = 'table'; end
                undoFromUI(keyAct,keyFrom);
        end
    end

    function deleteSelected()
        if isempty(S.sel), return; end
        stepsWas = undoCount;
        keep = true(size(S.el,1),1);
        keep(S.sel(S.sel<=numel(keep))) = false;
        S.el = S.el(keep,:); S.elRC = S.elRC(keep,:);
        S.sel = []; S.pending = [];
        checkLargeArrayWarning('');
        refreshAll();
        announceUndoable(stepsWas, sprintf('Deleted %d %s', sum(~keep), ...
            plural(sum(~keep),'element')));
    end

    % ------------------------------------------------------------ tapering
    function applyTaper()
        if isempty(S.el) || strcmp(S.taper,'Manual (table)'), return; end
        S.taperFailed = false;
        wx = taperVec(S.N); wy = taperVec(S.M);
        for n = 1:size(S.el,1)
            rr = min(max(S.elRC(n,1),1),S.M);
            cc = min(max(S.elRC(n,2),1),S.N);
            S.el(n,3) = wy(rr)*wx(cc);
        end
        % Say so, once, rather than letting the dropdown keep claiming a
        % taper the array does not actually have. Reverts to Uniform so
        % the control and the amplitudes agree afterwards -- leaving it
        % on "Chebyshev" would re-trigger this on every geometry change.
        if S.taperFailed
            failed = S.taper;
            S.taper = 'Uniform'; ddTap.Value = 'Uniform'; syncSLLControl();
            wx = taperVec(S.N); wy = taperVec(S.M);
            for n = 1:size(S.el,1)
                rr = min(max(S.elRC(n,1),1),S.M);
                cc = min(max(S.elRC(n,2),1),S.N);
                S.el(n,3) = wy(rr)*wx(cc);
            end
            uialert(fig,sprintf(['The %s taper could not produce usable weights for this array. ' ...
                'Its window function may be unavailable, or this array dimension may be too small. ' ...
                'The taper and all amplitudes have been set to Uniform.'],failed), ...
                'Taper unavailable','Icon','warning');
        end
    end

    function w = taperVec(L)
        if L < 2, w = 1; return; end
        try
            switch S.taper
                case 'Uniform',   w = ones(L,1);
                % Computed from their closed forms rather than by calling
                % the Signal Processing Toolbox's hamming()/hann().
                % Those calls sit inside this function's try/catch, so on
                % a machine WITHOUT that toolbox they threw and silently
                % fell back to Uniform -- the dropdown would still read
                % "Hamming" while the array was actually untapered,
                % turning a ~-43 dB design sidelobe into ~-13 dB with no
                % indication anything was wrong. That is a silent wrong
                % answer on someone else's machine, so the dependency is
                % removed outright. Verified identical to MATLAB's
                % symmetric hamming/hann to 5.6e-16 across L = 2..64.
                case 'Hamming',   w = 0.54 - 0.46*cos(2*pi*(0:L-1)'/(L-1));
                case 'Hanning',   w = 0.5  - 0.5 *cos(2*pi*(0:L-1)'/(L-1));
                case 'Chebyshev', w = chebwin(L,S.sll);
                case 'Taylor'
                    % nbar DERIVED from the sidelobe target, not fixed at
                    % 2. nbar is the number of near-in sidelobes held at
                    % the design level; the standard requirement is
                    % nbar >= 2*A^2 + 0.5 with A = acosh(10^(SLL/20))/pi.
                    % Below that the design simply cannot reach the
                    % requested level. With nbar hardcoded to 2 the taper
                    % saturated near -28 dB however high the Sidelobe
                    % level spinner was set -- measured on L=8:
                    %   target -20 -> -19.9 (ok, nbar 2 is enough there)
                    %   target -30 -> -24.1  (5.9 dB short)
                    %   target -40 -> -26.6 (13.4 dB short)
                    %   target -50 -> -28.2 (21.8 dB short)
                    % and with the derived nbar: -28.3, -38.3, -51.4.
                    % Chebyshev hits its target exactly from the same
                    % spinner, so there was nothing on screen to suggest
                    % Taylor was quietly ignoring most of the request.
                    % The derivation is taylorNbar(), shared with the
                    % n̄ note beside the Sidelobe level spinner.
                    w = taylorwin(L,taylorNbar(S.sll),-S.sll);
                case 'Binomial'
                    % coefficients of (1+x)^(L-1), i.e. row L of Pascal's
                    % triangle: w(1)=1, w(k+1)=w(k)*(L-k)/k. Gives the
                    % narrowest-sidelobe / widest-beamwidth extreme of
                    % the taper family -- for spacing <= 0.5 lambda the
                    % resulting pattern has NO sidelobes at all (the
                    % array factor collapses to cos^(L-1)(psi/2), a
                    % single-lobed function with no ripple).
                    w = zeros(L,1); w(1) = 1;
                    for k = 1:L-1
                        w(k+1) = w(k)*(L-k)/k;
                    end
                otherwise,        w = ones(L,1);
            end
        catch
            % Toolbox window unavailable (Chebyshev/Taylor still need the
            % Signal Processing Toolbox; Hamming/Hanning/Binomial were
            % reimplemented above precisely to avoid this). Falling back
            % to Uniform silently is the same class of silent wrong
            % answer documented above: the dropdown would keep reading
            % "Chebyshev" while the array ran untapered, turning a -30 dB
            % design into ~-13 dB with nothing on screen to say so. Flag
            % it so applyTaper can tell the user once.
            S.taperFailed = true;
            w = ones(L,1);
        end
        % Guarded normalization. hann(2) is exactly [0;0] -- with L=2
        % both samples land on the window's own zero endpoints -- so an
        % unguarded w/max(w) evaluates 0/0 = NaN. That NaN then went
        % straight into every element amplitude, made
        % ref = sum(abs(amp)) NaN, and surfaced as a misleading "Every
        % element amplitude is zero" alert with the pattern refusing to
        % compute at all. Reachable on any array with M=2 or N=2, which
        % is exactly the shape of a 2x2 sequential-rotation CP tile.
        % Any degenerate window falls back to Uniform rather than
        % poisoning S.el with NaN.
        mx = max(w);
        if any(~isfinite(w)) || ~isfinite(mx) || mx <= 0
            S.taperFailed = true;
            w = ones(L,1);
        else
            w = w(:)/mx;
        end
    end

    % ------------------------------------------------------ table handling
    function onTableEdit(~,evt)
        n = evt.Indices(1); col = evt.Indices(2);
        v = evt.NewData;
        % Column 6 (Feed phase) is ColumnEditable=false, so this
        % shouldn't fire for it -- guarded anyway since S.el only has 5
        % columns and a col==6 write would otherwise error.
        if col > 5, refreshTable(); return; end
        % isreal() added: MATLAB's numeric table-cell parsing accepts
        % complex literals like "3+4i" (confirmed: str2double('3+2i')
        % returns a complex double), and isnumeric/isfinite are both
        % true for complex values -- without this check a complex entry
        % would silently corrupt downstream real-valued math (taper
        % efficiency, dB conversions) rather than error or get rejected.
        if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || ~isreal(v)
            refreshTable(); return;
        end
        if col <= 2 && abs(v) > MAX_POS_LAMBDA
            uialert(fig,sprintf(['Element position must be within +/-%g wavelengths. ' ...
                'The entered coordinate cannot be plotted or phased reliably.'], ...
                MAX_POS_LAMBDA),'Position out of range');
            refreshTable(); return;
        end
        % Rotation wrapped to [0,360) like every other path that writes
        % it (setRot, stepRot, applySeqRot all use mod). Typing this
        % column was the one route that skipped the wrap. The pattern
        % maths is periodic either way, but mixedRot is computed as
        % numel(unique(S.el(:,5)))>1 -- so a typed 360 sitting alongside
        % a 0 counted as two distinct rotations and made the plot title
        % claim "mixed rot: Total ~= EF x AF" for an array whose
        % elements are all physically identical.
        if col == 4, v = wrapManualPhase(v); end
        if col == 5, v = mod(v,360); end
        % AMPLITUDE bound. "Finite and real" was not enough: the pattern
        % math squares and sums these, so a finite input can still produce
        % a non-finite intermediate. An amplitude of 1e200 passed the test
        % above, then overflowed sum(w.^2) to Inf (taper efficiency
        % Inf/Inf = NaN) and the radiated-power integral to Inf, whose
        % absolute offset is -Inf -- and the 3D colour limits then threw
        % outright, killing the recompute.
        %
        % These are RELATIVE excitation weights, so their absolute scale
        % carries no information; 1e6 of dynamic range is far past any
        % real taper and leaves every squared sum comfortably finite.
        if col == 3
            AMP_MAX = 1e6;
            if abs(v) > AMP_MAX
                uialert(fig, sprintf(['Amplitude %.3g is outside the supported range ' ...
                    '(+/-%g).' newline newline 'These are RELATIVE weights -- only their ' ...
                    'ratios matter -- and values this large overflow the power sums ' ...
                    'the pattern is built from. Scale your whole taper down instead; ' ...
                    'the pattern is unchanged by a uniform scale.'], v, AMP_MAX), ...
                    'Amplitude out of range');
                refreshTable(); return;
            end
        end
        S.el(n,col) = v;
        if col == 3
            S.taper = 'Manual (table)'; ddTap.Value = 'Manual (table)';
            syncSLLControl();
        end
        % full refresh: editing Amp / Phase / Rot changes the pattern, the
        % taper efficiency and the layout, not just the drawing
        refreshAll();
    end

    function sgn = seqRotSign()
        %SEQROTSIGN  Which way the sequential-rotation feed phase must go.
        %
        %   Rotating a circularly polarised element by alpha multiplies
        %   its field by exp(-j*alpha) for one hand and exp(+j*alpha) for
        %   the other. The compensating feed phase has to cancel that,
        %   so its sign is set by the ELEMENT's handedness -- it is not a
        %   free choice.
        %
        %   This used to be hard-wired to +rot. That is right for an
        %   RHCP-dominant element and exactly wrong for an LHCP one,
        %   where it DOUBLES the rotation phase instead of cancelling it:
        %   across a 0/90/180/270 block the four elements then carry
        %   1, -1, 1, -1 and cancel in pairs. Measured on two real CST
        %   exports of the same kind of dual-fed patch, differing only in
        %   which port got the +90 deg:
        %     RHCP-dominant file   21.16 -> 22.69 dBi with +rot  (correct)
        %     LHCP-dominant file   21.16 ->  2.36 dBi with +rot,
        %                                   21.11 dBi with -rot
        %   Flipping the sequential-rotation STEP sign does not help,
        %   because the feed phase follows the rotation and both flip
        %   together.
        %
        %   Handedness is read off the element's own pattern on a small
        %   theta ring, at the reference orientation. Only a CLEAR
        %   imbalance flips the sign: a linear element has |E_R| = |E_L|
        %   and either sign works (they just produce opposite hands of
        %   CP at the array level), so the historical +rot is kept there
        %   and every existing design is unaffected.
        sgn = 1;
        if isempty(S.el), return; end
        % The circular-component formulas below, E_R/E_L = (Eth -+ j*Eph),
        % are only meaningful when Eth/Eph ARE the theta/phi spherical
        % pair. A Ludwig-3 Copol/Cross export (or any other orthogonal
        % basis) substitutes into them perfectly happily and returns a
        % confident handedness that means nothing -- and this one feeds
        % the rotation compensation sign, so a wrong answer here flips
        % every feed phase in the array.
        %
        % Disabling the RHCP/LHCP readouts elsewhere does not protect
        % this: that gate is about what gets DISPLAYED, and this runs
        % regardless. On an unsupported basis, decline to guess and keep
        % the historical +rot, which is at least a documented default the
        % user can override by hand.
        % Requires POSITIVE evidence of a theta/phi pair, not merely the
        % absence of a recognised wrong one. impBasisOK() returns true
        % when the component names are empty -- it cannot tell an
        % unnamed theta/phi export from an unnamed Ludwig-3 one -- and
        % that is not a basis to feed E_R/E_L formulas whose answer sets
        % every feed phase in the array.
        if strcmp(S.efType,'Imported (CST far-field)')
            if isempty(S.impFF) || ~isfield(S.impFF,'compThetaPhi') || ...
                    ~S.impFF.compThetaPhi || S.impFF.noComponents
                return;   % keep the documented +rot default; do not guess
            end
        end
        if ~impBasisOK()
            return;
        end
        thRing = 10; phRing = (0:10:350)';
        try
            [eR_th, eR_ph] = elementFactor(thRing*ones(size(phRing)), phRing, 0);
        catch
            sgn=NaN; return;   % never export a sign from a failed formula
        end
        % Handedness is scale-independent; normalize before squaring.
        ringScale=max([abs(eR_th(:));abs(eR_ph(:))]);
        if ringScale==0, return; end
        eR_th=eR_th/ringScale; eR_ph=eR_ph/ringScale;
        pR = sum(abs((eR_th + 1j*eR_ph)/sqrt(2)).^2);
        pL = sum(abs((eR_th - 1j*eR_ph)/sqrt(2)).^2);
        HAND_TOL = 1.0;   % dB of dominance before the sign is allowed to flip
        % Compared as a product, never as a ratio in dB. The first form
        % was  pR > 0 && pL > 0 && 10*log10(pL/pR) > HAND_TOL,  whose
        % pR > 0 guard existed only to keep log10 away from a divide by
        % zero -- and which therefore refused to fire in exactly the
        % cleanest case it was written for: a PERFECT LHCP element has
        % pR identically 0, so the guard failed and the element fell
        % back to the +rot default that destroys it. Multiplying instead
        % of dividing needs no guard, and handles pure LHCP (pR = 0),
        % pure RHCP (pL = 0) and a dead field (both 0) correctly.
        if isfinite(pR) && isfinite(pL) && pL > 10^(HAND_TOL/10)*pR
            sgn = -1;
        end
    end

    function phTot = effectivePhaseDeg()
        % The ACTUAL total phase applied to each element's feed: manual
        % offset (S.el(:,4)) PLUS the steering-induced phase from the
        % current theta_s/phi_s -- i.e. exactly the phFeed term
        % computePattern uses internally, just exposed here for
        % display/export. This is the fix for "changing steering angle
        % doesn't change the phase column": the table's editable Phase
        % column was always ONLY the manual offset by design (so typing
        % a value there doesn't get clobbered by steering), but nothing
        % previously showed the RESULTING total -- so steering appeared
        % to do nothing to phase even though it was correctly being
        % applied internally the whole time. Returns degrees, wrapped to
        % [-180,180) -- e.g. an exact +180 deg phase wraps to -180, not
        % +180. Rotation's contribution here is gated by S.seqPhase (the
        % "rot angle -> feed phase" checkbox) -- OFF by default, so
        % rotation is purely element orientation unless you opt in.
        if isempty(S.el), phTot = zeros(0,1); return; end
        x = S.el(:,1); y = S.el(:,2); ph0 = S.el(:,4); rot = S.el(:,5);
        us = sind(S.theta_s)*cosd(S.phi_s);
        vs = sind(S.theta_s)*sind(S.phi_s);
        phFeed = -S.k*phaseFreqRatio()*(x*us + y*vs) + deg2rad(ph0);
        if S.seqPhase, phFeed = phFeed + seqRotSign()*deg2rad(rot); end
        phTot = mod(rad2deg(phFeed)+180,360)-180;
    end

    function onTableSelect(~,evt)
        if isempty(evt.Indices), S.sel = []; else, S.sel = unique(evt.Indices(:,1)); end
        refreshLayout();      % highlight the matching element in the layout
        refreshPick(); refreshInfo();
    end

    function tf = isAdditiveClick()
        % true while Ctrl (Win/Linux) or Cmd (Mac) is held during a click
        try
            mods = fig.CurrentModifier;
        catch
            mods = {};
        end
        tf = any(ismember(mods, {'control','command'}));
    end

    function sel = toggleSel(sel, n, additive)
        if additive
            if ismember(n, sel)
                sel = sel(sel ~= n);   % clicking a selected element deselects it
            else
                sel = [sel(:); n];
            end
        else
            sel = n;
        end
    end

    function [exportWeights,exportPhases] = elementExportValues()
        exportWeights=S.el(:,3);
        exportPhases=effectivePhaseDeg();
        if any(~isfinite(S.el(:))) || any(~isfinite(exportPhases(:)))
            error('PAD:ElementExport','Element data or feed phases are invalid. Correct the field formula or element values before exporting.');
        end
    end

    function exportCSV()
        if ~patternReady(true), return; end
        if isempty(S.el)
            uialert(fig,'No elements to export.','Nothing to export');
            return;
        end
        [f,p] = uiputfile('array_elements.csv','Export elements to CSV');
        figure(fig);   % see loadImportedFF's comment -- native file dialogs
                       % can leave the uifigure behind other windows on macOS
        if isequal(f,0), return; end
        % Feed_phase_deg is the TOTAL applied phase (steering + manual
        % offset, plus rotation when "rot angle -> feed phase" is on) --
        % the sixth column of the on-screen table, and the same quantity
        % the TSV writes as its "Phase". Without it the CSV carried only
        % Phase_deg, the MANUAL offset, which is zero on any array that
        % was steered purely with theta_s/phi_s: the file therefore
        % described an UNSTEERED array while the app showed a steered
        % one, under a column name that read like the phase. Both columns
        % are written now, so the CSV matches what the table displays and
        % the two exporters agree about what "phase" means.
        try
            [~,exportPhases]=elementExportValues();
            T = array2table([S.el, exportPhases, S.elRC], 'VariableNames', ...
            {'x_lambda','y_lambda','Amp','Phase_deg','Rot_deg','Feed_phase_deg', ...
             'lattice_row','lattice_col'});
            % Keep exported feed phases self-describing. In retuned mode the
            % phase reference is the operating frequency; in Beam-squint mode
            % it is the design frequency even though the pattern is evaluated
            % at the operating frequency.
            T.Phase_reference_GHz = repmat(phaseReferenceGHz(),height(T),1);
            T.Design_frequency_GHz = repmat(S.freqGHz,height(T),1);
            T.Operating_frequency_GHz = repmat(S.freqOpGHz,height(T),1);
            T.Beam_squint = repmat(~S.retunePhase,height(T),1);
            writetable(T, fullfile(p,f));
            setStatus(sprintf('Exported the element table to %s: %d elements', ...
                f,size(S.el,1)),'good');
        catch err
            uialert(fig, err.message, 'Export failed');
            setStatus(['CSV export failed: ' err.message],'bad');
        end
    end

    % ---- config-value validators, used only by loadConfig ----
    % Each takes the target widget and the raw value read from the .mat
    % file, and returns something the widget will definitely accept.
    % They deliberately fall back to the widget's CURRENT value rather
    % than to a hardcoded default, so a bad field leaves that one
    % control untouched instead of resetting it to something arbitrary.
    function validateConfigCandidate(cfgCandidate)
        if ~isstruct(cfgCandidate) || ~isscalar(cfgCandidate), error('PAD:Config','Expected one configuration structure.'); end
        numberFields = {'M','N','dx','dy','subOff','theta_s','phi_s','sll','efQ','efBeamAz','efBeamEl', ...
            'seqBlockM','seqBlockN','gridAngle','stagger','freqGHz','freqOpGHz', ...
            'cutFixedTheta','cutPhi','impTotEffPct','vDynRange'};
        boolFields = {'seqPhase','dualFeedCP','squintMode','retunePhase','fullSphere', ...
            'cutPhiFollow','impUnitCell','portMapSet','impNeedsPattern', ...
            'vAbsLevel','vARCut','vAuto'};
        textFields = {'mode','taper','efType','customFormula','customFormulaPh', ...
            'arrayShape','freqUnit','cutMode','angleConvention', ...
            'impFFName','impFFFreqFrom','vShow','vPol','vScale'};
        for cf = numberFields
            if isfield(cfgCandidate,cf{1})
                v=cfgCandidate.(cf{1});
                if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v))
                    error('PAD:Config','%s must be a real finite numeric scalar.',cf{1});
                end
            end
        end
        for cf = boolFields
            if isfield(cfgCandidate,cf{1}) && ~validBoolean(cfgCandidate.(cf{1}))
                error('PAD:Config','%s must be a scalar boolean (0 or 1).',cf{1});
            end
        end
        for cf = textFields
            if isfield(cfgCandidate,cf{1}) && ~validText(cfgCandidate.(cf{1}))
                error('PAD:Config','%s must contain one line of text.',cf{1});
            end
        end
        if isfield(cfgCandidate,'freqUnit') && ~any(strcmp(cfgCandidate.freqUnit,{'GHz','MHz'}))
            error('PAD:Config','Saved display frequency unit must be GHz or MHz.');
        end
        % NaN is a legal value here (the file stated no frequency), so
        % this one cannot sit in numberFields with the finite ones.
        if isfield(cfgCandidate,'impFFGHz') && ~validImportGHz(cfgCandidate.impFFGHz)
            error('PAD:Config','impFFGHz must be a positive frequency in GHz, or NaN.');
        end
        for cf={'el','elRC'}
            if isfield(cfgCandidate,cf{1}) && ~(isnumeric(cfgCandidate.(cf{1})) && ismatrix(cfgCandidate.(cf{1})))
                error('PAD:Config','%s must be a two-dimensional numeric matrix.',cf{1});
            end
        end
        if isfield(cfgCandidate,'targets') && isempty(sanitizeTargets(cfgCandidate.targets))
            error('PAD:Config',['targets must be one structure of numeric ' ...
                'scalars (NaN = no target) and a boolean noGrating.']);
        end
        if isfield(cfgCandidate,'impFF') && ~isempty(cfgCandidate.impFF) && ~validSavedFF(cfgCandidate.impFF)
            error('PAD:Config',['Saved imported pattern is invalid or predates the strict realized-gain schema. ' ...
                'Open the older app/config and re-export settings without the pattern, or re-import the original ' ...
                'realized-gain file in the new app. The current design has not been changed.']);
        end
    end
    function candidate=normalizeConfigNumbers(candidate)
        % Integer MAT arrays must not quantize subsequent fractional edits.
        numericNames={'M','N','dx','dy','subOff','theta_s','phi_s','sll','efQ','efBeamAz','efBeamEl', ...
            'seqBlockM','seqBlockN','gridAngle','stagger','freqGHz','freqOpGHz', ...
            'cutFixedTheta','cutPhi','impTotEffPct','vDynRange','el','elRC'};
        for numericName=numericNames
            if isfield(candidate,numericName{1}) && isnumeric(candidate.(numericName{1}))
                candidate.(numericName{1})=double(candidate.(numericName{1}));
            end
        end
        if isfield(candidate,'el') && isnumeric(candidate.el) && ...
                ismatrix(candidate.el) && size(candidate.el,2)==5 && isreal(candidate.el)
            % A phase offset is periodic. Reducing it before the first
            % degrees-to-radians conversion avoids large-angle argument
            % reduction error in every field and export path.
            candidate.el(:,4)=wrapManualPhase(candidate.el(:,4));
            candidate.el(:,5)=mod(candidate.el(:,5),360);
        end
    end

    function phase = wrapManualPhase(phase)
        % Reduce first, then choose the familiar signed representation.
        % mod(phase+180,360)-180 loses the 180-degree offset when phase is
        % large enough that adding 180 rounds away in double precision.
        phase=mod(phase,360);
        phase(phase>=180)=phase(phase>=180)-360;
    end

    function tf = validBoolean(v)
        tf=(islogical(v)||isnumeric(v)) && isscalar(v) && isreal(v) && ...
            isfinite(v) && (v==0 || v==1);
    end
    function tf = validText(v)
        tf=(ischar(v) && (isrow(v)||isempty(v))) || (isstring(v)&&isscalar(v)&&~ismissing(v));
        if tf, tf=~any(ismember(char(v),[newline sprintf('\r')])); end
    end
    function tf = validSavedFF(q)
        tf=false;
        required={'FReEth','FImEth','FReEph','FImEph','peak','peakMag','eff', ...
            'compThetaPhi','noComponents','compNames','thMin','thMax', ...
            'phiFull','thetaFull','quantity','schemaVersion'};
        if ~isstruct(q)||~isscalar(q)||~all(isfield(q,required)), return; end
        if ~validText(q.quantity)||~strcmp(q.quantity,'realized_gain') || ...
                ~isnumeric(q.schemaVersion)||~isscalar(q.schemaVersion)||q.schemaVersion~=2
            return;
        end
        for vf={'compThetaPhi','noComponents','phiFull','thetaFull'}
            if ~validBoolean(q.(vf{1})), return; end
        end
        if ~q.compThetaPhi||q.noComponents||~q.phiFull||~q.thetaFull, return; end
        if ~iscell(q.compNames)||numel(q.compNames)~=2|| ...
                ~all(cellfun(@validText,q.compNames))|| ...
                ~strcmpi(q.compNames{1},'theta')||~strcmpi(q.compNames{2},'phi')
            return;
        end
        for vf={'peakMag','thMin','thMax'}
            vv=q.(vf{1});
            if ~(isnumeric(vv)&&isscalar(vv)&&isreal(vv)&&isfinite(vv)), return; end
        end
        if q.peakMag<=0||q.thMin~=0||q.thMax~=180, return; end
        if ~(isnumeric(q.peak)&&isreal(q.peak)&&isvector(q.peak)&&numel(q.peak)==3&&all(isfinite(q.peak))), return; end
        if ~(isnumeric(q.eff)&&isreal(q.eff)&&isscalar(q.eff)&& ...
                (isnan(q.eff)||(isfinite(q.eff)&&q.eff>0&&q.eff<=1))), return; end
        if q.peak(2)<0||q.peak(2)>180||q.peak(3)<0||q.peak(3)>=360, return; end
        if abs(20*log10(q.peakMag)-q.peak(1))>1e-6, return; end
        try
            obj0=q.FReEth;
            if ~isa(obj0,'scatteredInterpolant')||size(obj0.Points,2)~=2, return; end
            basePts=unique(round(obj0.Points(obj0.Points(:,2)>=0 & obj0.Points(:,2)<360,:),6),'rows');
            ts=unique(basePts(:,1)); ps=unique(basePts(:,2));
            if numel(ts)<3||numel(ps)<8||ts(1)~=0||ts(end)~=180, return; end
            gs=diff([ps;ps(1)+360]);
            if max(gs)>45+1e-6||max(abs(gs-median(gs)))>1e-5|| ...
                    max(abs(diff(ts)-median(diff(ts))))>1e-5|| ...
                    size(basePts,1)~=numel(ts)*numel(ps)||size(obj0.Points,1)~=3*size(basePts,1)
                return;
            end
            expectedPts=sortrows([basePts+[0 -360];basePts;basePts+[0 360]]);
            if ~isequal(sortrows(round(obj0.Points,6)),expectedPts), return; end
            [inBase,baseIndex]=ismember(basePts,round(obj0.Points,6),'rows');
            wrappedPts=round([obj0.Points(:,1),mod(obj0.Points(:,2),360)],6);
            [inPeriod,periodIndex]=ismember(wrappedPts,basePts,'rows');
            if ~all(inBase)||~all(inPeriod), return; end
            for vf={'FReEth','FImEth','FReEph','FImEph'}
                obj=q.(vf{1});
                if ~isa(obj,'scatteredInterpolant') || size(obj.Points,2)~=2 || ...
                        ~isreal(obj.Values)||any(~isfinite(obj.Values)), return; end
                if ~strcmp(obj.ExtrapolationMethod,'none') || ~strcmp(obj.Method,'linear'), return; end
                vals=obj([0;45;90;135;180],[0;90;180;270;359.9]);
                if ~isreal(vals)||any(~isfinite(vals)), return; end
                if ~isequal(obj.Points,q.FReEth.Points), return; end
                expectedVals=obj.Values(baseIndex(periodIndex));
                if any(abs(obj.Values-expectedVals)>1e-10*q.peakMag), return; end
            end
            valuePeak=max(q.FReEth.Values.^2+q.FImEth.Values.^2+q.FReEph.Values.^2+q.FImEph.Values.^2);
            if ~isfinite(valuePeak)||abs(valuePeak/q.peakMag^2-1)>1e-6, return; end
            peakCheck=q.FReEth(q.peak(2),q.peak(3))^2+q.FImEth(q.peak(2),q.peak(3))^2+ ...
                q.FReEph(q.peak(2),q.peak(3))^2+q.FImEph(q.peak(2),q.peak(3))^2;
            if ~isfinite(peakCheck)||abs(peakCheck/q.peakMag^2-1)>0.05, return; end
        catch
            return;
        end
        tf=true;
    end

    function v = pickNum(cfg, fn, dflt)
        %PICKNUM  A real finite scalar from a config field, or a default.
        %   Used before any arithmetic touches the value, so a cell array
        %   or a string in the file cannot throw inside round().
        v = dflt;
        if isfield(cfg,fn)
            q = cfg.(fn);
            if isnumeric(q) && isscalar(q) && isreal(q) && isfinite(q)
                v = double(q);
            end
        end
    end

    function tf = okElements(el, elRC, M, N)
        %OKELEMENTS  Is a loaded element table safe to install?
        %   Checked before assignment, not after: once S.el is replaced
        %   there is no earlier copy to fall back to.
        tf = false;
        if ~isnumeric(el) || ~isnumeric(elRC) || ~ismatrix(el) || ~ismatrix(elRC), return; end
        % Empty is legal, but [] is not the same SHAPE as zeros(0,5): the
        % callers index columns of these arrays, so a 0x0 passes here and
        % fails later. Accepted, and normalised at the call site.
        if isempty(el) && isempty(elRC), tf = true; return; end
        if size(el,2) ~= 5 || size(elRC,2) ~= 2, return; end
        if size(el,1) ~= size(elRC,1), return; end
        if ~all(isfinite(el(:))) || ~all(isfinite(elRC(:))), return; end
        if ~isreal(el) || ~isreal(elRC), return; end
        if any(abs(el(:,1:2)) > MAX_POS_LAMBDA,'all'), return; end
        if any(abs(el(:,3)) > 1e6), return; end
        if any(elRC(:) < 1) || any(elRC(:) ~= round(elRC(:))), return; end
        % Indices have to fall inside the lattice being loaded with them.
        % A table saved from a 16x16 dropped into a 4x4 leaves row/column
        % numbers with no cell to belong to, and every path that maps an
        % element back to its lattice position then reads out of range.
        if nargin >= 4
            if any(elRC(:,1) > M) || any(elRC(:,2) > N), return; end
        end
        tf = true;
    end

    function v = cfgNum(h, v)
        if ~isnumeric(v) || ~isscalar(v) || ~isfinite(v) || ~isreal(v)
            v = h.Value; return;
        end
        lim = h.Limits;
        v = double(min(max(v, lim(1)), lim(2)));
    end
    function v = cfgItem(h, v)
        % Dropdowns throw on any string outside Items -- which is
        % exactly what a config from a version with different element
        % factors / shapes / tapers would carry.
        choices = h.Items;
        if ~isempty(h.ItemsData), choices = h.ItemsData; end
        if validText(v) && ismember(char(v), choices)
            v = char(v);
        else
            v = h.Value;
        end
    end
    function v = cfgBool(h, v)
        if isscalar(v) && (islogical(v) || (isnumeric(v) && isfinite(v) && isreal(v)))
            v = logical(v);
        else
            v = h.Value;
        end
    end
    function v = cfgText(h, v)
        if validText(v)
            v = char(v);
        else
            v = h.Value;
        end
    end

    function lam = lambdaMM()
        % c = 299792458 m/s -> lambda(mm) = 299.792458 / freq(GHz)
        lam = 299.792458 / S.freqGHz;
    end
    function propagationRatioValue = freqRatio(fOpGHz)
        % fOpGHz (default: the Operating frequency) lets an analysis that
        % sweeps frequency -- the band sweep -- use this same ratio at
        % each of its frequencies instead of a copy of it.
        if nargin < 1, fOpGHz = S.freqOpGHz; end
        propagationRatioValue = fOpGHz/S.freqGHz;
    end
    function k = freqScale()
        % Spinner value = frequency-in-GHz * k.
        if strcmp(S.freqUnit,'MHz'), k = 1000; else, k = 1; end
    end
    function assignFreqUnit(u)
        S.freqUnit = u;
        applyFreqUnit();
        % No maybeCompute: nothing about the DESIGN changed, only how it
        % is displayed. Recomputing here would be pure churn.
    end
    function applyFreqUnit()
        %APPLYFREQUNIT  Re-express both spinners in the current unit.
        %   Limits move with the unit as well as the value: leaving them
        %   at the GHz range would cap an MHz entry at 1e6 MHz = 1 THz on
        %   the top end and, worse, refuse anything below 0.001 MHz while
        %   claiming to accept 0.001 of the displayed unit.
        k = freqScale();
        % Order matters, and the obvious order is wrong. Setting Limits
        % first throws whenever the value still on screen falls outside
        % the NEW range -- switching 1e6 GHz to MHz makes the value
        % 1e9 MHz, and switching back then tries to impose a 1e6 ceiling
        % on a spinner still reading 1e9. Found by the fuzz section, not
        % by hand. Widening to infinite first makes every transition
        % legal regardless of which direction it goes.
        for sp = [spFreq spFreqOp]
            sp.Limits = [-Inf Inf];
            sp.Step   = 0.5*k;
        end
        spFreq.Value   = S.freqGHz*k;
        spFreqOp.Value = S.freqOpGHz*k;
        for sp = [spFreq spFreqOp]
            sp.Limits = [0.001 1e6]*k;
        end
        lblLambdaMM.Text = sprintf('Design λ = %.4g mm', lambdaMM());
        refreshSpacingMM();
    end
    function assignFreq(v)
        S.freqGHz = v/freqScale();
        lblLambdaMM.Text = sprintf('Design λ = %.4g mm', lambdaMM());
        % Design frequency now affects the actual pattern math whenever
        % squint mode is on (it's the frequency the phase shifters were
        % "set" at) -- previously this only fed the TSV export, so this
        % call wasn't needed before squint mode existed.
        refreshAll();
        warnImportFrequency();
    end
    function assignFreqOp(v)
        S.freqOpGHz = v/freqScale();
        refreshAll();
        warnImportFrequency();
    end
    function phaseRatioValue = phaseFreqRatio()
        %PHASEFREQRATIO  The ratio the STEERING PHASE is scaled by.
        %   Every feed-phase site in the app goes through this one
        %   function, which is what makes the squint switch a single
        %   decision rather than eight of them: the phase map, the TSV
        %   and CST exports, computePattern, the cut, the max-scan sweep
        %   and the scan-loss sweep all read it. (The band sweep is the
        %   one exception: it sweeps frequency itself, so it hands the
        %   same law its own ratio per frequency -- see runBandSweepCore.)
        %
        %   Retuned: the same ratio the PROPAGATION term uses
        %   (freqRatio()), so the two cancel and the stationary-phase
        %   point sits at U = us -- the commanded angle, at every
        %   frequency.
        %
        %   Squinting: 1. The phases are the ones computed at the design
        %   frequency and they do not move, while propagation still
        %   scales by freqRatio(). The beam therefore lands where
        %   k*fr*x*U = k*x*us, i.e. sin(theta_beam) = sin(theta_s)/fr.
        if S.retunePhase
            phaseRatioValue = S.freqOpGHz/S.freqGHz;
        else
            phaseRatioValue = 1;
        end
    end
    function fGHz = phaseReferenceGHz()
        %PHASEREFERENCEGHZ  Frequency at which the current steering phases are defined.
        %   Retuned mode -> operating frequency. Beam-squint mode -> design frequency.
        %   Keeping this as one helper prevents exports/readouts from labelling frozen
        %   design-frequency phases as operating-frequency phases.
        fGHz = S.freqGHz*phaseFreqRatio();
    end

    function drawSteeringPointer(ax)
        %DRAWSTEERINGPOINTER  Commanded direction plus squinted direction when phases are frozen.
        usP = sind(S.theta_s)*cosd(S.phi_s);
        vsP = sind(S.theta_s)*sind(S.phi_s);
        hold(ax,'on');
        if S.retunePhase || size(S.el,1) <= 1
            % For a single element there is no progressive AF phase gradient,
            % so a frozen phase cannot create beam squint.
            plot3(ax,[0 1.05*usP],[0 1.05*vsP],[0 1.05*cosd(S.theta_s)], ...
                'r-','LineWidth',2.5);
        else
            % Dashed red = commanded design-frequency direction. Solid green =
            % ideal array-factor squint prediction at the operating frequency.
            % The total-pattern peak may differ slightly because of element-factor
            % pull, manual phase offsets, taper, or a competing grating lobe.
            plot3(ax,[0 1.05*usP],[0 1.05*vsP],[0 1.05*cosd(S.theta_s)], ...
                'r--','LineWidth',1.3);
            sPtr = squintSinTheta(S.theta_s);
            if abs(sPtr) <= 1
                thPtr = asind(sPtr);
                uPtr = sind(thPtr)*cosd(S.phi_s);
                vPtr = sind(thPtr)*sind(S.phi_s);
                plot3(ax,[0 1.05*uPtr],[0 1.05*vPtr],[0 1.05*cosd(thPtr)], ...
                    'g-','LineWidth',2.5);
            end
        end
        hold(ax,'off');
    end

    function s = squintSinTheta(thNominalDeg, phRatio, fr)
        %SQUINTSINTHETA  sin(theta) predicted by the commanded steering phase term.
        %   sin(theta_s) when retuned; sin(theta_s)*f_design/f_op when
        %   the phases are frozen. phRatio and fr default to the main
        %   window's phaseFreqRatio() and freqRatio(); the band sweep
        %   passes its own for each swept frequency and steering mode.
        %   May exceed 1 in magnitude, which means
        %   the commanded steering term has no stationary direction in visible space.
        %   Do not clamp that case to grazing and call it a beam direction;
        %   a finite array can still have a visible maximum that must be
        %   found from the actual pattern.
        %
        %   One implementation, deliberately: this law was previously
        %   written out inline at the one call site that needed it, and a
        %   second reader of it is exactly how the two would drift.
        if nargin < 2, phRatio = phaseFreqRatio(); end
        if nargin < 3, fr = freqRatio(); end
        if ~isfinite(fr) || fr <= 0, s = sind(thNominalDeg); return; end
        s = sind(thNominalDeg)*phRatio/fr;
    end
    function assignSquint(v)
        % The box is the inverse of the field -- see S.retunePhase.
        S.retunePhase = ~logical(v);
        refreshAll();
    end
    function syncFrequencyMode()
        % squintMode is not user-switchable; retunePhase is, so this must
        % MIRROR it rather than force it. Forcing it here is what made a
        % loaded config silently discard the squint setting it carried.
        S.squintMode = true;
        % set() rather than cbSquint.Value = ...: this runs only from
        % loadConfig, so if the handle were ever not yet built, a dot
        % assignment would quietly turn cbSquint into a STRUCT and the
        % box would stop tracking the field with no error anywhere --
        % the same trap already documented at lblLambdaMM. set() reads
        % the handle, so it throws instead, and the read is also what
        % clears the Code Analyzer's STRNU warning here.
        set(cbSquint,'Value',~S.retunePhase);
        showFreqOp(true);
    end
    function showFreqOp(~)
        lblFreqOp.Visible='on'; spFreqOp.Visible='on'; lblLambdaMM.Visible='on';
    end

    function exportTSV()
        if ~patternReady(true), return; end
        % Matches a REAL CST-exported .tsv sample the user provided
        % (CST Studio Suite 2025, "unit: meters"), not a guessed format:
        %
        %   # Created by <...>
        %   # On <date>
        %   # unit: meters
        %   # design frequency: -
        %   # Element	X	Y	Z	Magnitude	Phase	Phi	Theta	Gamma
        %   1,1	0	0	0	<mag>	<phase>	<phi>	<theta>	<gamma>
        %   1,2	0	0.01	0	...
        %
        % Confirmed directly from that sample (high confidence):
        %  - tab-separated, "#"-prefixed comment/header block first
        %  - Element ID is the literal text "Xindex,Yindex" (verified
        %    from the sample: e.g. "2,1" has X=0.01,Y=0 and "1,2" has
        %    X=0,Y=0.01 -- first number tracks X, second tracks Y),
        %    not two separate numeric columns
        %  - positions are in METERS
        %  - Magnitude is the element amplitude. Phase is the app's ELEMENT-
        %    LEVEL steering/manual phase (plus the optional extra rotation term
        %    when enabled). This TSV does NOT contain the separate H/V port
        %    phases or L/R/B/T dual-feed topology; use the Phase Scheme Map /
        %    CST excitation macro for those.
        %
        % INFERRED, not confirmed by the sample (it only ever showed
        % Theta=0, Gamma=0, and a Phi that never varied in the rows
        % shown): Phi/Theta/Gamma are assumed to be per-element
        % orientation angles for a flat, untilted array, so this maps
        % Phi = this app's per-element Rotation (deg), Theta = 0,
        % Gamma = 0 for every element. If your elements need real 3D
        % tilt (theta != 0) this mapping will need revisiting -- this
        % app only ever models coplanar, in-plane-rotated elements.
        if isempty(S.el)
            uialert(fig,'No elements to export.','Nothing to export');
            return;
        end
        [f,p] = uiputfile('array_cst.tsv','Export TSV for CST array import');
        figure(fig);   % see loadImportedFF's comment -- native file dialogs
                       % can leave the uifigure behind other windows on macOS
        if isequal(f,0), return; end
        try
            [mag,phase]=elementExportValues();
            lambdaM=lambdaMM()/1000;
            Xm=S.el(:,1)*lambdaM; Ym=S.el(:,2)*lambdaM;
            if any(~isfinite([Xm;Ym]))
                error('PAD:ElementExport','Converted physical coordinates exceed numeric range.');
            end
            Zm=zeros(size(Xm));
            phase=mod(phase+180*(mag<0),360);
            mag=abs(mag);
        catch exportErr
            uialert(fig,exportErr.message,'Export failed'); return;
        end
        phi   = S.el(:,5);          % this app's per-element Rotation -> Phi
        theta = zeros(size(Xm));    % coplanar array: no out-of-plane tilt
        gamma = zeros(size(Xm));
        % fid tracked outside the try so the catch can close it: a
        % failure part-way through the fprintf loop (disk full, a
        % revoked permission, a network drive dropping) would otherwise
        % report the error while leaving the handle open, and every
        % retry would leak another one.
        % Element IDs are what CST keys the array on, so they have to be
        % UNIQUE. The ID is built from S.elRC, which identifies the
        % LATTICE CELL rather than the element -- and Sub-position mode
        % deliberately places up to 9 elements inside a single cell, all
        % of which inherit that cell's elRC. Those elements collided:
        % exporting a 4x4 grid with one sub-position element added wrote
        % 17 rows carrying only 16 distinct IDs (two rows both labelled
        % "1,1", at different physical positions), so the import kept one
        % and silently dropped the other.
        %
        % Cells with a single occupant -- every uniform and sparse array,
        % i.e. the overwhelmingly common case -- keep the exact "col,row"
        % form confirmed against the real CST sample. Only the 2nd and
        % later occupants of a shared cell take a "_2"/"_3" suffix, so
        % nothing that exported correctly before changes at all.
        [~,~,cellOf] = unique(S.elRC,'rows','stable');
        occ = ones(numel(cellOf),1);
        seen = zeros(max(cellOf),1);
        for n = 1:numel(cellOf)
            seen(cellOf(n)) = seen(cellOf(n)) + 1;
            occ(n) = seen(cellOf(n));
        end
        nShared = sum(occ > 1);
        fid = -1;
        try
            fid = fopen(fullfile(p,f),'w');
            if fid < 0, error('Could not open file for writing.'); end
            fprintf(fid,'# Created by MATLAB phasedArrayDesigner.m (not CST Studio Suite)\n');
            fprintf(fid,'# On %s\n', char(datetime('now')));
            fprintf(fid,'# unit: meters\n');
            fprintf(fid,'# design frequency: %.6g GHz\n', S.freqGHz);
            fprintf(fid,'# phase reference frequency: %.6g GHz\n',phaseReferenceGHz());
            fprintf(fid,'# operating frequency: %.6g GHz\n',S.freqOpGHz);
            fprintf(fid,'# beam squint: %s\n',ternStr(~S.retunePhase,'ON (phases frozen at design frequency)','OFF (phases retuned at operating frequency)'));
            fprintf(fid,'# Element\tX\tY\tZ\tMagnitude\tPhase\tPhi\tTheta\tGamma\n');
            for n = 1:size(S.el,1)
                if occ(n) > 1
                    eid = sprintf('%d,%d_%d', S.elRC(n,2), S.elRC(n,1), occ(n));
                else
                    eid = sprintf('%d,%d', S.elRC(n,2), S.elRC(n,1));
                end
                fprintf(fid,'%s\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\t%.10g\n', ...
                    eid, Xm(n), Ym(n), Zm(n), ...
                    mag(n), phase(n), phi(n), theta(n), gamma(n));
            end
            [writeMessage,writeCode]=ferror(fid);
            closeCode=fclose(fid); fid=-1;
            if writeCode~=0 || closeCode~=0
                error('PAD:TSVWrite','Could not finish writing TSV: %s',writeMessage);
            end
            % Physical element rotation and extra electrical rotation phase
            % are intentionally independent.  For this dual-feed L/R/B/T
            % topology the port-side/CP phase scheme is handled in the Phase
            % scheme map; do not warn that S.seqPhase must be enabled.
            rotTxt = '';
            if nShared > 0
                shareTxt = sprintf(['\n\n%d element(s) share a lattice cell ' ...
                    '(Sub-position mode). Their IDs carry a _2/_3 suffix so ' ...
                    'every row stays uniquely addressable -- without it CST ' ...
                    'would keep only one element per cell.'], nShared);
            else
                shareTxt = '';
            end
            % Steering reported through directionText, like every other
            % message, so it names the convention on screen. It used to be
            % a bare "(%.0f,%.0f) deg" of the canonical theta/phi: in az/el
            % mode az 210 / el 60 read as "(30,210)" -- apparently az 30,
            % el 210 -- and a design converted from theta -30 read "(-30,20)"
            % while its controls showed az 200, el 60.
            uialert(fig, [sprintf(['Wrote %d elements to %s (positions in meters, ' ...
                'freq=%.4g GHz -> λ=%.4g mm). Phase column includes your current ' ...
                'steering direction (%s), not just the manual per-element offset.' ...
                '\n\nColumn layout matches a real CST-exported .tsv sample. Phi/Theta/' ...
                'Gamma mapping (Phi=element rotation, Theta=Gamma=0) is inferred, not ' ...
                'confirmed by that sample -- double-check it imports correctly, ' ...
                'especially if any elements use CP/sequential rotation.' newline newline ...
                'IMPORTANT: TSV Phase is element-level only. It does NOT encode separate H/V ' ...
                'port IDs, L/R/B/T feed sides, or the Dual Feed LP/CP port-phase scheme. ' ...
                'For physical port excitation use EXPORT ▸ CST macro (also in the Phase scheme map).'], ...
                size(S.el,1), f, S.freqGHz, lambdaMM(), directionText(S.theta_s,S.phi_s)) ...
                sprintf('\nPhase reference: %.4g GHz. Operating frequency: %.4g GHz. Beam squint: %s.', ...
                    phaseReferenceGHz(),S.freqOpGHz,ternStr(~S.retunePhase,'ON','OFF')) shareTxt rotTxt], ...
                'TSV exported', 'Icon','info');
            setStatus(sprintf('Exported the CST array to %s: %d elements', ...
                f,size(S.el,1)),'good');
        catch err
            if fid >= 0, fclose(fid); end
            uialert(fig, err.message, 'Export failed');
            setStatus(['CST .tsv export failed: ' err.message],'bad');
        end
    end

    function ok = saveConfig()
        %SAVECONFIG  File > Save and the Save button.
        %   Writes straight to the open design's file; a design not saved
        %   yet asks for a name first, as Save As does. ok is false when
        %   nothing was written (the question was cancelled, or the write
        %   failed and said so) -- closing the window depends on it.
        if isempty(docFile)
            ok = saveConfigAs();
        else
            ok = writeConfig(docFile);
        end
    end

    function ok = writeConfig(cfgFile)
        %WRITECONFIG  Save the design to cfgFile; it becomes the open file.
        cfg = struct('M',S.M,'N',S.N,'dx',S.dx,'dy',S.dy,'subOff',S.subOff, ...
            'theta_s',S.theta_s,'phi_s',S.phi_s, ...
            'angleConvention',S.angleConvention,'mode',S.mode,'el',S.el, ...
            'elRC',S.elRC,'taper',S.taper,'sll',S.sll,'efType',S.efType, ...
            'efQ',S.efQ,'efBeamAz',S.efBeamAz,'efBeamEl',S.efBeamEl, ...
            'seqBlockM',S.seqBlockM,'seqBlockN',S.seqBlockN, ...
            ... % The imported far-field travels WITH the design. Saving only
            ... % efType = 'Imported (CST far-field)' and not the data meant a
            ... % reopened design silently fell back to isotropic, or -- worse --
            ... % picked up whatever unrelated pattern happened to be loaded in
            ... % that session, which also flips seqRotSign() and so the whole
            ... % rotation compensation. The struct is a few MB at worst.
            'impFF',S.impFF,'impFFName',S.impFFName, ...
            ... % ...with the frequency it states, so a reopened design
            ... % still warns when it is evaluated at another one.
            'impFFGHz',S.impFFGHz,'impFFFreqFrom',S.impFFFreqFrom, ...
            ... % Saved so a substitution stays visible across sessions. The
            ... % blocking path drops efType to Isotropic, and without this
            ... % the file records a deliberate isotropic design -- the
            ... % substitution's history erased by the act of saving it.
            'impNeedsPattern',S.impNeedsPattern, ...
            'seqPhase',S.seqPhase,'dualFeedCP',S.dualFeedCP, ...
            'portMap',S.portMap,'portMapSet',S.portMapSet, ...
            'customFormula',S.customFormula, ...
            'customFormulaPh',S.customFormulaPh,'arrayShape',S.arrayShape, ...
            'gridAngle',S.gridAngle,'stagger',S.stagger,'freqGHz',S.freqGHz, ...
            'steeringSchema',2,'squintMode',S.squintMode,'retunePhase',S.retunePhase,'freqOpGHz',S.freqOpGHz,'freqUnit',S.freqUnit, ...
            'cutMode',S.cutMode,'cutFixedTheta',S.cutFixedTheta, ...
            'fullSphere',S.fullSphere, ...
            'cutPhi',S.cutPhi,'cutPhiFollow',S.cutPhiFollow,'impUnitCell',S.impUnitCell, ...
            'impTotEffPct',S.impTotEffPct, ...
            ... % View settings. These are read straight off their widgets
            ... % rather than mirrored in S, so they are saved by widget
            ... % value and restored the same way (see loadConfig) instead
            ... % of going through the generic S.(field) loop. Three of
            ... % them -- surface source, polarization readout and
            ... % absolute level -- change the NUMBERS on screen, not just
            ... % the styling, so a reopened design that ignored them
            ... % displayed different figures than when it was saved. The
            ... % cut plane, full-sphere, unit-cell and efficiency
            ... % settings were already being saved; these complete the set.
            'vShow',ddShow.Value,'vPol',ddPol.Value,'vScale',ddScale.Value, ...
            'vDynRange',spDR.Value,'vAbsLevel',cbAbs.Value, ...
            'vARCut',cbARCut.Value,'vAuto',cbAuto.Value);
        % The design targets are the spec this design is judged against,
        % so they travel with it (loadConfigCore restores them).
        cfg.targets = S.targets;
        % Which version wrote the file, and when: loadConfig ignores both,
        % they are for the person who finds the file later. The names
        % contain "version" and "time" on purpose -- the test suites drop
        % such fields when they compare two saves of one design.
        cfg.appVersion = APP_VERSION;
        cfg.savedTime = char(datetime('now','Format','yyyy-MM-dd HH:mm:ss'));
        try
            save(cfgFile,'cfg');
        catch err
            ok = false;
            setStatus(['Save failed: ' err.message],'bad');
            uialert(fig, err.message, 'Save failed');
            return;
        end
        ok = true;
        docFile = cfgFile; docSaved = designSnapshot();
        rememberFile(cfgFile);
        refreshDocTitle();
        docStatus(sprintf('Saved %s at %s', fileLabel(cfgFile), ...
            char(datetime('now','Format','HH:mm'))), false);
    end

    % ------------------------------------------ the open design file (C3)
    % The window title names the open file and marks unsaved changes;
    % Save writes to that file, Save As picks a new one; File > Open Recent
    % keeps the last RECENT_MAX files; the close button asks before
    % unsaved changes are lost. "Unsaved" compares the design with the one
    % last saved or opened (docSaved), not a flag set on every edit, so
    % undoing back to the saved design clears the marker again.

    function ok = saveConfigAs()
        %SAVECONFIGAS  File > Save As…: a new file, which becomes the open one.
        [f,p] = uiputfile('*.mat','Save design as',suggestedDesignPath());
        figure(fig);   % see loadImportedFF's comment -- native file dialogs
                       % can leave the uifigure behind other windows on macOS
        ok = false;
        if isequal(f,0), return; end
        % save() adds .mat to a bare name; the title and the recent list
        % must name the file it actually wrote.
        [~,~,fExt] = fileparts(f);
        if isempty(fExt), f = [f '.mat']; end
        ok = writeConfig(fullfile(p,f));
    end

    function p = suggestedDesignPath()
        % Save As offers the open file itself; a new design gets a name
        % that says what it is (not the same array_config.mat every time),
        % in the folder last used.
        if ~isempty(docFile), p = docFile; return; end
        p = fullfile(lastFolder(), sprintf('array_%dx%d_%gGHz.mat', ...
            S.M, S.N, S.freqGHz));
    end

    function d = lastFolder()
        d = prefValue('lastDir', '');
        if ~(ischar(d) && isfolder(d)), d = pwd; end
    end

    function dirty = designDirty()
        dirty = ~isempty(docSaved) && ~isequaln(designSnapshot(), docSaved);
    end

    function refreshDocTitle()
        %REFRESHDOCTITLE  "App — file.mat *": the open file, * if unsaved.
        %   The bare app name while nothing has been saved or changed, so
        %   a freshly started window is still found by its plain name.
        dirty = designDirty();
        docTitle = APP_NAME;
        if ~isempty(docFile)
            docTitle = [APP_NAME ' — ' fileLabel(docFile)];
        elseif dirty
            docTitle = [APP_NAME ' — Untitled'];
        end
        if dirty, docTitle = [docTitle ' *']; end
        if isgraphics(fig) && ~strcmp(fig.Name, docTitle), fig.Name = docTitle; end
    end

    function nm = fileLabel(fPath)
        [~,fName,fExt] = fileparts(fPath);
        nm = [fName fExt];
    end

    function docStatus(msg, keepWarning)
        % A save or open confirmed in the status bar -- after the compute
        % an open runs. A warning still true stays readable after it:
        % plots out of date with Auto off, an empty design, a failed
        % compute, and for an open (keepWarning) whatever the load itself
        % warned about, e.g. an imported pattern at another frequency.
        warnNow = any(strcmp(statusTone, {'warn','bad'})) && ...
            (keepWarning || ~S.radiationValid);
        if warnNow
            why = statusMsg;
            if ~endsWith(why, {'.','!','?'}), why = [why '.']; end
            setStatus(sprintf('%s — %s', msg, why), statusTone);
        elseif ~S.radiationValid
            setStatus([msg ' · plots out of date: click Compute pattern'],'warn');
        else
            setStatus(msg,'good');
        end
    end

    function docWas = docSwap(docTo)
        % Undo / Redo of an Open: put back the file that was open on the
        % other side of it, and hand over the one being left.
        docWas = [];
        if isempty(docTo), return; end
        docWas = struct('file',docFile,'saved',docSaved);
        docFile = docTo.file; docSaved = docTo.saved;
    end

    % ---- recent files (MATLAB preferences) ----------------------------------
    function v = prefValue(name, default)
        % A preference that cannot be read -- none yet, or a damaged
        % preferences file -- is simply not there.
        v = default;
        try
            if ispref(PREF_GROUP, name), v = getpref(PREF_GROUP, name); end
        catch
        end
    end

    function setPrefValue(name, v)
        % Remembering is a convenience: a preferences folder that cannot
        % be written must not make a save or an open fail.
        try
            setpref(PREF_GROUP, name, v);
        catch
        end
    end

    function files = recentFiles()
        files = prefValue('recentFiles', {});
        if ~iscell(files) || ~all(cellfun(@(x) ischar(x) && isrow(x), files))
            files = {};
        end
        files = files(:)';
    end

    function rememberFile(fPath)
        % Most recent first, each file once, at most RECENT_MAX.
        files = recentFiles();
        files = [{fPath} files(~strcmp(files, fPath))];
        setPrefValue('recentFiles', files(1:min(end, RECENT_MAX)));
        setPrefValue('lastDir', fileparts(fPath));
        refreshRecentMenu();
    end

    function refreshRecentMenu()
        %REFRESHRECENTMENU  Rebuild File > Open Recent from the preferences.
        %   Files that no longer exist are dropped from the stored list as
        %   well, so a moved design does not come back on every start.
        if isempty(mRecent) || ~isgraphics(mRecent), return; end
        files = recentFiles();
        there = cellfun(@isfile, files);
        if ~all(there)
            files = files(there);
            setPrefValue('recentFiles', files);
        end
        % Rebuilt only when the list changed: this also runs as File opens.
        if isequal(mRecent.UserData, files) && ~isempty(mRecent.Children)
            return;
        end
        mRecent.UserData = files;
        delete(mRecent.Children);
        if isempty(files)
            uimenu(mRecent,'Text','No recent designs','Tag','menuRecentNone', ...
                'Enable','off', ...
    'Tooltip','No recent designs. Apply this setting or action to the current design.');
            return;
        end
        names = cellfun(@fileLabel, files, 'UniformOutput', false);
        for kRec = 1:numel(files)
            recTxt = names{kRec};
            % Two files of one name: say which folder each is in.
            if nnz(strcmp(names, recTxt)) > 1
                [~,recDir] = fileparts(fileparts(files{kRec}));
                recTxt = sprintf('%s — %s', recTxt, recDir);
            end
            uimenu(mRecent,'Text',recTxt,'Tag','menuRecentFile', ...
                'UserData',files{kRec},'Tooltip',files{kRec}, ...
                'MenuSelectedFcn',@(s,~)openRecent(s.UserData));
        end
        uimenu(mRecent,'Text','Clear Menu','Tag','menuRecentClear', ...
            'Separator','on','Tooltip','Forget the recent designs (the files stay)', ...
            'MenuSelectedFcn',@(s,e)clearRecent());
    end

    function openRecent(fPath)
        if ~isfile(fPath)
            refreshRecentMenu();
            setStatus(sprintf('%s no longer exists: removed from Open Recent.', ...
                fileLabel(fPath)),'warn');
            uialert(fig, sprintf(['%s\n\nThis file no longer exists (moved, ' ...
                'renamed or deleted). It has been removed from Open Recent.'], ...
                fPath), 'Design not found', 'Icon','warning');
            return;
        end
        loadConfig(fPath);
    end

    function clearRecent()
        setPrefValue('recentFiles', {});
        refreshRecentMenu();
        setStatus('Open Recent cleared.');
    end

    % ---- New and close share the Save / Discard / Cancel decision -----------
    function onCloseRequest(isNew)
        if nargin == 0, isNew = false; end
        if closeAsking, return; end
        try
            if ~designDirty()
                if isNew, resetNewProject(); else, delete(fig); end
                return;
            end
            closeAsking = true;
            if isNew, verb = 'starting a new design';
            else, verb = 'closing'; end
            if isempty(docFile)
                question = sprintf('This design has not been saved. Save it before %s?',verb);
            else
                question = sprintf('Save the changes to %s before %s?', ...
                    fileLabel(docFile),verb);
            end
            if isNew
                consequence = 'Discard starts a new design and loses the changes.';
            else
                consequence = 'Discard closes the window and loses the changes.';
            end
            uiconfirm(fig, [question newline newline consequence], ...
                'Unsaved changes', 'Options',{'Save','Discard','Cancel'}, ...
                'DefaultOption',1, 'CancelOption',3, 'Icon','warning', ...
                'CloseFcn',@(~,evt)leavingAnswer(isNew,evt.SelectedOption));
        catch err
            closeAsking = false;
            % If the confirmation UI fails, keep the unsaved design open.
            % Closing here would silently discard the user's changes.
            if isgraphics(fig)
                if isNew
                    action = 'New design';
                else
                    action = 'Close';
                end
                setStatus([action ' cancelled: ' err.message],'bad');
                uialert(fig,err.message,[action ' cancelled']);
            end
        end
    end

    function leavingAnswer(isNew,choice)
        % A failed save keeps the window open (the save said why); asking
        % again then offers Discard. A cancelled Save As also stops New.
        closeAsking = false;
        if ~isgraphics(fig), return; end
        if strcmp(choice,'Cancel')
            setStatus('Action cancelled: the design is still open.');
            return;
        end
        if strcmp(choice,'Save') && ~saveConfig(), return; end
        if isNew, resetNewProject(); else, delete(fig); end
    end

    function resetNewProject()
        % Close result panes first so no old result, phase map, or CST
        % comparison can remain visible over the new document.
        for kResult = 1:numel(resultTaskNames)
            if ~isempty(resultPanels{kResult}) && isgraphics(resultPanels{kResult})
                closeResult(resultTaskNames{kResult});
            end
        end
        if ~isempty(targetsWin) && isgraphics(targetsWin), delete(targetsWin); end
        if ~isempty(elementGalleryPopup) && isgraphics(elementGalleryPopup)
            elementGalleryPopup.Visible = 'off';
        end
        if ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup)
            shapeGalleryPopup.Visible = 'off';
        end
        S.refCut = [];
        S.cstOverlay = [];
        S.coverageCache = [];
        S.pattern2DBasis = [];
        S.cstScheme = 'single';
        cstComparisonClosed = false;
        syncRefControls();

        % View settings are deliberately outside Undo's design snapshot.
        % Restore them explicitly, then restore the launch design in one
        % recompute. The file identity and Undo history belong to the old
        % document and are cleared only after that restore succeeds.
        S.angleConvention = defaultView.angleConvention;
        S.freqUnit = defaultView.freqUnit;
        S.cutMode = defaultView.cutMode;
        S.cutFixedTheta = defaultView.cutFixedTheta;
        S.fullSphere = defaultView.fullSphere;
        S.cutPhi = defaultView.cutPhi;
        S.cutPhiFollow = defaultView.cutPhiFollow;
        S.targets = defaultView.targets;
        cbFullSphere.Value = S.fullSphere;
        cbCutFollow.Value = S.cutPhiFollow;
        spCutPhi.Enable = ternStr(S.cutPhiFollow,'off','on');
        ddShow.Value = defaultView.show;
        ddPol.Value = defaultView.pol;
        ddScale.Value = defaultView.scale;
        spDR.Value = defaultView.dynRange;
        cbAbs.Value = defaultView.absLevel;
        cbARCut.Value = defaultView.arCut;
        cbAuto.Value = defaultView.auto;
        layoutView = defaultView.layoutView;
        cbLayoutIndex.Value = layoutView.indices;
        cbLayoutTaper.Value = layoutView.taper;
        cbLayoutAxes.Value = layoutView.localAxes;
        cbLayoutAnnotation.Value = layoutView.annotation;
        shapePreviewOlderModel = 'Diamond';
        shapePreviewPreviousModel = 'Circle';
        restoreDesign(defaultDesign);
        docFile = '';
        docSaved = designSnapshot();
        docNextStep = [];
        undoStack = {};
        redoStack = {};
        undoBase = docSaved;
        lastUndoUI = [];
        updateUndoMenus();
        selectTask('Array');
        refreshDocTitle();
        setStatus('New design ready.');
    end

    % ---------------------------------------------- CST path (C1)
    % Import: what the CST export must look like, what was read from it,
    % and whether its frequency is the one the array is evaluated at.
    % Export: the CST macro, from the ribbon/menu and from the Phase
    % scheme map through one implementation (exportCstExcitation).

    function exportCstMacro()
        %EXPORTCSTMACRO  EXPORT > CST macro in the main window.
        %   The Phase scheme map's export, reached without opening the
        %   map; the form starts on the scheme last exported.
        exportCstExcitation(S.cstScheme, fig, false);
    end

    function [cmNames, cmKeys] = cstSchemes()
        % Display names as the Phase scheme map words them, and the keys
        % cstPortWeights takes.
        cmNames = {'Single feed','Dual feed LP','Dual feed CP'};
        cmKeys  = {'single','lp','cp'};
    end

    function openCstMacroDialog(scheme, cParent, fromMap)
        %OPENCSTMACRODIALOG  The CST macro settings, in one form window.
        %   Returns at once: Export (cmDoExport) finishes the job. The form
        %   is modal, so the design cannot change while it is open -- its
        %   default name and notes describe the design as it was when it
        %   opened, and the export then writes that same design.
        [cmNames, cmKeys] = cstSchemes();
        delete(findall(0,'Type','figure','Tag','cstMacroDlg'));
        cmPos = [0 0 580 460];
        try   % centred on the window that asked
            cmPos(1:2) = cParent.Position(1:2) + ...
                max((cParent.Position(3:4) - cmPos(3:4))/2, 0);
        catch
        end
        cmDlg = uifigure('Name','CST macro export','Tag','cstMacroDlg', ...
            'Position',cmPos,'WindowStyle','modal');
        matchTheme(cmDlg);
        cmG = uigridlayout(cmDlg,[8 2]);
        cmG.RowHeight = {22,22,22,22,22,22,'1x',28};
        cmG.ColumnWidth = {170,'1x'};

        uilabel(cmG,'Text','Feed scheme');
        cmScheme = uidropdown(cmG,'Items',cmNames, ...
            'Value',cmNames{strcmp(cmKeys,scheme)},'Tag','cstScheme', ...
            'Tooltip',['Ports driven per element: one (single feed), or H ' ...
            'and V (dual feed, linear or circular). Dual-feed ports follow ' ...
            'the feed layout declared in the Phase scheme map.']);
        % From the Phase scheme map, the map's own Feed scheme decides, as
        % it always has: the map exports the port phases it is showing,
        % and a second, different choice here would contradict them.
        if fromMap
            cmScheme.Enable = 'off';
            cmScheme.Tooltip = ['Set by the Phase scheme map''s Feed scheme: ' ...
                'the export writes the port phases the map shows.'];
        end
        uilabel(cmG,'Text','Phases');
        cmPhase = uidropdown(cmG,'Items',{'Constant phases','Parametric phases'}, ...
            'Tag','cstPhases', ...
            'Tooltip',['Constant: fixed phase values (°) for the current ' ...
            'beam. Parametric: steering θ/φ become CST parameters, so the ' ...
            'beam can be re-steered inside CST.']);
        uilabel(cmG,'Text','CST frequency unit');
        cmUnit = uidropdown(cmG,'Items',{'GHz','MHz','kHz','Hz'}, ...
            'Value',S.freqUnit,'Tag','cstUnit', ...
            'Tooltip',['Frequency unit of the CST project (its Units ' ...
            'setting); the reference frequency is written in it.']);
        uilabel(cmG,'Text','Waveguide port mode');
        cmPort = uispinner(cmG,'Limits',[1 2147483647],'Step',1, ...
            'RoundFractionalValues','on','ValueDisplayFormat','%.0f', ...
            'Value',1,'Tag','cstPortMode', ...
            'Tooltip','Mode number of the waveguide ports that carries the excitation (usually 1).');
        uilabel(cmG,'Text','Combination name');
        cmName = uieditfield(cmG,'text','Tag','cstName', ...
            'Tooltip',['Name of the excitation combination created in CST: ' ...
            'letters, digits, _ - . (at most 120 characters). Cleared, it ' ...
            'goes back to the name made from the design.']);
        uilabel(cmG,'Text','Parameter prefix');
        cmPrefix = uieditfield(cmG,'text','Value','Beam','Tag','cstPrefix', ...
            'Tooltip',['Parametric phases only: prefix of the CST parameters ' ...
            '(Beam_Theta_deg, ...). Combinations with the same prefix share ' ...
            'their steering parameters.']);
        cmNotes = uitextarea(cmG,'Editable','off','FontSize',11,'Tag','cstNotes', ...
            'Tooltip','Conventions this export uses; check them against your CST ports.');
        cmNotes.Layout.Column = [1 2];
        cmBtns = uigridlayout(cmG,[1 3]);
        cmBtns.Layout.Column = [1 2];
        cmBtns.ColumnWidth = {'1x',100,100}; cmBtns.Padding = [0 0 0 0];
        uilabel(cmBtns,'Text','');
        uibutton(cmBtns,'Text','Export…','Tag','cstExport', ...
            'Tooltip','Choose where to save the .bas macro, then write it', ...
            'ButtonPushedFcn',@(~,~)cmDoExport());
        uibutton(cmBtns,'Text','Cancel','Tag','cstCancel', ...
            'Tooltip','Close without exporting', ...
            'ButtonPushedFcn',@(~,~)delete(cmDlg));

        % The name follows the scheme and phase mode until it is typed in.
        cmNameEdited = false;
        cmScheme.ValueChangedFcn = @(~,~)cmSync();
        cmPhase.ValueChangedFcn = @(~,~)cmSync();
        cmName.ValueChangedFcn = @(~,~)cmNameChanged();
        cmSync();

        function cmSync()
            cmKey = cmKeys{strcmp(cmNames,cmScheme.Value)};
            cmParam = strcmp(cmPhase.Value,'Parametric phases');
            cmPrefix.Enable = ternStr(cmParam,'on','off');
            if ~cmNameEdited, cmName.Value = cstDefaultLabel(cmKey,cmParam); end
            cmNotes.Value = splitlines(cstMacroNote(cmKey));
        end
        function cmNameChanged()
            cmNameEdited = ~isempty(strtrim(cmName.Value));
            if ~cmNameEdited, cmSync(); end
        end
        function cmDoExport()
            cmKey = cmKeys{strcmp(cmNames,cmScheme.Value)};
            cmParam = strcmp(cmPhase.Value,'Parametric phases');
            cmLabel = strtrim(cmName.Value);
            cmPre = strtrim(cmPrefix.Value);
            % Wrong answers are sent back to the form, which stays open.
            if cmParam && isempty(regexp(cmPre,'^[A-Za-z][A-Za-z0-9_]{0,23}$','once'))
                uialert(cmDlg,['Use a parameter prefix of 1-24 letters, digits ' ...
                    'or underscores, starting with a letter.'],'Invalid prefix');
                return;
            end
            if isempty(regexp(cmLabel,'^[A-Za-z0-9][A-Za-z0-9_.-]*$','once')) || ...
                    numel(cmLabel) > 120
                uialert(cmDlg,['Use a name starting with a letter or number, ' ...
                    'at most 120 characters, containing only letters, numbers, ' ...
                    'underscores, hyphens and dots.'],'Invalid combination name');
                return;
            end
            cmMode = cmPort.Value;
            cmUnitTxt = cmUnit.Value;
            % Closed BEFORE the file dialog: a native save panel opened
            % from a modal window can come up behind it.
            delete(cmDlg);
            S.cstScheme = cmKey;
            writeCstMacro(cmKey,cParent,cmParam,cmUnitTxt,cmMode,cmLabel,cmPre);
        end
    end

    function writeCstMacro(scheme,cParent,cParametric,cUnitTxt,cMode,cLabel,cPrefix)
        %WRITECSTMACRO  Build the macro for the answers given and save it.
        %   The file is buildCstMacro's text, unchanged from the old
        %   dialog chain: the same answers still give the same bytes.
        if ~isgraphics(cParent) || isempty(S.el), return; end
        try
            [cPorts,cAmp,cPh,~,cUseRot]=cstPortWeights(scheme);
            cParamSnapshot=cstParametricData(cPorts,cPh,scheme,'Beam');
        catch exportErr
            uialert(cParent,exportErr.message,'Invalid excitation'); return;
        end
        [cTheta,cPhi,cFreq,cGeometryFreq] = cstSteering();
        cParam=[];
        if cParametric
            cParam=cParamSnapshot;
            cParam.prefix=cPrefix;
        end
        cUnits = {'ghz','mhz','khz','hz'};
        cScale = [1 1e3 1e6 1e9];
        cUnitIdx = find(strcmpi(cUnitTxt,cUnits),1);
        if isempty(cUnitIdx) || ~validPortMode(cMode)
            uialert(cParent,'Use GHz/MHz/kHz/Hz and a positive integer mode.', ...
                'Invalid CST settings'); return;
        end
        cRef = cFreq*cScale(cUnitIdx);
        figure(cParent);
        [cFile,cPath] = uiputfile('*.bas','Export CST VBA macro',[cLabel '.bas']);
        if ~isgraphics(cParent), return; end
        figure(cParent);
        if isequal(cFile,0), return; end
        try
            cText=buildCstMacro(cPorts,cMode,cAmp,cPh,cRef,cLabel,scheme,cUseRot, ...
                cTheta,cPhi,cFreq,cGeometryFreq,cUnits{cUnitIdx},cParam);
            cFid = fopen(fullfile(cPath,cFile),'wb');
            if cFid < 0, error('PAD:CSTWrite','Cannot open the macro file.'); end
            cCleanup = onCleanup(@()closeIfOpen(cFid));
            cWritten = fprintf(cFid,'%s', [cText sprintf('\r\n')]);
            if cWritten ~= numel(cText)+2
                error('PAD:CSTWrite','The macro file could not be completely written.');
            end
            if fclose(cFid)~=0, error('PAD:CSTWrite','Could not finish writing the macro.'); end
            clear cCleanup;
            uialert(cParent,sprintf(['Saved %s\n\nIn CST, create a VBA macro and paste the file contents.\n' ...
                'Run Main, then select %s in Combine Excitation.\n' ...
                'Check port assignments and amplitudes/phases before simulation.'], ...
                cFile,cLabel),'CST macro saved','Icon','info');
            [cmNames,cmKeys] = cstSchemes();
            setStatus(sprintf('Saved CST macro %s: %d ports, %s, %s phases', ...
                cFile,numel(cPorts),lower(cmNames{strcmp(cmKeys,scheme)}), ...
                ternStr(cParametric,'parametric','constant')),'good');
        catch cErr
            uialert(cParent,cErr.message,'CST export failed');
            setStatus(['CST macro export failed: ' cErr.message],'bad');
        end
    end

    function [cTheta,cPhi,cFreq,cGeometryFreq] = cstSteering()
        % Steering and frequencies as the macro states them. CST
        % Theta/Phi parameters use a nonnegative polar angle. A negative
        % signed-theta command has the same phase gradient at positive
        % theta and phi+180, including before any Az/El edit.
        cGeometryFreq = S.freqGHz; cFreq = phaseReferenceGHz();
        cTheta = S.theta_s; cPhi = S.phi_s;
        if cTheta < 0
            cTheta = -cTheta;
            cPhi = mod(cPhi+180,360);
        end
    end

    function cLabel = cstDefaultLabel(scheme, cParametric)
        % The combination name the form offers: the design in one word.
        try
            [~,~,~,cN,cUseRot] = cstPortWeights(scheme);
        catch
            cLabel = 'CST_beam'; return;
        end
        [cTheta,cPhi,cFreq,cGeometryFreq] = cstSteering();
        if cParametric
            cLabel = sprintf('%dx%d_%del_%s_parametric_geom%gGHz_rot%s', ...
                S.M,S.N,cN,upper(scheme),cGeometryFreq,ternStr(cUseRot,'ON','OFF'));
        else
            cLabel = sprintf('%dx%d_%del_%s_th%g_ph%g_%gGHz_geom%gGHz_rot%s', ...
                S.M,S.N,cN,upper(scheme),cTheta,cPhi,cFreq,cGeometryFreq, ...
                ternStr(cUseRot,'ON','OFF'));
        end
        cLabel = strrep(cLabel,'e+','e');
    end

    function cNote = cstMacroNote(scheme)
        % The conventions the macro relies on, shown in the form.
        try
            [cPorts,~,~,~,cUseRot] = cstPortWeights(scheme);
        catch noteErr
            cNote = noteErr.message; return;
        end
        if S.retunePhase
            cFreqModeNote = ['Beam squint is OFF: parametric frequency changes ' ...
                'RECALCULATE the steering phases, so the commanded angle is maintained.'];
        else
            cFreqModeNote = ['Beam squint is ON: phase frequency is FROZEN at the ' ...
                'design frequency. Parametric export varies steering θ/φ only; change ' ...
                'the CST operating frequency separately to observe squint.'];
        end
        cNote = sprintf(['Exports %d waveguide ports, scheme %s.\n' ...
            'Extra rotation-angle feed phase applied in this export: %s.\n' ...
            'Port IDs follow the element export order and the declared H/V layout ' ...
            '(Phase scheme map ▸ Edit feed layout).\n' ...
            'Amplitudes are relative wave amplitudes (largest = 1), not watts.\n' ...
            'Dual feed uses equal H/V amplitudes; CP adds +90° to V.\n' ...
            'Opposite feed sides use the existing, uncalibrated 180° preset.\n' ...
            'Check these conventions against your CST ports.\n\n%s\n' ...
            'CST macro parameters use θ/φ; azimuth/elevation steering is converted for export.\n' ...
            'Parametric mode keeps geometry and feed/rotation offsets fixed.\n' ...
            'The macro creates a combination only; select it in CST afterwards.'], ...
            numel(cPorts),upper(scheme),ternStr(cUseRot,'ON','OFF'),cFreqModeNote);
    end

    function txt = importGuideText()
        % What loadCstFarfieldASCII accepts, in one short paragraph, so
        % the format is known before the first refusal rather than
        % learned from it. Ludwig-3 is named because CST offers it next
        % to Theta/Phi and the loader refuses it: its components cannot
        % be turned into Theta/Phi without the reference polarization.
        txt = ['From CST: Farfield ▸ Export ▸ ASCII of the realized gain ' ...
            '(dBi) with complex Theta/Phi components (Ludwig-3 is refused), ' ...
            'θ 0–180° and φ 0–360° on a regular grid, element boresight on ' ...
            '+z (θ = 0°). A frequency in the file name, as in CST''s ' ...
            '"farfield (f=20)", is checked against the operating frequency.'];
    end

    function [fGHz, fromWhat] = cstFileFrequency(fileName, filePath)
        %CSTFILEFREQUENCY  The frequency an exported far-field states, in GHz.
        %   CST names a far-field monitor "farfield (f=19.45)" and its
        %   ASCII export after it, so the file name usually carries the
        %   frequency; a comment line above the column titles may state
        %   it instead. NaN when neither states exactly one value.
        %   fromWhat says where it came from, and when the unit was
        %   assumed (see freqInText).
        fromWhat = '';
        [~,impStem] = fileparts(char(fileName));
        [fGHz,impAssumed,impMany] = freqInText(impStem);
        impWhere = 'file name';
        if isnan(fGHz) && ~impMany && nargin > 1
            % Only the lines ABOVE the column titles: the loader has just
            % read the whole file, and a 1° export runs to megabytes.
            impHead = '';
            impFid = fopen(filePath,'r');
            if impFid >= 0
                impClose = onCleanup(@()closeIfOpen(impFid));
                for kHeadLine = 1:50
                    impLine = fgetl(impFid);
                    if ~ischar(impLine) || (contains(impLine,'Theta','IgnoreCase',true) ...
                            && contains(impLine,'Abs(','IgnoreCase',true))
                        break;
                    end
                    impHead = [impHead ' ' impLine]; %#ok<AGROW>
                end
                clear impClose;
            end
            [fGHz,impAssumed,impMany] = freqInText(impHead);
            impWhere = 'file header';
        end
        if impMany
            fromWhat = [impWhere ' names several frequencies'];
        elseif ~isnan(fGHz)
            fromWhat = impWhere;
            if ~isempty(impAssumed)
                fromWhat = sprintf('%s, %s assumed',impWhere,impAssumed);
            end
        end
    end

    function [fGHz, assumedUnit, many] = freqInText(txt)
        % "f=19.45", "f = 20 GHz", "Frequency: 20 GHz" first; otherwise a
        % number with an explicit unit ("19.45GHz"). Several different
        % values -> none (many = true).
        %
        % CST writes "f=20" in the PROJECT's frequency unit, which the
        % file does not name (and the app's own display unit says nothing
        % about). CST's units are a factor 1000 apart, so the one that
        % puts the number nearest the operating frequency is taken:
        % 20 -> 20 GHz, 20000 -> 20 GHz, 900 -> 0.9 GHz in a 20 GHz
        % design. A mismatch up to ~30x still shows as one; assumedUnit
        % names the unit so the readout can say it was not stated.
        fGHz = NaN; assumedUnit = ''; many = false;
        if isempty(txt), return; end
        num = '(\d+(?:\.\d*)?(?:[eE][-+]?\d+)?)';
        tok = regexpi(txt,['(?<![A-Za-z0-9])(?:f|freq|frequency)\s*[=:]\s*' ...
            num '\s*(THz|GHz|MHz|kHz|Hz)?(?![A-Za-z])'],'tokens');
        if isempty(tok)
            tok = regexpi(txt,['(?<![A-Za-z0-9.])' num '\s*(THz|GHz|MHz|kHz)(?![A-Za-z])'],'tokens');
        end
        if isempty(tok), return; end
        unitNames = {'THz','GHz','MHz','kHz','Hz'};
        unitGHz = [1e3 1 1e-3 1e-6 1e-9];
        vals = zeros(1,numel(tok)); guessed = repmat({''},1,numel(tok));
        for kTok = 1:numel(tok)
            tokVal = str2double(tok{kTok}{1});
            kUnit = find(strcmpi(tok{kTok}{2},unitNames),1);
            if isempty(kUnit)
                [~,kUnit] = min(abs(log(tokVal*unitGHz/S.freqOpGHz)));
                guessed{kTok} = unitNames{kUnit};
            end
            vals(kTok) = tokVal*unitGHz(kUnit);
        end
        ok = isfinite(vals) & vals > 0 & vals <= 1e6;
        vals = vals(ok); guessed = guessed(ok);
        if isempty(vals), return; end
        if any(abs(vals-vals(1)) > 1e-9*vals(1))
            many = true; return;
        end
        fGHz = vals(1); assumedUnit = guessed{1};
    end

    function tf = validImportGHz(v)
        tf = isnumeric(v) && isscalar(v) && isreal(v) && ...
            (isnan(v) || (isfinite(v) && v > 0));
    end

    function tf = importFreqMismatch()
        % 1 %: the pattern of a real element barely moves within that,
        % while the deck's band edges (17.7 vs 21.2 GHz) are far outside.
        tf = ~isempty(S.impFF) && isfinite(S.impFFGHz) && ...
            abs(S.impFFGHz - S.freqOpGHz) > 0.01*S.freqOpGHz;
    end

    function txt = importFreqWarning(inFull)
        impAt = sprintf('%.4g GHz',S.impFFGHz);
        opAt = sprintf('%.4g GHz',S.freqOpGHz);
        if ~inFull
            txt = sprintf(['pattern is for %s, but the array is evaluated at ' ...
                '%s (operating frequency)'],impAt,opAt);
            return;
        end
        txt = sprintf(['This far-field was exported at %s (%s), but the array ' ...
            'is evaluated at the operating frequency, %s.\n\nAn imported element ' ...
            'pattern is used exactly as exported: it is not rescaled with ' ...
            'frequency, so its gain, beamwidth and polarization stay those of ' ...
            '%s.\n\nRe-export the far-field from CST at %s, or set Operating ' ...
            'frequency to %s if that is the frequency you mean to study.'], ...
            impAt,S.impFFFreqFrom,opAt,impAt,opAt,impAt);
    end

    function txt = importReadText()
        % Grid, peak and frequency of the loaded pattern, for the readout
        % and the status bar. The grid is read off the interpolant's own
        % nodes (the copy in 0 <= phi < 360), so a pattern restored from a
        % saved design reports the same as a fresh import.
        try
            impPts = S.impFF.FReEth.Points;
            impBase = impPts(:,2) >= 0 & impPts(:,2) < 360;
            impTh = unique(impPts(impBase,1)); impPh = unique(impPts(impBase,2));
            gridTxt = sprintf(['θ %g–%g° × φ %g–%g° in %g° × %g° steps ' ...
                '(%d × %d points)'],impTh(1),impTh(end),impPh(1),impPh(end), ...
                median(diff(impTh)),median(diff(impPh)),numel(impTh),numel(impPh));
        catch
            gridTxt = 'grid not available';
        end
        impPk = S.impFF.peak;
        if impPk(2) == 0
            pkTxt = sprintf('peak realized gain %.2f dBi at θ = 0°',impPk(1));
        else
            pkTxt = sprintf('peak realized gain %.2f dBi at θ = %g°, φ = %g°', ...
                impPk(1),impPk(2),impPk(3));
        end
        if isfinite(S.impFFGHz)
            fTxt = sprintf('%.4g GHz (from the %s)',S.impFFGHz,S.impFFFreqFrom);
        elseif ~isempty(S.impFFFreqFrom)
            fTxt = sprintf('frequency unclear: the %s',S.impFFFreqFrom);
        else
            fTxt = 'no frequency stated in the file name or header';
        end
        txt = sprintf('%s; %s; %s',gridTxt,pkTxt,fTxt);
    end

    function refreshImportInfo()
        %REFRESHIMPORTINFO  The readout under the import button, and the
        %   two import buttons' tooltips. Runs with every redraw and theme
        %   change (S.portMapRefresh), so a frequency edit or a new theme
        %   is reflected at once.
        if isempty(lblImportInfo) || ~isgraphics(lblImportInfo), return; end
        P = pal();
        guide = importGuideText();
        if isempty(S.impFF)
            impTxt = guide; impCol = P.muted;
            impTip = guide;
        else
            impTxt = ['Read ' importReadText() '.'];
            impCol = P.muted;
            if importFreqMismatch()
                impWarn = importFreqWarning(false);
                impTxt = sprintf('%s\n⚠ %s%s.',impTxt,upper(impWarn(1)),impWarn(2:end));
                impCol = P.warn;
            end
            impTip = sprintf('Loaded: %s\n%s\n\nClick to load another.\n%s', ...
                S.impFFName,impTxt,guide);
        end
        lblImportInfo.Text = impTxt;
        lblImportInfo.FontColor = impCol;
        lblImportInfo.Tooltip = impTip;
        bImportFF.Tooltip = impTip;
        bImportRibbon.Tooltip = impTip;
    end

    function importStatus(msg, tone)
        % Posted after the compute an import or frequency change runs. A
        % failed compute keeps its own line; plots left out of date (Auto
        % off) are said as well rather than hidden by this message.
        if ~S.radiationValid
            if strcmp(statusTone,'bad'), return; end
            msg = [msg ' · plots out of date: click Compute pattern'];
            tone = 'warn';
        end
        setStatus(msg,tone);
    end

    function warnImportFrequency()
        % The compute that a frequency change runs replaces the status
        % line; while the imported element is in use at another frequency,
        % the line says that instead.
        if strcmp(S.efType,'Imported (CST far-field)') && importFreqMismatch()
            impWarn = importFreqWarning(false);
            importStatus(['Imported ' impWarn],'warn');
        end
    end

    function restoreImportFreq(cfg)
        % Follows whichever pattern loadConfigCore left in S.impFF: none ->
        % no frequency; the session's own (the file carries no impFF) ->
        % its frequency stays; the file's -> the frequency saved with it,
        % or for an older file, what its saved file name states.
        if isempty(S.impFF)
            S.impFFGHz = NaN; S.impFFFreqFrom = ''; return;
        end
        if ~isfield(cfg,'impFF'), return; end
        if isfield(cfg,'impFFGHz')
            S.impFFGHz = double(cfg.impFFGHz);   % validated on the way in
            S.impFFFreqFrom = '';
            if isfinite(S.impFFGHz), S.impFFFreqFrom = 'saved design'; end
            if isfield(cfg,'impFFFreqFrom') && validText(cfg.impFFFreqFrom) && ...
                    ~isempty(strtrim(char(cfg.impFFFreqFrom)))
                S.impFFFreqFrom = char(cfg.impFFFreqFrom);
            end
        else
            [S.impFFGHz,S.impFFFreqFrom] = cstFileFrequency(S.impFFName);
        end
    end

    % ------------------------------------------------ undo / redo (C2)
    % Snapshot-based: a step holds the whole design (the UNDO_FIELDS of S)
    % as it was before a change, so Undo lands exactly on what was on
    % screen instead of replaying an inverse that could drift. Steps are
    % taken in maybeCompute, which every design change reaches -- one hook
    % rather than a call in each of thirty callbacks, and a callback added
    % later is covered without anyone remembering to. View-only changes go
    % through safeCompute and never get there; and since only UNDO_FIELDS
    % are compared and restored, the cut plane, the 3D display and Auto are
    % neither recorded nor reverted. A snapshot shares its arrays with S
    % until S changes (copy on write), so a step costs only what differs;
    % an imported pattern is shared the same way, not copied.

    function uSnap = designSnapshot()
        uSnap = struct();
        for kUndoF = 1:numel(UNDO_FIELDS)
            uSnap.(UNDO_FIELDS{kUndoF}) = S.(UNDO_FIELDS{kUndoF});
        end
    end

    function undoCheckpoint(uWhat)
        %UNDOCHECKPOINT  Make the change since the last checkpoint a step.
        %   The step holds the design BEFORE the change (undoBase). It is
        %   named by comparing the two designs unless the caller names it.
        %   isequaln, not isequal: an unstated import frequency is NaN, and
        %   NaN ~= NaN would turn every compute into a step.
        if undoHold > 0, return; end
        uNow = designSnapshot();
        if isempty(undoBase)
            undoBase = uNow;              % start-up: the first design
        elseif ~isequaln(uNow, undoBase)
            if nargin < 1 || isempty(uWhat)
                uWhat = changeName(undoBase, uNow);
            end
            % 'doc' is the file an Open left (docNextStep), [] otherwise.
            undoStack{end+1} = struct('design',undoBase,'what',uWhat, ...
                'doc',docNextStep);
            docNextStep = [];
            if numel(undoStack) > UNDO_MAX, undoStack(1) = []; end
            redoStack = {};               % a new change ends the redo chain
            undoBase = uNow;
            undoCount = undoCount + 1;
            refreshDocTitle();            % the unsaved marker
        end
        updateUndoMenus();
    end

    function uTok = holdUndo(uWhat)
        %HOLDUNDO  No steps until uTok dies; then at most one, named uWhat.
        %   An onCleanup token, like beginBusy's, so a load that throws
        %   still releases the hold.
        undoHold = undoHold + 1;
        uTok = onCleanup(@() releaseUndo(uWhat));
    end

    function releaseUndo(uWhat)
        undoHold = max(undoHold - 1, 0);
        if undoHold == 0 && isgraphics(fig), undoCheckpoint(uWhat); end
        % An Open that changed nothing made no step to carry its file.
        if undoHold == 0, docNextStep = []; end
    end

    function undoFromUI(uAct, uFrom)
        %UNDOFROMUI  Undo or Redo asked for by the Edit menu or a key.
        %   uAct 'undo'|'redo'; uFrom 'menu', 'window' or 'table'. One
        %   keystroke can reach the app twice -- as the menu Accelerator
        %   and as a KeyPressFcn event (of the table, too, when it has
        %   focus) -- and which of those a uifigure delivers is not
        %   documented. So a request from ANOTHER source within half a
        %   second of the last one is taken as the same keystroke: the
        %   same action is not repeated. And as an Accelerator cannot see
        %   Shift, Undo's may also fire for Shift+Cmd+Z: a menu Undo right
        %   after a key Redo is dropped, and a key Redo right after a menu
        %   Undo first takes that Undo back.
        takeBack = false;
        if ~isempty(lastUndoUI) && ~strcmp(lastUndoUI.from, uFrom) && ...
                toc(lastUndoUI.at) < 0.5
            if strcmp(lastUndoUI.act, uAct), return; end
            if strcmp(uAct,'undo') && strcmp(uFrom,'menu'), return; end
            takeBack = strcmp(lastUndoUI.from,'menu') && lastUndoUI.did;
        end
        try
            if takeBack, redoDesign(); end
            if strcmp(uAct,'undo')
                uDid = undoDesign();
            else
                uDid = redoDesign();
            end
        catch uErr
            % Not rethrown: a menu or key callback has nobody to report to.
            uDid = false;
            setStatus(sprintf('%s failed: %s', [upper(uAct(1)) uAct(2:end)], ...
                uErr.message),'bad');
            uialert(fig, uErr.message, 'Undo / Redo');
        end
        lastUndoUI = struct('act',uAct,'from',uFrom,'at',tic,'did',uDid);
    end

    function uDid = undoDesign()
        %UNDODESIGN  Edit > Undo (Cmd+Z): the design before the last step.
        undoCheckpoint();   % anything not yet recorded is a step of its own
        uDid = ~isempty(undoStack);
        if ~uDid, setStatus('Nothing to undo.'); return; end
        uEntry = undoStack{end}; undoStack(end) = [];
        % Undoing an Open also goes back to the file open before it.
        redoStack{end+1} = struct('design',undoBase,'what',uEntry.what, ...
            'doc',docSwap(uEntry.doc));
        restoreDesign(uEntry.design);
        updateUndoMenus();
        undoStatus(['Undid ' uEntry.what], ['Redo with ' undoKey(true)]);
    end

    function uDid = redoDesign()
        %REDODESIGN  Edit > Redo (Shift+Cmd+Z, Cmd+Y): re-apply an undone step.
        undoCheckpoint();   % a change made since the undo clears redo
        uDid = ~isempty(redoStack);
        if ~uDid, setStatus('Nothing to redo.'); return; end
        uEntry = redoStack{end}; redoStack(end) = [];
        undoStack{end+1} = struct('design',undoBase,'what',uEntry.what, ...
            'doc',docSwap(uEntry.doc));
        if numel(undoStack) > UNDO_MAX, undoStack(1) = []; end
        restoreDesign(uEntry.design);
        updateUndoMenus();
        undoStatus(['Redid ' uEntry.what], ['Undo with ' undoKey(false)]);
    end

    function restoreDesign(uSnap)
        %RESTOREDESIGN  Put a snapshot back: S, its controls, one compute.
        %   The controls are set from S as loadConfigCore sets them after a
        %   load, and refreshAll recomputes once. No step may be taken
        %   meanwhile, and the base afterwards is what the restore actually
        %   produced, so the next change is compared with what is on screen.
        undoHold = undoHold + 1;
        try
            for kUndoF = 1:numel(UNDO_FIELDS)
                S.(UNDO_FIELDS{kUndoF}) = uSnap.(UNDO_FIELDS{kUndoF});
            end
            S.sel = []; S.pending = [];   % indices into another element table
            syncDesignControls();
            refreshAll();
            warnImportFrequency();        % as after a load or a frequency edit
        catch restoreErr
            undoHold = max(undoHold - 1, 0);
            undoBase = designSnapshot();
            refreshDocTitle();
            rethrow(restoreErr);
        end
        undoHold = max(undoHold - 1, 0);
        undoBase = designSnapshot();
        refreshDocTitle();
    end

    function syncDesignControls()
        %SYNCDESIGNCONTROLS  Show S's design in every control that edits it.
        %   Setting Value fires no callback. Each value was a valid control
        %   value when its snapshot was taken, so nothing needs clamping;
        %   steering is re-expressed by syncAngleControls because its
        %   spinners follow the angle convention, a view setting that may
        %   have changed since.
        spM.Value = S.M;            spN.Value = S.N;
        spDx.Value = S.dx;          spDy.Value = S.dy;
        ddMode.Value = S.mode;
        ddShape.Value = S.arrayShape;
        syncShapeGallery();
        spGridAngle.Value = S.gridAngle;
        spStagger.Value = S.stagger;
        spStaggerAngle.Value = clampedStaggerAngleDeg(S.stagger, S.dy);
        ddTap.Value = S.taper;
        if ~contains(ddTap.Tooltip, sprintf(' %g dB',S.sll))   % as a load does
            ddTap.Tooltip = sprintf('Chebyshev and Taylor preset: %g dB.',S.sll);
        end
        ddEF.Value = S.efType;      spQ.Value = S.efQ;
        spEfBeamAz.Value = S.efBeamAz; spEfBeamEl.Value = S.efBeamEl;
        syncElementGallery();
        efCustom.Value = S.customFormula;
        efCustomPh.Value = S.customFormulaPh;
        cbSeq.Value = S.seqPhase;
        spSeqBlkM.Value = S.seqBlockM; spSeqBlkN.Value = S.seqBlockN;
        cbUnitCell.Value = S.impUnitCell;
        spTotEff.Value = S.impTotEffPct;
        syncFrequencyMode();        % the Beam squint box
        applyFreqUnit();            % both frequency spinners, the λ readout
        syncAngleControls();        % steering, and a cut plane that follows it
        refreshPolAvailability();   % readouts the element type supports
        refreshImportLabel();       % import button and readout
    end

    function uWhat = changeName(uA, uB)
        %CHANGENAME  What changed between two designs, for "Undo <what>".
        %   The first group that differs names it. Causes come before their
        %   consequences: an import also switches the element type, a new M
        %   or a shape also rebuilds the element table.
        groups = { ...
            {'impFF','impFFName','impFFGHz','impFFFreqFrom'},  'import'
            {'M','N'},                                         'grid size'
            {'dx','dy','subOff'},                              'spacing'
            {'arrayShape'},                                    'shape'
            {'gridAngle','stagger'},                           'lattice'
            {'mode'},                                          'placement mode'
            {'taper','sll'},                                   'taper'
            {'efType','efQ','efBeamAz','efBeamEl','customFormula','customFormulaPh', ...
             'dualFeedCP','impNeedsPattern'},                  'element type'
            {'impUnitCell','impTotEffPct'},                    'import settings'
            {'theta_s','phi_s'},                               'steering'
            {'freqGHz','freqOpGHz','squintMode','retunePhase'}, 'frequency'
            {'seqPhase','seqBlockM','seqBlockN'},              'sequential rotation'
            {'portMap','portMapSet'},                          'feed layout'};
        nA = size(uA.el,1); nB = size(uB.el,1);
        for kGroup = 1:size(groups,1)
            if any(cellfun(@(fld) ~isequaln(uA.(fld),uB.(fld)), groups{kGroup,1}))
                uWhat = groups{kGroup,2};
                % Typing an amplitude switches the taper to Manual (table):
                % that is the table edit, not a taper choice.
                if strcmp(uWhat,'taper') && strcmp(uB.taper,'Manual (table)') ...
                        && nA == nB && ~isequal(uA.el,uB.el)
                    uWhat = 'table edit';
                end
                return;
            end
        end
        % Only the element table changed.
        if nB == 0,                             uWhat = 'clear';
        elseif nB < nA,                         uWhat = 'delete';
        elseif nB > nA,                         uWhat = 'add element';
        elseif ~isequal(uA.el(:,5),uB.el(:,5)), uWhat = 'rotation';
        elseif ~isequal(uA.el(:,1:2),uB.el(:,1:2)), uWhat = 'element move';
        else,                                   uWhat = 'table edit';
        end
    end

    function updateUndoMenus()
        % Enabled only with a step to take, and named after it. Set only on
        % a change: this runs with every compute.
        if isempty(mUndo) || ~isgraphics(mUndo), return; end
        setUndoMenu(mUndo, 'Undo', undoStack);
        setUndoMenu(mRedo, 'Redo', redoStack);
    end

    function setUndoMenu(hMenu, verb, uStack)
        if isempty(uStack)
            txt = verb; onOff = 'off';
        else
            txt = [verb ' ' uStack{end}.what]; onOff = 'on';
        end
        if ~strcmp(hMenu.Text, txt), hMenu.Text = txt; end
        if ~strcmp(char(hMenu.Enable), onOff), hMenu.Enable = onOff; end
    end

    function undoStatus(msg, hint)
        % msg and the key hint on the status line. A warning the compute
        % just left (no active elements, plots out of date with Auto off,
        % an imported pattern missing) stays readable between the two: it
        % is why the plots are empty.
        if any(strcmp(statusTone, {'warn','bad'}))
            why = statusMsg;
            if ~endsWith(why, {'.','!','?'}), why = [why '.']; end
            setStatus(sprintf('%s — %s %s', msg, why, hint), statusTone);
        else
            setStatus(sprintf('%s — %s', msg, hint));
        end
    end

    function announceUndoable(stepsWas, msg)
        %ANNOUNCEUNDOABLE  "<msg> — Undo with ⌘Z", if the action made a step.
        if undoCount > stepsWas
            undoStatus(msg, ['Undo with ' undoKey(false)]);
        end
    end

    function announceAllElements(stepsWas, fmt)
        % "Apply to selected" and Move act on EVERY element when none is
        % selected. Legitimate (it is how a whole array is rotated), but
        % not what the caption says, so the status line spells it out.
        % fmt holds one %s for "all N elements".
        if isempty(S.sel)
            nAll = size(S.el,1);
            announceUndoable(stepsWas, sprintf([fmt ' (none were selected)'], ...
                sprintf('all %d %s', nAll, plural(nAll,'element'))));
        end
    end

    function keyTxt = undoKey(isRedo)
        % The shortcut as this platform writes it.
        if ismac
            if isRedo, keyTxt = '⇧⌘Z'; else, keyTxt = '⌘Z'; end
        elseif isRedo
            keyTxt = 'Ctrl+Y';
        else
            keyTxt = 'Ctrl+Z';
        end
    end

    function wordOut = plural(nCount, wordOut)
        if nCount ~= 1, wordOut = [wordOut 's']; end
    end

    function loadConfig(cfgIn)
        % Roll back both state and editable controls if a load fails after
        % validation. A half-loaded design must never remain live.
        % The whole load is ONE undo step ('design load'), taken when this
        % returns: loadConfigCore computes part-way through, and a
        % rolled-back load must leave no step at all.
        % cfgIn: the file to open (Open Recent); omitted, a file is asked for.
        if nargin < 1, cfgIn = ''; end
        undoTok = holdUndo('design load'); %#ok<NASGU>
        beforeState=S;
        beforeHandles=findall(fig,'-property','Value');
        beforeProps=cell(numel(beforeHandles),1);
        for lc=1:numel(beforeHandles)
            % struct('Value',[]) creates an empty struct array in MATLAB.
            % Some open result controls legitimately have an empty Value.
            beforeProps{lc}=struct();
            beforeProps{lc}.Value=beforeHandles(lc).Value;
            if isprop(beforeHandles(lc),'Limits'), beforeProps{lc}.Limits=beforeHandles(lc).Limits; end
            if isprop(beforeHandles(lc),'Items'), beforeProps{lc}.Items=beforeHandles(lc).Items; end
            if isprop(beforeHandles(lc),'Enable'), beforeProps{lc}.Enable=beforeHandles(lc).Enable; end
        end
        try
            [cfgLoaded, cfgFile, cfgNote] = loadConfigCore(cfgIn);
        catch loadErr
            S=beforeState;
            for lc=1:numel(beforeHandles)
                hRestore=beforeHandles(lc); oldProps=beforeProps{lc};
                if ~isgraphics(hRestore), continue; end
                try
                    if isfield(oldProps,'Limits'), hRestore.Limits=[-Inf Inf]; end
                    if isfield(oldProps,'Items'), hRestore.Items=oldProps.Items; end
                    hRestore.Value=oldProps.Value;
                    if isfield(oldProps,'Limits'), hRestore.Limits=oldProps.Limits; end
                    if isfield(oldProps,'Enable'), hRestore.Enable=oldProps.Enable; end
                catch
                end
            end
            try
                applyFreqUnit(); syncFrequencyMode(); refreshPolAvailability();
                refreshImportLabel(); refreshAll();
            catch
                invalidateRadiation('Configuration load failed; previous settings restored.');
            end
            uialert(fig,loadErr.message,'Configuration rejected - previous settings restored');
            setStatus(['Design rejected, previous settings restored: ' ...
                loadErr.message],'bad');
            return;
        end
        if ~cfgLoaded, return; end
        % The opened file becomes the open design. Set before this function
        % returns, which is when the load becomes an undo step: the step
        % carries the file that was open before (docNextStep), so Undo
        % returns to it rather than leaving this file's name on the old
        % design -- where Save would then overwrite this file with it.
        docNextStep = struct('file',docFile,'saved',docSaved);
        docFile = cfgFile; docSaved = designSnapshot();
        rememberFile(cfgFile);
        refreshDocTitle();
        docStatus(['Opened ' fileLabel(cfgFile) cfgNote], true);
    end

    function [loaded, cfgFile, cfgNote] = loadConfigCore(cfgIn)
        % loaded is true only when a design was actually applied; a
        % cancelled question or a refused file leaves the open file alone.
        loaded = false; cfgNote = '';
        if ~isempty(cfgIn)
            cfgFile = cfgIn;
        else
            if isempty(docFile), openDir = lastFolder(); else, openDir = fileparts(docFile); end
            [f,p] = uigetfile('*.mat','Open design',[openDir filesep]);
            figure(fig);   % see loadImportedFF's comment -- native file dialogs
                           % can leave the uifigure behind other windows on macOS
            if isequal(f,0), cfgFile = ''; return; end
            cfgFile = fullfile(p,f);
        end
        try
            L = load(cfgFile,'cfg');
            cfg = L.cfg;
        catch err
            setStatus(['Open failed: ' err.message],'bad');
            uialert(fig, err.message, 'Load failed');
            return;
        end
        % Older configs retain geometry/settings but must re-import their
        % pattern: they did not record a verified quantity/coverage schema.
        if isstruct(cfg)&&isscalar(cfg)&&isfield(cfg,'impFF')&& ...
                isstruct(cfg.impFF)&&isscalar(cfg.impFF)&&~isempty(cfg.impFF)&& ...
                ~isfield(cfg.impFF,'schemaVersion')
            cfg.impFF=[];
            if isfield(cfg,'efType') && strcmp(cfg.efType,'Imported (CST far-field)')
                cfg.impNeedsPattern=true;
            end
        end
        validateConfigCandidate(cfg);
        cfg=normalizeConfigNumbers(cfg);
        % In old files an inactive operating field did not affect the design.
        % Do not suddenly activate its stale hidden value during migration.
        if ~isfield(cfg,'steeringSchema') && isfield(cfg,'squintMode') && ...
                ~cfg.squintMode && isfield(cfg,'freqGHz')
            cfg.freqOpGHz=cfg.freqGHz;
        end
        if ~isfield(cfg,'retunePhase'), cfg.retunePhase=true; end
        % SHAPE first. Everything below does isfield(cfg,...) and cfg.(fn),
        % which need a scalar struct: a .mat carrying a struct ARRAY, a
        % cell, or a plain number under the name 'cfg' reached those lines
        % and threw partway through -- after some of S had already been
        % overwritten, leaving a half-loaded design with no way back.
        if ~isstruct(cfg) || ~isscalar(cfg)
            uialert(fig, ['This file does not contain a design.' newline newline ...
                'It stores a variable named "cfg", but not the single settings ' ...
                'structure this app saves. Nothing has been changed.'], ...
                'Load failed');
            return;
        end
        if isfield(cfg,'el') && isnumeric(cfg.el) && ismatrix(cfg.el) && ...
                size(cfg.el,2) == 5 && isreal(cfg.el) && ...
                any(abs(cfg.el(:,3)) > 1e6)
            uialert(fig,'Saved amplitudes exceed +/-1e6. Scale the weights before loading. Nothing has been changed.', ...
                'Amplitude out of range');
            return;
        end
        % Field-by-field with isfield checks: a config saved by an older
        % or newer version of this app, missing a field, must not crash
        % the whole load -- missing fields simply keep their current value.
        missing = {};
        fields = {'M','N','dx','dy','subOff','theta_s','phi_s','mode', ...
            'el','elRC','taper','sll','efType','efQ','efBeamAz','efBeamEl', ...
            'seqBlockM','seqBlockN','seqPhase', ...
            'dualFeedCP','customFormula','customFormulaPh', ...
            'arrayShape','gridAngle','stagger','freqGHz','squintMode','retunePhase','freqOpGHz','freqUnit', ...
            'cutMode','cutFixedTheta','fullSphere','cutPhi','cutPhiFollow','impUnitCell', ...
            'angleConvention', ...
            'impTotEffPct'};
        % ---- validate the CANDIDATE before anything reaches S ----------
        % Everything here is checked against the geometry that will be
        % loaded WITH it, and el/elRC are treated as one inseparable
        % object. The earlier version checked them only when both were
        % present -- so a file carrying just one could replace half the
        % pair -- and it dropped the element table while still loading
        % the new M/N, leaving indices that pointed outside the lattice
        % they now belonged to.
        %
        % M and N are also settled here, before any rounding: round() on
        % a cell array or a string throws, and that throw would land
        % after S had already been partly overwritten.
        rebuildFromSize = false;
        % Snapshot of everything rebuildUniform() reads to place elements,
        % taken BEFORE the field loop overwrites it. Compared against the
        % same quantities afterwards to decide whether a geometry-only
        % config has moved the lattice -- see the geometry check below.
        % A snapshot is used rather than candidate values pulled out of
        % cfg because these fields pass through several validators on
        % their way into S, and what matters is where they LANDED.
        geoWas = [S.M S.N S.dx S.dy S.gridAngle S.stagger];
        modeWas = S.mode;
        shapeWas = S.arrayShape;
        taperWas = S.taper;
        sllWas = S.sll;
        cM = pickNum(cfg,'M',S.M); cN = pickNum(cfg,'N',S.N);
        cM = max(1,min(64,round(cM))); cN = max(1,min(64,round(cN)));
        hasEl = isfield(cfg,'el'); hasRC = isfield(cfg,'elRC');
        if ~hasEl && ~hasRC && (cM ~= S.M || cN ~= S.N)
            % Geometry-only file. Keeping the current elements next to a
            % different M/N leaves lattice indices pointing at cells that
            % no longer exist -- the same inconsistency the element-table
            % check prevents, arriving by the one route that skipped it.
            % Rebuilding from the loaded size is the only reading of a
            % geometry-only config that stays self-consistent.
            rebuildFromSize = true;
        end
        % Normalised to the column counts every reader assumes, so an
        % empty table stored as [] cannot become a 0x0 that fails the
        % first time something indexes column 5.
        if isfield(cfg,'el') && isnumeric(cfg.el) && isempty(cfg.el), cfg.el = zeros(0,5); end
        if isfield(cfg,'elRC') && isnumeric(cfg.elRC) && isempty(cfg.elRC), cfg.elRC = zeros(0,2); end
        if hasEl || hasRC
            bothThere = hasEl && hasRC;
            if ~bothThere || ~okElements(cfg.el, cfg.elRC, cM, cN)
                if ~bothThere
                    why = 'the file carries only one half of the element table';
                else
                    why = sprintf(['the element table has the wrong shape, ' ...
                           'non-finite values, positions outside +/-%g wavelengths, ' ...
                           'mismatched row counts, or row/column indices outside ' ...
                           'its saved array size'],MAX_POS_LAMBDA);
                end
                uialert(fig, sprintf(['This configuration''s elements could not ' ...
                    'be loaded: %s.\n\nThe rest of the file has been applied ' ...
                    'and the array has been rebuilt from its saved size, so ' ...
                    'the elements on screen match the geometry rather than ' ...
                    'being left over from the previous design.'], why), ...
                    'Element table rejected','Icon','warning');
                if hasEl,  cfg = rmfield(cfg,'el');   end
                if hasRC,  cfg = rmfield(cfg,'elRC'); end
                rebuildFromSize = true;
            end
        end
        cfg.M = cM; cfg.N = cN;
        % Frequencies reach applyFreqUnit and the spinners directly, so
        % they are screened here rather than by the widget validators.
        for fq = {'freqGHz','freqOpGHz'}
            if isfield(cfg,fq{1})
                v = cfg.(fq{1});
                if ~(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) ...
                        && v >= 1e-3 && v <= 1e6)
                    % Removing the field is enough: the fields loop below
                    % then finds it absent and reports it. Appending here as
                    % well listed it twice -- "freqGHz, freqGHz" in the
                    % load dialog.
                    cfg = rmfield(cfg,fq{1});
                end
            end
        end
        for i = 1:numel(fields)
            fn = fields{i};
            if isfield(cfg,fn)
                S.(fn) = cfg.(fn);
            else
                missing{end+1} = fn; %#ok<AGROW>
            end
        end
        % M/N were settled in the candidate check above, so nothing more
        % is needed here.
        S.sel = []; S.pending = [];
        % Design targets: a file without them (older, or never given any)
        % specifies none. Keeping this session's targets would judge the
        % loaded design against another design's spec without saying so.
        % Validated above; sanitizeTargets clamps them to the window's
        % ranges. refreshAll below re-judges the cards.
        if isfield(cfg,'targets')
            S.targets = sanitizeTargets(cfg.targets);
        else
            S.targets = sanitizeTargets(struct());
        end
        % The aim memory belongs to the design that was on screen, not to
        % the one being loaded. Left set, a loaded steer angle that
        % happened to equal the old command would be read as "already
        % aimed" and silently re-aim at the PREVIOUS design's target.
        % This has to sit here, in loadConfig: the first attempt put it
        % in rebuildUniform, which a config load never calls -- it
        % installs the saved element table directly and calls refreshAll.

        % reflect the loaded values back into every widget -- S alone
        % doesn't refresh what's on screen, the controls must be told too
        % Every value here came straight out of a .mat file, so nothing
        % guarantees it sits inside the current widgets' Limits/Items --
        % a config written by a different version of this app, or one
        % hand-edited, can carry anything. Assigning an out-of-range
        % Value to a uispinner (or a string that isn't in a dropdown's
        % Items) THROWS, and only the load() call above is inside a
        % try/catch: the throw would land partway through this block,
        % AFTER S was already overwritten, leaving S describing one
        % array while the widgets showed a half-updated different one
        % and refreshAll() below never running at all. The field-by-
        % field isfield() guard above already protects against MISSING
        % fields for exactly this version-skew reason; this protects
        % against INVALID ones. Each value is clamped/rejected first and
        % the corrected value written BACK into S, so state and UI
        % cannot drift apart.
        S.M = cfgNum(spM,S.M);                     spM.Value = S.M;
        S.N = cfgNum(spN,S.N);                     spN.Value = S.N;
        S.dx = cfgNum(spDx,S.dx);                  spDx.Value = S.dx;
        S.dy = cfgNum(spDy,S.dy);                  spDy.Value = S.dy;
        S.subOff = min(max(S.subOff,0.01),5);
        if ~(isnumeric(S.theta_s) && isscalar(S.theta_s) && ...
                isreal(S.theta_s) && isfinite(S.theta_s)), S.theta_s = 0; end
        if ~(isnumeric(S.phi_s) && isscalar(S.phi_s) && ...
                isreal(S.phi_s) && isfinite(S.phi_s)), S.phi_s = 0; end
        S.theta_s = min(max(double(S.theta_s),-90),90);
        S.phi_s = min(max(double(S.phi_s),-180),360);
        S.angleConvention = cfgItem(ddAngleConvention,S.angleConvention);
        S.mode = cfgItem(ddMode,S.mode);           ddMode.Value = S.mode;
        S.taper = cfgItem(ddTap,S.taper);          ddTap.Value = S.taper;
        % Clamped to the spinner's range, and the spinner set here --
        % before applyTaper below builds the amplitudes from S.sll -- so
        % the control and the weights cannot disagree. An older file
        % without 'sll' keeps the current level.
        S.sll = min(max(double(S.sll),13),80);
        syncSLLControl();
        S.efType = cfgItem(ddEF,S.efType);         ddEF.Value = S.efType;
        S.efQ = cfgNum(spQ,S.efQ);                 spQ.Value = S.efQ;
        S.efBeamAz = cfgNum(spEfBeamAz,S.efBeamAz); spEfBeamAz.Value = S.efBeamAz;
        S.efBeamEl = cfgNum(spEfBeamEl,S.efBeamEl); spEfBeamEl.Value = S.efBeamEl;
        syncElementGallery();
        S.seqPhase = cfgBool(cbSeq,S.seqPhase);    cbSeq.Value = S.seqPhase;
        % Cleared before the restore decides. Left over from an earlier
        % failed load, the flag survived a subsequent GOOD load and kept
        % stamping "NEEDS AN IMPORTED PATTERN" on a design that had one.
        S.impNeedsPattern = false;

        % ---- imported far-field -------------------------------------
        % Restored BEFORE the element type is acted on, so an
        % 'Imported (CST far-field)' design finds its own pattern rather
        % than whatever happened to be loaded in this session. Validated
        % on the fields every downstream reader indexes; anything short
        % of that is rejected whole and reported, never half-installed.
        % Every field the app later READS, with its type and shape --
        % not just the ones that came to mind. The first list omitted
        % eff (refreshImportLabel uses it immediately), compNames
        % (impBasisOK indexes compNames{1}), and thMin/thMax (the field
        % evaluator clamps to them), so a struct could pass validation,
        % replace S.impFF, and then throw while being displayed.
        ffOK = isfield(cfg,'impFF') && validSavedFF(cfg.impFF);
        if isfield(cfg,'impFF')
            if ffOK
                S.impFF = refreshImportedPowerMetadata(cfg.impFF);
                if isfield(cfg,'impFFName') && (ischar(cfg.impFFName) || isstring(cfg.impFFName))
                    S.impFFName = char(cfg.impFFName);
                else
                    S.impFFName = '(from saved config)';
                end
                refreshImportLabel();
            elseif isempty(cfg.impFF)
                % The file says, explicitly, that this design carries NO
                % imported pattern -- so this session's pattern, loaded by
                % hand for some OTHER design, must not survive the load.
                % It used to: nothing assigned S.impFF on this path, so the
                % import button went on naming the previous file, switching
                % the element dropdown to Imported silently used a stranger's
                % pattern, and -- because saveConfig persists S.impFF -- the
                % next save wrote that stranger's pattern INTO this design's
                % file. That is the "keyed on what the FILE carries" rule
                % stated below, which this one case did not follow.
                % Also reached by every pre-schema config that held a
                % pattern, since the migration above empties it on purpose.
                %
                % (The branch this replaces, for a non-empty malformed
                % pattern, could never run: validateConfigCandidate already
                % rejects the whole config in that case, deliberately --
                % see its "strict realized-gain schema" message.)
                S.impFF = []; S.impFFName = '';
                refreshImportLabel();
            end
        end
        % A design that SAYS it uses an imported pattern but carries none
        % must not quietly become isotropic, and must not adopt a
        % stranger's pattern either -- both change the answer silently,
        % and the second also flips seqRotSign() and the whole rotation
        % compensation with it.
        % Keyed on what the FILE carries, not on what this session happens
        % to have loaded. Testing isempty(S.impFF) meant a config that
        % named an imported element but stored none passed silently
        % whenever any unrelated pattern was already open -- which is
        % precisely the "reopen while another pattern is loaded" failure,
        % the one that also flips seqRotSign() underneath the design.
        % A file that recorded an unresolved substitution keeps saying so.
        if isfield(cfg,'impNeedsPattern') && isscalar(cfg.impNeedsPattern) && ...
                (islogical(cfg.impNeedsPattern) || isnumeric(cfg.impNeedsPattern)) && ...
                isreal(cfg.impNeedsPattern) && isfinite(cfg.impNeedsPattern) && ...
                cfg.impNeedsPattern == 1
            S.impNeedsPattern = true;
        end

        % Keyed on whether restoration actually SUCCEEDED (ffOK), not on
        % the field merely being present and non-empty -- a malformed
        % struct satisfied the latter while having been rejected above.
        if isfield(cfg,'efType') && strcmp(S.efType,'Imported (CST far-field)') ...
                && ~ffOK
            missing{end+1} = ['efType (this design uses an imported pattern, but the ' ...
                              'file carries none - import it before computing)'];
            % BLOCKED, not merely warned. Leaving efType on Imported let
            % the design compute against whatever pattern this session
            % happened to have open -- a stranger's -- or fall back to
            % isotropic in a fresh one, and both silently change the
            % answer AND seqRotSign() with it. Clearing the stale pattern
            % and dropping to a named built-in makes the substitution
            % visible in the Element factor box instead of invisible.
            S.impFF = []; S.impFFName = '';
            S.efType = 'Isotropic';
            S.impNeedsPattern = true;
            % The widgets are synced EARLIER in this function than this
            % branch runs, so they have to be corrected here as well --
            % otherwise S says Isotropic while the dropdown still reads
            % "Imported (CST far-field)" and the button still names the
            % stranger's file, which is the same lie in a new place.
            ddEF.Value = S.efType;
            refreshPolAvailability();
            refreshImportLabel();
        end

        % Restored only when structurally intact -- a hand-edited or
        % older config must not install a half-built layout that the map
        % would then index into.
        % CONTENTS and TYPES, not only field names and sizes. The first
        % version checked isfield + size, which a hand-edited or foreign
        % .mat passes while carrying a cell array for hOdd, a side code
        % of 7, or a 1x1 cell for portMapSet that then reaches logical()
        % and throws. Everything here is indexed into on every redraw, so
        % a malformed unit has to be rejected whole rather than installed
        % and discovered later.
        if isfield(cfg,'portMap') && isstruct(cfg.portMap) && isscalar(cfg.portMap) && ...
                all(isfield(cfg.portMap,{'hOdd','hSide','vSide'}))
            pmC = cfg.portMap;
            okPM = (islogical(pmC.hOdd) || (isnumeric(pmC.hOdd) && all(ismember(pmC.hOdd(:),[0 1])))) ...
                && isnumeric(pmC.hSide) && isnumeric(pmC.vSide) ...
                && ~isempty(pmC.hOdd) && ismatrix(pmC.hOdd) ...
                && isequal(size(pmC.hOdd),size(pmC.hSide)) ...
                && isequal(size(pmC.hOdd),size(pmC.vSide)) ...
                && size(pmC.hOdd,1) <= 8 && size(pmC.hOdd,2) <= 8 ...
                && all(ismember(pmC.hSide(:),[1 2])) && all(ismember(pmC.vSide(:),[1 2]));
            if okPM
                % A finite, real 0 or 1 -- not "any numeric scalar then
                % logical() it", which accepted 7, NaN and 2+3i.
                pmsOK = isfield(cfg,'portMapSet') && isscalar(cfg.portMapSet) && ...
                        (islogical(cfg.portMapSet) || ...
                         (isnumeric(cfg.portMapSet) && isreal(cfg.portMapSet) && ...
                          isfinite(cfg.portMapSet) && ismember(double(cfg.portMapSet),[0 1])));
                if pmsOK && logical(cfg.portMapSet)
                    % Explicit user declaration: preserve it exactly.
                    S.portMap = struct('hOdd',logical(pmC.hOdd), ...
                                       'hSide',double(pmC.hSide), ...
                                       'vSide',double(pmC.vSide));
                    S.portMapSet = true;
                else
                    % The saved map was merely that older app version's
                    % built-in default. A config with portMapSet=false must
                    % follow THIS version's default rather than silently
                    % resurrecting a historical default and still labelling
                    % it "Default" in the phase-map window.
                    S.portMap = DEFAULT_PORT_MAP;
                    S.portMapSet = false;
                end
            else
                % Do not inherit the feed map from the previously open design.
                % A malformed saved map cannot safely identify the CST ports,
                % so fall back to this version's explicit built-in default.
                S.portMap = DEFAULT_PORT_MAP;
                S.portMapSet = false;
                missing{end+1} = 'portMap (malformed - reset to the current default feed layout)';
            end
        elseif isfield(cfg,'portMap')
            % Present but not even shaped like a portMap (not a struct, a
            % struct ARRAY, or missing fields). Do not silently reuse the
            % previous design's map: that can export the wrong physical ports.
            S.portMap = DEFAULT_PORT_MAP;
            S.portMapSet = false;
            missing{end+1} = 'portMap (unrecognised - reset to the current default feed layout)';
        else
            % Older configurations predate the declared physical feed map.
            % They must start from the current documented default, not inherit
            % a custom map left in memory by the design loaded before them.
            S.portMap = DEFAULT_PORT_MAP;
            S.portMapSet = false;
            missing{end+1} = 'portMap (not stored by this older config - using the current default feed layout)';
        end
        if S.dualFeedCP && strcmp(S.efType,'Patch (cos^q x lin pol)')
            missing{end+1}='Built-in patch CP removed: patch is now linear; import a CP field if required';
        end
        S.dualFeedCP = false;
        S.seqBlockM = round(cfgNum(spSeqBlkM,S.seqBlockM)); spSeqBlkM.Value = S.seqBlockM;
        S.seqBlockN = round(cfgNum(spSeqBlkN,S.seqBlockN)); spSeqBlkN.Value = S.seqBlockN;
        S.customFormula = cfgText(efCustom,S.customFormula);        efCustom.Value = S.customFormula;
        S.customFormulaPh = cfgText(efCustomPh,S.customFormulaPh);  efCustomPh.Value = S.customFormulaPh;
        S.arrayShape = cfgItem(ddShape,S.arrayShape);      ddShape.Value = S.arrayShape;
        S.gridAngle = cfgNum(spGridAngle,S.gridAngle);     spGridAngle.Value = S.gridAngle;
        S.stagger = cfgNum(spStagger,S.stagger);           spStagger.Value = S.stagger;
        spStaggerAngle.Value = clampedStaggerAngleDeg(S.stagger,S.dy);
        % Read the raw GHz value, not the spinner: at this point the
        % spinner is still showing the OLD unit, so cfgNum's fallback
        % would fold a stale scale factor into the stored frequency.
        if isfield(cfg,'freqGHz') && isnumeric(cfg.freqGHz) && isscalar(cfg.freqGHz)
            S.freqGHz = cfg.freqGHz;
        end
        S.squintMode = logical(S.squintMode);
        S.retunePhase = logical(S.retunePhase);
        syncFrequencyMode();
        if isfield(cfg,'freqOpGHz') && isnumeric(cfg.freqOpGHz) && isscalar(cfg.freqOpGHz)
            S.freqOpGHz = cfg.freqOpGHz;
        end
        if isfield(cfg,'freqUnit') && any(strcmp(cfg.freqUnit,{'GHz','MHz'}))
            S.freqUnit = cfg.freqUnit;
        end
        ddFreqUnit.Value = S.freqUnit;
        applyFreqUnit();   % rescales both spinners and the lambda readout
        showFreqOp(S.squintMode);
        % Limits BEFORE Value, and the value clamped into range: a config
        % saved in full-sphere mode can carry cutFixedTheta > 90, which
        % MATLAB rejects outright against the hemisphere-mode [0 90]
        % limits. Setting limits first (and clamping) makes either
        % direction of transition safe.
        S.fullSphere = cfgBool(cbFullSphere,S.fullSphere);
        cbFullSphere.Value = S.fullSphere;
        S.cutMode = cfgItem(ddCutMode,S.cutMode);
        % Screened before the arithmetic: min/max on a cell array or a
        % vector throws, and that throw would land after S was already
        % partly overwritten.
        if ~(isnumeric(S.cutFixedTheta) && isscalar(S.cutFixedTheta) && ...
                isreal(S.cutFixedTheta) && isfinite(S.cutFixedTheta))
            S.cutFixedTheta = 0;
        end
        cutMax = 90;
        if S.fullSphere, cutMax = 180; end
        S.cutFixedTheta = min(max(S.cutFixedTheta,0),cutMax);
        % Cut plane. Restored through the same cfg* helpers, so a config
        % predating these fields keeps the current (linked) behaviour
        % rather than erroring or silently unlinking.
        S.cutPhiFollow = cfgBool(cbCutFollow,S.cutPhiFollow);
        S.cutPhi = cfgNum(spCutPhi,S.cutPhi);   % cfgNum already clamps to Limits
        S.impUnitCell = cfgBool(cbUnitCell,S.impUnitCell); cbUnitCell.Value = S.impUnitCell;
        S.impTotEffPct = cfgNum(spTotEff,S.impTotEffPct); spTotEff.Value = S.impTotEffPct;
        cbCutFollow.Value = S.cutPhiFollow;
        if S.cutPhiFollow, spCutPhi.Enable = 'off'; else, spCutPhi.Enable = 'on'; end
        syncAngleControls();

        % The ribbon's shape gallery is a SECOND view of S.arrayShape, and
        % this path restores the shape by writing ddShape.Value directly
        % rather than going through assignShape -- which is where the
        % gallery highlight is kept in step. Without this, loading a
        % config left the dropdown reading the loaded shape while the
        % ribbon still lit the previous one. (R10 did not catch it: that
        % test drives the dropdown through its callback, which config
        % load deliberately bypasses.)
        % An element table that was rejected leaves S.el describing the
        % PREVIOUS design while M/N now describe this one. Rebuilding
        % from the loaded size puts the two back in step instead of
        % leaving a mismatched pair live.
        % A geometry-only config (no el, no elRC) that moved ANY of the
        % quantities the lattice is built from has to rebuild, not just
        % one that changed M or N. The M/N test alone left a config with
        % the same 4x4 size but a different dx loading the new spacing
        % into the spinner while every element kept its old coordinate --
        % the panel then said 1.25 lambda while the array it plotted,
        % exported and sent to CST was still on 0.5. Wrong in the same
        % way the M/N case was, and quieter, because nothing about the
        % element count looks out of place.
        %
        % Guarded on ~hasEl && ~hasRC, which is the whole point: a config
        % that DOES carry an element table is describing hand-placed
        % positions, and rebuilding would throw them away.
        if ~hasEl && ~hasRC && ...
                (~isequal([S.M S.N S.dx S.dy S.gridAngle S.stagger], geoWas) || ...
                 ~isequal(S.arrayShape, shapeWas) || ...
                 xor(strcmp(S.mode,'Uniform grid'),strcmp(modeWas,'Uniform grid')))
            rebuildFromSize = true;
        end
        if rebuildFromSize
            rebuildUniform();
        elseif ~hasEl && ~hasRC && ...
                (~strcmp(S.taper,taperWas) || S.sll ~= sllWas)
            % A settings-only file can change the taper without changing
            % geometry. Keep the excitation table in step with its label.
            applyTaper();
        end

        syncShapeGallery();

        % Run ONCE here and again after the view block. Two different
        % jobs, and both are needed:
        %
        %   here  -- rebuild the list of CHOICES for the element type the
        %            config just loaded, so cfgItem() below validates the
        %            saved view against the right list;
        %   after -- ENFORCE, so a saved view the loaded state cannot
        %            support is cleared rather than restored.
        %
        % With only the second call, loading a perfectly valid saved
        % Patch design while a magnitude-only import happened to be live
        % lost the saved Axial Ratio view outright: the AR entry was
        % still trimmed out of ddShow.Items when cfgItem() ran, so it
        % rejected 'Axial Ratio (dB)' as an unknown item and fell back to
        % the current value; the later call then put the entry back into
        % the list, but nothing put the selection back. The user's own
        % saved view vanished silently on load.
        refreshPolAvailability();
        refreshTable();

        % ---- view settings ----------------------------------------------
        % Restored last, through the same cfg* validators as everything
        % else, and each guarded by isfield: a config written before these
        % were saved simply leaves the current view untouched rather than
        % erroring or snapping back to defaults. Assigned directly to the
        % widgets because these controls have no mirror in S -- the rest
        % of the app reads them off the widget. No callback is fired here;
        % refreshAll() below recomputes once for the whole load, which is
        % also why this cannot sit in the generic field loop above.
        if isfield(cfg,'vShow'),     ddShow.Value  = cfgItem(ddShow,cfg.vShow);      end
        if isfield(cfg,'vPol'),      ddPol.Value   = cfgItem(ddPol,cfg.vPol);        end
        if isfield(cfg,'vScale'),    ddScale.Value = cfgItem(ddScale,cfg.vScale);    end
        if isfield(cfg,'vDynRange'), spDR.Value    = cfgNum(spDR,cfg.vDynRange);     end
        if isfield(cfg,'vAbsLevel'), cbAbs.Value   = cfgBool(cbAbs,cfg.vAbsLevel);   end
        if isfield(cfg,'vARCut'),    cbARCut.Value = cfgBool(cbARCut,cfg.vARCut);    end
        if isfield(cfg,'vAuto'),     cbAuto.Value  = cfgBool(cbAuto,cfg.vAuto);      end

        % The imported pattern's frequency goes with whichever pattern the
        % steps above left in place, and the readout compares it with the
        % frequencies loaded above, so both come after them.
        restoreImportFreq(cfg);
        refreshImportLabel();

        % A loaded design can name an element type or restore a view
        % selection the current import cannot support, and it reaches
        % neither the element-type callback nor the import path.
        %
        % LAST, after the view settings above -- not before them, where
        % this used to sit. Run first, it cleared the unsupported
        % selections and then watched the restore block put them straight
        % back: a config saved with vPol = 'RHCP component' and vARCut =
        % true came back with the polarization dropdown reading RHCP and
        % the AR checkbox ticked while both controls were greyed out, so
        % the plots showed polarization readouts the imported file cannot
        % support and the disabled controls gave the user no way to see
        % why. A guard has to be the last word on the state it guards.
        refreshPolAvailability();

        refreshAll();
        warnImportFrequency();
        % Both notices merged into one dialog. rebuildUniform() is what
        % normally fires the large-array notice; a loaded config sets
        % S.el directly and bypasses it, so it is checked again here --
        % and a config can legitimately raise BOTH (version skew leaves
        % fields missing while the array itself is large). Since uialert
        % does not block, showing them separately would let the second
        % overwrite the first.
        loadMsgs = {};
        if ~isempty(missing)
            loadMsgs{end+1} = sprintf(['Config loaded, but these missing fields kept ' ...
                'their current values: %s'], strjoin(missing,', ')); 
        end
        bigMsg = largeArrayNotice('Loaded ');
        if ~isempty(bigMsg), loadMsgs{end+1} = bigMsg; end 
        if isscalar(loadMsgs)
            uialert(fig, loadMsgs{1}, 'Config loaded', 'Icon','warning');
        elseif numel(loadMsgs) > 1
            uialert(fig, strjoin(loadMsgs, sprintf('\n\n')), 'Config loaded', 'Icon','warning');
        end
        % When and by which version the file was saved, if it says; to the
        % minute (the status line is one line).
        if isfield(cfg,'savedTime') && validText(cfg.savedTime)
            savedWhen = char(cfg.savedTime);
            cfgNote = [' (saved ' savedWhen(1:min(end,16))];
            if isfield(cfg,'appVersion') && validText(cfg.appVersion)
                cfgNote = [cfgNote ', v' char(cfg.appVersion)];
            end
            cfgNote = [cfgNote ')'];
        end
        loaded = true;
    end

    function result = showScanLossCore(coverageSamples)
        result = [];
        coverage = nargin>0;
        if ~patternReady(true), return; end
        normalizedAmp = calculationWeights();
        scanPhi = scanPhiVal();
        % FIRST statement in the function, ahead of every element
        % evaluation this analysis performs -- the rotation cache below,
        % seqRotSign(), and the per-angle point evaluations alike.
        %
        % It sat after the cache in the previous version, with a comment
        % claiming otherwise, and that ordering silently ERASED real
        % failures rather than merely inheriting stale ones. The failing
        % case is exact: E_theta = 1./(th-89) is singular at 89 deg,
        % which is on the cache's full grid but not on the handedness
        % ring (10 deg) nor the sweep's own points (0:2:80). So the cache
        % substituted ones(), the reset wiped the flag, the point
        % evaluations succeeded, and the plot mixed substituted cached
        % fields with valid ones under no warning at all.
        %
        % Reset here and NOWHERE else in this function.
        S.efFallback = false; S.efFallbackWhich = [false false];
        % Sweeps theta_s from 0 to 80 deg, holding the CURRENT design
        % (geometry, taper, element factor, phi_s) fixed, and tracks the
        % INTENDED-DIRECTION DIRECTIVITY at each steering angle, via full
        % solid-angle integration -- not just peak field magnitude along
        % a cut. This distinction matters: peak field magnitude only
        % captures the element-factor contributor to scan loss. It is
        % structurally blind to the array-side foreshortening
        % contributor, since the array factor's own peak amplitude never
        % drops with steering (AF(0)=N always) -- foreshortening only
        % shows up via beamwidth, which requires a solid-angle integral
        % to detect. Verified: isotropic elements steered to 60 deg gave
        % exactly 0 dB loss under the old peak-magnitude method, vs a
        % true -3.65 dB from directivity (matching the cos(60deg)
        % foreshortening prediction of -3.01 dB) -- the old method was
        % completely missing this contributor for every array.
        %
        % Labelled DIRECTIVITY throughout, not "gain": no radiation
        % efficiency, mismatch loss, or mutual coupling is modelled here
        % (a genuine, correctly-flagged review finding) -- calling this
        % "true gain" overstated what's actually computed.
        if isempty(S.el)
            uialert(fig,'Place at least one element first.','No elements');
            return;
        end
        if sum(abs(normalizedAmp)) <= 0
            uialert(fig,['Every element amplitude is zero - nothing to plot. ' ...
                'Set at least one amplitude to a nonzero value.'],'Zero amplitude');
            return;
        end
        sweepAngles = 0:2:80;   % stop short of grazing incidence; the separate
                                 % max-scan tool extends farther (to 88 deg). This
                                 % avoids the uninformative flat endfire plateau here.
        % peakDb / uAtDb are no longer preallocated: both are now built in
        % one shot from DatRaw / UatRaw after the sweep, so the loop never
        % indexes into them.
        % Unit-cell mode scores this plot on a different quantity --
        % see the note where the curve is assembled below.
        isUnitCellMode = ~isempty(unitCellFigures());
        % Array-factor-only curve removed for now (to be reintroduced later).

        % Seconds on a large array, so a progress bar with Cancel. Being
        % modal, it also stops the controls changing the design halfway
        % through the sweep. Indeterminate while the element patterns are
        % cached, then one step per steering angle.
        if coverage, sweepAngles = coverageSamples(:,1).'; end
        gratingFlags = false(size(sweepAngles));
        tScan = tic;
        scanTitle = 'Scan loss'; if coverage, scanTitle = 'Scan coverage'; end
        [slProg, slProgDone] = startProgress(scanTitle, ...
            'Evaluating the element pattern…', true); %#ok<ASGLU>

        % fixed observation grid, built once -- only the feed phase
        % changes per swept angle, not what directions are evaluated
        % Same integration-domain rule as computePattern's calcFull: an
        % imported CST pattern has real back-hemisphere power, and
        % leaving it out of Prad distorts directivity at every swept
        % angle.
        %
        % This curve is plotted RELATIVE to broadside, which cancels any
        % error that is CONSTANT across the sweep -- and the omitted
        % back-hemisphere power is exactly constant in the special case
        % where the element's back lobe is a scaled MIRROR of its front
        % lobe (the array factor is mirror-symmetric, so the whole
        % integrand mirrors and the ratio is scan-independent). Measured
        % error there: 0.000 dB, so that case genuinely needs nothing.
        %
        % Real back lobes are not mirrors of the front -- diffraction
        % round a finite ground plane gives a broad, often horizon-
        % tilted back pattern -- and then the omitted fraction DOES move
        % with scan angle and does not cancel. Measured against known
        % truth on an 8x8 lambda/2 array, all with a -10 dB back lobe:
        % broad (cos^0.3) back lobe 1.03 dB max error, near-isotropic
        % back 1.69 dB, horizon-tilted back 2.60 dB -- and in every case
        % the error GROWS with steering angle (worst at 80 deg), which
        % is precisely where a scan-loss plot is being read. Integrating
        % the full sphere removes it outright.
        %
        % Keeping the rule identical in both places also stops this tool
        % and the main pattern view from disagreeing -- the same class
        % of drift that already bit the squint scaling here once.
        % Full sphere also when the USER has asked for it, not only for an
        % imported pattern. A Custom formula written in terms of cosd(th)
        % rather than the clamped ct radiates behind the array, and so
        % does the Dipole; integrating only the front hemisphere for
        % those left this tool disagreeing with the main pattern view,
        % which does honour S.fullSphere -- the same one-array-two-answers
        % drift the paragraph above exists to prevent.
        thetaG = 0:1:180;  % same full-sphere power basis as the main calculation
        phiG = 0:2:360;
        [THg,PHg] = meshgrid(thetaG,phiG);
        Ug = sind(THg).*cosd(PHg); Vg = sind(THg).*sind(PHg);

        % Element-factor cache, same rotation-keyed idea and same memory
        % guard as computePattern -- but hoisted OUTSIDE the sweep loop
        % as well, because the grid-evaluated element pattern depends
        % only on rotation and not at all on the swept steering angle.
        % Without it the grid call runs numel(sweepAngles) x nElements
        % times (41 x N here); with it, once per distinct rotation for
        % the whole sweep. On a 1000-element array with four rotations
        % that is 41,000 evaluations down to 4.
        rotAllS = S.el(:,5);
        [uniqRotS,~,rotIdxS] = unique(rotAllS);
        nURS = numel(uniqRotS);
        bytesPerRotS = 2*numel(THg)*16;
        % NOTE the missing "nURS < numel(rotAllS)" term that
        % computePattern's equivalent guard carries. That term says "do
        % not bother caching when every element has its own rotation,
        % since the degenerate all-unique case does the same number of
        % evaluations either way". True for computePattern, which makes
        % ONE pass over the elements. False here: this function sweeps
        % ~45 steering angles, so even an all-unique array reuses each
        % element's pattern 45 times -- the element factor depends on
        % (theta,phi,rotation) and never on the steering angle.
        %
        % Copying that term across cost 89.7 s for an 8x8 import with
        % distinct rotations, against 0.8 s cached; a 16x16 would have
        % been about six minutes. Only the memory budget belongs here.
        useEFCacheS  = (nURS*bytesPerRotS <= 256e6);
        if useEFCacheS
            EthCacheS = cell(nURS,1); EphCacheS = cell(nURS,1);
            for i = 1:nURS
                [EthCacheS{i}, EphCacheS{i}] = elementFactor(THg,PHg,uniqRotS(i));
            end
        end

        % Beam squint: apply the SAME propagation-frequency scaling as
        % computePattern (fr = operating/design when squint mode is on).
        % This was previously missing here entirely, so squint mode had
        % no effect on this plot -- a real inconsistency with the rest
        % of the app, now fixed.
        fr = freqRatio();
        % Hoisted: seqRotSign() evaluates the element factor on a ring,
        % which must not happen once per element inside the sweep.
        seqSgn = seqRotSign();

        DatRaw = nan(size(sweepAngles));   % directivity at the steer angle
        UatRaw = zeros(size(sweepAngles)); % raw intensity there
        powerKernelBlocks = {}; % reused across scan angles for analytic element models
        for i = 1:numel(sweepAngles)
            if (coverage && ~isgraphics(S.coverageWin)) || ...
                    stepProgress(slProg, (i-1)/numel(sweepAngles), sprintf( ...
                    'Steering angle %d of %d', i, numel(sweepAngles)))
                reportAnalysis('Scan loss cancelled: no curve drawn');
                return;
            end
            ts = sweepAngles(i);
            if coverage, scanPhi = coverageSamples(i,2); end
            AFcoverage = zeros(size(THg));
            us_ = sind(ts)*cosd(scanPhi); vs_ = sind(ts)*sind(scanPhi);
            Eth_tot = zeros(size(THg)); Eph_tot = zeros(size(THg));
            % Exact single-point accumulators, evaluated directly AT
            % (ts,S.phi_s) -- NOT looked up from the grid. This is the
            % fix for a real, confirmed accuracy issue: the old code
            % found U_at via nearest-grid-point lookup (thetaG step 1
            % deg, phiG step 2 deg), which is exact for theta (every
            % swept ts already lands exactly on the 1-deg grid) but NOT
            % for phi -- phi_s can be any value (the Steer phi spinner's
            % own step is 5 deg, not a multiple of phiG's 2-deg grid),
            % so a real quantization error was possible. Verified
            % numerically: a 16x16 array at phi_s=37 deg (1 deg off the
            % nearest 36-deg grid point) showed a genuine 0.115 dB
            % error -- small here, but grows for larger/narrower-beam
            % arrays, exactly as flagged. Evaluating a single point
            % directly costs almost nothing extra (one more per-element
            % sum, not a finer grid), so there's no reason to keep the
            % approximation now that it's identified.
            Eth_at = 0; Eph_at = 0;
            % The single-POINT element factor is evaluated at (ts,phi_s),
            % which moves with the sweep, so it cannot be hoisted out
            % like the grid one -- but within this iteration it still
            % only depends on rotation, so it is computed once per
            % distinct rotation here rather than once per element.
            Eth0U = zeros(nURS,1); Eph0U = zeros(nURS,1);
            % Index named iu, NOT i: this loop is nested inside the sweep
            % loop, which is itself indexed by i and still needs it after
            % this point for "peakDb(i) = ..." at the end of the body.
            % Reusing i left it holding nURS there, so every swept angle
            % wrote its result into the SAME slot of peakDb and the other
            % 40 entries kept their initialised zeros. The plotted scan
            % loss was a two-level step (0 dB at broadside, one constant
            % value everywhere else), not a curve -- confirmed on the
            % default 8x8: 2 distinct values across 41 sweep points.
            for iu = 1:nURS
                [Eth0U(iu), Eph0U(iu)] = elementFactor(ts,scanPhi,uniqRotS(iu));
            end
            for n = 1:size(S.el,1)
                xn=S.el(n,1); yn=S.el(n,2); amp=normalizedAmp(n); ph0=S.el(n,4); rot=S.el(n,5);
                phFeed = -S.k*phaseFreqRatio()*(xn*us_+yn*vs_) + deg2rad(ph0);
                if S.seqPhase, phFeed = phFeed + seqSgn*deg2rad(rot); end
                steer = amp*exp(1j*phFeed)*exp(1j*S.k*fr*(xn*Ug + yn*Vg));
                if coverage, AFcoverage = AFcoverage + steer; end
                if useEFCacheS
                    Eth = EthCacheS{rotIdxS(n)};  Eph = EphCacheS{rotIdxS(n)};
                else
                    [Eth,Eph] = elementFactor(THg,PHg,rot);
                end
                Eth_tot = Eth_tot + steer.*Eth;
                Eph_tot = Eph_tot + steer.*Eph;

                steer_at = amp*exp(1j*phFeed)*exp(1j*S.k*fr*(xn*us_ + yn*vs_));
                Eth_at = Eth_at + steer_at*Eth0U(rotIdxS(n));
                Eph_at = Eph_at + steer_at*Eph0U(rotIdxS(n));
            end
            % ALWAYS total power (both components), never routed through
            % polCombine/ddPol.Value: this integral is compared against
            % the L_AF+L_EF reference below, which is itself a TOTAL-
            % power quantity. If the Polarization readout dropdown were
            % left on "RHCP component" or "LHCP component", using
            % polCombine here would silently integrate only PART of the
            % radiated power, inflating the computed directivity.
            P_sweep = fieldPower(Eth_tot,Eph_tot);

            if ismember(S.efType,{'Isotropic','cos^q(theta)','Short dipole (z-axis)', ...
                    'Dipole (linear pol)','Patch (cos^q x lin pol)'})
                feedPhase = -S.k*phaseFreqRatio()* ...
                    (S.el(:,1)*us_+S.el(:,2)*vs_) + deg2rad(S.el(:,4));
                if S.seqPhase, feedPhase = feedPhase + seqSgn*deg2rad(S.el(:,5)); end
                if strcmp(S.efType,'Isotropic')
                    Prad = exactIsotropicPrad(normalizedAmp.*exp(1j*feedPhase));
                elseif strcmp(S.efType,'cos^q(theta)')
                    [Prad,powerKernelBlocks] = exactCosQPrad( ...
                        normalizedAmp.*exp(1j*feedPhase),powerKernelBlocks);
                elseif strcmp(S.efType,'Short dipole (z-axis)')
                    [Prad,powerKernelBlocks] = exactZDipolePrad( ...
                        normalizedAmp.*exp(1j*feedPhase),powerKernelBlocks);
                else
                    [Prad,powerKernelBlocks] = exactVectorPrad( ...
                        normalizedAmp.*exp(1j*feedPhase),S.el(:,5),powerKernelBlocks);
                end
            elseif needsFinePowerIntegral()
                feedPhase = -S.k*phaseFreqRatio()* ...
                    (S.el(:,1)*us_+S.el(:,2)*vs_) + deg2rad(S.el(:,4));
                if S.seqPhase, feedPhase = feedPhase + seqSgn*deg2rad(S.el(:,5)); end
                Prad = quadraturePrad(normalizedAmp.*exp(1j*feedPhase));
            else
                [~, Prad] = fieldDirectivity(P_sweep, thetaG, phiG, THg);
            end

            % directivity specifically AT the intended (theta_s,phi_s)
            % direction, referenced to the SAME total radiated power --
            % this is what correctly captures both contributors
            % together. U_at now comes from the EXACT single-point
            % evaluation above, not a grid lookup.
            U_at = fieldPower(Eth_at,Eph_at);
            if Prad > 0
                Dat = 4*pi*U_at/Prad;
            else
                Dat = NaN;
            end
            % Kept RAW here and converted to dB after the sweep. The eps
            % floors that used to sit on U_at and Dat at this point were
            % absolute, and U_at is an absolute intensity the Custom
            % element formula scales freely: below about 1e-10 on an 8x8
            % every angle's U_at fell under eps together, so the
            % numerator froze at eps while Prad went on shrinking with
            % steer angle, and Dat = const/Prad ROSE. The window then
            % drew +4.39 dB of scan GAIN at 80 degrees where the true
            % curve reads -10.82 dB of loss -- not merely wrong, but the
            % wrong sign, for a planar array that cannot have it.
            DatRaw(i) = Dat;
            UatRaw(i) = U_at;
            if coverage
                gratingFlags(i) = coverageGrating(ts,scanPhi,AFcoverage,thetaG,phiG);
            end
            % Unit-cell mode is scored on RADIATED INTENSITY instead.
            % Dat divides by Prad from the pattern integral, and for a
            % periodic unit cell that integral is not the radiated power
            % -- the same invalid basis already corrected in the title,
            % the 3D surface and the cut plot. Here it was worse than a
            % constant offset: Prad drifts with steering, so the curve
            % carried a scan-dependent error and this plot disagreed
            % with the gain in the main window (-2.50 dB against the
            % title's -1.44 dB, both labelled directivity).
            %
            % U_at needs no basis at all: input power is fixed, so gain
            % scales with intensity at the steer direction, and its
            % ratio to broadside is exactly the embedded element's own
            % rolloff -- which is what the title reports.
        end
        closeProgress(slProg);   % before any refusal alert below
        if ~isgraphics(fig), return; end   % app closed during the last step
        % dB conversion AFTER the sweep, so each floor can be set relative
        % to the sweep's own maximum instead of to an absolute constant.
        % 1e-12 of the maximum is 120 dB down -- far below any scan loss
        % worth plotting (the axis shows tens of dB), so the floor exists
        % only to stop a true zero becoming -Inf, which is what the eps
        % floors were for. Unlike eps it moves with the data, so the
        % curve is now identical at every field scale.
        peakDb = relDb(DatRaw);
        uAtDb  = relDb(UatRaw);
        % The BROADSIDE entry is the reference the whole curve is quoted
        % against, so it has to be real radiation and not the floor.
        % relDb's floor keeps a true zero from becoming -Inf, which is
        % right for a single dead sweep point -- but when the DEAD point
        % is element 1, subtracting it turns the floor's own 120 dB span
        % into the plotted curve, and the window drew a smooth rise to
        % exactly +120.000 dB: scan GAIN, for a planar array, its value
        % set entirely by the floor constant. An element with a genuine
        % broadside null (a custom 'sind(th)', a monopole-like pattern)
        % reaches this through the UI.
        if isUnitCellMode, refRaw = UatRaw(1); else, refRaw = DatRaw(1); end
        if ~(isfinite(refRaw) && refRaw > 0)
            uialert(fig, ['This element radiates nothing at broadside, so there ' ...
                'is no broadside level to measure scan loss against -- every ' ...
                'point on the curve would be quoted against zero.' newline newline ...
                'Scan loss is defined as a RATIO to the broadside beam. Use an ' ...
                'element with a broadside beam, or read the absolute levels in ' ...
                'the main window instead.'], 'Scan loss undefined');
            return;
        end
        if isUnitCellMode
            peakDb = uAtDb - uAtDb(1);
        else
            peakDb = peakDb - peakDb(1);   % relative to broadside (theta_s=0)
        end

        % REFUSES, rather than plotting a warned curve. Preserving the
        % warning was the previous step and it was not enough: every
        % point on that curve is computed from a substituted element, so
        % the plot answers a question about a different antenna. Same
        % decision as Max scan, for the same reason -- a caveated wrong
        % answer is still a wrong answer, and this one is a picture the
        % user will read levels off.
        %
        % This analysis has its own grid, so it needs its own check:
        % computePattern's runs on a different one and cannot speak for
        % it, in either direction.
        if S.efFallback
            uialert(fig, ['No valid scan-loss curve.' newline newline ...
                'The Custom element formula could not be evaluated across ' ...
                'this analysis'' grid, so an element was substituted for it. ' ...
                'Every point on the curve would describe that substitute ' ...
                'rather than your design, so none is plotted.' newline newline ...
                'Note the formula can be perfectly valid at the angles the ' ...
                'sweep itself samples and still fail on the wider grid the ' ...
                'directivity integral needs -- 1./(th-89) is exactly that.'], ...
                'Scan loss', 'Icon','warning');
            return;
        end
        if coverage
            result = struct('samples',coverageSamples,'lossDb',peakDb(:), ...
                'grating',gratingFlags(:),'unitCell',isUnitCellMode);
            return;
        end
        if useAzEl(), sweepDisplay = 90-sweepAngles;
        else, sweepDisplay = scanThetaSign()*sweepAngles; end
        plotScanLoss(sweepAngles, sweepDisplay, peakDb, isUnitCellMode);
        reportAnalysis(sprintf('Scan loss: %d steering angles in %.1f s', ...
            numel(sweepAngles), toc(tScan)));
    end


    function y = relDb(v)
        %RELDB  10*log10 with a floor 120 dB below the data's own maximum.
        %   Absolute floors (eps, 1e-12) cannot be used on a quantity whose
        %   scale a Custom element formula sets: once every sample sits
        %   under the floor the whole vector flattens, and any ratio taken
        %   against it afterwards is meaningless. Scaling the floor with
        %   the data keeps it doing its only real job -- keeping a true
        %   zero from becoming -Inf.
        m = max(v(isfinite(v)));
        if isempty(m) || ~(m > 0)
            fl = realmin;
        else
            fl = m * 1e-12;
        end
        y = 10*log10(max(v, fl));
    end
    function findMaxScanAngleCore()
        if ~patternReady(true), return; end
        normalizedAmp = calculationWeights();
        scanPhi = scanPhiVal();
        % Complements the live grating-lobe readout in the info panel
        % (which reports the grating lobe AT the CURRENT theta_s/phi_s).
        % This instead sweeps theta_s (holding phi_s fixed at its
        % current value) and finds the first angle where a grating lobe
        % enters visible space -- the actual safe-scan threshold for
        % THIS specific array (whatever its current spacing/stagger/
        % grid-angle/shape happen to be), found the same
        % geometry-agnostic way as the live readout: computing the real
        % array factor and running the full-hemisphere flood-fill search
        % at each swept angle, rather than a closed-form formula.
        if isempty(S.el)
            uialert(fig,'Place at least one element first.','No elements');
            return;
        end
        ref = sum(abs(normalizedAmp));
        if ~(ref > 0)
            uialert(fig,'Every element amplitude is zero.','Zero amplitude');
            return;
        end
        theta = 0:1:90; phi = 0:2:360;
        [TH,PH] = meshgrid(theta,phi);
        U = sind(TH).*cosd(PH); V = sind(TH).*sind(PH);
        x = S.el(:,1); y = S.el(:,2); amp = normalizedAmp; ph0 = S.el(:,4); rot = S.el(:,5);

        % Swept to 88, not 80: the whole point is to answer "can I steer
        % there", and the steering control itself allows 90.
        sweepAngles = 0:2:88;
        % A TRUE grating lobe is a repeat of the main beam -- every
        % element back in phase -- so it arrives at essentially the SAME
        % level as the main beam, taper included, because the taper
        % repeats with it. That is what separates it from an ordinary
        % sidelobe or the main lobe's own shoulder.
        %
        % The old threshold was -6 dB, which is far too loose for that
        % test: at 0.5 lambda spacing the main lobe's shoulder reaches
        % -4.5 dB at the horizon around 60 deg scan and tripped it, so
        % the tool announced a grating lobe for a lattice that cannot
        % produce one at any scan angle. Checked against the closed form
        % (sin(theta_s) = lambda/d - 1) across 0.5/0.6/0.75/1.0 lambda.
        % Level alone CANNOT make this call. As the scan angle grows the
        % grating lobe's shoulder climbs continuously toward 0 dB at the
        % horizon, so any level threshold just picks an arbitrary point
        % on that ramp: -6 dB fired at 60 deg for a 0.5 lambda lattice,
        % -3 dB merely moved it to 64 deg. Neither is a grating lobe.
        %
        % The real test is geometric: a grating lobe has ENTERED visible
        % space when its PEAK sits at theta < 90, rather than the horizon
        % slicing through its skirt. Requiring an interior peak
        % reproduces the closed form sin(theta_s) = lambda/d - 1.
        % Swept against that closed form at 0.50 (never), 0.52, 0.5137,
        % 0.55, 0.5644, 0.58, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
        % 0.95 and 1.00 lambda, plus eight rectangular dx/dy pairs: every
        % onset lands within one 2 deg sweep step of theory, always on
        % the late (conservative) side.
        %
        % An earlier version of this comment claimed the same agreement
        % at 1.0 lambda and it did not hold -- the tool reported 6 deg
        % against a theoretical 0. That was the tied-lobe masking fixed
        % in findGratingLobe, not a limitation of the interior test.
        %
        % The level gate stays only to rule out calling an ordinary
        % interior sidelobe a grating lobe; a uniform array's first
        % sidelobe is -13 dB, well clear of it.
        GL_REL_DB = -6;
        onsetScanTheta = NaN; onsetLobeTheta = NaN; onsetPhi = NaN; onsetDb = NaN;
        worstDb = -Inf; worstTs = NaN;
        grazeTs = NaN;   % first angle a strong lobe reaches the horizon
        % Same operating/design-frequency propagation scaling as
        % showScanLoss and computePattern -- this was previously missing
        % here, so squint mode had no effect on this search (a real,
        % confirmed inconsistency: the two scan-related tools could
        % disagree about safe scan range under squint).
        fr = freqRatio();

        % Progress with Cancel, as for Scan loss: indeterminate while the
        % element patterns are cached, then one step per steering angle.
        tScanMS = tic;
        [msProg, msProgDone] = startProgress('Max scan', ...
            'Evaluating the element pattern…', true); %#ok<ASGLU>

        % ---- element factors, hoisted out of BOTH loops ----------------
        % The element pattern is a function of (theta,phi,rotation) only
        % -- it does NOT depend on the steering angle -- so evaluating it
        % per element per swept angle repeated identical work 45 times
        % over. For an IMPORTED pattern each call is four
        % scatteredInterpolant evaluations across the whole grid, and the
        % cost was brutal: measured 36 s for an 8x8 import, scaling
        % linearly with element count (a 16x16 would have taken ~2.5 min,
        % a 32x32 about ten). Caching by distinct rotation turns 45 x N
        % calls into one per rotation -- for the usual uniform-rotation
        % array, ONE.
        %
        % Guarded on memory the same way computePattern's cache is: the
        % number of distinct rotations is not bounded (sequential
        % rotation with block = M x N gives every element its own), and
        % one grid pair per rotation would otherwise run to gigabytes.
        % Past the budget it falls back to evaluating inside the loop,
        % which is slow but bounded.
        [uRotMS,~,rIdxMS] = unique(rot);
        nRotMS = numel(uRotMS);
        % Cleared BEFORE the element-factor cache below, not after it.
        % Not inherited (computePattern resets these too, but with
        % auto-recompute off it may not have run since the user fixed the
        % formula), and not cleared late either: the cache is where
        % elementFactor actually gets called for this analysis, so a reset
        % placed after it wiped the very flag the cache had just set and
        % the sweep -- reading the cache -- never set it again.
        S.efFallback = false; S.efFallbackWhich = [false false];
        cacheMS = (nRotMS*2*numel(TH)*16) <= 256e6;
        if cacheMS
            EthRot = cell(nRotMS,1); EphRot = cell(nRotMS,1);
            for q = 1:nRotMS
                [EthRot{q}, EphRot{q}] = elementFactor(TH,PH,uRotMS(q));
            end
        end

        % Hoisted out of the per-element loop below: seqRotSign()
        % evaluates the element factor on a ring.
        seqSgn = seqRotSign();
        deadSweep = false;   % set if any swept angle radiates nothing
        nValidTs  = 0;       % angles that actually had a pattern to score
        nEvalTs   = 0;       % angles the loop actually reached
        nDeadTs   = 0;       % of those, how many radiated nothing
        try
            for i = 1:numel(sweepAngles)
                if stepProgress(msProg, (i-1)/numel(sweepAngles), sprintf( ...
                        'Steering angle %d of %d', i, numel(sweepAngles)))
                    reportAnalysis('Max scan cancelled: no result');
                    return;
                end
                ts = sweepAngles(i);
                nEvalTs = nEvalTs + 1;   % reached, whatever happens below
                us = sind(ts)*cosd(scanPhi); vs = sind(ts)*sind(scanPhi);
                % TOTAL field (element factor x array factor), not the
                % array factor alone. Sweeping the AF was a real defect:
                % the element pattern falls off toward the horizon, so it
                % pulls a grating lobe's apparent PEAK inward. An AF lobe
                % that merely grazes at theta = 90 -- which the interior
                % test correctly rejects -- becomes a genuine interior
                % lobe once the element pattern multiplies it in.
                %
                % Measured on a 4x4 unit-cell import at 0.5 lambda: this
                % tool reported "no grating lobe, practical limit 44 deg"
                % while steering to 60 deg put a lobe at theta = 71.6 deg,
                % 3.6 dB below the peak, and the cut plot flagged it. The
                % cut plot analyses the Total curve; this now does too, so
                % all three readouts describe the same pattern.
                Eaf = zeros(size(TH));
                Eth = zeros(size(TH)); Eph = zeros(size(TH));
                for n = 1:numel(x)
                    phFeed = -S.k*phaseFreqRatio()*(x(n)*us+y(n)*vs) + deg2rad(ph0(n));
                    if S.seqPhase, phFeed = phFeed + seqSgn*deg2rad(rot(n)); end
                    prop = amp(n)*exp(1j*phFeed)*exp(1j*S.k*fr*(x(n)*U+y(n)*V));
                    Eaf = Eaf + prop;
                    if cacheMS
                        q = rIdxMS(n);
                        Eth = Eth + prop.*EthRot{q};
                        Eph = Eph + prop.*EphRot{q};
                    else
                        [et,ep] = elementFactor(TH,PH,rot(n));
                        Eth = Eth + prop.*et;
                        Eph = Eph + prop.*ep;
                    end
                end
                % Two grids, deliberately. The ARRAY FACTOR supplies the
                % null structure that bounds the main lobe; the TOTAL
                % supplies the level, because the element pattern is what
                % turns a grating lobe grazing the horizon into a lobe
                % you can actually see well inside visible space.
                afDb  = 20*log10(max(abs(Eaf)/ref,1e-4));
                % Floored RELATIVE to its own peak, not at an absolute
                % 1e-8. ref = sum|amp| divides out the taper scale but
                % NOT the element-factor scale, so a Custom formula of
                % magnitude A puts this grid's peak at A^2 -- and below
                % A = 1e-4 every cell clamped to the same -80 dB. On a
                % flat grid findGratingLobe's level comparisons all read
                % 0.0 dB and its interior falloff test never fires, so
                % the dialog reported "NO GRATING LOBE up to 88 deg,
                % treat theta_s = 0 deg as the practical limit" for a
                % 1.0 lambda lattice whose real onset is 14 deg. Only
                % differences from the peak are ever used downstream, so
                % a peak-relative floor changes nothing where the old one
                % did not already bite.
                pwTot = fieldPower(Eth,Eph)/ref^2;
                pwMax = max(pwTot(:));
                if ~(pwMax > 0)
                    % Nothing radiates at this steer angle. ACTUALLY skip
                    % it: the first version set totDb = -Inf, flagged
                    % deadSweep and then fell through to findGratingLobe
                    % anyway, so an all-dead sweep still produced a
                    % verdict -- with a warning bolted on top of a
                    % "safe to steer across the whole range" conclusion
                    % drawn from a field that does not exist.
                    deadSweep = true; nDeadTs = nDeadTs + 1;
                    continue
                end
                totDb = 10*log10(max(pwTot, pwMax*1e-12));
                nValidTs = nValidTs + 1;
                % inTh/inP are the lobe's OWN location, and the message
                % below reports them as such. It once reported the
                % STEERING angle (ts) instead, which was wrong and
                % confirmed so: steering to 40 deg with a lobe actually
                % appearing at 90 deg read "appears at theta=40 deg",
                % falsely implying the lobe sits on the main beam. Only
                % the LEVEL is taken from the overall strongest lobe
                % (glDb, for the worst-case and grazing readouts); the
                % onset decision and the reported location both come from
                % the strongest INTERIOR lobe.
                [~,~,glDb,inTh,inP,inDb] = findGratingLobe(afDb,theta,phi,totDb);
                % Track the worst secondary lobe across the whole sweep,
                % so a "no grating lobe" answer can still say how close
                % the array actually came instead of just "none".
                if glDb > worstDb, worstDb = glDb; worstTs = ts; end
                % A strong lobe can sit ON the horizon without its peak
                % ever entering visible space -- a half-wave lattice does
                % exactly this near grazing, reaching full strength at
                % theta = 90. That is not a grating lobe by the interior
                % test below, but it is emphatically not "safe" either,
                % so it is tracked separately and reported.
                if glDb > GL_REL_DB && isnan(grazeTs), grazeTs = ts; end
                % Interior test by FALLOFF, not by a fixed cutoff angle.
                % A grating lobe whose peak is genuinely inside visible
                % space has lower levels just beyond it; the horizon
                % slicing through a rising skirt does not. A fixed
                % "theta < 89" gate could not express that -- it still
                % fired at 88 deg for a 0.5 lambda lattice, where the
                % repeat's true peak is at u = -1.0006, fractionally
                % outside visible space.
                %
                % The falloff test now lives inside findGratingLobe and is
                % applied to EVERY candidate, which is why this reads the
                % in* outputs rather than re-testing gl* here: gl* is the
                % strongest lobe anywhere and can be a grazing one that
                % ties with a real interior lobe, in which case testing
                % gl* alone declared the array safe. See the interior
                % block in findGratingLobe.
                if inDb > GL_REL_DB
                    onsetScanTheta = ts; onsetLobeTheta = inTh;
                    onsetPhi = inP; onsetDb = inDb;
                    break;
                end
            end
        catch err
            closeProgress(msProg);
            setStatus(['Max scan failed: ' err.message],'bad');
            if isgraphics(fig), uialert(fig, err.message, 'Scan failed'); end
            return;
        end
        closeProgress(msProg);   % before the result or refusal alert
        if ~isgraphics(fig), return; end   % app closed during the last step

        if isnan(onsetScanTheta) && isnan(grazeTs)
            msg = sprintf(['NO GRATING LOBE at any steering direction up to %s ' ...
                '(%s).\n\nNo secondary lobe above the sampled -6 dB threshold was detected ' ...
                'in the evaluated range; this is not a general safe-scan guarantee.\n\nThe strongest secondary lobe found ' ...
                'anywhere in the sweep was %.1f dB below the peak, at %s -- an ordinary sidelobe.'], ...
                scanDirectionText(sweepAngles(end)),scanPlaneText(), ...
                abs(worstDb),scanDirectionText(worstTs));
        elseif isnan(onsetScanTheta)
            msg = sprintf(['NO GRATING LOBE enters visible space up to %s ' ...
                '(%s) -- but treat %s as ' ...
                'the practical limit.\n\nFrom %s onward a near-full-strength ' ...
                'lobe GRAZES the horizon: its peak stays just outside visible ' ...
                'space, so it is not a grating lobe by the strict test, but what ' ...
                'you see at %s reaches %.1f dB below your main beam ' ...
                'by %s.\n\nSwept in 2 deg steps.'], ...
                scanDirectionText(sweepAngles(end)),scanPlaneText(), ...
                scanDirectionText(grazeTs),scanDirectionText(grazeTs), ...
                scanDirectionText(90),abs(worstDb),scanDirectionText(worstTs));
        else
            % "from X onward" claimed something the sweep never checked:
            % it stops at the FIRST detection, so every angle past X is
            % unexamined. What is known is where it starts.
            msg = sprintf(['STRONG SECONDARY LOBE detected at %s ' ...
                '(%s).\n\nThe first sampled threshold crossing is %s; earlier sweep points tested ' ...
                'below the threshold.\n\nAt %s a repeat of the main beam sits at ' ...
                '%s, only %.1f dB below the peak ' ...
                '-- above the chosen detection threshold ' ...
                'under this heuristic classification.\n\nSwept in 2 deg steps, so the true ' ...
                'onset could be up to ~2 deg before this sampled direction. The sweep STOPPED here, ' ...
                'so later steering directions were not examined -- they are ' ...
                'not implied to be worse or better.'], ...
                scanSteerText(onsetScanTheta),scanPlaneText(), ...
                scanDirectionText(onsetScanTheta),scanDirectionText(onsetScanTheta), ...
                directionText(onsetLobeTheta,onsetPhi),abs(onsetDb));
        end
        if isnan(onsetScanTheta) && isnan(grazeTs) && ~isfinite(worstDb)
            msg=sprintf('No distinct secondary lobe was resolved at the valid sampled steering directions up to %s. This sampled heuristic does not establish a safe scan limit.', ...
                scanDirectionText(sweepAngles(nEvalTs)));
        end
        % Every result above is read off a pattern the element factor
        % produced, so a substituted formula or a dead sweep angle makes
        % the verdict describe something other than the design. Said
        % here because this analysis calls elementFactor on its own grid
        % and never passes through computePattern's check.
        % No valid angle means there is no analysis, not a permissive one.
        % Replacing the verdict rather than prefixing it: a "safe to steer
        % across the whole range" sentence under a warning still reads as
        % a result, and this one was computed from nothing.
        if nValidTs == 0
            showMaxScanResult(['No valid radiation data.' newline newline ...
                'Every steering angle in the sweep produced a pattern that ' ...
                'radiates nothing, so there is no grating-lobe analysis to ' ...
                'report -- not a clean result. Check the element factor: a ' ...
                'Custom formula that evaluates to zero, or an imported pattern ' ...
                'with no power in the swept range, will do this.']);
            return;
        end
        % Refuses rather than advises. A recommendation computed from a
        % substituted element is not a weaker recommendation, it is a
        % recommendation about a different antenna.
        if S.efFallback
            showMaxScanResult(['No valid scan analysis.' newline newline ...
                'The Custom element formula could not be evaluated on this ' ...
                'analysis grid, so an element was substituted for it. Any ' ...
                'maximum-scan verdict would describe that substitute, not your ' ...
                'design, so none is given.' newline newline ...
                'Fix the formula and run this again.']);
            return;
        end
        if deadSweep
            showMaxScanResult(sprintf(['Incomplete scan analysis: %d of %d evaluated angles had no radiation. ' ...
                'No continuous scan limit can be determined across those missing angles.'], ...
                nDeadTs,nEvalTs));
            return;
        end
        reportAnalysis(sprintf('Max scan: %d steering angles checked in %.1f s', ...
            nEvalTs, toc(tScanMS)));
        showMaxScanResult(msg);
    end

    function idx = selIdx()
        % selected elements, or all of them when nothing is selected
        if isempty(S.el)
            idx = [];
        elseif isempty(S.sel)
            idx = (1:size(S.el,1))';
        else
            idx = S.sel(S.sel <= size(S.el,1));
        end
    end

    function setRot(deg)
        idx = selIdx(); stepsWas = undoCount;
        if isempty(idx), return; end
        S.el(idx,5) = mod(deg,360);
        refreshAll();
        announceAllElements(stepsWas, ...
            sprintf('Set the rotation of %%s to %g°', mod(deg,360)));
    end

    function stepRot(delta)
        idx = selIdx(); stepsWas = undoCount;
        if isempty(idx), return; end
        S.el(idx,5) = mod(S.el(idx,5) + delta, 360);
        refreshAll();
        announceAllElements(stepsWas, sprintf('Rotated %%s by %+g°', delta));
    end

    function applySeqRot(step,blockM,blockN)
        % Rotation only -- sets each element's physical/pattern
        % orientation (S.el(:,5)) for polarization purposes (e.g. a
        % sequentially-rotated CP subarray of patches). This function itself
        % never edits the manual Phase column. If the separate ADVANCED
        % 'extra rot angle -> feed phase' checkbox is enabled, downstream
        % pattern math may additionally use S.el(:,5) as an electrical phase.
        if isempty(S.el), return; end
        % Block-relative indexing, NOT a running index across the whole
        % array: rotation depends only on (row,col) mod (blockM,blockN),
        % so the SAME small 0/step/2*step/... pattern repeats identically
        % in every blockM x blockN tile, keeping elRC (each element's
        % fixed [row col] lattice identity, stable across moveSelected)
        % as the position reference. Set blockM=M,
        % blockN=N to instead have one continuously-incrementing pattern
        % across the whole array.
        blockM = max(1,round(blockM)); blockN = max(1,round(blockN));
        rowInBlock = mod(S.elRC(:,1)-1, blockM);
        colInBlock = mod(S.elRC(:,2)-1, blockN);
        % Boustrophedon (snake) traversal starting at the BOTTOM-LEFT
        % corner of each block. This follows the array's element-numbering
        % convention: element 1 is lattice row 1 / column 1, i.e. the
        % physical bottom-left element, and it is the 0-degree reference.
        %
        % The sequence first moves left-to-right along the BOTTOM row,
        % then moves up one row and comes back right-to-left, continuing
        % through physically adjacent elements. For the classic 2x2 tile:
        %
        %       TL = 3*step       TR = 2*step
        %       BL = 0            BR = 1*step
        %
        % so with the default +90-degree step:
        %
        %       TL = 270 deg      TR = 180 deg
        %       BL =   0 deg      BR =  90 deg
        %
        % This keeps element 1 in its default orientation and starts the
        % sequential rotation from that physical reference element.
        rowFromBottom = rowInBlock;
        isOddRowFromBottom = mod(rowFromBottom,2)==1;
        colSnake = colInBlock;
        colSnake(isOddRowFromBottom) = blockN-1 - colInBlock(isOddRowFromBottom);
        idxInBlock = rowFromBottom*blockN + colSnake;
        S.el(:,5) = mod(idxInBlock*step, 360);
        refreshAll();
        % Said HERE, at the moment the mixed-rotation state comes into
        % existence. The import warns about this too, but that dialog is
        % dismissed long before anyone presses this button, and the
        % failure it describes is silent: the plot redraws looking
        % perfectly ordinary, just tens of dB low.
        if magOnlyMixedRot()
            uialert(fig, ['The loaded pattern is a magnitude-only CST "Abs" export: it ' ...
                'carries no split between the two polarisation components, so this app ' ...
                'has had to put the whole magnitude into E_theta with E_phi = 0.' newline newline ...
                'That is fine while every element shares one orientation. Sequential ' ...
                'rotation has just given them different ones, and rotating an element ' ...
                'rotates its polarisation -- so the array sum is now combining ' ...
                'components that do not exist in the file. The gain and pattern shown ' ...
                'from here on are WRONG, not approximate: the same element exported ' ...
                'both ways differs by 8.6 dB on an 8x8, and by far more in unit-cell ' ...
                'mode.' newline newline ...
                'Re-export from CST with a component pair (Theta/Phi) selected rather ' ...
                'than Abs, and the rotation will be handled correctly. Until then the ' ...
                'plot titles carry a NOT VALID note.'], ...
                'Rotation needs polarisation data');
        end
    end

    function moveSelected(dx,dy)
        % shift the selected element(s) by an exact, arbitrary distance -
        % this is the perturbation delta_n from the sidelobe-decorrelation
        % study: x_n -> x_n + dx, y_n -> y_n + dy.
        idx = selIdx(); stepsWas = undoCount;
        if isempty(idx), return; end
        S.el(idx,1) = S.el(idx,1) + dx;
        S.el(idx,2) = S.el(idx,2) + dy;
        S.pending = [];   % cancel any armed 9-pt cluster -- it could be
                           % stale/invalid after this move, and would
                           % otherwise silently hijack the user's next click
        refreshAll();
        announceAllElements(stepsWas, sprintf('Moved %%s by (%g, %g) λ', dx, dy));
    end


    function resetSelectedToLattice()
        % undo any perturbation: snap selected element(s) back to their
        % own cell's ideal, perfectly periodic lattice position.
        % Uses latticePosAt (skewed + staggered), NOT the plain
        % orthogonal (col-1)*dx,(row-1)*dy formula this used to have --
        % that was a real bug: with a non-90 Grid angle or a nonzero Row
        % stagger active, Reset was snapping elements to the WRONG
        % position (the plain rectangular grid, not the actual current
        % lattice), silently undoing the skew/stagger for whichever
        % elements got reset. Found on a full pass through every place
        % that still referenced (row,col)->(x,y) directly, rather than
        % through latticePosAt.
        idx = selIdx(); stepsWas = undoCount;
        if isempty(idx), return; end
        % Sparse and Sub-position modes PLACE elements on the plain
        % orthogonal dx,dy grid -- onGridClick and refreshLayout both
        % ignore grid angle and stagger there. Resetting through the
        % skewed lattice therefore snapped an element to a position it
        % could never have been placed at, moving it instead of undoing
        % a perturbation.
        if strcmp(S.mode,'Uniform grid')
            [gx,gy] = latticePosAt(S.elRC(idx,1),S.elRC(idx,2),S.dx,S.dy,S.gridAngle,S.stagger);
        else
            [gx,gy] = latticePosAt(S.elRC(idx,1),S.elRC(idx,2),S.dx,S.dy,90,0);
        end
        S.el(idx,1) = gx;
        S.el(idx,2) = gy;
        S.pending = [];
        refreshAll();
        announceAllElements(stepsWas, 'Reset %s to the lattice');
    end

    % --------------------------------------------------------- refreshing
    function refreshAll()
        % refreshInfo() runs here ONLY when nothing recomputed.
        % computePattern ends by calling refreshInfo itself, so calling it
        % twice would repeat the O(n^2) nearest-neighbour scan for no gain.
        % When Auto is off (or radiation is blocked), maybeCompute() now
        % INVALIDATES old pattern metrics rather than leaving stale D/SLL/
        % grating-lobe numbers next to newly changed controls; refreshInfo()
        % then rebuilds the panel from that explicitly invalid state.
        refreshLayout(); refreshTable(); refreshPick();
        if ~maybeCompute()
            refreshInfo();
        end
    end

    function refreshTable()
        if patternReady(false)
            % A phase that is a whole number of turns in exact arithmetic
            % can wrap to +-1e-13 instead of 0 (theta 45, phi 45, element
            % at y = 2 lambda), and the table printed '-1.1369e-13' in a
            % column that otherwise reads -90 / 180. Display only: the
            % exports and the pattern call effectivePhaseDeg themselves.
            phShown = effectivePhaseDeg();
            phShown(abs(phShown) < 1e-9) = 0;
            tbl.Data = [S.el, phShown];
        else
            tbl.Data = [S.el, nan(size(S.el,1),1)];
            for blockedMap = 1:numel(S.portMapRefresh)
                if ~isempty(S.portMapRefresh{blockedMap})
                    try S.portMapRefresh{blockedMap}(); catch, end
                end
            end
        end
        % Selection AFTER Data, and here rather than in the callers.
        % Every add sets S.sel to the new element's index while tbl.Data
        % is still one row short, so the assignment threw "Selection
        % indices are out of data boundary" and the newly placed element
        % ended up not selected at all. Data and Selection have to move
        % together, and this is the only place that sees both.
        try
            tbl.Selection = S.sel(:).';
        catch
        end
        refreshPhaseMaps();
    end

    function refreshControlVisibility()
        % Disabled values remain intact so switching models is reversible.
        imported = strcmp(S.efType,'Imported (CST far-field)');
        custom = strcmp(S.efType,'Custom (formula)');
        setEnabled(spQ,any(strcmp(S.efType,{'cos^q(theta)','Patch (cos^q x lin pol)'})));
        usesBeamwidth = any(strcmp(S.efType,{'Gaussian','Sinc','3GPP TR 38.901 shape'}));
        setEnabled([spEfBeamAz spEfBeamEl],usesBeamwidth);
        setEnabled([efCustom efCustomPh],custom);
        setEnabled([bImportFF bImportRibbon cbUnitCell],imported);
        setEnabled(spTotEff,imported && S.impUnitCell);
        uniform = strcmp(S.mode,'Uniform grid');
        setEnabled([ddShape spGridAngle spStagger spStaggerAngle],uniform);
        setEnabled(shapeBtns,uniform);
        setEnabled([shapePreviewOlder shapePreviewPrevious ...
            shapePreviewCurrent btnShapeGallery],uniform);
        if ~uniform && ~isempty(shapeGalleryPopup) && isgraphics(shapeGalleryPopup)
            shapeGalleryPopup.Visible = 'off';
        end
        phiCut = strcmp(S.cutMode,'Phi cut (fixed theta)');
        setEnabled(spCutTheta,phiCut);
        setEnabled(spCutPhi,~phiCut && ~S.cutPhiFollow);
        setEnabled([bCut0 bCut90 findall(fig,'Tag','btnCutFollow')],~phiCut);
        syncSLLControl();
    end

    function setEnabled(h,yes)
        if yes, set(h,'Enable','on'); else, set(h,'Enable','off'); end
    end

    function refreshSpacingMM()
        if isempty(lblSpacingMM) || ~isgraphics(lblSpacingMM), return; end
        lblSpacingMM.Text = sprintf('dx × dy = %.4f × %.4f mm\nDesign frequency: %g %s', ...
            S.dx*lambdaMM(),S.dy*lambdaMM(),S.freqGHz*freqScale(),S.freqUnit);
        lblSpacingMM.UserData = [S.dx S.dy S.subOff]*lambdaMM();
        if strcmp(S.mode,'Sub-position (9-pt)')
            lblSpacingMM.Text = sprintf('%s · offset %.4f mm', ...
                lblSpacingMM.Text,S.subOff*lambdaMM());
        end
    end

    function showHelp(which)
        switch which
            case 'guide'
                web('https://github.com/msgokdol/antenna-array-designers/blob/main/docs/phased-array-guide.md','-browser');
                return;
            case 'about'
                heading = 'About Phased Array Designer';
                txt = sprintf(['%s\nVersion %s\n\nMuhammed Said Gökdöl\n' ...
                    '<saidgokdol@gmail.com>\n\nMATLAB %s (%s)'], ...
                    APP_NAME,APP_VERSION,version,version('-release'));
            otherwise
                heading = 'Quick start';
                txt = ['1. Array: choose a template or set rows, columns and spacing (λ). ' ...
                    'The mm readout uses the design frequency.' newline newline ...
                    '2. Element: choose a model or import a full-sphere CST realized-gain ' ...
                    'ASCII export with complex Theta/Phi components. Select unit cell only ' ...
                    'for an embedded-element export.' newline newline ...
                    '3. Beam: set steering angles and taper. Angle convention changes labels ' ...
                    'and coordinates, not the physical direction. Beam squint fixes the ' ...
                    'phases at the design frequency.' newline newline ...
                    '4. View: compute, choose a cut and its polarization. Pin a reference ' ...
                    'or load a whole-array CST overlay to compare. 2D pattern opens ' ...
                    'polar azimuth/elevation views and a signed-U cut. Analysis offers scan loss, ' ...
                    'band sweep and scan coverage.' newline newline ...
                    '5. Export: save the design; use CST .tsv for element excitations or ' ...
                    'CST macro for physical feed ports. The Phase scheme map also exports macros.'];
        end
        uialert(fig,txt,heading,'Icon','info');
    end

    function applyTemplate(which)
        % Deck slides 2 and 30 specify 17.7–21.2 GHz, design 21.4 GHz,
        % 0.45λ pitch and uniform amplitude. Slide 32 projects a 16×16.
        % The analytic patch is explicitly a placeholder, not the CST model.
        undoCheckpoint();
        docNextStep = struct('file',docFile,'saved',docSaved);
        candidate = defaultDesign;
        if strcmp(which,'ka')
            candidate.M = 16; candidate.N = 16;
            candidate.dx = 0.45; candidate.dy = 0.45;
            candidate.freqGHz = 21.4; candidate.freqOpGHz = 19.45;
            candidate.efType = 'Patch (cos^q x lin pol)';
        elseif strcmp(which,'blank')
            candidate.mode = 'Sparse (click cells)';
        end
        for kt = 1:numel(UNDO_FIELDS)
            S.(UNDO_FIELDS{kt}) = candidate.(UNDO_FIELDS{kt});
        end
        S.sel = []; S.pending = [];
        if strcmp(which,'blank'), S.el = zeros(0,5); S.elRC = zeros(0,2);
        else, rebuildUniform(); end
        syncDesignControls(); refreshAll();
        % A preset starts an untitled document. Its predecessor remains in Undo.
        docFile = ''; docSaved = defaultDesign; docNextStep = []; refreshDocTitle();
        if strcmp(which,'ka')
            reportAnalysis(['Ka-band template: 17.7–21.2 GHz band; analytic patch ' ...
                'approximation — import your CST element for calibrated results']);
        else
            reportAnalysis('Template applied — Undo restores the previous design');
        end
    end

    function styleElementTable()
        %STYLEELEMENTTABLE  Show which table columns can be edited.
        %   The computed columns (position and total feed phase) get a
        %   shaded background, x and y also muted text; the editable ones
        %   stay on the theme's plain cells, like an input field. Column
        %   styles survive every tbl.Data update (rows come and go with
        %   the design), but their colours are explicit, so paintChrome
        %   re-runs this on a theme change.
        %
        %   Values are CENTRED. The table's scrollbar is an overlay on
        %   macOS: it takes no width and sits on top of the last column's
        %   right edge, where right-aligned digits disappeared under it.
        %   uistyle has no padding, and centring every column keeps the
        %   grid consistent.
        %
        %   Only this function's own styles are replaced, so a style
        %   another feature puts on the table survives a theme change.
        %   uistyle is a value object with no handle identity: they are
        %   found again by comparing against the copies kept in
        %   elemTblStyles.
        P = pal();
        cfg = tbl.StyleConfigurations;
        if ~isempty(cfg) && ~isempty(elemTblStyles)
            mine = false(height(cfg),1);
            for kStyle = 1:height(cfg)
                mine(kStyle) = string(cfg.Target(kStyle)) == "column" ...
                    && any(arrayfun(@(s) isequal(s, cfg.Style(kStyle)), ...
                    elemTblStyles));
            end
            if any(mine), removeStyle(tbl, find(mine)); end
        end
        % The total feed phase is the number sent to the hardware, so it
        % keeps full-strength text; only the positions are muted.
        elemTblStyles = [uistyle('HorizontalAlignment','center'), ...
            uistyle('BackgroundColor',P.readOnlyBg), ...
            uistyle('FontColor',P.muted)];
        addStyle(tbl, elemTblStyles(1), 'column', 1:numel(tbl.ColumnEditable));
        addStyle(tbl, elemTblStyles(2), 'column', find(~tbl.ColumnEditable));
        addStyle(tbl, elemTblStyles(3), 'column', [1 2]);
    end

    function [ax, note] = scanLossWindow()
        %SCANLOSSWINDOW  Reuse the main workspace pane for each sweep.
        win = S.scanLossWin;
        if isempty(win) || ~isgraphics(win)
            win = resultPane('Scan Loss','scanLossWindow');
            g = uigridlayout(win,[2 1]);
            g.RowHeight = {'1x','fit'};
            uiaxes(g,'Tag','scanAxes');
            note = uilabel(g,'Tag','scanNote','Text','','WordWrap','on');
            % Muted caption: paintChrome recolours it on a theme change.
            mutedLbls = [mutedLbls(isgraphics(mutedLbls)), note];
            S.scanLossWin = win;
        else
            selectTask('Scan Loss');
        end
        ax = findall(win,'Tag','scanAxes');
        note = findall(win,'Tag','scanNote');
    end

    function plotScanLoss(sweepAngles, sweepDisplay, curveDb, isUnitCellMode)
        %PLOTSCANLOSS  Draw a finished scan-loss sweep in its window.
        %   sweepAngles are the polar steering angles swept (deg), shown
        %   at sweepDisplay on the x axis (signed theta, or elevation in
        %   the azimuth/elevation convention); curveDb is relative to
        %   broadside. A short title names the scan plane, the subtitle
        %   the design, and the caveats go in a wrapped note under the
        %   plot -- the old single-line title ran off both window edges.
        [ax, note] = scanLossWindow();
        P = pal();
        cla(ax);
        % The idealized aperture projection is a comparison only; the
        % blue trace below is the computed scan loss for this design.
        cosCurve = 10*log10(max(cosd(sweepAngles),1e-6));
        hold(ax,'on');
        if useAzEl()
            refName = 'sin(elevation)  [idealized aperture]';
        else
            refName = 'cos(\theta_s)  [idealized aperture]';
        end
        hRef = plot(ax, sweepDisplay, cosCurve, '-.', 'Color',P.traceAR, ...
            'LineWidth',1.4, 'Tag','scanCosRef', 'DisplayName',refName);
        hTotal = plot(ax, sweepDisplay, curveDb, '-', 'Color',P.traceTotal, ...
            'LineWidth',2, 'Tag','scanTotal', ...
            'DisplayName','Total  [vector-field model]');
        grid(ax,'on');
        legend(ax,[hRef hTotal],'Location','southwest','AutoUpdate','off');
        % Limits always take in 0 dB and the -3 dB guide, so a shallow
        % curve is not stretched to fill the axes and read as a steep one.
        yl = [min([curveDb(:); cosCurve(:); -3]), ...
            max([curveDb(:); cosCurve(:); 0])];
        yPad = 0.08*diff(yl);
        ylim(ax, [yl(1)-yPad, yl(2)+yPad]);
        xlim(ax, sort(sweepDisplay([1 end])));
        yline(ax, -3, '--', '-3 dB', 'Tag','scanRef3dB', 'Color',P.muted, ...
            'LineWidth',1, 'LabelHorizontalAlignment','left');
        % Mark the design's own steering angle on the curve, when the
        % sweep reaches it. Its value is interpolated between the 2 deg
        % sweep points, which is why it is quoted to 0.1 dB only.
        thNow = abs(S.theta_s);
        if thNow <= sweepAngles(end)
            if useAzEl()
                xNow = 90 - thNow;
                nowTxt = sprintf('elevation %.4g°', xNow);
            else
                xNow = S.theta_s;
                nowTxt = sprintf('\\theta_s = %.4g°', xNow);
            end
            dbNow = interp1(sweepAngles, curveDb, thNow);
            % The label goes on the side where the curve falls away, at
            % the top, clear of the trace; near the far end there is no
            % room there, so it moves toward broadside, at the bottom.
            awayRight = sweepDisplay(end) > sweepDisplay(1);
            reach = abs(xNow - sweepDisplay(1)) / ...
                abs(sweepDisplay(end) - sweepDisplay(1));
            if reach < 0.6
                vAlign = 'top'; toRight = awayRight;
            else
                vAlign = 'bottom'; toRight = ~awayRight;
            end
            xline(ax, xNow, ':', sprintf('Current steering: %s, %.1f dB', ...
                nowTxt, dbNow), 'Tag','scanNow', 'Color',P.ink, ...
                'LineWidth',1.5, 'LabelOrientation','horizontal', ...
                'LabelVerticalAlignment',vAlign, ...
                'LabelHorizontalAlignment',ternStr(toRight,'right','left'));
        end
        hold(ax,'off');

        % Title: what and in which plane. The subtitle: the design.
        if useAzEl()
            planeTxt = sprintf('azimuth held at %.0f°', scanPhiVal());
        else
            planeTxt = sprintf('\\phi_s = %.0f°', S.phi_s);
            if S.theta_s < 0, planeTxt = [planeTxt ', negative-\theta side']; end
        end
        nEl = size(S.el,1);
        % Sampled AT THE COMMANDED direction, not at the moving peak, so
        % with frozen phases the curve includes the pointing loss from
        % beam squint. Said on the plot: otherwise it can look
        % inconsistent with a main-window figure taken at the squinted
        % beam, though the two answer different questions.
        if S.retunePhase
            freqTxt = sprintf('phases retuned at %.4g GHz', S.freqOpGHz);
        elseif nEl <= 1
            freqTxt = 'single element, no array squint';
        else
            freqTxt = sprintf(['phases frozen at %.4g GHz, includes squint ' ...
                'pointing loss'], phaseReferenceGHz());
        end
        title(ax, ['Scan loss · ' planeTxt], sprintf('%d %s (%s) · %s', ...
            nEl, ternStr(nEl == 1,'element','elements'), S.efType, freqTxt));
        ax.Subtitle.Interpreter = 'none';   % element names hold ^ and _
        if useAzEl()
            xlabel(ax, 'Commanded elevation (°)');
        else
            xlabel(ax, 'Commanded steering angle \theta_s (°)');
        end
        % Named for what is actually plotted, which differs by mode.
        % Labelling both "directivity" is how the disagreement with the
        % main window stayed invisible.
        if isUnitCellMode
            ylabel(ax, 'Gain relative to broadside (dB)');
            note.Text = ['Gain of the total field toward the commanded ' ...
                'direction, relative to the broadside beam, on the ' ...
                'unit-cell basis: scored on radiated intensity, because ' ...
                'the pattern integral of a periodic unit cell is not its ' ...
                'radiated power.'];
        else
            ylabel(ax, 'Directivity relative to broadside (dB)');
            note.Text = ['Directivity of the total field toward the ' ...
                'commanded direction, relative to the broadside beam, ' ...
                'from a full-sphere power integral. Radiation efficiency, ' ...
                'mismatch and mutual coupling are not modelled.'];
        end
        note.FontColor = P.muted;
        note.Tooltip = sprintf(['What the curve measures and what it ' ...
            'leaves out. Steering angles %g-%g° in %g° steps.'], ...
            sweepAngles(1), sweepAngles(end), sweepAngles(2)-sweepAngles(1));
        selectTask('Scan Loss');
    end

    function refreshInfo()
        %REFRESHINFO  Fill the result cards and the Details table from S.
        %   Reads only state (the last compute's results and the design),
        %   so it is cheap to call after any change; computePattern calls
        %   it last, once every figure it reports has been written.
        if ~patternReady(false)
            showDetailsMessage(['Imported pattern required. Radiation ' ...
                'results and phase exports are unavailable.']);
            return;
        end
        nEl = size(S.el,1);
        if nEl == 0
            showDetailsMessage('No elements placed.');
            return;
        end
        % taper efficiency from |w|: measures amplitude-taper aperture
        % utilisation only, deliberately ignoring phase reversals (so a
        % difference / monopulse design reads its amplitude efficiency,
        % not 0). The directivity readout captures the phase cost.
        ap = sum(abs(S.el(:,3)));
        metricW=calculationWeights();
        den = nEl*sum(metricW.^2);
        if den > 0, eff = sum(abs(metricW))^2/den; else, eff = NaN; end
        if isempty(S.sel)
            selTxt = 'none (Move / Rotate / Reset apply to ALL elements)';
        else
            % sprintf, not [char string] concatenation: mixing a char
            % array with a MATLAB `string` via square brackets silently
            % creates a 2-element string ARRAY instead of one merged
            % string, which desynced the sprintf below and crashed it.
            selTxt = sprintf('#%s', strjoin(string(S.sel(:)'),', '));
        end
        if isnan(eff)
            % Reachable ONLY when every amplitude is zero: den is
            % nEl*sum(w.^2), a sum of squares, so it can never be
            % negative and is zero only in that one case. The label used
            % to read "mixed-sign amplitudes", which this branch cannot
            % actually detect -- the numerator sum(abs(w)) is
            % deliberately sign-insensitive precisely so a difference/
            % monopulse design still reports its amplitude efficiency
            % instead of collapsing to 0, so mixed signs always land in
            % the normal branch below.
            effTxt = 'N/A (all amplitudes are zero)';
        else
            % eff <= 1 (Cauchy-Schwarz), so the loss is |10 log10 eff|;
            % written so, a uniform taper reads 0.00 dB rather than -0.00.
            effTxt = sprintf('%.3f (%.2f dB loss)', eff, abs(10*log10(eff)));
        end
        % Nearest-neighbour centre-to-centre distance: the actual thing
        % that determines whether physical elements (patches etc.) will
        % overlap after Row stagger / Alt-row rotation bring diagonal
        % neighbours closer together -- dx,dy alone no longer tell you
        % this once a stagger offset is in play. Brute-force O(n^2) pair
        % distances; capped at 2000 elements so it doesn't slow down
        % auto-recompute on a large array. Deliberately a LOOSER cap
        % than checkLargeArrayWarning's 1500 (this comment used to claim
        % they were the same threshold, which they are not): that one
        % warns about computePattern's full-grid recompute, while this
        % is a cheap O(n^2) distance matrix that stays comfortable a bit
        % further up.
        if nEl > 1 && nEl <= 2000
            ddx = S.el(:,1) - S.el(:,1).';
            ddy = S.el(:,2) - S.el(:,2).';
            distAll = hypot(ddx,ddy);
            distAll(1:nEl+1:end) = Inf;   % exclude each element's distance to itself
            nnDist = min(distAll(:));
            nnTxt = sprintf('%.4f λ (design)', nnDist);
        elseif nEl > 2000
            nnTxt = 'skipped (> 2000 elements)';
        else
            nnTxt = 'N/A (single element)';
        end
        % Grating-lobe readout, from the full-hemisphere search done in
        % computePattern (findGratingLobe) -- watch THIS number as you
        % adjust Row stagger / Grid angle to see the effect directly:
        % staggering pushes the grating lobe to a larger theta and/or
        % pushes glRelDb further below the main peak (weaker), since it
        % changes the array from a single periodic lattice into a
        % coarser lattice with a 2-point basis whose structure factor
        % suppresses some of what would otherwise be closer-in lobes.
        % glState drives the Grating lobes card: '' not computed, 'none',
        % 'grazing' or 'possible'.
        if isnan(S.glTheta)
            if isinf(S.glRelDb)
                glState = 'none';
                glRows = {'Secondary lobe', ...
                    'none (one main lobe fills the whole hemisphere)', ''};
            else
                glState = '';
                glRows = {'Secondary lobe', ...
                    'not yet computed (click Compute pattern)', ''};
            end
        else
            tag = 'highest sidelobe'; glState = 'none'; glTone = '';
            if S.glRelDb > -6
                if S.glInterior
                    tag = 'possible grating lobe';
                    glState = 'possible'; glTone = 'bad';
                else
                    % Strong, but its peak is outside visible space --
                    % the horizon is cutting through the skirt. Naming it
                    % separately keeps this panel and the Max scan tool
                    % telling the same story.
                    tag = 'strong lobe grazing the horizon';
                    glState = 'grazing'; glTone = 'warn';
                end
            end
            % One row, so the verdict and the lobe it is about stay
            % together however the value column wraps.
            glRows = {'Secondary lobe', sprintf('%s at %s, %.1f dB below peak', ...
                tag, angleText(S.glTheta,S.glPhi), abs(S.glRelDb)), glTone};
        end
        % ---- Metrics: Area, Ideal Directivity, Calculated Directivity,
        % Aperture Efficiency -- standard aperture-antenna figures of
        % merit. Area = nElements * dx * dy * sin(gridAngle): each
        % element "owns" one unit-cell parallelogram of area dx*dy*sin
        % (gridAngle) -- the sin(gridAngle) is the general parallelogram-
        % area correction for a skewed (non-90 deg) lattice, reducing to
        % the familiar dx*dy at gridAngle=90. Shape-masked arrays are
        % handled correctly for free since nElements already reflects
        % however many cells the Diamond/Hexagon/Circle/etc mask kept.
        % Ideal Directivity is the theoretical max for a UNIFORMLY
        % illuminated aperture of that area: D_ideal = 4*pi*Area/lambda^2
        % (Area already in units of lambda^2, so just 4*pi*Area).
        % Calculated Directivity is this app's own computed TOTAL
        % directivity (S.DpkTot, persisted from computePattern) --
        % includes taper loss, scan/foreshortening loss, and element
        % pattern effects, so it's normally somewhat below ideal.
        % Aperture Efficiency = Calculated/Ideal (the gap between them).
        % sind(gridAngle) only applies in Uniform grid mode -- Sparse and
        % Sub-position modes always place elements on the plain
        % orthogonal dx,dy lattice regardless of Grid angle (see
        % onGridClick/onGeom), so applying the skew factor there would
        % misreport Area/Ideal directivity if Grid angle was left
        % non-90 from an earlier Uniform-grid session.
        % Counted per LATTICE CELL, not per element. Sub-position mode can
        % place more than one element inside a single dx-by-dy cell, and
        % nEl*dx*dy then charged a full cell to each occupant: adding a
        % second element inside one cell doubled the reported area while
        % the physical boundary did not move.
        if isempty(S.elRC)
            nCell = nEl;
        else
            nCell = size(unique(S.elRC,'rows'),1);
        end
        if strcmp(S.mode,'Uniform grid')
            areaLam2 = nCell * S.dx * S.dy * sind(S.gridAngle);
        else
            areaLam2 = nCell * S.dx * S.dy;
        end
        % dx/dy are DESIGN wavelengths, so the area above is in design
        % lambda^2 -- a statement about geometry. The ideal-directivity
        % bound 4*pi*A/lambda^2 needs the aperture in OPERATING
        % wavelengths, which differ by freqRatio() whenever operating and design frequencies differ. Without
        % the scaling the bound stayed put while the real one moves by
        % 20*log10(fr) per dimension: doubling the operating frequency on
        % fixed geometry is +6.02 dB that the readout never showed.
        frMet = freqRatio();
        areaOpLam2 = areaLam2 * frMet^2;
        idealDirDb = 10*log10(max(4*pi*areaOpLam2, eps));
        % Unit-cell mode gets its OWN metrics block rather than the generic
        % one. "Calculated directivity" and the aperture efficiency derived
        % from it are both built on S.DpkTot -- the directivity of the
        % pattern shape -- which is not the embedded element's directivity
        % and would sit in the panel contradicting the figures below it.
        % The first card follows the same split: array gain in unit-cell
        % mode, realized gain for an imported gain export (Details lists
        % both D and G there), directivity otherwise.
        uc = unitCellFigures();
        isUC = ~isempty(uc);
        cm = blankCardModel();
        if isUC
            % The file's peak IS the embedded element's realized gain, and
            % its directivity follows from the efficiency you typed:
            %   D = RG - 10log10(eff)
            % That is the same route CST uses internally (its own Dir.
            % readout equals its realized gain divided by Tot. Effic.).
            % 4pi*Area/lambda^2 is listed alongside the array figures as
            % the geometric ceiling they should approach, not as a
            % competing answer: it knows only the aperture size, while the
            % array gain/directivity above it carry the element's actual
            % efficiency and the amplitude taper.
            metricRows = {'Area', sprintf('%.2f λ² (design)',areaLam2), ''
                'Unit-cell realized gain', sprintf('%.3f dBi (file peak)',uc.rg), ''
                'Efficiency (entered)', sprintf('%.1f %%',uc.effPct), ''
                'Unit-cell directivity', sprintf('%.3f dBi (RG − 10 log10 eff)',uc.d), ''
                'Array', sprintf('%d elements, steered %s',uc.n, ...
                    angleText(S.theta_s,S.phi_s)), ''
                'Effective elements', sprintf('%.2f (taper)',uc.nEff), ''};
            cm.gainCap = 'Array gain (dBi)';
            if uc.beamVisible
                if uc.n <= 1
                    evalKind = 'single element';
                elseif S.retunePhase
                    evalKind = 'commanded';
                else
                    evalKind = 'squint prediction';
                end
                evalTxt = sprintf('%s (%s)', ...
                    angleText(uc.evalThetaDeg,uc.evalPhiDeg), evalKind);
                metricRows = [metricRows
                    {'Evaluated toward', evalTxt, ''
                    'Element gain there', sprintf('%.2f dBi (%.2f dB below peak)', ...
                        uc.rgSteer, uc.scanLossDb), ''
                    'Est. array gain there', sprintf('%.2f dBi',uc.gArr), ''
                    'Est. array D there', sprintf('%.2f dBi',uc.dArr), ''}];
                if isfinite(uc.gArr)
                    cm.gainDb = uc.gArr;
                    cm.gainNote = 'estimate, embedded';
                else
                    metricRows(end+1,:) = {'Note', ['no radiation toward the ' ...
                        'evaluation direction: gain and directivity there are unavailable'], 'warn'};
                    cm.gainNote = 'no radiation there';
                end
                cm.gainTip = sprintf(['Estimated realized gain of the whole array ' ...
                    'toward %s (dBi): the embedded element''s gain there plus the ' ...
                    'coherent array build-up.'], evalTxt);
            else
                metricRows = [metricRows
                    {'Steering-term prediction', 'outside visible space', 'warn'
                    'Directional gain', 'N/A (not clamped to grazing)', ''}];
                cm.gainNote = 'beam not visible';
                cm.gainTip = ['The steering-term prediction lies outside visible ' ...
                    'space, so there is no direction to estimate the array gain toward.'];
            end
            metricRows(end+1,:) = {'Aperture ceiling', ...
                sprintf('%.2f dBi (4πA/λ², broadside)',idealDirDb), ''};
            if uc.beamVisible && isfinite(uc.dArr)
                % Same definition as the generic row (directivity over the
                % 4*pi*A/lambda^2 ceiling), on the array's estimated
                % directivity: for an ideal unit cell the two meet.
                metricRows(end+1,:) = {'Aperture efficiency', sprintf( ...
                    '%.1f%% (est. array D vs 4πA/λ²)', ...
                    10^((uc.dArr - idealDirDb)/10)*100), ''};
            end
        elseif isnan(S.DpkTot) || S.DpkTot <= 0
            metricRows = {'Area', sprintf('%.2f λ² (design)',areaLam2), ''
                'Ideal directivity', sprintf('%.2f dBi',idealDirDb), ''
                'Calculated directivity', 'N/A (click Compute pattern)', ''};
            cm.gainNote = 'not computed';
        else
            calcDirDb = 10*log10(S.DpkTot);
            apEffDb = calcDirDb - idealDirDb;
            apEffPct = 10^(apEffDb/10)*100;
            metricRows = {'Area', sprintf('%.2f λ² (design)',areaLam2), ''
                'Ideal directivity', sprintf('%.2f dBi',idealDirDb), ''
                'Calculated directivity', sprintf('%.2f dBi',calcDirDb), ''};
            cm.gainDb = calcDirDb;
            cm.gainNote = 'total-field peak';
            cm.gainTip = ['Directivity of the whole array at its beam peak (dBi), ' ...
                'from the full-sphere power integral.'];
            effImp = impEffLin();
            if strcmp(S.efType,'Imported (CST far-field)') && effImp < EFF_TOL
                % The imported file carries its own losses, so the number
                % that matters is gain: the card shows it, and the target
                % "minimum directivity / gain" is checked against it.
                gDb = 10*log10(S.DpkTot*effImp);
                metricRows(end+1,:) = {'Realized gain', ...
                    sprintf('%.2f dBi (efficiency %.1f %%)',gDb,100*effImp), ''};
                cm.gainCap = 'Gain (dBi)';
                cm.gainDb = gDb;
                cm.gainNote = sprintf('eff %.1f %%',100*effImp);
                cm.gainTip = sprintf(['Realized gain at the beam peak (dBi): ' ...
                    'directivity %.2f dBi times the imported element''s ' ...
                    'efficiency (%.1f %%).'], calcDirDb, 100*effImp);
            end
            metricRows(end+1,:) = {'Aperture efficiency', ...
                sprintf('%.1f%% (%.2f dB)',apEffPct,apEffDb), ''};
            if ~strcmp(S.mode,'Uniform grid')
                metricRows(end+1,:) = {'Note', ['aperture efficiency uses the ' ...
                    'occupied-cell area as its reference, not a physical footprint'], ''};
            elseif apEffPct > 100.5
                metricRows(end+1,:) = {'Note', ['above 100 %: the filled-cell ' ...
                    'area is not a physical bound for this configuration'], 'warn'};
            end
            if S.efFallback
                cm.gainNote = 'element substituted'; cm.gainTone = 'warn';
            end
        end
        % Beam figures measured on the principal cut; HPBW and SLL are on
        % the cards, the rest of the cut's story is here.
        cutRows = cell(0,3);
        cmC = S.cutMetrics;
        if cmC.measured
            planeTxt = cmC.plane;
            if cmC.offBeam, planeTxt = [planeTxt ' (off the beam)']; end
            if isfinite(cmC.fnbw)
                fnbwTxt = sprintf('%.2f°',cmC.fnbw);
            else
                fnbwTxt = 'N/A (no null on both sides)';
            end
            cutRows = {'Cut plane', planeTxt, ''
                'FNBW', fnbwTxt, ''};
        end
        % This is center-to-center span, deliberately not called aperture:
        % the Area metric uses occupied unit-cell area, which includes the
        % half-cell extent beyond the outer element centers on a uniform
        % lattice. The author credit that used to close this panel lives
        % in Help > About now: this table is engineering results only.
        rows = [{'Elements', sprintf('%d',nEl), ''
            'Sum of amplitudes', sprintf('%.2f',ap), ''
            'Amplitude taper', taperDetailText(), ''
            'Taper efficiency', effTxt, ''
            'Element-center span', sprintf('%.2f × %.2f λ (design)', ...
                range0(S.el(:,1)), range0(S.el(:,2))), ''
            'Nearest spacing', nnTxt, ''}
            metricRows; cutRows; glRows; frequencyRows()
            {'Selected', selTxt, ''}];
        setDetails(rows);
        refreshCards(cm, glState, S.cutMetrics);
    end

    % ------------------------------------------- sidelobe level (B3)
    function assignSLL(v)
        % A new design sidelobe level goes the way a taper change does
        % (applyTaper, then refreshAll), so the table, the result cards
        % and Auto-compute follow it. With any other taper the value is
        % only kept -- the spinner is disabled then -- and reshaping
        % those amplitudes would be a pointless recompute.
        S.sll = min(max(double(v),13),80);
        syncSLLControl();
        if any(strcmp(S.taper,{'Chebyshev','Taylor'}))
            applyTaper(); refreshAll();
        end
    end

    function syncSLLControl()
        %SYNCSLLCONTROL  Spinner value, enable state and n̄ note from S.
        %   Called wherever S.taper or S.sll changes (taper dropdown,
        %   table edit -> Manual, taper fallback, loadConfig), so the
        %   control never offers a level the active taper ignores.
        if ~isgraphics(spSLL), return; end
        designed = any(strcmp(S.taper,{'Chebyshev','Taylor'}));
        spSLL.Value = S.sll;
        if designed, spSLL.Enable = 'on'; else, spSLL.Enable = 'off'; end
        if strcmp(S.taper,'Taylor')
            lblNbar.Text = sprintf('n̄ = %d (auto)', taylorNbar(S.sll));
        else
            lblNbar.Text = '';
        end
    end

    function nbar = taylorNbar(sllDb)
        % nbar = the number of near-in sidelobes a Taylor taper holds at
        % the design level. The standard requirement for reaching a
        % level at all is nbar >= 2*A^2 + 0.5, A = acosh(10^(SLL/20))/pi;
        % the smallest such integer (at least 2) is used. Why it must be
        % derived rather than fixed: see the Taylor case in taperVec.
        A = acosh(10^(sllDb/20))/pi;
        nbar = max(2, ceil(2*A^2 + 0.5));
    end

    function txt = taperDetailText()
        % The Details table's taper row: which taper, and for the two
        % designed ones the level (and Taylor's n̄) it was built for --
        % the target to hold the SLL card against. (No "lobe" in this
        % text: the grating-lobe row is found by that word.)
        switch S.taper
            case 'Chebyshev'
                txt = sprintf('Chebyshev, designed for -%g dB', S.sll);
            case 'Taylor'
                txt = sprintf('Taylor, designed for -%g dB, n̄ = %d', ...
                    S.sll, taylorNbar(S.sll));
            otherwise
                txt = S.taper;
        end
    end

    function loadCstOverlay()
        [file,folder] = uigetfile({'*.txt;*.dat;*.csv','CST far-field ASCII'}, ...
            'Load CST result of the WHOLE array');
        if isequal(file,0), return; end
        try
            data = loadCstFarfieldASCII(fullfile(folder,file));
            [ghz,~] = cstFileFrequency(file,fullfile(folder,file));
            S.cstOverlay = struct('ff',data,'name',file,'GHz',ghz);
            cstComparisonClosed = false;
            redrawCstOverlay();
            reportAnalysis(['CST array comparison loaded: ' file]);
        catch err
            setStatus(['CST overlay failed: ' err.message],'bad');
            uialert(fig,err.message,'CST overlay');
        end
    end

    function clearCstOverlay()
        S.cstOverlay = [];
        delete(findall(axCut,'Tag','cutCst'));
        if isgraphics(overlayWin)
            set(findall(overlayWin,'Tag','overlayName'),'Text','No CST array loaded');
            set(findall(overlayWin,'Tag','overlayNote'),'Text', ...
                'Load a whole-array CST far-field file on the View tab.');
            set(findall(overlayWin,'Tag','overlayTable'),'Data',cell(0,4));
            cla(findall(overlayWin,'Tag','overlayAxes'));
            overlayWin.UserData = [];
        end
        redrawRefTrace();
        caveats=cellstr(string(axCut.Subtitle.String));
        caveats=caveats(~startsWith(caveats,'CST gain overlay hidden'));
        setPlotTitle(axCut,axCut.Title.String,caveats);
        setStatus('CST overlay cleared.','good');
    end

    function redrawCstOverlay()
        delete(findall(axCut,'Tag','cutCst'));
        if isempty(S.cstOverlay), return; end
        if ~S.radiationValid || ~strcmp(S.cutView.kind,'gain') || isempty(S.cutX)
            if isgraphics(overlayWin)
                set(findall(overlayWin,'Tag','overlayNote'),'Text', ...
                    'A gain cut is required. Compute the pattern and show gain to compare it with CST.');
            end
            if strcmp(S.cutView.kind,'ar')
                addCaveat(axCut,'CST gain overlay hidden while axial ratio is shown');
            end
            return;
        end
        [et,ep] = evalImportedFF(S.cstOverlay.ff,S.cutTheta,S.cutPhiSamples);
        cdb = 20*log10(abs(polCombine(et,ep))) + S.cstOverlay.ff.peak(1);
        mdb = S.cutRaw;
        if ~any(isfinite(cdb)) || ~any(isfinite(mdb)), return; end
        % Analytic models report directivity, not realized gain. Compare
        % their shapes and say so; never silently subtract G from D.
        absolute = S.cutView.abs && strcmp(S.efType,'Imported (CST far-field)');
        offset = 0;
        if ~absolute
            offset = max(mdb);
            mdb = mdb-max(mdb); cdb = cdb-max(cdb);
        end
        P = pal(); wasHeld = ishold(axCut); hold(axCut,'on');
        plot(axCut,S.cutX,max(cdb+offset,min(S.cutDb)), '--', ...
            'Color',P.overlayTrace,'LineWidth',1.6,'Tag','cutCst', ...
            'DisplayName','CST array','HitTest','off');
        if ~wasHeld, hold(axCut,'off'); end
        axCut.YLim(2)=max(axCut.YLim(2),max(cdb+offset)+2);
        redrawRefTrace();
        if cstComparisonClosed, return; end
        periodic = S.cutView.isPhi || S.fullSphere;
        mat = comparisonMetrics(S.cutX,mdb,periodic);
        cst = comparisonMetrics(S.cutX,cdb,periodic);
        delta = cst-mat;
        delta(1) = mod(delta(1)+180,360)-180;
        % The app uses signed theta coordinates for an elevation cut.
        % Apply exactly the same display conversion to the table.
        azRows = cell(0,4);
        angleName = 'Peak signed θ (°)';
        if S.cutView.isPhi, angleName = 'Peak azimuth φ (°)';
        elseif useAzEl()
            maz=mod(S.cutView.plane+180*(mat(1)<0),360);
            caz=mod(S.cutView.plane+180*(cst(1)<0),360);
            if abs(mat(1))<1e-9, maz=NaN; end
            if abs(cst(1))<1e-9, caz=NaN; end
            azRows={'Peak azimuth (°)',maz,caz,mod(caz-maz+180,360)-180};
            mat(1) = 90-abs(mat(1)); cst(1) = 90-abs(cst(1));
            delta(1) = cst(1)-mat(1); angleName = 'Peak elevation (°)';
        end
        levelName = 'Peak relative level (dB)';
        note = 'Peak-normalized cut shapes; Δ = CST − MATLAB. Relative peaks are 0 dB.';
        if absolute
            levelName = 'Peak realized gain (dBi)';
            note = 'Absolute realized gain; Δ = CST − MATLAB.';
        elseif S.cutView.abs
            note = [note ' CST is aligned to the MATLAB peak on the plot; MATLAB shows directivity.'];
        end
        if isfinite(S.cstOverlay.GHz) && abs(S.cstOverlay.GHz-S.freqOpGHz)>1e-6
            note = sprintf('%s Frequency mismatch: CST %g GHz, MATLAB %g GHz.', ...
                note,S.cstOverlay.GHz,S.freqOpGHz);
        end
        if isempty(overlayWin) || ~isgraphics(overlayWin)
            overlayWin = resultPane('CST Compare','cstOverlayWin');
            og = uigridlayout(overlayWin,[4 1]);
            og.RowHeight = {25,'1x','0.55x',65};
            uilabel(og,'Text','','Tag','overlayName','Interpreter','none');
            uiaxes(og,'Tag','overlayAxes');
            uitable(og,'Tag','overlayTable','ColumnName', ...
                {'Metric','MATLAB','CST','Δ CST − MATLAB'}, ...
                'ColumnWidth',{'1x',115,115,130},'RowName',{}, ...
                'Tooltip','Same-cut metrics; NaN means no defined crossing or sidelobe');
            uilabel(og,'Text','','Tag','overlayNote','WordWrap','on');
        end
        set(findall(overlayWin,'Tag','overlayName'),'Text',S.cstOverlay.name);
        set(findall(overlayWin,'Tag','overlayNote'),'Text',note);
        overlayAxes = findall(overlayWin,'Tag','overlayAxes');
        cla(overlayAxes);
        plot(overlayAxes,S.cutX,mdb,'LineWidth',1.8, ...
            'Color',P.traceTotal,'DisplayName','MATLAB', ...
            'Tag','overlayMatlab');
        hold(overlayAxes,'on');
        plot(overlayAxes,S.cutX,cdb,'--','LineWidth',1.7, ...
            'Color',P.overlayTrace,'DisplayName','CST', ...
            'Tag','overlayCst');
        hold(overlayAxes,'off');
        grid(overlayAxes,'on');
        legend(overlayAxes,'Location','best');
        if S.cutView.isPhi
            xlabel(overlayAxes,'Azimuth (°)');
        elseif useAzEl()
            xlabel(overlayAxes,'Signed elevation (°)');
        else
            xlabel(overlayAxes,'Signed theta (°)');
        end
        ylabel(overlayAxes,ternStr(absolute,'Realized gain (dBi)','Relative level (dB)'));
        set(findall(overlayWin,'Tag','overlayTable'),'Data', ...
            [azRows; [{angleName;levelName;'HPBW (°)';'SLL (dB)'} num2cell([mat(:) cst(:) delta(:)])]]);
        overlayWin.UserData = struct('matlab',mat,'cst',cst,'difference',delta, ...
            'absolute',absolute);
    end

    function vals = comparisonMetrics(x,db,periodic)
        % Share the cut's null walks and crossing interpolation, including
        % a lobe spanning the azimuth seam and cuts without two crossings.
        [pk,ip] = max(db); n = numel(db)-double(periodic);
        if ip>n, ip=1; end
        [rn,nr,wr] = walkDown(db,ip,1,n,periodic);
        [ln,nl,wl] = walkDown(db,ip,-1,n,periodic);
        [okR,ra] = crossOut(db,x,ip,1,n,periodic,pk-3);
        [okL,la] = crossOut(db,x,ip,-1,n,periodic,pk-3);
        hp=NaN; sl=NaN;
        if okR && okL
            hp=ra-la; if periodic, hp=mod(hp,360); end
        end
        full = strcmp(wr,'circle') || strcmp(wl,'circle') || nr+nl+1>=n;
        mask = true(1,n);
        if periodic
            mask(mod((ip-nl:ip+nr)-1,n)+1)=false;
        else
            mask(ln:rn)=false;
        end
        if ~full && any(mask), sl=max(db(find(mask)))-pk; end %#ok<FNDSB>
        vals=[x(ip) pk hp sl];
    end

    % ------------------------------------------ pinned reference trace (B4)
    function pinReference()
        %PINREFERENCE  Keep the principal cut shown now as the reference.
        %   Copies the plotted Total curve (S.cutX / S.cutDb, the curve the
        %   click probe reads), its beam figures and the cutView it was
        %   drawn under, and draws it at once without a recompute -- until
        %   the design changes it lies on the Total curve. Only a current
        %   gain cut can be pinned: an out-of-date or substituted plot is
        %   not the design, and an axial ratio is not a level to compare.
        if ~S.radiationValid || isempty(S.cutX)
            setStatus(['Nothing to pin: the cut is out of date or not ' ...
                'valid. Compute the pattern first.'],'warn');
            return;
        end
        if ~strcmp(S.cutView.kind,'gain')
            setStatus(['Pin as reference takes the gain cut: untick "Cut ' ...
                'plot shows Axial Ratio" first.'],'warn');
            return;
        end
        replaced = ~isempty(S.refCut);
        S.refCut = struct('x',S.cutX,'db',S.cutDb, ...
            'time',char(datetime('now','Format','HH:mm')), ...
            'desc',refDesignText(), ...
            'planeText',cutPlaneText(S.cutView.isPhi,S.cutFixedTheta), ...
            'hpbw',S.cutMetrics.hpbw,'sll',S.cutMetrics.sll, ...
            'peak',max(S.cutDb),'view',S.cutView);
        redrawRefTrace();
        syncRefControls();
        msg = sprintf('Pinned the %s as reference %s: %s.', ...
            S.refCut.planeText, S.refCut.time, S.refCut.desc);
        if replaced, msg = ['Replaced the reference. ' msg]; end
        setStatus(msg,'good');
    end

    function clearReference()
        if isempty(S.refCut)
            setStatus('No reference trace to clear.');
            return;
        end
        S.refCut = [];
        redrawRefTrace();
        syncRefControls();
        setStatus('Reference trace cleared.','good');
    end

    function syncRefControls()
        % Clear is live only while something is pinned; its tooltip says
        % what, since the legend entry carries only the time.
        pinned = ~isempty(S.refCut);
        onOff = {'off','on'};
        btnClearRef.Enable = onOff{pinned+1};
        menuClearRef.Enable = onOff{pinned+1};
        if pinned
            btnClearRef.Tooltip = sprintf(['Remove the reference trace ' ...
                'pinned at %s (%s; %s)'], S.refCut.time, S.refCut.planeText, ...
                S.refCut.desc);
        else
            btnClearRef.Tooltip = 'Remove the pinned reference trace from the cut';
        end
    end

    function s = refDesignText()
        % The design a reference comes from, in a few words: enough to
        % tell two pins apart in the status line and the Clear tooltip.
        tap = S.taper;
        if any(strcmp(tap,{'Chebyshev','Taylor'}))
            tap = sprintf('%s %g dB', tap, S.sll);
        end
        s = sprintf('%d elements, %s, %g × %g λ, %g GHz, steered %s', ...
            size(S.el,1), tap, S.dx, S.dy, S.freqOpGHz, ...
            angleText(S.theta_s,S.phi_s));
    end

    function h = plotRefTrace()
        %PLOTREFTRACE  The pinned reference on the gain cut, if it applies.
        %   Replaces any earlier copy. Returns the line, or an empty handle
        %   when nothing is pinned, no gain cut is drawn, or the reference
        %   was pinned under a different view (refCaveat says why).
        delete(findall(axCut,'Tag','cutRef'));
        h = gobjects(0);
        if isempty(S.refCut) || isempty(findall(axCut,'Tag','cutTotal')) || ...
                ~isempty(refMismatch())
            return;
        end
        P = pal();
        wasHeld = ishold(axCut);
        hold(axCut,'on');
        h = plot(axCut, S.refCut.x, S.refCut.db, '--', 'LineWidth',1.5, ...
            'Color',P.refTrace, 'DisplayName',['Reference ' S.refCut.time], ...
            'Tag','cutRef', 'HitTest','off');
        if ~wasHeld, hold(axCut,'off'); end
    end

    function why = refMismatch()
        %REFMISMATCH  How the cut differs from the one the reference was
        %   pinned on: '' when the reference applies. The plane, the angle
        %   convention and the level basis each change what a point of the
        %   curve means (another direction, another label, dBi against dB
        %   below the peak). Worded to follow "pinned ".
        v = S.cutView; r = S.refCut.view;
        if ~strcmp(v.kind,'gain')
            why = 'as a gain cut; the axial ratio is shown';
            return;
        end
        parts = {};
        if v.isPhi ~= r.isPhi || ...
                abs(mod(v.plane - r.plane + 180, 360) - 180) > 1e-6
            parts{end+1} = ['on the ' S.refCut.planeText];
        end
        if ~strcmp(v.conv, r.conv)
            if strcmp(r.conv,'Azimuth / elevation')
                parts{end+1} = 'in the az/el convention';
            else
                parts{end+1} = 'in the θ/φ convention';
            end
        end
        if ~strcmp(v.pol,r.pol), parts{end+1} = 'with a different polarization'; end
        if v.full ~= r.full, parts{end+1} = 'with a different angular range'; end
        if v.abs && r.abs && v.gain~=r.gain
            parts{end+1} = 'with a different gain/directivity basis';
        end
        if v.abs ~= r.abs
            if r.abs, parts{end+1} = 'with absolute levels (dBi)';
            else, parts{end+1} = 'with relative levels'; end
        end
        why = strjoin(parts, ', ');
    end

    function txt = refCaveat()
        %REFCAVEAT  The cut's caveat line about the pinned reference.
        %   Shown: how the beam figures moved since the pin (now minus
        %   reference; SLL is relative to each curve's own peak, so a
        %   negative change means lower sidelobes). Hidden: why. '' with
        %   nothing pinned. refLineRe below recognises both forms.
        txt = '';
        if isempty(S.refCut), return; end
        why = refMismatch();
        if ~isempty(why)
            txt = sprintf('Reference %s hidden: pinned %s', S.refCut.time, why);
            return;
        end
        m = S.cutMetrics; r = S.refCut;
        d = {};
        % "+ 0" turns a -0 from rounding into 0, so no change reads +0.
        if isfinite(m.hpbw) && isfinite(r.hpbw)
            d{end+1} = sprintf('HPBW %+.2f°', round(m.hpbw - r.hpbw, 2) + 0);
        end
        if isfinite(m.sll) && isfinite(r.sll)
            d{end+1} = sprintf('SLL %+.1f dB', round(m.sll - r.sll, 1) + 0);
        end
        pk = max(S.cutDb);
        if isfinite(pk) && isfinite(r.peak)
            d{end+1} = sprintf('peak %+.2f dB', round(pk - r.peak, 2) + 0);
        end
        if ~isempty(d)
            txt = sprintf('Now vs reference %s: %s', r.time, strjoin(d, ', '));
        end
    end

    function redrawRefTrace()
        %REDRAWREFTRACE  Pin or Clear applied to the cut as drawn.
        %   Only the reference line, the legend and its caveat line change
        %   -- the design and its curves do not, so a recompute (seconds on
        %   a large array, and none at all with Auto off) would be pure
        %   cost. Blank (out-of-date) plots get nothing: the next compute
        %   draws the reference with the curves.
        refLine = plotRefTrace();
        hTotal = findall(axCut,'Tag','cutTotal');
        if ~isempty(hTotal)
            legend(axCut, [findall(axCut,'Tag','cutAF'); ...
                findall(axCut,'Tag','cutEF'); hTotal; refLine; findall(axCut,'Tag','cutCst')], ...
                'Location','eastoutside','Box','off','AutoUpdate','off');
        elseif isempty(findall(axCut,'Tag','cutAR'))
            return;
        end
        refLineRe = '^(Now vs r|R)eference \d{2}:\d{2}';
        cur = cellstr(string(axCut.Subtitle.String));
        cur = cur(~cellfun('isempty', cur));
        cur = cur(cellfun('isempty', regexp(cur, refLineRe, 'once')));
        note = refCaveat();
        if ~isempty(note), cur{end+1} = note; end
        setPlotTitle(axCut, axCut.Title.String, cur(:)');
    end

    % --------------------------------------------- plot titles and caveats
    function setPlotTitle(ax, head, caveats)
        %SETPLOTTITLE  A pattern plot's one-line title and its caveats.
        %   The title says what is drawn; the numbers are on the result
        %   cards. caveats (cell of char, {} for none) go one per line in
        %   a smaller muted subtitle, so a warning never lengthens the
        %   title into the plot above it. The subtitle carries the role
        %   Tag plotCaveat, which recolorPlots repaints on a theme change.
        P = pal();
        title(ax, head, 'FontSize', 10);
        if isempty(caveats), caveats = ''; end
        subtitle(ax, caveats, 'FontSize', 9, 'FontWeight', 'normal', ...
            'Color', P.muted, 'Tag', 'plotCaveat');
    end

    function addCaveat(ax, txt)
        % One more caveat line under a plot's title, added once.
        if ~isgraphics(ax), return; end
        cur = cellstr(string(ax.Subtitle.String));
        cur = cur(~cellfun('isempty', cur));
        if any(contains(cur, txt)), return; end
        setPlotTitle(ax, ax.Title.String, [cur(:)' {txt}]);
    end

    function blankPlots(msg3D, msgCut)
        % Both pattern plots emptied, for a state with no current result.
        % cla leaves the cut's legend and both subtitles behind: an empty
        % legend box and a caveat about a pattern that is no longer drawn.
        cla(ax3D); cla(axCut); legend(axCut, 'off');
        setPlotTitle(ax3D, msg3D, {}); setPlotTitle(axCut, msgCut, {});
    end

    function s = cutTitleText(isPhiCut, thFix)
        % Names the principal cut the way the Cut dropdown does, with the
        % angle held fixed along it.
        if isPhiCut
            if useAzEl(), s = sprintf('Azimuth cut at el = %g°', 90-thFix);
            else, s = sprintf('Phi cut at θ = %g°', thFix); end
        else
            if useAzEl(), s = sprintf('Elevation cut at az = %g°', cutPhiVal());
            else, s = sprintf('Theta cut at φ = %g°', cutPhiVal()); end
        end
    end

    function side = xlineLabelSide(ax, x, xOther)
        %XLINELABELSIDE  Which side of a vertical marker at x its label goes.
        %   Away from another marker at xOther (NaN: none) while that side
        %   still has a fifth of the axes width; otherwise, and with no
        %   other marker, toward the wider side, so the label is not cut
        %   off at the axes edge.
        xl = ax.XLim; minRoom = 0.2*diff(xl);
        roomL = x - xl(1); roomR = xl(2) - x;
        if ~isnan(xOther) && xOther > x && roomL >= minRoom
            side = 'left';
        elseif ~isnan(xOther) && xOther < x && roomR >= minRoom
            side = 'right';
        elseif roomR >= roomL
            side = 'right';
        else
            side = 'left';
        end
    end

    % ------------------------------------------ result cards and targets
    function clearCutMetrics()
        % No beam figures until a compute measures them again. Called
        % wherever the radiation result is discarded, so a card never
        % shows the HPBW of a pattern that is no longer on screen.
        S.cutMetrics = noCutMetrics();
    end

    function m = noCutMetrics()
        m = struct('hpbw',NaN,'fnbw',NaN,'sll',NaN, ...
            'measured',false,'plane','','offBeam',false,'why','');
    end

    function cm = blankCardModel()
        % What refreshInfo tells refreshCards about the first card; a NaN
        % value shows as a dash.
        cm = struct('gainCap','Directivity (dBi)','gainDb',NaN, ...
            'gainNote','','gainTone','','gainTip','');
    end

    function s = cutPlaneText(isPhiCut, thFix)
        % Names the principal cut the beam figures were measured on:
        % HPBW is per plane, so a card showing it has to say which one.
        if isPhiCut
            if useAzEl(), s = sprintf('el = %g° cone', 90-thFix);
            else, s = sprintf('θ = %g° cone', thFix); end
        else
            if useAzEl(), s = sprintf('az = %g° cut', cutPhiVal());
            else, s = sprintf('φ = %g° cut', cutPhiVal()); end
        end
    end

    function s = angleText(th, ph)
        % directionText's twin for the new readouts: same convention
        % switch, written with symbols and "°" to fit a narrow column.
        if useAzEl()
            [az,el] = azEl(th,ph);
            s = sprintf('az %.1f°, el %.1f°',az,el);
        else
            s = sprintf('θ %.1f°, φ %.1f°',th,ph);
        end
    end

    function rows = frequencyRows()
        %FREQUENCYROWS  The frequency / phase mode, one fact per row.
        %   This reports the GENERATED STEERING TERM, not a measured
        %   AF/total peak: manual/rotation phase offsets, a directive
        %   element factor, taper or competing lobes can move the actual
        %   maximum away from this prediction.
        if S.retunePhase
            modeTxt = 'retuned phases (steering follows f_op)';
        elseif size(S.el,1) <= 1
            modeTxt = 'phases fixed (single element: no squint)';
        else
            modeTxt = 'beam squint (phases fixed)';
        end
        rows = {'Frequency mode', modeTxt, ''
            'Design / operating', sprintf('%.4g / %.4g GHz',S.freqGHz,S.freqOpGHz), ''
            'Phase reference', sprintf('%.4g GHz',phaseReferenceGHz()), ''};
        if ~S.retunePhase && size(S.el,1) > 1
            sBeam = squintSinTheta(S.theta_s);
            if abs(sBeam) <= 1
                rows(end+1,:) = {'Steering-term squint', sprintf('%s → %s', ...
                    angleText(S.theta_s,S.phi_s),angleText(asind(sBeam),S.phi_s)), ''};
            else
                rows(end+1,:) = {'Steering-term squint', ...
                    'prediction outside visible space', 'warn'};
            end
        end
        if strcmp(S.efType,'Imported (CST far-field)') && ...
                abs(S.freqOpGHz-S.freqGHz) > 1e-12*max([1 abs(S.freqGHz) abs(S.freqOpGHz)])
            rows(end+1,:) = {'Imported pattern', ['not frequency-scaled: ' ...
                'load data for the operating frequency'], 'warn'};
        end
    end

    function showDetailsMessage(msg)
        % A state with nothing to report (no elements, pattern missing,
        % results discarded): one line saying why, and every card blank.
        setDetails({'Status', msg, 'warn'});
        refreshCards(blankCardModel(), '', noCutMetrics());
    end

    function setDetails(rows)
        %SETDETAILS  Show rows {quantity, value, tone} in the Details table.
        %   tone '' is plain; 'warn'/'bad' colour the value and prefix it
        %   with a symbol, so a warning does not rest on colour alone.
        %   Label pairs are created only when the table grows past what
        %   it has held before; spare pairs are hidden, not deleted.
        P = pal();
        nRow = size(rows,1);
        for kRow = size(detailLbls,1)+1:nRow
            qL = uilabel(detailsGrid,'Text','','FontSize',11,'FontColor',P.muted, ...
                'VerticalAlignment','top','Tag','detailName');
            qL.Layout.Row = kRow; qL.Layout.Column = 1;
            vL = uilabel(detailsGrid,'Text','','FontSize',11,'WordWrap','on', ...
                'VerticalAlignment','top','Tag','detailValue');
            vL.Layout.Row = kRow; vL.Layout.Column = 2;
            mutedLbls(end+1) = qL; %#ok<AGROW>
            detailLbls(kRow,:) = [qL vL];
        end
        heights = repmat({'fit'},1,size(detailLbls,1));
        for kRow = 1:size(detailLbls,1)
            inUse = kRow <= nRow;
            if inUse
                tone = rows{kRow,3};
                glyph = '';
                if any(strcmp(tone,{'warn','bad'})), glyph = '⚠ '; end
                setText(detailLbls(kRow,1), rows{kRow,1});
                setText(detailLbls(kRow,2), [glyph rows{kRow,2}]);
                setTone(detailLbls(kRow,2), tone, P);
            else
                setText(detailLbls(kRow,1), ''); setText(detailLbls(kRow,2), '');
                heights{kRow} = 0;
            end
            detailLbls(kRow,1).Visible = inUse; detailLbls(kRow,2).Visible = inUse;
        end
        detailsGrid.RowHeight = heights;
    end

    function setText(h, txt)
        % Skips the write when nothing changed: refreshInfo runs on every
        % selection click, and most of its text is the same each time.
        if ~strcmp(h.Text, txt), h.Text = txt; end
    end

    function setTone(h, tone, P)
        %SETTONE  Colour a result label, card frame or lamp by its role.
        %   tone '' leaves a label on the theme's own text colour and a
        %   frame on the divider colour; otherwise it is a pal() role
        %   ('good','warn','bad','muted'). The role is kept in UserData so
        %   repaintResultTones can re-apply it for a new theme.
        h.UserData = tone;
        if isa(h,'matlab.ui.control.Lamp')
            h.Color = P.(tone);
        elseif isa(h,'matlab.ui.container.Panel')
            if isempty(tone), h.BorderColor = P.rule; else, h.BorderColor = P.(tone); end
        elseif isempty(tone)
            if isprop(h,'FontColorMode'), h.FontColorMode = 'auto';
            else, h.FontColor = P.ink; end
        else
            h.FontColor = P.(tone);
        end
    end

    function repaintResultTones()
        % After a theme change: every card and Details value takes its
        % role's colour from the new palette.
        P = pal();
        objs = [findall(cardH.strip,'-property','UserData'); ...
            findall(detailsGrid,'Tag','detailValue')];
        for kObj = 1:numel(objs)
            if ischar(objs(kObj).UserData)
                setTone(objs(kObj), objs(kObj).UserData, P);
            end
        end
    end

    function H = buildResultCards(parent)
        %BUILDRESULTCARDS  The strip of result cards above the pattern plots.
        %   Each card is a caption, the value in large type, and a note
        %   line that carries either context (which cut, which efficiency)
        %   or the verdict against the design target. The value labels
        %   carry the Tags scripts and tests read: cardDirectivity,
        %   cardHPBW, cardSLL and cardGrating. Aperture efficiency is a
        %   Details row, not a card (the user's choice).
        P = pal();
        H.strip = uigridlayout(parent,[1 5],'Tag','resultCards');
        H.strip.Layout.Row = 1; H.strip.Layout.Column = 1;
        H.strip.ColumnWidth = {'1x','1x','1x','1x',72};
        H.strip.Padding = [0 0 0 0]; H.strip.ColumnSpacing = 6;
        spec = { ...  key    Tag                caption              tooltip
            'dir',   'cardDirectivity', 'Directivity (dBi)', ...
                'Directivity of the whole array at its beam peak (dBi).'
            'hpbw',  'cardHPBW',  'HPBW (°)', ['Half-power (−3 dB) beamwidth ' ...
                'on the principal cut shown below, in degrees.']
            'sll',   'cardSLL',   'SLL (dB)', ['Highest sidelobe on the principal ' ...
                'cut, relative to the beam peak (dB). More negative is better.']
            'gl',    'cardGrating', 'Grating lobes', ['Strongest secondary lobe ' ...
                'anywhere in visible space: none (green), a strong lobe grazing ' ...
                'the horizon (amber), or a possible grating lobe within 6 dB ' ...
                'of the peak (red).']};
        for kCard = 1:size(spec,1)
            key = spec{kCard,1}; tag = spec{kCard,2};
            pn = uipanel(H.strip,'BorderType','line','BorderColor',P.rule, ...
                'Tag',[tag 'Panel'],'UserData','');
            g = uigridlayout(pn,[3 1]);
            g.RowHeight = {14,'1x',14}; g.Padding = [8 2 6 3]; g.RowSpacing = 0;
            cap = uilabel(g,'Text',spec{kCard,3},'FontSize',10, ...
                'FontColor',P.muted,'Tag',[tag 'Caption'],'Tooltip',spec{kCard,4});
            mutedLbls(end+1) = cap; %#ok<AGROW>
            if strcmp(key,'gl')
                % Lamp beside the text; the lamp sits in its own 14 px
                % cell so it stays a small round light at any card height.
                gv = uigridlayout(g,[1 2]);
                gv.ColumnWidth = {14,'1x'}; gv.Padding = [0 0 0 0];
                gv.ColumnSpacing = 6;
                gl = uigridlayout(gv,[3 1]);
                gl.RowHeight = {'1x',14,'1x'}; gl.Padding = [0 0 0 0];
                H.lamp = uilamp(gl,'Color',P.muted,'Tag','cardGratingLamp', ...
                    'Tooltip',spec{kCard,4},'UserData','muted');
                H.lamp.Layout.Row = 2;
                % A word beside the lamp rather than a number, so a size
                % smaller: "Possible" still fits a laptop-width window.
                parentV = gv; valSize = 18;
            else
                parentV = g; valSize = 20;
            end
            val = uilabel(parentV,'Text','—','FontSize',valSize, ...
                'FontWeight','bold','Tag',tag,'Tooltip',spec{kCard,4},'UserData','');
            note = uilabel(g,'Text','','FontSize',10,'Tag',[tag 'Note'], ...
                'FontColor',P.muted,'UserData','muted');
            note.Layout.Row = 3;
            H.pan.(key) = pn; H.cap.(key) = cap;
            H.val.(key) = val; H.note.(key) = note; H.tip.(key) = spec{kCard,4};
        end
        gb = uigridlayout(H.strip,[3 1]);
        gb.RowHeight = {'1x',30,'1x'}; gb.Padding = [0 0 0 0];
        H.btn = uibutton(gb,'Text','Targets…','Tag','btnTargets', ...
            'Tooltip',['Set design targets: minimum directivity or gain (dBi), ' ...
            'maximum HPBW (°), maximum sidelobe level (dB), no grating lobes. ' ...
            'Cards turn green when a target is met and red when it is missed.'], ...
            'ButtonPushedFcn',@(s,e)editTargets());
        H.btn.Layout.Row = 2;
    end

    function refreshCards(cm, glState, C)
        %REFRESHCARDS  Fill the result cards and judge them against S.targets.
        %   cm is the card model refreshInfo built (the first card), C the
        %   cut metrics HPBW and SLL come from, and
        %   glState the grating state ('' / none / grazing / possible).
        P = pal();
        T = S.targets;
        verdicts = {};
        % Directivity / gain
        cardH.cap.dir.Text = cm.gainCap;
        tip = cm.gainTip; if isempty(tip), tip = cardH.tip.dir; end
        cardH.val.dir.Tooltip = tip;
        showCard('dir', cm.gainDb, '%.2f', cm.gainDb, T.minGain, false, ...
            '%.2f dBi', cm.gainNote, cm.gainTone);
        % HPBW and SLL, measured on the principal cut. A beam that never
        % falls 3 dB on the cut is wider than any HPBW target; a cut with
        % no sidelobe at all meets any SLL target.
        % A cut pinned away from the steering plane misses the beam, so
        % both figures describe something else: said in amber.
        if C.measured
            if C.offBeam
                hpNote = 'off-beam cut'; cutTone = 'warn';
            else
                hpNote = C.plane; cutTone = '';
            end
            sllNote = hpNote;
            if isfinite(C.hpbw)
                hpShow = C.hpbw; hpJudge = C.hpbw;
            else
                hpShow = NaN; hpJudge = Inf; hpNote = 'wider than cut';
            end
            if isfinite(C.sll)
                sllShow = C.sll; sllJudge = C.sll;
            else
                sllShow = -Inf; sllJudge = -Inf; sllNote = 'no sidelobes';
            end
        else
            % C.why names an expected blank (axial-ratio view, no field
            % on the cut); an empty one means not computed yet.
            hpShow = NaN; hpJudge = NaN; hpNote = C.why; cutTone = '';
            sllShow = NaN; sllJudge = NaN; sllNote = C.why;
        end
        cutTip = sprintf(' Measured on the %s.', C.plane);
        if ~C.measured, cutTip = ''; end
        cardH.val.hpbw.Tooltip = [cardH.tip.hpbw cutTip];
        cardH.val.sll.Tooltip = [cardH.tip.sll cutTip];
        % HPBW to two decimals: the cut title used to carry that
        % precision, and this card is where it is read now.
        showCard('hpbw', hpShow, '%.2f', hpJudge, T.maxHPBW, true, ...
            '%.1f°', hpNote, cutTone);
        showCard('sll', sllShow, '%.1f', sllJudge, T.maxSLL, true, ...
            '%.1f dB', sllNote, cutTone);
        % Grating lobes: the lamp carries the state, the value names it.
        switch glState
            case 'none'
                glTxt = 'None'; lampTone = 'good';
                if isfinite(S.glRelDb)
                    glNote = sprintf('next lobe %.1f dB', S.glRelDb);
                else
                    glNote = 'one main lobe';
                end
            case 'grazing'
                glTxt = 'Grazing'; lampTone = 'warn';
                glNote = sprintf('horizon, %.1f dB', S.glRelDb);
            case 'possible'
                glTxt = 'Possible'; lampTone = 'bad';
                if useAzEl()
                    [glAz,glEl] = azEl(S.glTheta,S.glPhi);
                    glNote = sprintf('az %.0f°, el %.0f°', glAz, glEl);
                else
                    glNote = sprintf('θ %.0f°, φ %.0f°', S.glTheta, S.glPhi);
                end
            otherwise
                glTxt = '—'; lampTone = 'muted'; glNote = '';
        end
        setTone(cardH.lamp, lampTone, P);
        setText(cardH.val.gl, glTxt);
        % "Not allowed" follows the app's own classification, the one Max
        % scan uses too: a lobe grazing the horizon peaks outside visible
        % space, so it meets the target -- the amber lamp stays as the
        % caution.
        glVerdict = '';
        if T.noGrating
            switch glState
                case 'possible'
                    glVerdict = 'bad'; glNote = '✕ not allowed';
                case 'grazing'
                    glVerdict = 'good'; glNote = '✓ met (grazing)';
                case 'none'
                    glVerdict = 'good'; glNote = '✓ target met';
                otherwise
                    glNote = 'target: none';
            end
        end
        glNoteTone = glVerdict;
        if isempty(glNoteTone)
            if strcmp(glState,'grazing'), glNoteTone = 'warn'; else, glNoteTone = 'muted'; end
        end
        setText(cardH.note.gl, glNote); setTone(cardH.note.gl, glNoteTone, P);
        setTone(cardH.val.gl, glVerdict, P); setTone(cardH.pan.gl, glVerdict, P);
        if ~isempty(glVerdict), verdicts{end+1} = glVerdict; end
        cardVerdicts = verdicts;

        function showCard(key, shown, fmt, judged, target, isMax, goalFmt, ctxNote, ctxTone)
            % One numeric card: its value, and either its context note or
            % its verdict against the target (NaN target = none set).
            if isnan(shown)
                vTxt = '—';
            elseif isinf(shown)
                vTxt = 'none';
            else
                % A level that rounds to zero reads 0.0, not -0.0.
                vTxt = regexprep(sprintf(fmt, shown), '^-(0\.?0*)$', '$1');
            end
            setText(cardH.val.(key), vTxt);
            verdict = ''; noteTxt = ctxNote;
            if isempty(ctxTone), noteTone = 'muted'; else, noteTone = ctxTone; end
            noteTip = '';
            if ~isnan(target)
                if isMax, op = '≤'; else, op = '≥'; end
                goal = sprintf(['%s ' goalFmt], op, target);
                if isnan(judged)
                    noteTxt = ['target ' goal]; noteTone = 'muted';
                    noteTip = sprintf('Design target %s: no result to check yet.', goal);
                elseif (isMax && judged <= target) || (~isMax && judged >= target)
                    verdict = 'good'; noteTxt = ['✓ ' goal];
                    noteTip = sprintf('Meets the design target (%s).', goal);
                else
                    verdict = 'bad'; noteTxt = ['✕ ' goal];
                    noteTip = sprintf('Misses the design target (%s).', goal);
                end
                if ~isempty(verdict), noteTone = verdict; verdicts{end+1} = verdict; end
            end
            setText(cardH.note.(key), noteTxt);
            cardH.note.(key).Tooltip = noteTip;
            setTone(cardH.note.(key), noteTone, P);
            setTone(cardH.val.(key), verdict, P);
            setTone(cardH.pan.(key), verdict, P);
        end
    end

    function editTargets()
        %EDITTARGETS  The Design targets window (one at a time).
        %   Modal but non-blocking: nothing here waits for it, the
        %   buttons do the work. An empty field means "no target".
        if ~isempty(targetsWin) && isgraphics(targetsWin)
            fillTargetsWindow(); figure(targetsWin); return;
        end
        P = pal();
        fp = fig.Position; w = 420; h = 262;
        targetsWin = uifigure('Name','Design targets','Tag','targetsWindow', ...
            'Position',[fp(1)+(fp(3)-w)/2, fp(2)+(fp(4)-h)/2, w, h], ...
            'Resize','off','WindowStyle','modal');
        matchTheme(targetsWin);
        % Goes when the app does, instead of lingering over a closed
        % design. The listener lives in the window, so it goes with it.
        winNow = targetsWin;
        targetsWin.UserData = listener(fig,'ObjectBeingDestroyed', ...
            @(~,~) delete(winNow(isgraphics(winNow))));
        g = uigridlayout(targetsWin,[6 2]);
        g.ColumnWidth = {'1x',120}; g.RowHeight = {24,24,24,24,'1x',28};
        g.Padding = [14 12 14 12]; g.RowSpacing = 8;
        % One row per entry of TARGET_NUMS, in that order.
        fields = { ...
            'Minimum directivity / gain (dBi)', 'tgtMinGain', ...
                ['Lowest acceptable value of the first card (dBi): directivity, ' ...
                'or realized gain for an imported gain pattern, or estimated ' ...
                'array gain in unit-cell mode. Empty = no target.']
            'Maximum HPBW (°)', 'tgtMaxHPBW', ...
                ['Widest acceptable half-power beamwidth on the principal cut ' ...
                '(degrees). Empty = no target.']
            'Maximum sidelobe level (dB)', 'tgtMaxSLL', ...
                ['Highest acceptable sidelobe relative to the peak (dB, so ' ...
                'negative): -20 means every sidelobe at least 20 dB down. ' ...
                'Empty = no target.']};
        tgtCtl = struct();
        for kField = 1:size(fields,1)
            lb = uilabel(g,'Text',fields{kField,1},'Tooltip',fields{kField,3});
            lb.Layout.Row = kField; lb.Layout.Column = 1;
            ed = uieditfield(g,'numeric','AllowEmpty','on','Value',[], ...
                'Limits',TARGET_LIMITS{kField},'Placeholder','none', ...
                'Tag',fields{kField,2},'Tooltip',fields{kField,3});
            ed.Layout.Row = kField; ed.Layout.Column = 2;
            tgtCtl.(TARGET_NUMS{kField}) = ed;
        end
        cb = uicheckbox(g,'Text','Grating lobes not allowed','Tag','tgtNoGrating', ...
            'Tooltip',['Fail the Grating lobes card when a possible grating ' ...
            'lobe (within 6 dB of the peak) is found anywhere in visible space.']);
        cb.Layout.Row = 4; cb.Layout.Column = [1 2];
        tgtCtl.noGrating = cb;
        hint = uilabel(g,'Text',['Leave a field empty for no target. Cards ' ...
            'turn green when their target is met and red when it is missed. ' ...
            'Targets are saved with the design.'], ...
            'WordWrap','on','FontSize',11,'FontColor',P.muted, ...
            'VerticalAlignment','top');
        hint.Layout.Row = 5; hint.Layout.Column = [1 2];
        mutedLbls(end+1) = hint;
        gb = uigridlayout(g,[1 4]);
        gb.Layout.Row = 6; gb.Layout.Column = [1 2];
        gb.ColumnWidth = {84,'1x',84,84}; gb.Padding = [0 0 0 0];
        uibutton(gb,'Text','Clear all','Tag','tgtClear', ...
            'Tooltip','Empty every field (applied only with Apply)', ...
            'ButtonPushedFcn',@(s,e)clearTargetsWindow());
        bCancel = uibutton(gb,'Text','Cancel','Tag','tgtCancel', ...
            'Tooltip','Close without changing the targets', ...
            'ButtonPushedFcn',@(s,e)delete(targetsWin));
        bCancel.Layout.Column = 3;
        bApply = uibutton(gb,'Text','Apply','Tag','tgtApply', ...
            'Tooltip','Use these targets and close', ...
            'ButtonPushedFcn',@(s,e)applyTargetsWindow());
        bApply.Layout.Column = 4;
        fillTargetsWindow();
    end

    function fillTargetsWindow()
        % Shows S.targets in the open window (NaN = empty field).
        if isempty(targetsWin) || ~isgraphics(targetsWin), return; end
        for kTgt = 1:numel(TARGET_NUMS)
            v = S.targets.(TARGET_NUMS{kTgt});
            if isnan(v), v = []; end
            tgtCtl.(TARGET_NUMS{kTgt}).Value = v;
        end
        tgtCtl.noGrating.Value = S.targets.noGrating;
    end

    function clearTargetsWindow()
        for kTgt = 1:numel(TARGET_NUMS)
            tgtCtl.(TARGET_NUMS{kTgt}).Value = [];
        end
        tgtCtl.noGrating.Value = false;
    end

    function applyTargetsWindow()
        % Read the window into S.targets, close it, and re-judge the cards.
        T = S.targets;
        for kTgt = 1:numel(TARGET_NUMS)
            v = tgtCtl.(TARGET_NUMS{kTgt}).Value;
            if isempty(v), v = NaN; end
            T.(TARGET_NUMS{kTgt}) = double(v);
        end
        T.noGrating = logical(tgtCtl.noGrating.Value);
        delete(targetsWin);
        assignTargets(T);
    end

    function assignTargets(T)
        % Targets are part of the design but change no result, so this
        % re-judges the cards (refreshInfo) instead of recomputing.
        S.targets = T;
        refreshInfo();
        nSet = ~isnan(T.minGain) + ~isnan(T.maxHPBW) + ~isnan(T.maxSLL) + T.noGrating;
        nMet = nnz(strcmp(cardVerdicts,'good'));
        nJudged = numel(cardVerdicts);
        if nSet == 0
            setStatus('Design targets cleared.');
        elseif nJudged < nSet
            setStatus(sprintf(['Design targets set: %d of %d met, %d not ' ...
                'checked yet (no result).'], nMet, nSet, nSet-nJudged), 'warn');
        elseif nMet == nSet
            setStatus(sprintf('Design targets: all %d met.', nSet), 'good');
        else
            setStatus(sprintf('Design targets: %d of %d met.', nMet, nSet), 'warn');
        end
    end

    function T = sanitizeTargets(raw)
        % A saved targets struct, checked field by field and clamped to
        % the ranges the Targets window accepts. [] when it is not one,
        % so loadConfig can report it rather than install it.
        T = [];
        if ~isstruct(raw) || ~isscalar(raw), return; end
        out = struct('minGain',NaN,'maxHPBW',NaN,'maxSLL',NaN,'noGrating',false);
        for kLim = 1:numel(TARGET_NUMS)
            fn = TARGET_NUMS{kLim}; lim = TARGET_LIMITS{kLim};
            if ~isfield(raw,fn), continue; end
            v = raw.(fn);
            if ~(isnumeric(v) && isscalar(v) && isreal(v) && ~isinf(v)), return; end
            if ~isnan(v), v = min(max(double(v),lim(1)),lim(2)); end
            out.(fn) = double(v);
        end
        if isfield(raw,'noGrating')
            if ~validBoolean(raw.noGrating), return; end
            out.noGrating = logical(raw.noGrating);
        end
        T = out;
    end
    function v = range0(x)
        if isempty(x), v = 0; else, v = max(x)-min(x); end
    end

    function u = unitCellFigures()
        %UNITCELLFIGURES  Embedded-element and array-level figures, unit-cell mode.
        %   Returns [] unless a unit-cell import is active with elements
        %   placed. Single source of truth for these numbers: the pattern
        %   title and the Metrics panel both read it, so the two cannot
        %   disagree the way they did when each computed its own.
        u = [];
        if ~(S.impUnitCell && strcmp(S.efType,'Imported (CST far-field)') ...
                && ~isempty(S.impFF) && isfield(S.impFF,'peak') && ~isempty(S.el))
            return;
        end
        u.n      = size(S.el,1);
        u.rg     = S.impFF.peak(1);          % embedded-element realized gain
        u.effPct = S.impTotEffPct;
        % Directivity from the entered efficiency, NOT from the pattern
        % integral: for a periodic unit cell that integral does not equal
        % the radiated power, so it cannot supply a directivity here.
        u.d      = u.rg - 10*log10(u.effPct/100);
        % Array build-up uses the EFFECTIVE element count
        % (sum|w|)^2 / sum(w^2), not the raw count. An amplitude taper
        % spends aperture to buy sidelobe suppression, and a flat
        % 10*log10(N) would credit the array with gain the taper has
        % already given away (~1 dB on a -30 dB Taylor). Reduces exactly
        % to N for a uniform array (N^2/N), the common case here.
        w  = calculationWeights();
        % Weights normalised before squaring. nEff is scale-invariant by
        % construction, but computing it on raw values overflows: two
        % amplitudes of 1e200 give sum(|w|)^2 = Inf and sum(w^2) = Inf, so
        % Inf/Inf = NaN and the taper efficiency -- and everything derived
        % from it -- became NaN for a design the amplitude validator had
        % accepted. Dividing by the largest weight first keeps every
        % intermediate below 1 without changing the answer.
        wmx = max(abs(w));
        if wmx > 0 && isfinite(wmx)
            wn = abs(w)/wmx;
            swn = sum(wn.^2);
            if swn > 0, u.nEff = sum(wn)^2 / swn; else, u.nEff = 0; end
        else
            u.nEff = 0;
        end
        % Embedded-element realized gain toward the STEER direction, not
        % at its peak. This is the whole scan-loss term, and leaving it
        % out was a real bug: the reported array gain used the file's
        % peak, which is a constant, so steering to 60 deg changed the
        % plotted pattern while the gain and directivity readouts sat
        % perfectly still. An array that loses nothing off broadside is
        % not a phased array.
        %
        % For an element embedded in a periodic lattice the element
        % pattern already carries the aperture foreshortening (ideally
        % G_emb(th) = 4*pi*A_cell*cos(th)/lambda^2), so evaluating it at
        % the steer angle captures the scan loss in full. The array
        % factor contributes nothing here: its peak value at the steer
        % direction is sum|w| regardless of where that direction is.
        %
        % Rotations are power-averaged, weighted by |w|^2. All elements
        % normally share a rotation, in which case this is exact; under
        % sequential rotation they genuinely see different element gains
        % toward the steer angle, and their mean is the right scaling
        % for the summed array.
        % Canonicalise the steer direction FIRST. The steering spinner
        % allows theta down to -90, but an imported pattern is only
        % defined over theta 0..180 -- so a negative steer angle landed
        % outside the interpolant, evaluated to zero, and hit the floor
        % below: -30 deg reported -96.93 dBi instead of 22.47. Negative
        % theta at azimuth phi is the same physical direction as
        % +|theta| at phi+180 (both give u = -sin|th|cos(phi),
        % v = -sin|th|sin(phi)), so reflecting it is exact, not an
        % approximation.
        thS = S.theta_s; phS = S.phi_s;
        % Under SQUINT an ARRAY's ideal progressive-phase stationary
        % direction moves away from the nominal command. A single element is
        % different: it has no progressive phase gradient at all, so changing
        % the frequency cannot create array-factor squint; only the element
        % pattern itself may change with frequency.
        if u.n <= 1
            sBeam = sind(thS);
        else
            sBeam = squintSinTheta(thS);
        end
        u.beamVisible = abs(sBeam) <= 1;
        if ~u.beamVisible
            % Do NOT clamp an out-of-visible stationary point to grazing and
            % call the grazing value the beam result. There is no real
            % commanded steering-term stationary direction in visible space in this case; the
            % finite array can still have a boundary/interior maximum, but
            % that must be found from the plotted pattern, not fabricated by
            % clipping the analytic squint equation to 90 degrees.
            u.evalThetaDeg = NaN; u.evalPhiDeg = NaN;
            u.scanLossDb = NaN;  u.rgSteer = NaN;
            u.buildUpDb = NaN;   u.gArr = NaN; u.dArr = NaN;
            return;
        end
        thS = asind(sBeam);
        if thS < 0, thS = -thS; phS = phS + 180; end
        phS = mod(phS,360);
        u.evalThetaDeg = thS;
        u.evalPhiDeg = phS;
        rotAll = S.el(:,5); wAmp = abs(calculationWeights());
        [uRot,~,rIdx] = unique(rotAll);
        pwSum = 0; wSum = 0;
        for q = 1:numel(uRot)
            [etS,epS] = elementFactor(thS, phS, uRot(q));
            wq = sum(wAmp(rIdx == q).^2);
            pwSum = pwSum + wq*fieldPower(etS,epS);
            wSum  = wSum + wq;
        end
        if wSum > 0, relPw = pwSum/wSum; else, relPw = 1; end
        % elementFactor hands back a UNIT-PEAK pattern for an import, so
        % relPw is already peak-relative. A true null is deliberately NOT
        % floored: it becomes infinite scan loss and therefore -Inf gain,
        % which is more truthful than inventing a finite floor value.
        u.scanLossDb = -10*log10(relPw); % true null -> Inf loss
        u.rgSteer    = u.rg - u.scanLossDb;

        % Array build-up from the COHERENT VECTOR SUM, not from an
        % amplitude-only effective element count.
        %
        % nEff = (sum|w|)^2 / sum(w^2) is a TAPER efficiency. It assumes
        % every element arrives in phase, so it cannot see a manual phase
        % offset, a negative amplitude, or the different vector fields
        % that mixed rotations produce. Two elements driven 0 and 180 deg
        % cancel at broadside, and the old form still reported the full
        % +3.01 dB build-up while the plotted pattern correctly moved --
        % a reported gain that did not depend on half the excitation.
        %
        % This instead evaluates the standard embedded-element relation
        %
        %   G_arr(u) = G_emb(u) * |sum_n w_n E_n(u) e^{jk r_n.u}|^2
        %                       / sum_n |w_n|^2
        %
        % with w_n the complex feed weight, E_n the element's own vector
        % field at that rotation, and the sum taken per rotation group so
        % each group's field is applied once. E_n is unit-peak, so the
        % element rolloff toward the steer angle is already inside the
        % sum -- which is why this replaces rgSteer + 10log10(nEff)
        % entirely rather than being added to it. It reduces exactly to
        % that old expression for the uniform, co-phased, single-rotation
        % case it silently assumed.
        us_ = sind(thS)*cosd(phS); vs_ = sind(thS)*sind(phS);
        % Propagation is evaluated at the operating frequency while the
        % element coordinates remain in design-wavelength units. This
        % local ratio is required here just as it is in computePatternCore;
        % without it the unit-cell metrics branch referenced an undefined
        % variable `fr` and failed as soon as it evaluated array build-up.
        fr = freqRatio();
        phF = deg2rad(effectivePhaseDeg());
        ampAll = calculationWeights();
        EthS = 0; EphS = 0;
        for q = 1:numel(uRot)
            [etq,epq] = elementFactor(thS, phS, uRot(q));
            sel = (rIdx == q);
            wq = ampAll(sel).*exp(1j*phF(sel)) .* ...
                 exp(1j*2*pi*fr*(S.el(sel,1)*us_ + S.el(sel,2)*vs_));
            EthS = EthS + sum(wq)*etq;
            EphS = EphS + sum(wq)*epq;
        end
        pIn = sum(ampAll.^2);
        if pIn > 0
            buildUp = fieldPower(EthS,EphS)/pIn;
        else
            buildUp = 0;
        end
        % A fully cancelling excitation gives exactly zero build-up. Keep
        % that as -Inf dB rather than hiding the cancellation behind a floor.
        u.buildUpDb = -Inf;
        if buildUp>0, u.buildUpDb=10*log10(buildUp); end
        u.gArr = u.rg + u.buildUpDb;
        u.dArr = u.gArr - 10*log10(u.effPct/100);
    end

    function refreshLayout()
        P = pal();
        if isempty(S.el), amp = []; else, amp = abs(calculationWeights()); end
        drawElementLayout(axLay,S,layoutView,P,amp,@latticePosAt);
        updateLayoutPanelTitle(pLay,S,layoutView);
    end

    % --------------------------------------------------- element factor
    function E = polCombine(E_th, E_ph)
        % Combine the theta/phi field components into one magnitude, in
        % whichever polarization basis the user asked to view.
        % IEEE convention (e^{jwt}), verified independently: a pure RHCP
        % field has (E_theta,E_phi) = E0*(1,-j), so extracting its RHCP
        % AMPLITUDE via the Hermitian inner product with the RHCP unit
        % vector e_R=(1,-j)/sqrt2 gives c_R = (E_th + j*E_ph)/sqrt2 (the
        % PLUS sign) -- not minus. Sign flip vs. an earlier version of
        % this function, which mapped a pure-RHCP field to zero.
        switch ddPol.Value
            case 'RHCP component'
                E = abs(E_th + 1j*E_ph)/sqrt(2);
            case 'LHCP component'
                E = abs(E_th - 1j*E_ph)/sqrt(2);
            otherwise
                E = hypot(abs(E_th),abs(E_ph));
        end
    end

    function g = safeOff(effLin, Prad)
        %SAFEOFF  Absolute-dBi offset from the radiated-power integral.
        %   G(th,ph) = 4*pi*eff*P(th,ph)/Prad, and the plotted curves are
        %   already P on the same scale Prad was integrated over, so the
        %   whole conversion is this one additive constant. Zero radiated
        %   power (an all-zero excitation) has no absolute level to
        %   report, so the display stays peak-relative rather than
        %   returning -Inf.
        if isfinite(Prad) && Prad > 0
            g = 10*log10(4*pi*effLin/Prad);
        else
            % No radiated power: directivity is 0/0, undefined -- not
            % zero. Returning 0 here left the floored -80 dB curve on
            % screen still labelled as absolute dBi, so a field that does
            % not exist was being reported as a gain. NaN propagates into
            % the title and the axis label instead, which is what
            % "undefined" should look like.
            g = NaN;
        end
    end

    function powerValue = fieldPower(eth,eph)
        magnitude=hypot(abs(eth),abs(eph));
        powerValue=magnitude.^2;
        if any(~isfinite(powerValue(:)))
            error('PAD:PowerRange','Field power exceeds numeric range. Rescale both custom field formulas by the same factor.');
        end
        if any(magnitude(:)>0) && ~any(powerValue(:)>0)
            error('PAD:PowerRange','Nonzero field power underflowed to zero. Increase both custom field formulas by the same factor.');
        end
    end

    function [D, Prad] = fieldDirectivity(P_base, theta, phi, TH)
        % D = 4*pi*P_max / integral(P dOmega). Calculation callers use
        % the full sphere regardless of the display hemisphere setting.
        % Unlike a raw sum, trapz needs BOTH ends of
        % a periodic domain to correctly close the integral -- dropping
        % the phi=360 "duplicate" here would leave a sliver of the loop
        % unintegrated and inflate D. (Verified: keeping both endpoints
        % gives exactly the same Prad as an independent periodic
        % uniform-weight integration, to 5 decimal places -- dropping the
        % endpoint does not.)
        %
        % P_base (not the possibly polarization-projected display value)
        % is what gets integrated: when a polarization readout (RHCP/
        % LHCP) is selected, the displayed field only carries ONE
        % polarization component, and integrating that alone would
        % silently exclude cross-pol power from Prad -- inflating
        % directivity by exactly the amount of leakage, which defeats
        % the purpose of a polarization readout on a CP array. P_base
        % always carries the TRUE total (both-pol) power for whichever
        % field it belongs to, so cross-pol loss shows up in the
        % reported number instead of disappearing from it. Verified
        % numerically: a deliberately imperfect 4-element sequential-
        % rotation array showed a real 0.35 dB gap between "RHCP
        % directivity vs. its own RHCP power" and "RHCP partial
        % directivity vs. true total power" -- confirms this is a real,
        % non-negligible effect, not just a formality.
        %
        % Pmax is now taken DIRECTLY from P_base itself (no longer
        % accepted as a separate peak_db argument) -- this is the fix
        % for a real bug: computePattern used to pass max(totDb(:))/
        % max(efDb(:)) as the peak, which reflect whatever the
        % Polarization dropdown (RHCP/LHCP) had selected, while P_base
        % always carries TOTAL power -- pairing an RHCP-only peak with
        % a total-power Prad is internally inconsistent (numerator and
        % denominator from different fields). Deriving Pmax from P_base
        % itself makes that mismatch structurally impossible.
        P_int = P_base.*sind(TH);
        inner = trapz(deg2rad(theta), P_int, 2);
        Prad  = trapz(deg2rad(phi), inner);
        if Prad > 0
            D = 4*pi*max(P_base(:))/Prad;
        else
            D = NaN;
        end
    end

    function Prad = exactIsotropicPrad(feedWeights, fr)
        % Each isotropic point-source pair integrates over the full sphere
        % to 4*pi*sin(k*r)/(k*r). Summing pairs avoids missing narrow beams
        % on the display grid. Blocks bound memory for arbitrary layouts.
        % fr: propagation ratio, default freqRatio() (the band sweep
        % passes its own, as it does to the three integrals below).
        if nargin < 2, fr = freqRatio(); end
        feedWeights = feedWeights(:);
        n = numel(feedWeights);
        powerSum = 0;
        kOp = S.k*fr;
        for first = 1:128:n
            block = first:min(first+127,n);
            dx = S.el(block,1) - S.el(:,1).';
            dy = S.el(block,2) - S.el(:,2).';
            kd = kOp*hypot(dx,dy);
            kernel = ones(size(kd));
            nonzero = kd ~= 0;
            kernel(nonzero) = sin(kd(nonzero))./kd(nonzero);
            powerSum = powerSum + real(sum(conj(feedWeights(block)).* ...
                (kernel*feedWeights)));
        end
        Prad = 4*pi*max(powerSum,0);
    end

    function [Prad, kernelBlocks] = exactCosQPrad(feedWeights, kernelBlocks, fr)
        % A cos^q element radiates cos(theta)^q power into the front
        % hemisphere. The azimuthal integral of each planar source pair is
        % a Bessel function; summing those pairs resolves narrow beams
        % without making the display grid impractically dense.
        % Cached kernelBlocks belong to ONE fr: pass {} for a new one.
        if nargin < 2, kernelBlocks = {}; end
        if nargin < 3, fr = freqRatio(); end
        feedWeights = feedWeights(:);
        n = numel(feedWeights);
        nu = (S.efQ+1)/2;
        powerSum = 0;
        kOp = S.k*fr;
        keepBlocks = nargout > 1 && n*n*8 <= 160e6;
        if keepBlocks && isempty(kernelBlocks)
            kernelBlocks = cell(ceil(n/128),1);
        end
        for first = 1:128:n
            block = first:min(first+127,n);
            slot = ceil(first/128);
            if slot <= numel(kernelBlocks) && ~isempty(kernelBlocks{slot})
                kernel = kernelBlocks{slot};
            else
                dx = S.el(block,1) - S.el(:,1).';
                dy = S.el(block,2) - S.el(:,2).';
                kd = kOp*hypot(dx,dy);
                kernel = ones(size(kd))/(S.efQ+1);
                nonzero = kd > 1e-8;
                if S.efQ == 0
                    kernel(nonzero) = sin(kd(nonzero))./kd(nonzero);
                else
                    b = kd(nonzero);
                    kernel(nonzero) = 2^(nu-1)*gamma(nu)* ...
                        besselj(nu,b)./b.^nu;
                end
                if keepBlocks, kernelBlocks{slot} = kernel; end
            end
            powerSum = powerSum + real(sum(conj(feedWeights(block)).* ...
                (kernel*feedWeights)));
        end
        Prad = 2*pi*max(powerSum,0);
    end

    function [Prad, kernelBlocks] = exactZDipolePrad(feedWeights,kernelBlocks,fr)
        % Exact full-sphere integral for planar arrays of z-directed short
        % dipoles: each pair has sin^2(theta) element power. With b = k*d,
        % the angular kernel is 4*pi*(sin(b)/b + cos(b)/b^2 - sin(b)/b^3).
        % The Taylor branch gives the correct 8*pi/3 self-pair limit and
        % avoids subtracting large terms for closely spaced elements.
        % Cached blocks belong to one propagation ratio fr.
        if nargin < 2, kernelBlocks = {}; end
        if nargin < 3, fr = freqRatio(); end
        feedWeights = feedWeights(:);
        n = numel(feedWeights);
        powerSum = 0;
        kOp = S.k*fr;
        keepBlocks = nargout > 1 && n*n*8 <= 160e6;
        if keepBlocks && isempty(kernelBlocks)
            kernelBlocks = cell(ceil(n/128),1);
        end
        for first = 1:128:n
            block = first:min(first+127,n);
            slot = ceil(first/128);
            if slot <= numel(kernelBlocks) && ~isempty(kernelBlocks{slot})
                kernel = kernelBlocks{slot};
            else
                dx = S.el(block,1) - S.el(:,1).';
                dy = S.el(block,2) - S.el(:,2).';
                b = kOp*hypot(dx,dy);
                kernel = zeros(size(b));
                small = abs(b) < 0.01;
                z = b(small);
                kernel(small) = 4*pi*(2/3 - 2*z.^2/15 + z.^4/140);
                z = b(~small);
                kernel(~small) = 4*pi*(sin(z)./z + cos(z)./z.^2 - sin(z)./z.^3);
                if keepBlocks, kernelBlocks{slot} = kernel; end
            end
            powerSum = powerSum + real(sum(conj(feedWeights(block)).* ...
                (kernel*feedWeights)));
        end
        Prad = max(powerSum,0);
    end

    function [Prad, kernelBlocks] = exactVectorPrad(feedWeights,rotations,kernelBlocks,fr)
        % Dipole and patch fields contain only azimuthal orders 0 and 2.
        % Integrating each ordered element pair analytically avoids the
        % display mesh's narrow-beam power error, including when elements
        % have different physical rotations.
        % Cached kernelBlocks belong to ONE fr: pass {} for a new one.
        if nargin < 3, kernelBlocks = {}; end
        if nargin < 4, fr = freqRatio(); end
        feedWeights = feedWeights(:);
        rotations = rotations(:);
        n = numel(feedWeights);
        isPatch = strcmp(S.efType,'Patch (cos^q x lin pol)');
        if isPatch, q = S.efQ; else, q = 0; end
        nu = (q+1)/2;
        scale = 2^(nu-1)*gamma(nu);
        scale2 = 2^nu*gamma(nu+1);
        az = deg2rad(90-rotations);
        kOp = S.k*fr;
        powerSum = 0;
        keepBlocks = nargout > 1 && n*n*8 <= 160e6;
        if keepBlocks && isempty(kernelBlocks)
            kernelBlocks = cell(ceil(n/128),1);
        end
        for first = 1:128:n
            block = first:min(first+127,n);
            slot = ceil(first/128);
            if slot <= numel(kernelBlocks) && ~isempty(kernelBlocks{slot})
                kernel = kernelBlocks{slot};
            else
                dx = S.el(block,1) - S.el(:,1).';
                dy = S.el(block,2) - S.el(:,2).';
                b = kOp*hypot(dx,dy);
                c0 = ones(size(b))/(q+1);
                c2 = ones(size(b))/(q+3);
                s2 = zeros(size(b));
                nonzero = b > 1e-8;
                z = b(nonzero);
                if q == 0
                    c0(nonzero) = sin(z)./z;
                else
                    c0(nonzero) = scale*besselj(nu,z)./z.^nu;
                end
                c2(nonzero) = scale2*besselj(nu+1,z)./z.^(nu+1);
                s2(nonzero) = scale*besselj(nu+2,z)./z.^nu;
                cosDelta = cos(az(block)-az.');
                cosSum = cos(2*atan2(dy,dx)-az(block)-az.');
                if isPatch
                    kernel = pi*(cosDelta.*(c0+c2)-cosSum.*s2);
                else
                    kernel = 2*pi*(cosDelta.*(c0+c2)+cosSum.*s2);
                end
                if keepBlocks, kernelBlocks{slot} = kernel; end
            end
            powerSum = powerSum + real(sum(conj(feedWeights(block)).* ...
                (kernel*feedWeights)));
        end
        Prad = max(powerSum,0);
    end

    function tf = needsFinePowerIntegral()
        % True for models without an exact pair-sum power integral. The
        % built-ins above get exact radiated power, while these
        % asymmetric or user-supplied patterns use fine quadrature at
        % every array size. Custom and Imported
        % used to reach the quadrature only above a 12-wavelength span
        % and otherwise fell back to trapz on the 1 x 2 degree display
        % grid -- the very mechanism that under-resolves a narrow
        % broadside beam in the DENOMINATOR. Measured just under the old
        % threshold: +0.11 dB on D at 25x25 (0.5 lambda); a 20x20 at
        % 0.6 lambda reading ABOVE the 4*pi*A ideal, i.e. aperture
        % efficiency over 100%; a spurious 0.11 dB jump between 25x25
        % and 26x26 as the array crossed the cutoff; and the same
        % physical antenna reading 0.11 dB apart depending on whether it
        % was entered as the built-in cos^q or as the identical Custom
        % formula. The extra cost is modest: at the minimum n = 96 the
        % quadrature adds 18432 directions, against the 32761-point grid
        % the display computes anyway.
        tf = ismember(S.efType,{'Custom (formula)','Imported (CST far-field)', ...
            'Cardioid','Gaussian','Sinc','Crossed dipole (RHCP)', ...
            '3GPP TR 38.901 shape'});
    end

    function [nodes,weights] = legendreNodes(n)
        % Gauss-Legendre nodes on [-1,1], using Newton's method on P_n.
        j = (1:n).';
        nodes = cos(pi*(j-0.25)/(n+0.5));
        for iteration = 1:12
            p0 = ones(n,1); p1 = nodes;
            for degree = 2:n
                p2 = ((2*degree-1)*nodes.*p1-(degree-1)*p0)/degree;
                p0 = p1; p1 = p2;
            end
            derivative = n*(nodes.*p1-p0)./(nodes.^2-1);
            step = p1./derivative;
            nodes = nodes-step;
            if max(abs(step)) < 1e-14, break; end
        end
        weights = 2./((1-nodes.^2).*derivative.^2);
    end

    function Prad = quadraturePrad(feedWeights, fr)
        % For an arbitrary element factor, integrate over the visible
        % direction-cosine disk. Let u=x and v=sqrt(1-u^2)*y. The solid-
        % angle Jacobian becomes du*dy/sqrt(1-y^2), so Gauss-Legendre in
        % u and Gauss-Chebyshev in y resolve the rapidly oscillating array
        % factor with far fewer samples than a fine theta/phi rectangle.
        % fr: propagation ratio, default freqRatio().
        if nargin < 2, fr = freqRatio(); end
        feedWeights = feedWeights(:);
        % Sized from the bounding-box DIAGONAL, not the larger of the x
        % and y extents. The fastest oscillation in |AF|^2 comes from the
        % largest separation between any two elements, and along the
        % Gauss-Legendre direction a separation (dx,dy) oscillates at
        % k*hypot(dx,dy) -- up to sqrt(2) faster than max(dx,dy). The
        % 4-nodes-per-wavelength rule has about 27% margin over the
        % pi*L nodes needed, so a diagonal layout overran it: a 46-
        % element line along x=y read D 0.068 dB high at n=64, and still
        % 0.035 dB high at n=208. The diagonal bounds the largest
        % separation for any layout (and equals it for a filled
        % rectangle), which removes the error entirely.
        span = hypot(max(S.el(:,1))-min(S.el(:,1)), ...
            max(S.el(:,2))-min(S.el(:,2)))*fr;
        % n=64 aliases a theta oscillation that the 1-degree display grid
        % can resolve: E_theta=cosd(64*th) on a single custom element
        % reports 3.79 dBi instead of its closed-form 3.01 dBi. n=96
        % resolves every integer cosd(m*th) case through m=90, the
        % display grid's theta Nyquist limit, while keeping the small-
        % array quadrature cheaper than the display calculation.
        n = max(96,16*ceil(4*span/16));
        % Even finite, editor-valid coordinates can imply millions of
        % nodes here (and legendreNodes takes O(n^2) work). Refuse a
        % element-field integral we cannot resolve in a bounded
        % time rather than allocating an effectively unbounded grid or
        % silently using too few nodes and reporting a wrong gain.
        if n > 2048
            error('PAD:QuadratureSize',[ ...
                'This element pattern needs %g integration nodes for its ' ...
                '%.4g-wavelength operating span, above the 2048-node limit. ' ...
                'Reduce the element separation or operating/design frequency ratio, ' ...
                'or use an analytic built-in element model.'],n,span);
        end
        [u,wu] = legendreNodes(n);
        y = cos(pi*((1:n)-0.5)/n).';
        kOp = S.k*fr;
        [uniqueRot,~,rotIdx] = unique(S.el(:,5));
        Prad = 0;
        for first = 1:16:n
            cols = first:min(first+15,n);
            [U,Y] = meshgrid(u(cols),y);
            V = sqrt(max(1-U.^2,0)).*Y;
            thFront = asind(min(hypot(U,V),1));
            thBack = 180-thFront;
            ph = mod(atan2d(V,U),360);
            frontTh = zeros(size(U)); frontPh = frontTh;
            backTh = frontTh; backPh = frontTh;
            cache = numel(uniqueRot)<numel(feedWeights) && ...
                numel(uniqueRot)*4*numel(U)*16 <= 128e6;
            if cache
                ef = cell(numel(uniqueRot),4);
                for rotationIndex = 1:numel(uniqueRot)
                    [ef{rotationIndex,1},ef{rotationIndex,2}] = ...
                        elementFactor(thFront,ph,uniqueRot(rotationIndex));
                    [ef{rotationIndex,3},ef{rotationIndex,4}] = ...
                        elementFactor(thBack,ph,uniqueRot(rotationIndex));
                end
            end
            for element = 1:numel(feedWeights)
                if cache
                    rotationIndex = rotIdx(element);
                    etF = ef{rotationIndex,1}; epF = ef{rotationIndex,2};
                    etB = ef{rotationIndex,3}; epB = ef{rotationIndex,4};
                else
                    [etF,epF] = elementFactor(thFront,ph,S.el(element,5));
                    [etB,epB] = elementFactor(thBack,ph,S.el(element,5));
                end
                wave = feedWeights(element)*exp(1j*kOp* ...
                    (S.el(element,1)*U+S.el(element,2)*V));
                frontTh = frontTh+wave.*etF; frontPh = frontPh+wave.*epF;
                backTh = backTh+wave.*etB; backPh = backPh+wave.*epB;
            end
            power = fieldPower(frontTh,frontPh)+fieldPower(backTh,backPh);
            Prad = Prad+(pi/n)*sum(power.*wu(cols).','all');
        end
    end

    function peaks = refinedPatternPeaks(grids,theta,phi)
        % Refine AF, reference-element and TOTAL-power peaks independently
        % of the display mesh. Never infer a peak from the command alone:
        % manual phases, element patterns and competing lobes can move it.
        peaks = cellfun(@(p)max(p(:)),grids);
        seeds = zeros(0,3); % theta, phi, quantity index
        for kind = 1:3
            p = grids{kind}(1:end-1,:); % periodic phi endpoint
            candidates = find(localMaxMask(p));
            % Every theta = 0 cell is the SAME direction (+z), and every
            % theta = 180 cell is -z, so their values tie exactly and all
            % of them pass the local-maximum test. Sorted and cut to four,
            % they filled every seed slot with copies of one direction --
            % leaving nothing to refine a beam lying just off the pole.
            % Keep a single representative of each pole.
            [~,cc] = ind2sub(size(p),candidates);
            atPole = theta(cc(:)) == 0 | theta(cc(:)) == 180;
            keep = ~atPole(:);
            for poleTheta = [0 180]
                firstAt = find(theta(cc(:)) == poleTheta, 1);
                if ~isempty(firstAt), keep(firstAt) = true; end
            end
            candidates = candidates(keep);
            [~,order] = sort(p(candidates),'descend');
            candidates = candidates(order(1:min(4,numel(order))));
            [ir,ic] = ind2sub(size(p),candidates);
            seeds = [seeds; theta(ic(:)).', phi(ir(:)).', ...
                repmat(kind,numel(ic),1)]; %#ok<AGROW>
        end
        sBeam = squintSinTheta(S.theta_s);
        if abs(sBeam)<=1
            th = asind(sBeam); ph = S.phi_s;
            if th<0, th=-th; ph=ph+180; end
            for kind = [1 3]
                seeds = [seeds; th mod(ph,360) kind; ...
                    180-th mod(ph,360) kind]; %#ok<AGROW>
            end
        end
        seeds = unique(seeds,'rows');
        % The horizon is shared by both hemispheres. A grid peak at
        % theta=90 can lie just inside either side, but a fixed front
        % hemisphere flag cannot step through it. Search that seed from
        % both sides so a rearward peak near 90 degrees is not lost.
        horizonSeeds=seeds(seeds(:,1)==90,:);
        seeds=[seeds;horizonSeeds];
        weights = calculationWeights();
        weights = weights/sum(abs(weights));
        feed = weights.*exp(1j*deg2rad(effectivePhaseDeg()));
        [rotations,~,groups] = unique(S.el(:,5));
        kOp = S.k*freqRatio();

        % ---- refinement, in DIRECTION COSINES (u,v) plus a hemisphere
        % flag, with a step that belongs to each seed ----
        % The previous search stepped in (theta,phi) with ONE step shared
        % by every seed and halved on every iteration, and that failed in
        % two separate ways:
        %  - theta = 0 is a coordinate singularity, and theta was CLAMPED
        %    there rather than passed through. A seed on the pole could
        %    only drift towards the phi it already carried, so a beam
        %    lying 0.45 deg off broadside on the far side got no
        %    refinement at all: -1.14 dB on D at 64x64.
        %  - forced halving caps total travel at 2 grid cells (2 deg in
        %    theta, 4 in phi). On an elongated array a directive element
        %    pulls the TOTAL peak many degrees along the fan beam's wide
        %    axis, and that ridge runs diagonally through (theta,phi), so
        %    every seed stalled part-way up it: -0.64 dB at 64x8.
        % In (u,v) the pole is just the ordinary point (0,0), and a
        % rectangular array's fan-beam ridge lies along the u or v axis.
        % Each seed's step now GROWS while it keeps finding higher ground
        % and halves only when it cannot, so travel is unbounded; its
        % starting size follows the beamwidth separately along each axis
        % (half-width ~0.5/L), the same geometry as the grating-lobe
        % test's elliptical separation.
        nSeed = size(seeds,1);
        uS = sind(seeds(:,1)).*cosd(seeds(:,2));
        vS = sind(seeds(:,1)).*sind(seeds(:,2));
        hemi = sign(90 - seeds(:,1)); hemi(hemi==0) = 1;   % +1 front, -1 back
        if ~isempty(horizonSeeds)
            hemi(end-size(horizonSeeds,1)+1:end) = -1;
        end
        kindS = seeds(:,3);
        apX = (max(S.el(:,1))-min(S.el(:,1)))*freqRatio();
        apY = (max(S.el(:,2))-min(S.el(:,2)))*freqRatio();
        stepU = repmat(min(0.25, 0.5/max(apX,eps)), nSeed, 1);
        stepV = repmat(min(0.25, 0.5/max(apY,eps)), nSeed, 1);
        % The lone reference element has no aperture of its own: its
        % pattern is broad, so it starts from one moderate step.
        stepU(kindS==2) = 0.1; stepV(kindS==2) = 0.1;
        capU = 4*stepU; capV = 4*stepV;
        [dt,dp] = ndgrid(-1:1,-1:1);
        dt = dt(:).'; dp = dp(:).';   % 9 offsets; index 5 is (0,0), the centre
        for iteration = 1:100
            act = find(max(stepU,stepV) > 1e-5);
            if isempty(act), break; end
            nA = numel(act);
            uq = uS(act) + stepU(act).*dt;   % nA x 9
            vq = vS(act) + stepV(act).*dp;
            rq = hypot(uq,vq);
            inDisk = rq <= 1;                % outside it is not a direction
            tq = asind(min(rq,1));
            back = hemi(act) < 0;
            tq(back,:) = 180 - tq(back,:);
            pq = mod(atan2d(vq,uq),360);
            tq = tq(:).'; pq = pq(:).';
            uqf = uq(:).'; vqf = vq(:).';
            af = zeros(size(tq)); et = af; ep = af;
            for group = 1:numel(rotations)
                ix = find(groups==group);
                factor = zeros(size(tq));
                % Bound temporary phase matrices even for large imports.
                for first = 1:128:numel(ix)
                    block = ix(first:min(first+127,numel(ix)));
                    phase = kOp*(S.el(block,1)*uqf+S.el(block,2)*vqf);
                    factor = factor + sum(feed(block).*exp(1j*phase),1);
                end
                [eth,eph] = elementFactor(tq,pq,rotations(group));
                af = af+factor; et = et+factor.*eth; ep = ep+factor.*eph;
            end
            [eth,eph] = elementFactor(tq,pq,S.el(1,5));
            powers = {abs(af).^2,fieldPower(eth,eph),fieldPower(et,ep)};
            vals = -Inf(nA,9);
            for kind = 1:3
                rows = kindS(act) == kind;
                if ~any(rows), continue; end
                pw = reshape(powers{kind},nA,9);
                vals(rows,:) = pw(rows,:);
            end
            vals(~inDisk) = -Inf;
            for kind = 1:3
                rows = kindS(act) == kind;
                if any(rows)
                    peaks(kind) = max(peaks(kind), max(vals(rows,:),[],'all'));
                end
            end
            [best,col] = max(vals,[],2);
            moved = col ~= 5 & best > vals(:,5);
            pick = sub2ind([nA 9],(1:nA).',col);
            mv = act(moved);
            uS(mv) = uq(pick(moved)); vS(mv) = vq(pick(moved));
            stepU(mv) = min(2*stepU(mv),capU(mv));
            stepV(mv) = min(2*stepV(mv),capV(mv));
            st = act(~moved);
            stepU(st) = stepU(st)/2; stepV(st) = stepV(st)/2;
        end
    end

    function [E_th, E_ph] = elementFactor(TH,PH,rot)
        % Vector element pattern, split into orthogonal polarisation
        % components. Summing E_theta and E_phi SEPARATELY across elements
        % (and only then taking the magnitude) is what makes mixed-rotation
        % (sequential rotation / CP) arrays come out right - a scalar
        % magnitude cannot represent polarisation and produces an
        % artificial broadside null for a 0/90/180/270 CP array.
        % Rotation is measured from the +y axis (0 deg = up on the layout),
        % positive = clockwise (towards +x); az = 90 - rot converts to the
        % +x-referenced azimuth the formulas use.
        az = 90 - rot;
        ct = max(cosd(TH),0);
        switch S.efType
            case 'Isotropic'
                E_th = ones(size(TH)); E_ph = zeros(size(TH));
            case 'cos^q(theta)'
                % q is the conventional POWER-pattern exponent (G(theta) =
                % cos^q(theta), Balanis/Mailloux convention) -- so the FIELD
                % gets the half exponent, ct.^(q/2), giving |E|^2 = ct.^q.
                % Applying q directly to the field (ct.^q) would give a
                % cos^(2q) power pattern, twice the intended falloff.
                % Masked to the front hemisphere explicitly -- see
                % patchEF for why ct.^(q/2) cannot be trusted to do it
                % at q = 0, where 0^0 = 1 turns this one-sided model
                % omnidirectional without warning.
                E_th = (ct > 0).*(ct.^(S.efQ/2));  E_ph = zeros(size(TH));
            case 'Cardioid'
                % Quarter-wavelength two-source cardioid used by the
                % gallery reference, aimed along +z with a rear null.
                % It reaches half power at the horizon (theta = 90 deg).
                E_th = cosd(45*(1-cosd(TH)));
                E_ph = zeros(size(TH));
            case {'Gaussian','Sinc','3GPP TR 38.901 shape'}
                % Local signed azimuth/elevation offsets about +z. The
                % element rotation turns an asymmetric pattern with its
                % physical feed, as it does for Patch and imported CST.
                dAz = atan2d(sind(TH).*cosd(PH-az),cosd(TH));
                % Elevation uses the full horizontal projection, not just
                % the boresight component; atan2(y,z) narrows diagonal
                % cuts even though the two principal cuts look correct.
                dEl = asind(sind(TH).*sind(PH-az));
                if strcmp(S.efType,'Gaussian')
                    amp = exp(-2*log(2)*((dAz/S.efBeamAz).^2 + ...
                        (dEl/S.efBeamEl).^2));
                elseif strcmp(S.efType,'Sinc')
                    xHalf = 1.39155737825151; % sin(x)/x = 1/sqrt(2)
                    uAz = xHalf*sind(dAz)/sind(S.efBeamAz/2);
                    uEl = xHalf*sind(dEl)/sind(S.efBeamEl/2);
                    sAz = ones(size(uAz)); sEl = ones(size(uEl));
                    idxAz = abs(uAz)>1e-9; idxEl = abs(uEl)>1e-9;
                    sAz(idxAz) = sin(uAz(idxAz))./uAz(idxAz);
                    sEl(idxEl) = sin(uEl(idxEl))./uEl(idxEl);
                    amp = sAz.*sEl.*(cosd(TH)>0); % front-facing aperture
                else
                    % 3GPP TR 38.901 single-element attenuation SHAPE,
                    % capped at 30 dB. Absolute gain still follows this
                    % app's radiated-power normalization, not a fixed 8 dBi.
                    atten = min(12*(dAz/S.efBeamAz).^2 + ...
                        12*(dEl/S.efBeamEl).^2,30);
                    amp = 10.^(-atten/20);
                end
                if strcmp(S.efType,'3GPP TR 38.901 shape')
                    E_th = amp.*cosd(PH-az);
                    E_ph = -amp.*sind(PH-az);
                else
                    E_th = amp; E_ph = zeros(size(TH));
                end
            case 'Short dipole (z-axis)'
                % A free-space Hertzian dipole along the array normal.
                % Its power is sin^2(theta), including exact axial nulls
                % at theta = 0 and 180 degrees.
                E_th = sind(TH);
                E_ph = zeros(size(TH));
            case 'Dipole (linear pol)'
                % SIGNED cosine, not the clamped ct. ct = max(cos,0)
                % exists so a ground-plane-backed element vanishes behind
                % the array, and it is right for cos^q and Patch. A
                % dipole is a free-space element that radiates both ways,
                % and clamping only E_theta while leaving E_phi intact
                % was not even self-consistent: it deleted half of one
                % component and none of the other, giving a full-sphere
                % directivity of 12/7 (2.34 dBi) where this field's true
                % value is 3/2 (1.76 dBi).
                %
                % Identical over the upper hemisphere -- cosd(TH) IS ct
                % for theta <= 90 -- so the default hemisphere directivity
                % of 4.771 dBi is unchanged. Only full-sphere mode moves.
                E_th = cosd(TH).*cosd(PH-az);
                E_ph = -sind(PH-az);
            case 'Crossed dipole (RHCP)'
                % Two orthogonal short dipoles in the array plane with
                % equal amplitudes and -90-degree relative feed phase.
                % The sign produces RHCP under this app's E_theta/E_phi
                % decomposition at +z broadside.
                E_th = (cosd(TH).*cosd(PH-az) - ...
                    1i*cosd(TH).*cosd(PH-az-90))/sqrt(2);
                E_ph = (-sind(PH-az) + ...
                    1i*sind(PH-az-90))/sqrt(2);
            case 'Patch (cos^q x lin pol)'
                % E_phi is tangential to the ground plane at theta=90 (phi_hat
                % lies in the xy-plane there, theta_hat is normal to it) --
                % the PEC boundary condition requires the tangential field to
                % vanish at grazing, independent of q. The extra ct factor
                % therefore belongs on E_phi, not E_theta. Same field/power
                % exponent convention as the cos^q case above: ct.^(q/2).
                [E_th,E_ph] = patchEF(TH,PH,az,S.efQ);
            case 'Imported (CST far-field)'
                % Rotates the QUERY direction, not the data itself: for a
                % pure z-axis (azimuthal) rotation, phi_hat/theta_hat
                % rotate identically to the pattern, so a plain phi shift
                % of the query point correctly captures a physical
                % rotation by rot degrees -- same substitution the
                % Patch/Dipole analytic formulas make via az=90-rot
                % inside cosd(PH-az): expanding that shows the unrotated
                % reference pattern is evaluated at (PH+rot), which is
                % the same shift applied here.
                % mod(...,360) is REQUIRED here in a way it is not for
                % any analytic element factor above: those evaluate
                % cosd/sind, which are inherently 360-periodic, so a
                % query at 405 deg is automatically the same as 45 deg.
                % This branch instead samples a scatteredInterpolant
                % whose data spans phi 0..360 only -- a query past that
                % falls outside the convex hull and is 'nearest'-CLAMPED
                % to the boundary rather than wrapped round, freezing a
                % whole wedge of the pattern at one value. The wedge is
                % as wide as rot itself, so it scales with rotation:
                % measured 8.8% of the azimuth wrong at rot=45, 21.5% at
                % 90, 42.0% at 180 and 63.5% at 270 -- and 90/180/270
                % are exactly the values sequential rotation assigns, so
                % a CP array built from an imported element hit this on
                % most of its elements.
                PHq = mod(PH + rot, 360);
                if isempty(S.impFF)
                    % Missing data must never become a substitute antenna.
                    error('phasedArrayDesigner:missingPattern','Import the required far-field pattern first.');
                else
                    [E_th,E_ph] = evalImportedFF(S.impFF, TH, PHq);
                end
            case 'Custom (formula)'
                % User-typed E_theta and E_phi expressions, in the
                % variables th, ph (degrees), st=sind(th), and ct --
                % which is max(cosd(th),0), CLAMPED at zero, not plain
                % cosd(th). The clamp is what makes the built-in
                % ground-plane-backed elements vanish below the horizon,
                % and it is shared with this branch. It matters in
                % full-sphere mode: a custom formula written in terms of
                % ct simply cannot produce back radiation, since ct is
                % already 0 for every theta past 90. Use cosd(th)
                % directly for a signed cosine that goes negative there.
                % Invalid formulas throw and the calculation wrapper clears stale results.
                % ph carries the element ROTATION, exactly as the
                % built-in factors do through az = 90 - rot: their
                % formulas all depend on (PH - az) = PH - 90 + rot, so a
                % rotation of +rot is a +rot shift of the observation
                % azimuth. A custom formula had no such term at all, so
                % rotating a custom element left its pattern completely
                % unchanged -- measured 0.0000 dB difference between
                % rot = 0 and rot = 90 for "ct.*cosd(ph)", against
                % 3.0103 dB for the built-in Dipole under the same test.
                % Sequential rotation still moved its FEED phase, so a
                % custom element silently behaved as if it were
                % omnidirectional in azimuth while everything else in the
                % array rotated properly.
                th = TH; ph = mod(PH + rot, 360); st = sind(TH); %#ok<NASGU>
                theta = th; phi = ph; %#ok<NASGU>  -- aliases of the same values
                % A SCALAR result is broadcast to the grid -- that is what
                % makes "0" or "1" usable as a formula. Anything else has
                % to match the grid EXACTLY.
                %
                % This used to be `V = V.*ones(size(TH))` for any
                % mismatch, which is implicit expansion, not validation.
                % A result that merely happened to be broadcast-
                % compatible was silently reshaped instead of rejected: a
                % stray transpose -- "ct'", an easy typo -- turned a
                % 1x721 cut grid into 721x721. The wrong-shaped field
                % then reached the cut plot, where plot() draws one line
                % per ROW, so legend([hAF hEF hTot]) threw "Dimensions of
                % arrays being concatenated are not consistent". That
                % error escaped from the formula field's OWN callback,
                % which left the app throwing on every recompute until
                % the text was corrected by hand.
                % Non-finite results are rejected too. eval does not throw
                % on Inf or NaN, so the try/catch alone never saw them:
                % "1./ct" is a natural thing to try and is Inf at every
                % theta >= 90, where ct is clamped to zero, and "ct./ct"
                % is 0/0 there. Those flowed straight into the pattern and
                % turned the whole readout -- directivity, cut, 3D surface
                % -- into NaN with nothing said. A radiated field is
                % finite everywhere by definition, so this can only ever
                % reject a formula that was already wrong.
                try
                    E_th = eval(S.customFormula);
                    E_th = formulaFieldValue(E_th,size(TH));
                catch formulaErr
                    S.efFallback=true; S.efFallbackWhich(1)=true;
                    error('PAD:CustomFormula','E_theta formula failed: %s',formulaErr.message);
                end
                try
                    E_ph = eval(S.customFormulaPh);
                    E_ph = formulaFieldValue(E_ph,size(TH));
                catch formulaErr
                    S.efFallback=true; S.efFallbackWhich(2)=true;
                    error('PAD:CustomFormula','E_phi formula failed: %s',formulaErr.message);
                end
            otherwise
                E_th = ones(size(TH)); E_ph = zeros(size(TH));
        end
    end

    % ------------------------------------------------------ pattern calc
    function msk = localMaxMask(gridDb)
        %LOCALMAXMASK  True where a cell is >= all 8 neighbours.
        %   Phi wraps; theta does not. Plateaus qualify (>=), so a lobe
        %   sampled flat across two cells is still found, and a lobe
        %   peaking exactly on the horizon still qualifies because the
        %   edge simply has fewer neighbours to beat.
        [nR,nC] = size(gridDb);
        msk = true(nR,nC);
        for dr = -1:1
            for dc = -1:1
                if dr == 0 && dc == 0, continue; end
                shifted = circshift(gridDb, dr, 1);        % phi wraps
                if dc ~= 0
                    shifted = circshift(shifted, dc, 2);
                    if dc > 0
                        shifted(:,1:dc) = -Inf;            % theta does not
                    else
                        shifted(:,end+dc+1:end) = -Inf;
                    end
                end
                msk = msk & (gridDb >= shifted);
            end
        end
    end

    function idx = afClimb(gridDb, startIdx)
        %AFCLIMB  Walk uphill from startIdx to a local maximum.
        %   Used to find which lobe of the array factor contains a given
        %   point, so exactly that lobe -- and no other -- is excluded as
        %   the main beam. 8-connected, wrapping in phi, same neighbour
        %   rule as the downhill fill so the two agree on what a lobe is.
        [nR,nC] = size(gridDb);
        idx = startIdx;
        for step = 1:(nR*nC)          % bounded: cannot loop forever
            [pr,pc] = ind2sub([nR nC], idx);
            best = gridDb(pr,pc); bestIdx = idx;
            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0, continue; end
                    rr = mod(pr-1+dr, nR) + 1;
                    cc = pc + dc;
                    if cc < 1 || cc > nC, continue; end
                    if gridDb(rr,cc) > best
                        best = gridDb(rr,cc); bestIdx = sub2ind([nR nC], rr, cc);
                    end
                end
            end
            if bestIdx == idx, return; end   % local maximum reached
            idx = bestIdx;
        end
    end

    function msk = downhillFill(gridDb, seedIdx)
        %DOWNHILLFILL  8-connected strictly-downhill flood fill from one
        %   seed, wrapping in phi. Same rules as the main-lobe fill in
        %   findGratingLobe -- factored out so the level grid's own main
        %   lobe is traced by exactly the same walk rather than a second,
        %   subtly different one.
        [nR,nC] = size(gridDb);
        msk = false(nR,nC); seen = false(nR,nC);
        [sr,sc] = ind2sub([nR nC], seedIdx);
        qq = zeros(nR*nC,2); qh = 1; qt = 2;
        qq(1,:) = [sr,sc]; seen(sr,sc) = true; msk(sr,sc) = true;
        while qh < qt
            pr = qq(qh,1); pc = qq(qh,2); qh = qh + 1;
            cur = gridDb(pr,pc);
            for dr = -1:1
                for dc = -1:1
                    if dr == 0 && dc == 0, continue; end
                    rr = mod(pr-1+dr, nR) + 1;    % phi wraps
                    cc = pc + dc;                 % theta does not
                    if cc < 1 || cc > nC, continue; end
                    if ~seen(rr,cc) && gridDb(rr,cc) <= cur + 1e-9
                        seen(rr,cc) = true; msk(rr,cc) = true;
                        qq(qt,:) = [rr,cc]; qt = qt + 1;
                    end
                end
            end
        end
    end

    function [glTheta,glPhi,glRelDb,inTheta,inPhi,inRelDb] = ...
            findGratingLobe(afDbIn,theta,phi,levelDbIn)
        % Returns TWO lobes: the strongest secondary lobe anywhere
        % (gl*, used for the "how close did it come" / grazing readout),
        % and the strongest one whose peak is genuinely INSIDE visible
        % space (in*, which is what decides whether to call it a grating
        % lobe). They differ whenever a grazing lobe ties with an
        % interior one -- see the interior block at the end.
        % levelDbIn (optional) is the grid the LEVEL is read from, while
        % afDbIn is the grid whose shape defines the main lobe. They
        % differ when the caller wants the null structure of the ARRAY
        % FACTOR -- which is where nulls actually come from -- but the
        % level of the TOTAL pattern, which is what the user sees.
        %
        % Splitting them is necessary, not cosmetic: the flood-fill walks
        % strictly downhill, and an imported element pattern carries
        % enough interpolation ripple to stall that walk one cell from
        % the peak. Feeding the total in as BOTH left the main beam
        % outside its own mask and reported a point 5 deg off boresight,
        % 0.43 dB down, as a "grating lobe" at every steering angle.
        % General, geometry-agnostic first-grating-lobe finder: works
        % for ANY current lattice (plain rectangular, grid-angle-skewed,
        % staggered/brick, arbitrary shape-masked, even randomly
        % perturbed) without needing a closed-form lattice model --
        % staggering in particular turns the array into a proper 2-point
        % basis on a coarser primitive lattice, whose grating-lobe
        % structure factor isn't the simple single-lattice formula used
        % earlier in this conversation, so searching the ACTUAL computed
        % array factor directly is the robust way to get this right.
        %
        % Method: flood-fill outward from afDbIn's own GLOBAL MAXIMUM
        % (NOT the nominal steering direction theta_s/phi_s -- see
        % below), following only the DOWNHILL direction (a neighbour
        % joins the main-lobe region only if its level is <= the cell
        % that reached it). This traces out the main lobe's true
        % footprint out to its first null "moat" -- same idea as the
        % existing 1D-cut null-search elsewhere in this file, just
        % generalized to the full 2D hemisphere via BFS, so it can catch
        % a grating lobe that a skewed/staggered lattice has pushed OFF
        % the current principal cut (phi_s) entirely, which the older
        % cut-only search structurally cannot see. Then reports the
        % strongest peak found OUTSIDE that footprint -- the true first
        % grating lobe (or highest ordinary sidelobe, if no true grating
        % lobe exists in visible space).
        %
        % BUG FIX: this used to seed at (theta_s,phi_s) instead of the
        % true peak. Since the flood-fill can only walk DOWNHILL from
        % its seed, seeding anywhere off the actual peak means the
        % peak itself is "uphill" and UNREACHABLE -- it stays outside
        % mainMask and gets reported as a full-strength (0 dB) false
        % "grating lobe" AT THE MAIN BEAM'S OWN LOCATION. This is real
        % and reachable any time a manual per-element phase offset (the
        % editable Phase column) shifts the true peak away from the
        % nominal steering angle -- verified numerically: an 8x8 array
        % with a 15 deg/element manual phase ramp shifted the true peak
        % to (5,180) while theta_s/phi_s stayed (0,0); seeding at the
        % old (theta_s,phi_s) reported a false 0 dB "grating lobe"
        % exactly at that shifted true peak. Seeding at the actual
        % global max instead fixes this by construction: the seed IS
        % the peak, so it can never be excluded from its own mainMask.
        if nargin < 4 || isempty(levelDbIn), levelDbIn = afDbIn; end
        phiU = phi(1:end-1);              % drop the duplicate phi=360==phi=0 point
        afU  = afDbIn(1:end-1,:);         % size [numel(phiU) x numel(theta)]
        lvU  = levelDbIn(1:end-1,:);
        nPhiU = numel(phiU); nTh = numel(theta);

        [~,seedIdx] = max(afU(:));
        [seedPr,seedPc] = ind2sub(size(afU),seedIdx);

        visited = false(nPhiU,nTh);
        mainMask = false(nPhiU,nTh);
        q = zeros(nPhiU*nTh,2); qh = 1; qt = 2;
        q(1,:) = [seedPr,seedPc];
        visited(seedPr,seedPc) = true; mainMask(seedPr,seedPc) = true;
        while qh < qt
            pr = q(qh,1); pc = q(qh,2); qh = qh+1;
            curVal = afU(pr,pc);
            prm = mod(pr-2,nPhiU)+1;   % phi neighbour, wrapped
            prp = mod(pr,nPhiU)+1;
            % Fixed-size preallocation (max 4 neighbours: phi-wrap x2,
            % theta x2) instead of growing cand by concatenation each
            % BFS step -- same result, no repeated reallocation.
            % EIGHT neighbours, not four. The main lobe is an ellipse in
            % (u,v); on a SKEWED lattice it is also inclined, so on the
            % (theta,phi) sampling grid its ridge runs diagonally. A
            % 4-connected strictly-downhill walk cannot follow a diagonal
            % ridge -- stepping in theta alone or phi alone goes UP before
            % it comes down -- so the fill stalled part-way and left a
            % piece of the MAIN BEAM outside its own mask. That piece was
            % then reported as a full-height repeat: at dx=dy=0.58 lambda,
            % grid angle 45 deg, steering to 36 deg, the tool announced
            % "GRATING LOBE ... theta = 37 deg, phi = 358 deg, 0.3 dB
            % below peak" -- one grid cell from a main beam at
            % (36 deg, 0 deg). Verified against the reciprocal lattice
            % (b1 = (1.724,-1.724), b2 = (0,2.438); the nearest points are
            % |g| = 1.87, needing sin(theta_s) > 1, so this lattice has NO
            % grating lobe at any scan angle) and against a direct array
            % factor, whose strongest lobe outside the beam is -13.9 dB at
            % every angle from 30 to 50 deg -- an ordinary sidelobe.
            %
            % Adding the diagonals cannot leak the fill into a real
            % grating lobe: the walk only ever moves downhill, and a
            % grating lobe is a local MAXIMUM, so the level rises on
            % approach and blocks it. It can only reach further around the
            % null moat that already bounds the main lobe.
            cand = zeros(8,2); nCand = 0;
            nCand = nCand+1; cand(nCand,:) = [prm,pc];
            nCand = nCand+1; cand(nCand,:) = [prp,pc];
            if pc>1
                nCand = nCand+1; cand(nCand,:) = [pr, pc-1];
                nCand = nCand+1; cand(nCand,:) = [prm,pc-1];
                nCand = nCand+1; cand(nCand,:) = [prp,pc-1];
            end
            if pc<nTh
                nCand = nCand+1; cand(nCand,:) = [pr, pc+1];
                nCand = nCand+1; cand(nCand,:) = [prm,pc+1];
                nCand = nCand+1; cand(nCand,:) = [prp,pc+1];
            end
            cand = cand(1:nCand,:);
            for kk = 1:size(cand,1)
                rr = cand(kk,1); cc = cand(kk,2);
                if ~visited(rr,cc) && afU(rr,cc) <= curVal + 1e-9
                    visited(rr,cc) = true; mainMask(rr,cc) = true;
                    q(qt,:) = [rr,cc]; qt = qt+1;
                end
            end
        end

        % ---- the LEVEL grid's own main lobe must also be masked -------
        % mainMask above is flood-filled on afDbIn (the array factor),
        % whose shape gives the reliable null structure. But peakDb below
        % is read from lvU (the TOTAL), and those two peak in DIFFERENT
        % places once elements carry different rotations: the rotation
        % feed phase is not a linear ramp, so the AF peak walks off the
        % steer direction while the total still peaks on it. The total's
        % own peak then sat OUTSIDE the AF-seeded mask and was reported
        % as a grating lobe 0.0 dB below the peak -- i.e. the main beam,
        % announced as its own grating lobe. Measured on a healthy 8x8
        % (23.19 dBi, 103.8%% aperture efficiency) with sequential
        % rotation and feed phase on: "GRATING LOBE: theta=0 deg,
        % phi=22 deg, 0.0 dB below peak", at the steer direction itself.
        %
        % So the level grid gets its own downhill fill from its own
        % maximum, and the two masks are unioned. Nothing that is the
        % top of either pattern's main lobe can be called a lobe beside
        % it. A genuine grating lobe is a local maximum in both and is
        % blocked from both fills exactly as before.
        [~,lvSeed] = max(lvU(:));
        if ~mainMask(lvSeed)
            % ONE region, not a union. The first attempt at this added a
            % second fill seeded on the level grid and OR-ed the two
            % masks, which fixed the false positive but created a worse
            % failure: when the AF and the total peak on GENUINELY
            % DIFFERENT lobes, the union swallows both, and a real second
            % beam only 3 dB down is masked out of existence with nothing
            % left to report.
            %
            % The main beam is the lobe the user is looking at, which is
            % the one carrying the TOTAL's peak -- that is also where
            % peakDb is measured, so it is the one that must be excluded.
            % Its extent is still traced on the ARRAY FACTOR, whose null
            % structure is clean (the total's interpolation ripple stalls
            % a strictly-downhill walk, which is why the shape came from
            % the AF in the first place).
            %
            % Climb the AF to the local maximum that OWNS the total's
            % peak, then fill downhill from there. That yields exactly the
            % AF main-lobe region containing the total's peak: one beam
            % excluded, every other lobe left as a candidate.
            mainMask = downhillFill(afU, afClimb(afU, lvSeed));
        end
        % ---- pole degeneracy ------------------------------------------
        % Every phi column at theta = 0 is the SAME physical direction,
        % so if one of them is in the main lobe they all are. Without
        % this the fill can mask the pole cell it happened to seed and
        % leave its duplicates outside, which is the other half of the
        % theta=0 false lobe above.
        if any(mainMask(:,1)), mainMask(:,1) = true; end
        if abs(theta(end) - 180) < 1e-9 && any(mainMask(:,end))
            mainMask(:,end) = true;
        end

        peakDb = max(lvU(:));
        % ---- a lobe is a LOCAL MAXIMUM, separated by a null ------------
        % The mask alone cannot decide this, and two attempts proved it.
        % Filling from the AF peak alone left the total's peak outside its
        % own mask (main beam reported as its own grating lobe at
        % theta=0). Unioning a second fill from the total's peak fixed
        % that but swallowed BOTH lobes whenever the two grids peak on
        % genuinely different ones, hiding a real second beam 3 dB down.
        % Filling one region from the total's peak fixed THAT and brought
        % back the first failure in its original form: "theta=1 deg,
        % 0.1 dB below peak" on a healthy 23.19 dBi array -- a cell one
        % degree off the beam, on its shoulder.
        %
        % None of those are lobes. A grating lobe is a LOCAL MAXIMUM of
        % the pattern being measured, and a shoulder cell 0.1 dB below the
        % peak is not one however the mask is drawn. Requiring local
        % maximality is a property of the candidate itself, so it holds
        % regardless of which region the fill happened to cover -- which
        % is what makes it robust where the masks were not.
        isPk = localMaxMask(lvU);
        % ...and SEPARATED from the main beam by at least half a
        % first-null spacing. Local maximality alone was not enough: near
        % theta = 0 every phi column is nearly the same direction, so
        % ripple across that ring manufactures local maxima a tenth of a
        % dB down, and the healthy array was still accused at
        % "theta=1 deg, 0.1 dB below peak".
        %
        % An aperture L wavelengths across puts its first null at
        % du ~ 1/L, so anything closer than 0.5/L in direction cosine is
        % inside the main beam by construction, whatever the sampling
        % does. For the 8x8 at 0.5 lambda that is du = 0.14, while
        % theta = 1 deg is du = 0.017 -- an order of magnitude inside.
        % A genuine grating lobe sits a full du = 1 away at 1.0 lambda
        % spacing, so this cannot suppress one.
        [pkR,pkC] = ind2sub(size(lvU), find(lvU == peakDb, 1));
        u0 = sind(theta(pkC))*cosd(phiU(pkR));
        v0 = sind(theta(pkC))*sind(phiU(pkR));
        [PHc,THc] = ndgrid(phiU, theta);
        % Aperture in OPERATING wavelengths, not design ones. S.el holds
        % positions in design wavelengths while the propagation term
        % carries freqRatio(), so in squint mode the two disagree by
        % exactly that factor -- and the disagreement suppresses real
        % lobes rather than false ones. Measured: a 2x2 isotropic array
        % at dx = dy = 0.5 design wavelengths, design 10 GHz, operating
        % 40 GHz, is physically on a 2 lambda lattice and has a true 0 dB
        % replica at theta = 30 deg (du = 0.5). The unscaled radius asked
        % for du > 1.0 and threw it away; scaled, it asks for du > 0.25
        % and keeps it.
        %
        % HEURISTIC, and worth naming as one: 0.5/L is the half-first-null
        % spacing of a uniformly excited aperture. It is not a guaranteed
        % main-beam boundary for every taper, lattice and element pattern
        % this app supports -- a heavily tapered array has a wider beam
        % and could in principle still place a shoulder sample outside it.
        % It is paired with the local-maximum test for that reason;
        % neither is sufficient alone.
        % ELLIPTICAL, one semi-axis per aperture dimension -- not a
        % single circular radius taken from the LARGER aperture. The main
        % beam's half-width is ~0.5/L separately along u (set by the x
        % aperture) and along v (set by the y aperture). A circle sized
        % from max(apX,apY) is the NARROW axis's width applied in every
        % direction, so on an elongated array the fan beam's wide axis
        % reached far outside it. On 64x8 at 0.5 lambda the circle was
        % du = 0.016 while the beam extends to dv ~ 0.29, so grid ripple
        % along the main beam's own crest was reported as a "possible
        % grating lobe: theta=21, phi=48, 0.1 dB below peak", and Max
        % scan claimed a grating lobe onset at theta_s = 2 deg on a
        % lattice that has none before 90.
        %
        % THRESHOLD 0.8, measured rather than chosen. Normalising by the
        % bounding-box extent (M-1)*d puts an M-element uniform main
        % lobe's -6 dB contour at ~0.53-0.59, just PAST the old 0.5, so a
        % sliver of the main beam above the reporting level still leaked
        % through: 8/60 random steers on 64x8 at 0.5 lambda (12/60
        % tapered) still read "possible grating lobe" on a lattice with
        % none. At 0.8, over 200 steers: 1/200 uniform, 1/200 moderately
        % tapered, 4/200 with a heavy near-Hann taper.
        %
        % It cannot suppress a genuine grating lobe. Along an axis holding
        % M elements at pitch d the replica sits at du = lambda/d while the
        % aperture is (M-1)*d, so its normalised distance is M-1 >= 1; the
        % tightest case, 2 elements at 1 lambda, lands at exactly 1.0 and
        % 0.8 leaves 20% headroom for grid sampling. Checked: 140 genuine-
        % lobe steers across 2-row, 2-column, 3-row, 1.0 and 1.5 lambda
        % lattices, none missed.
        %
        % Known residual, named rather than hidden: a HEAVY taper widens
        % the main lobe past any fixed threshold below 1.0, while the
        % nearest genuine lobe sits exactly at 1.0, so no single value
        % separates them perfectly -- hence the remaining 4/200. The
        % label says "possible" for that reason.
        %
        % Normalising each axis by its own aperture reduces EXACTLY to a
        % circular test when apX = apY. A single row (apY = 0) correctly
        % puts no limit along v, where a line array's pattern does not
        % vary; a single element (both 0) correctly reports no lattice
        % lobe.
        apX = (max(S.el(:,1)) - min(S.el(:,1))) * freqRatio();
        apY = (max(S.el(:,2)) - min(S.el(:,2))) * freqRatio();
        uG = sind(THc).*cosd(PHc); vG = sind(THc).*sind(PHc);
        farEnough = hypot((uG - u0)*apX, (vG - v0)*apY) > 0.8;
        outside = ~mainMask & isPk & farEnough;
        if any(outside(:))
            % Masked by the full candidate test, not by mainMask alone.
            % The two had drifted apart: `outside` gained the local-max
            % and separation requirements while this line still cleared
            % only the main-lobe region, so the filters gated whether a
            % lobe was reported but not WHICH -- and the shoulder cell at
            % theta = 1 deg went on winning max().
            afOut = lvU; afOut(~outside) = -Inf;
            [glDb,glIdx] = max(afOut(:));
            [glPr,glPc] = ind2sub(size(lvU),glIdx);
            glTheta = theta(glPc); glPhi = phiU(glPr);
            glRelDb = glDb - peakDb;
        else
            % the whole hemisphere is one connected downhill region from
            % the main beam -- no secondary lobe anywhere in visible space
            glTheta = NaN; glPhi = NaN; glRelDb = -Inf;
        end

        % ---- TIE-BREAK: an interior lobe among the equal-strongest -----
        % A lobe has entered visible space when its peak is interior in
        % theta: the level one step further out is LOWER. The horizon
        % slicing through a still-rising skirt is not. That is the same
        % falloff test the callers used to apply themselves; it lives
        % here now so it can be applied to more than one candidate.
        %
        % Compare local peak estimates, not the raw mesh samples. A narrow
        % 64x64 grating lobe at theta=68.47 deg can sample 0.2 dB below a
        % grazing lobe even though both have the same true height. A fixed
        % 0.05 dB raw-sample tie then hides the interior lobe. The small
        % parabolic correction below removes that sampling preference
        % without opening the search to unrelated weak sidelobes.
        TIE_DB = 0.05;
        falls = false(nPhiU,nTh);
        falls(:,1:nTh-1) = lvU(:,2:nTh) < lvU(:,1:nTh-1);
        if any(outside(:))
            inMask = outside & falls;
        else
            inMask = false(nPhiU,nTh);
        end
        if any(inMask(:))
            candidateIdx = find(inMask);
            [candidateR,candidateC] = ind2sub(size(lvU),candidateIdx);
            candidateDb = lvU(candidateIdx);
            refinedDb = candidateDb;
            insideTheta = candidateC > 1 & candidateC < nTh;
            if any(insideTheta)
                j = find(insideTheta);
                lo = lvU(sub2ind(size(lvU),candidateR(j),candidateC(j)-1));
                hi = lvU(sub2ind(size(lvU),candidateR(j),candidateC(j)+1));
                curvature = lo + hi - 2*candidateDb(j);
                offset = (lo-hi)./(2*curvature);
                good = curvature < -1e-9 & abs(offset) <= 0.5;
                refinedDb(j(good)) = refinedDb(j(good)) - ...
                    (hi(good)-lo(good)).^2./(8*curvature(good));
            end
            loR = mod(candidateR-2,nPhiU)+1;
            hiR = mod(candidateR,nPhiU)+1;
            lo = lvU(sub2ind(size(lvU),loR,candidateC));
            hi = lvU(sub2ind(size(lvU),hiR,candidateC));
            curvature = lo + hi - 2*candidateDb;
            offset = (lo-hi)./(2*curvature);
            good = curvature < -1e-9 & abs(offset) <= 0.5;
            refinedDb(good) = refinedDb(good) - ...
                (hi(good)-lo(good)).^2./(8*curvature(good));
            refinedDb(refinedDb < glDb - TIE_DB) = -Inf;
            [inDb,localIdx] = max(refinedDb);
            inIdx = candidateIdx(localIdx);
            if ~isfinite(inDb)
                inTheta = NaN; inPhi = NaN; inRelDb = -Inf;
                return;
            end
            inDb = refinedDb(localIdx); % use the estimated peak for the -6 dB decision
            [inPr,inPc] = ind2sub(size(lvU),inIdx);
            inTheta = theta(inPc); inPhi = phiU(inPr);
            inRelDb = inDb - peakDb;
        else
            inTheta = NaN; inPhi = NaN; inRelDb = -Inf;
        end
    end

    function computePatternCore()
        if ~patternReady(true), return; end
        normalizedAmp = calculationWeights();
        if isempty(S.el)
            uialert(fig,'Place at least one element first.','No elements');
            return;
        end
        % Cleared here and checked at the end: elementFactor sets it if a
        % Custom formula blew up on the real grid and had to be replaced.
        S.efFallback = false; S.efFallbackWhich = [false false];
        % Display clipping must never change radiated power or directivity.
        calcFull = true;
        theta = 0:1:180;
        phi   = 0:2:360;
        [TH,PH] = meshgrid(theta,phi);
        nThUp = numel(0:1:90);   % column count of the upper hemisphere --
                                 % TH is [numel(phi) x numel(theta)], so
                                 % columns 1:nThUp are theta = 0..90 in
                                 % both grid modes. Used to keep grating-
                                 % lobe detection on the hemisphere.

        us = sind(S.theta_s)*cosd(S.phi_s);
        vs = sind(S.theta_s)*sind(S.phi_s);

        U = sind(TH).*cosd(PH);
        V = sind(TH).*sind(PH);

        % Steering and propagation both use operating/design frequency;
        % physical positions stay referenced to the design wavelength.
        fr = freqRatio();

        Etot_th = zeros(size(TH));  % E-theta component of the total field
        Etot_ph = zeros(size(TH));  % E-phi component
        Eaf     = zeros(size(TH));  % array factor alone (elements isotropic)

        % ---- element-factor cache, keyed on ROTATION ----
        % The element pattern depends only on the element's rotation, so
        % evaluating it once per element re-does identical work for
        % every element sharing a rotation -- and most arrays use very
        % few distinct ones (a 2x2 sequential-rotation tile repeats just
        % 0/90/180/270 no matter how large the array gets). For an
        % Imported element that inner call is four scatteredInterpolant
        % evaluations over the whole angle grid, so on a 1000-element
        % array this was ~250x more work than needed.
        %
        % GUARDED, because the number of distinct rotations is NOT
        % bounded: applySeqRot with blockM=M/blockN=N assigns a
        % different rotation to every element, as does randomising them.
        % Caching blindly would allocate one grid pair per rotation --
        % 1.05 MB each on the full-sphere grid, so ~4.3 GB for a 64x64
        % array of unique rotations -- turning a slow UI into an
        % out-of-memory failure on exactly the large arrays this is
        % meant to speed up. Two conditions must hold:
        %   - fewer distinct rotations than elements, or there is no
        %     work to save in the first place (the degenerate all-unique
        %     case does exactly as many evaluations either way);
        %   - the cache fits a fixed byte budget.
        % unique's third output gives the per-element index directly, so
        % the lookup inside the loop is O(1) rather than a search.
        % Sign of the sequential-rotation feed phase. Hoisted above both
        % element loops in this function: seqRotSign() costs an element-
        % factor evaluation on a ring and must not run per element.
        seqSgn = seqRotSign();
        rotAll = S.el(:,5);
        [uniqRot,~,rotIdx] = unique(rotAll);
        nUR = numel(uniqRot);
        bytesPerRot = 2*numel(TH)*16;          % Eth + Eph, complex double
        useEFCache  = (nUR < numel(rotAll)) && (nUR*bytesPerRot <= 256e6);
        if useEFCache
            EthCache = cell(nUR,1); EphCache = cell(nUR,1);
            for i = 1:nUR
                [EthCache{i}, EphCache{i}] = elementFactor(TH,PH,uniqRot(i));
            end
        end

        for n = 1:size(S.el,1)
            xn = S.el(n,1); yn = S.el(n,2);
            amp = normalizedAmp(n); ph0 = S.el(n,4); rot = S.el(n,5);

            % steering phase: -k * r_n . u_hat(theta_s,phi_s). Its frequency
            % reference is operating frequency in retuned mode and design
            % frequency in Beam-squint mode (phaseFreqRatio()).
            phFeed = -S.k*phaseFreqRatio()*(xn*us + yn*vs) + deg2rad(ph0);
            if S.seqPhase, phFeed = phFeed + seqSgn*deg2rad(rot); end

            steer = amp*exp(1j*phFeed)*exp(1j*S.k*fr*(xn*U + yn*V));

            Eaf  = Eaf  + steer;                                % pure array factor

            % vector components summed SEPARATELY across elements; the
            % magnitude is taken only after the sum. This is what makes
            % mixed-rotation (sequential rotation / CP) arrays correct.
            if useEFCache
                Eth = EthCache{rotIdx(n)};  Eph = EphCache{rotIdx(n)};
            else
                [Eth, Eph] = elementFactor(TH,PH,rot);
            end
            Etot_th = Etot_th + steer.*Eth;
            Etot_ph = Etot_ph + steer.*Eph;
        end

        % combine polarisation components AFTER the sum
        if any(~isfinite(Etot_th(:)))||any(~isfinite(Etot_ph(:))) || ...
                ~any(abs(Etot_th(:))>0 | abs(Etot_ph(:))>0)
            error('PAD:Radiation','No finite nonzero total radiation is available.');
        end
        Etot = polCombine(Etot_th, Etot_ph);

        % single-element factor, evaluated at a reference rotation
        refRot = S.el(1,5);
        [Eef_th, Eef_ph] = elementFactor(TH,PH,refRot);
        Eef = polCombine(Eef_th, Eef_ph);
        mixedRot = numel(unique(S.el(:,5))) > 1;

        % common 0 dB reference: the max possible coherent sum. Uses |amp| so
        % mixed-sign designs (difference / monopulse beams) work correctly.
        ref = sum(abs(normalizedAmp));
        if ~(ref > 0)
            uialert(fig,['Every element amplitude is zero - nothing to plot. ' ...
                'Set at least one amplitude to a nonzero value.'],'Zero amplitude');
            return;
        end
        % Display floor: 80 dB below THIS pattern's own peak, not a fixed
        % 1e-4 on the ratio.
        %
        % The fixed floor silently assumed the ratio peaks near 1, which
        % it does for every built-in element -- and does not for a Custom
        % formula, whose amplitude the user sets freely. The floor is
        % applied before the absolute-level offset is added, and that
        % offset comes from the radiated-power integral of the UNFLOORED
        % field, so the two stopped describing the same pattern: a
        % hemispherical constant field with E_theta = 1e-7 has exactly
        % the directivity of one with E_theta = 1 (a directivity is a
        % ratio; a uniform rescale cannot move it), yet the floored curve
        % read 63.01 dBi against a true 3.01, and 103.01 dBi at 1e-9 --
        % 20 dB of pure artefact per decade, while the title alongside
        % went on reporting the correct 3.01.
        %
        % Scaling the floor with the peak makes it what it was always
        % meant to be: a fixed amount of DYNAMIC RANGE. For any pattern
        % whose ratio peaks near 1 the two are identical, so nothing
        % about the built-in elements changes.
        %
        % THREE references, one per curve, each taken from its OWN
        % quantity. The first version of this fix derived a single floor
        % from the total field and applied it to all three, on the
        % reasoning that a shared floor preserves the curves' positions
        % relative to each other. That only holds while all three peak
        % near 1 -- the very assumption the fixed 1e-4 floor was wrong
        % to make. An element amplitude of 1e7 puts the shared floor at
        % 1e3, which is above the ENTIRE array factor (|AF|/ref <= 1),
        % so the AF curve flattened and rose 60 dB. The array factor is
        % a property of positions, amplitudes and phases alone; nothing
        % about the element pattern's amplitude may reach it.
        %
        % Each reference is the BOTH-POLARIZATION magnitude over the
        % FULL grid:
        %   - both-polarization, so selecting an RHCP/LHCP readout that
        %     happens to be identically zero floors against the field
        %     that is actually there, and the polarization loss stays
        %     visible instead of being normalised away;
        %   - full-grid, so the cut block below can reuse these rather
        %     than deriving a floor from the cut it is about to draw.
        dispRef = [max(abs(Eaf(:)))/ref, ...
                   sqrt(max(abs(Eef_th(:)).^2  + abs(Eef_ph(:)).^2)), ...
                   sqrt(max(abs(Etot_th(:)).^2 + abs(Etot_ph(:)).^2))/ref];
        dispFloor = dispRef * 1e-4;
        % A reference of zero means that quantity is identically zero
        % everywhere, which is the one case with no peak to be relative
        % to; the old flat floor is right there and invents nothing,
        % because a globally dead pattern has no absolute offset either
        % (safeOff returns NaN and the caller zeroes it).
        dispFloor(~(dispFloor > 0)) = 1e-4;
        afFloor = dispFloor(1); efFloor = dispFloor(2); totFloor = dispFloor(3);
        afDb  = 20*log10(max(abs(Eaf)/ref,  afFloor));
        totDb = 20*log10(max(abs(Etot)/ref, totFloor));
        efDb  = 20*log10(max(abs(Eef),      efFloor));

        % ---- Axial ratio (dB), from the TOTAL field's own RHCP/LHCP
        % split -- AR = (|E_R|+|E_L|) / ||E_R|-|E_L||, the standard
        % relation between axial ratio and the same circular-basis
        % decomposition polCombine already uses (1 / 0 dB for pure CP,
        % -> Inf for pure linear). ARCEIL_DB caps the display: AR is
        % genuinely unbounded (blows up wherever the field is purely
        % linearly polarized), so an uncapped scale would only be
        % meaningful right at the CP sweet spot and useless everywhere
        % else. 30 dB is the standard "might as well call it linear"
        % practical ceiling.
        ARCEIL_DB = 30;
        % How far below the peak the field has to fall before AR is
        % treated as meaningless and pinned to the ceiling. Shared by the
        % 3D surface and the cut plot, which each apply it to their own
        % grid -- it used to be a bare 20 written out separately in both
        % places, ~480 lines apart, so the two could silently disagree
        % about what counts as a null.
        %
        % 40 dB, not 20. At 20 the mask was swallowing perfectly good
        % data: a single dual-fed CP patch has NO nulls at all in the
        % forward hemisphere and an exactly known AR of 1/cos(theta), yet
        % at q=2 the app reported the 30 dB "essentially linear" ceiling
        % from about theta = 82 deg onward, where the field is only 17 dB
        % down and the true AR is 17.1 dB. That is a smooth roll-off, not
        % a null, and 14% of peak amplitude is not ill-conditioned in
        % double precision. Checking element AR against angle is exactly
        % what a CP element design needs, and the ceiling hid it.
        %
        % 40 dB still covers what the mask exists for. Verified on arrays
        % of LINEARLY polarised elements, whose true AR is infinite so any
        % dip below the ceiling is the null ill-conditioning itself:
        % 8x8 broadside, 8x8 steered 40 deg, 16x16 steered 30 deg, 16x16
        % Chebyshev and 32x1 all returned a flat 30.00 dB with zero
        % samples below 29 dB.
        ARNULL_DB = 40;
        E_R_tot = (Etot_th + 1j*Etot_ph)/sqrt(2);
        E_L_tot = (Etot_th - 1j*Etot_ph)/sqrt(2);
        magSumTot = abs(E_R_tot) + abs(E_L_tot);
        magDifTot = abs(abs(E_R_tot) - abs(E_L_tot));
        % Formed as a RATIO in [0,1] rather than as sum/difference with
        % the difference floored at eps. That floor was an absolute one
        % standing next to a user-scalable field, and it inverted the
        % answer for small amplitudes: a purely LINEAR field has an
        % identically zero difference, so eps took over as the
        % denominator, and once the sum itself fell below eps the ratio
        % dropped under 1, clamped to 1, and reported 0 dB -- perfect
        % circular polarisation for a linear field. Measured at
        % E_theta = 1e-12 and 1e-18.
        %
        % r = ||E_R|-|E_L|| / (|E_R|+|E_L|) is bounded by 1 by the
        % triangle inequality, is 1 for pure CP and 0 for pure linear,
        % and is invariant under any uniform rescale of the field
        % because numerator and denominator scale together. AR = -20
        % log10(r) then gives the same 0 dB for CP and the ceiling for
        % linear, with no clamp needed to keep it non-negative.
        % Clamped to [0, ARCEIL_DB]. The lower clamp is not cosmetic: r is
        % bounded by 1 in exact arithmetic, but round-off can put it a few
        % ulp above, and -20*log10 of that is a small NEGATIVE axial
        % ratio -- a quantity that does not exist.
        rTot = magDifTot ./ max(magSumTot, realmin);
        arDb = min(max(-20*log10(max(rTot, 10^(-ARCEIL_DB/20))), 0), ARCEIL_DB);
        % Near/at a pattern NULL the total field is negligible, so AR
        % becomes numerically ill-conditioned (a ratio of two tiny,
        % floating-point-noise-dominated magnitudes) -- not a genuine
        % "perfect CP" result. Confirmed this was producing spurious
        % streaks/spikes right at every sidelobe null. Masked by GAIN
        % level relative to the pattern's own peak, not a tiny absolute
        % field-magnitude threshold (which was too narrow to cover the
        % whole ill-conditioned neighbourhood around each null): more
        % than 20 dB below peak, there's essentially no radiated power
        % there for a polarization state to be meaningful about.
        %
        % Uses a TRUE both-polarization power reference (from Etot_th/
        % Etot_ph directly), NOT totDb -- totDb goes through polCombine,
        % which depends on the Polarization readout dropdown. If that's
        % set to RHCP/LHCP component instead of Total, totDb would read
        % as a false "null" anywhere the field is purely the OPPOSITE
        % circular sense, even though real power (and a perfectly
        % meaningful AR) exists there -- same class of bug documented in
        % fieldDirectivity's own P_base comments, reintroduced here by
        % reusing totDb instead of an independent total-power quantity.
        % Two separate conditions, and only one of them was handled.
        %
        % RELATIVE nulls -- a direction far below the pattern peak -- are
        % ill-conditioned and get the ceiling, as before.
        %
        % ABSOLUTE zeros are different: where the field is identically
        % zero there is no polarisation state at all. AR came out as
        % 20*log10(max(0/eps,1)) = 0 dB there, i.e. PERFECT circular
        % polarisation reported for a direction that radiates nothing. A
        % Patch at theta = 90 does exactly this -- both components carry
        % a cos(theta) factor -- so an AR cut taken there was uniformly
        % 0 dB. The relative test cannot catch it either: with every
        % sample at the floor, none of them is 40 dB below the maximum.
        % Read from the UNFLOORED power, since the floor is what hides it.
        pTrue = fieldPower(Etot_th,Etot_ph)/ref^2;
        % RELATIVE to this pattern's own maximum, not an absolute floor.
        % The first version tested pTrue <= 1e-12 outright, but ref is
        % sum|amp| -- the element amplitudes -- not anything tied to the
        % field's own scale, so pTrue carries whatever absolute scale the
        % element factor happens to have. A Custom formula is the clear
        % case: E_theta = 1, E_phi = -1j and E_theta = 1e-7, E_phi =
        % -1j*1e-7 are the SAME polarisation state (uniform scaling
        % cannot change an axial ratio), yet the second tripped the
        % absolute floor and was reported as "no field" at 30 dB. There
        % is no calibrated absolute scale here to compare against.
        %
        % 1e-20 is -200 dB relative. Chosen to sit far below any null a
        % real pattern produces, yet far above the round-off floor of the
        % coherent sum (~(eps*N)^2, about 1e-25 for a 1024-element
        % array), so it separates "identically zero" from "very small"
        % without either misfiring. A true null that does fall below it
        % is already caught by nullMask and pinned to the same ceiling,
        % so nothing changes for those directions.
        %
        % pRefMax is deliberately the FULL-pattern maximum and is reused
        % by the cut below: a cut can be entirely dead (a Patch along
        % theta = 90) and its own maximum is then zero, which would make
        % a cut-local reference divide the test by nothing.
        pRefMax = max(pTrue(:));
        if pRefMax > 0
            noFieldMask = pTrue <= pRefMax * 1e-20;
        else
            noFieldMask = true(size(pTrue));   % nothing radiates anywhere
        end
        % Compared as a power RATIO against the pattern's own maximum.
        % The dB form this replaces floored the power at an absolute 1e-8
        % first, so once a rescaled custom field sat entirely below that
        % floor every sample flattened to -80 dB, the spread vanished and
        % the mask selected nothing -- the mask moved with the field's
        % scale, which is exactly what a relative test is supposed to
        % rule out. pRefMax is the unfloored maximum computed above.
        nullMask = pTrue < pRefMax * 10^(-ARNULL_DB/10);
        arDb(nullMask) = ARCEIL_DB;
        arDb = min(arDb, ARCEIL_DB);
        % Pinned to the ceiling, like any other null. NaN would be more
        % literally "no data", but it puts non-finite values into plotted
        % output, and the suite's blanket no-NaN check exists to catch
        % genuine computation failures -- weakening it to admit a
        % deliberate NaN would blind it to the accidental kind. The
        % ceiling reads as "not meaningfully circular", which is true of
        % a direction that radiates nothing, and the title says outright
        % when the whole cut is dead.
        arDb(noFieldMask) = ARCEIL_DB;

        % Full-hemisphere grating-lobe search on the ARRAY FACTOR (not
        % the total pattern) -- grating lobes are fundamentally an
        % array-factor/periodicity phenomenon, independent of element
        % pattern, matching the "Array factor only" display option.
        % Grating-lobe search stays on the UPPER HEMISPHERE even in full-
        % sphere mode, deliberately. A planar array in the z=0 plane has
        % an array factor that depends only on (u,v) = (sin(th)cos(ph),
        % sin(th)sin(ph)), and sin(th) = sin(180-th) -- so the AF is
        % exactly mirror-symmetric about the array plane and every lobe
        % appears TWICE on a full sphere. Searching the whole sphere
        % would therefore report the main beam's own mirror image at
        % 180-theta_s as a "grating lobe" at 0 dB down: a guaranteed
        % false positive on every array. Visible space is counted once,
        % here.
        [gThA,gPhA,gDbA,gThI,gPhI,gDbI] = findGratingLobe( ...
            afDb(:,1:nThUp), theta(1:nThUp), phi, totDb(:,1:nThUp));
        % Same interior-peak test the Max scan tool uses. Without it this
        % panel labelled any lobe above -6 dB a GRATING LOBE, so a
        % half-wave lattice steered to 60 deg read "GRATING LOBE" here
        % while Max scan correctly reported none -- two readouts in one
        % window, disagreeing, and this one wrong. A grating lobe's peak
        % lies INSIDE visible space; a skirt clipped by the horizon is
        % still rising at theta = 90.
        % Prefer a STRONG interior lobe when one exists -- that is the
        % lobe the "GRATING LOBE" verdict is actually about. The falloff
        % test itself now lives in findGratingLobe and runs over every
        % candidate, so a grazing lobe that ties in level with a real
        % interior one can no longer win max() and suppress the verdict.
        % -6 dB here is the same gate the message below applies; below it
        % the label reads "highest sidelobe" either way, so falling back
        % to the overall strongest lobe leaves the grazing and sidelobe
        % readouts behaving exactly as before.
        if gDbI > -6
            S.glTheta = gThI; S.glPhi = gPhI; S.glRelDb = gDbI;
            S.glInterior = true;
        else
            S.glTheta = gThA; S.glPhi = gPhA; S.glRelDb = gDbA;
            S.glInterior = false;
        end

        % Each field's OWN P_base (true total, both-pol power -- see the
        % note below on why this must NOT be the possibly polarization-
        % projected db) and its OWN directivity, computed independently.
        % Previously only ONE field (whichever the 3D dropdown happened
        % to have selected) got its directivity computed, and that SAME
        % offset was applied to all three cut-plot curves -- meaning the
        % other two curves' "absolute gain" values were wrong by however
        % much their true directivity actually differed. Verified this
        % was not a rounding-level issue: with default settings (Total
        % selected), the Element-factor-only curve was off by 8.6 dB
        % from its own true directivity, and Array-factor-only by 1.1 dB.
        P_base_af  = abs(Eaf).^2 / ref^2;
        P_base_ef  = fieldPower(Eef_th,Eef_ph);
        P_base_tot = fieldPower(Etot_th,Etot_ph) / ref^2;

        % AF represents isotropic radiators; include both hemispheres too.
        % Its power integral is evaluated exactly below.
        Dpk_af = NaN;
        [Dpk_ef,  Prad_ef]  = fieldDirectivity(P_base_ef,  theta, phi, TH);
        [Dpk_tot, Prad_tot] = fieldDirectivity(P_base_tot, theta, phi, TH);
        % The display mesh can miss narrow array-factor lobes at large
        % apertures. Its isotropic power integral has an exact pair-sum.
        feedPhase = -S.k*phaseFreqRatio()* ...
            (S.el(:,1)*us + S.el(:,2)*vs) + deg2rad(S.el(:,4));
        if S.seqPhase, feedPhase = feedPhase + seqSgn*deg2rad(S.el(:,5)); end
        Prad_af = exactIsotropicPrad(normalizedAmp.*exp(1j*feedPhase))/ref^2;
        if strcmp(S.efType,'Isotropic'), Prad_tot = Prad_af; end
        if strcmp(S.efType,'cos^q(theta)')
            Prad_ef = 2*pi/(S.efQ+1);
            Prad_tot = exactCosQPrad(normalizedAmp.*exp(1j*feedPhase))/ref^2;
        elseif strcmp(S.efType,'Short dipole (z-axis)')
            Prad_ef = 8*pi/3;
            Prad_tot = exactZDipolePrad(normalizedAmp.*exp(1j*feedPhase))/ref^2;
        elseif ismember(S.efType,{'Dipole (linear pol)','Patch (cos^q x lin pol)'})
            Prad_tot = exactVectorPrad(normalizedAmp.*exp(1j*feedPhase), ...
                S.el(:,5))/ref^2;
            if strcmp(S.efType,'Dipole (linear pol)')
                Prad_ef = 8*pi/3;
            else
                Prad_ef = pi*(1/(S.efQ+1)+1/(S.efQ+3));
            end
        elseif needsFinePowerIntegral()
            Prad_tot = quadraturePrad(normalizedAmp.*exp(1j*feedPhase))/ref^2;
        end
        peakPower = refinedPatternPeaks({P_base_af,P_base_ef,P_base_tot},theta,phi);
        if Prad_af>0, Dpk_af = 4*pi*peakPower(1)/Prad_af; end
        if Prad_ef>0, Dpk_ef = 4*pi*peakPower(2)/Prad_ef; end
        if Prad_tot>0, Dpk_tot = 4*pi*peakPower(3)/Prad_tot; end
        S.DpkTot = Dpk_tot;   % for the Metrics readout in refreshInfo (Aperture Efficiency etc.)
        % Reuse this full-sphere normalization for the optional 2D cuts.
        % A cut-local maximum or power integral would give a different
        % answer whenever the selected plane misses the actual beam.
        S.pattern2DBasis = struct('ref',ref,'floor',totFloor, ...
            'floorAF',afFloor,'floorEF',efFloor, ...
            'Prad',Prad_tot,'PradAF',Prad_af,'PradEF',Prad_ef, ...
            'seqSign',seqSgn, ...
            'peakDb',10*log10(max(peakPower(3),realmin)), ...
            'peakDbAF',10*log10(max(peakPower(1),realmin)), ...
            'peakDbEF',10*log10(max(peakPower(2),realmin)));

        % Trim back to the displayed hemisphere when the full sphere was
        % computed purely for the integrations above. Every array the
        % display path indexes is sliced together in one place, so
        % TH/PH/U/V stay dimensionally consistent with the pattern data
        % they pair with. Done AFTER findGratingLobe (which takes its own
        % hemisphere slice) and after all Dpk/arDb work, so none of those
        % see the truncated grid.
        %
        % `theta` itself is deliberately NOT sliced here. Its last read is
        % the fieldDirectivity block just above; below this point the cut
        % section reassigns it wholesale (theta = cutX), so a slice would
        % be dead code -- MATLAB's Code Analyzer flags it as such. Left
        % out on purpose; don't add it back for symmetry with the arrays
        % below, which ARE all read again by the display path.
        if calcFull && ~S.fullSphere
            keep  = 1:nThUp;
            TH    = TH(:,keep);    PH    = PH(:,keep);
            afDb  = afDb(:,keep);  totDb = totDb(:,keep);
            efDb  = efDb(:,keep);  arDb  = arDb(:,keep);
        end

        DR = spDR.Value;   % also used below by the AR branch's view(...) reuse
        isARView = strcmp(ddShow.Value, 'Axial Ratio (dB)');
        if isARView
            % ---- 3D surface: axial ratio, NOT gain -- unit-radius sphere
            % (radius carries no meaning here; AR is a per-direction
            % QUALITY metric, not a magnitude, so overloading radius with
            % it the way gain plots do would be misleading), colour = AR
            % (dB), capped at ARCEIL_DB. Colormap used in its NATURAL
            % (un-flipped) orientation: turbo/parula already run
            % blue(low)->red(high), and low AR is GOOD here, so that
            % lines up directly with the usual blue=cool/good,
            % red=hot/bad convention without needing to invert anything
            % -- flipping it would put good CP under red, which is
            % backwards from what a reader expects at a glance.
            R = ones(size(TH));
            X = R.*sind(TH).*cosd(PH);
            Y = R.*sind(TH).*sind(PH);
            Z = R.*cosd(TH);
            prevView = get(ax3D,'View');
            cla(ax3D);
            surf(ax3D,X,Y,Z,arDb,'EdgeColor','none');
            try
                colormap(ax3D,turbo);
            catch
                colormap(ax3D,parula);
            end
            clim(ax3D,[0 ARCEIL_DB]);
            axis(ax3D,'equal'); grid(ax3D,'on');
            xlabel(ax3D,'x'); ylabel(ax3D,'y'); zlabel(ax3D,'z');
            if isequal(prevView,[0 90]) || all(prevView==0)
                view(ax3D,135,25);
            else
                view(ax3D,prevView);
            end
            cb = colorbar(ax3D); cb.Label.String = 'AR (dB, low=good CP)';
            drawSteeringPointer(ax3D);
            % Evaluate the commanded direction itself. A narrow or custom
            % polarization pattern can change substantially between the
            % two-degree azimuth samples of the display mesh.
            steerTh = abs(S.theta_s);
            steerPh = mod(S.phi_s+180*(S.theta_s<0),360);
            steerEth = 0; steerEph = 0;
            for group = 1:nUR
                members = rotIdx==group;
                [ethPoint,ephPoint] = elementFactor(steerTh,steerPh,uniqRot(group));
                phasePoint = feedPhase(members)+S.k*fr* ...
                    (S.el(members,1)*us+S.el(members,2)*vs);
                contribution = sum(normalizedAmp(members).*exp(1j*phasePoint));
                steerEth = steerEth+contribution*ethPoint;
                steerEph = steerEph+contribution*ephPoint;
            end
            steerPower = fieldPower(steerEth,steerEph)/ref^2;
            if steerPower < pRefMax*10^(-ARNULL_DB/10)
                arAtSteer = ARCEIL_DB;
            else
                steerR = abs((steerEth+1j*steerEph)/sqrt(2));
                steerL = abs((steerEth-1j*steerEph)/sqrt(2));
                ratio = abs(steerR-steerL)/max(steerR+steerL,realmin);
                arAtSteer = min(max(-20*log10(max(ratio, ...
                    10^(-ARCEIL_DB/20))),0),ARCEIL_DB);
            end
            % Axial ratio is derived from the theta/phi components, so
            % it is only meaningful when the element pattern actually
            % supplies that basis. An imported CST file stored in some
            % other orthogonal basis (Ludwig-3 Copol/Cross is the common
            % one) still yields finite, plausible-looking AR numbers
            % that mean nothing. The load-time warning covers it once,
            % but the user can dismiss that and switch to this view much
            % later -- so the plot says so itself.
            % One line, as for the gain surface; the cap moves to the
            % colour bar it limits.
            cb.Label.String = sprintf('AR (dB, capped at %g; low = good CP)',ARCEIL_DB);
            arTitle = sprintf(['Axial Ratio (Total field), unit sphere  ·  ' ...
                'steered %s  ·  AR there = %.1f dB'], ...
                angleText(S.theta_s,S.phi_s),arAtSteer);
            arCaveats = {};
            if ~impBasisOK()
                arCaveats{end+1} = 'NOT VALID: imported pattern is not in a Theta/Phi basis';
            end
            setPlotTitle(ax3D, arTitle, arCaveats);
        else
            % ---- 3D surface: selectable radius, colour = dB ----
            % ---- select which field the 3D surface / title describe ----
            switch ddShow.Value
                case 'Array factor only'
                    db = afDb;  ttl = 'Array factor (isotropic elements)';
                    Dpk = Dpk_af;   effThis = 1;
                case 'Element factor only'
                    db = efDb;  ttl = sprintf('Element factor (rot %g°)',refRot);
                    Dpk = Dpk_ef;   effThis = impEffLin();
                otherwise
                    db = totDb; ttl = 'Total pattern (EF × AF)';
                    Dpk = Dpk_tot;  effThis = impEffLin();
            end
            caveats3 = {};   % muted lines under the title (setPlotTitle)

            % absolute gain option, for the 3D surface only: shift the
            % SELECTED field's curve up by its own peak directivity, so the
            % 3D display reads real dBi instead of "dB below its own peak".
            % (The 2D cut-plot curves below get their OWN per-field offsets,
            % computed from Dpk_af/Dpk_ef/Dpk_tot above -- not this gOff.)
            % effThis folds in the imported pattern's own efficiency, so
            % the absolute readout is real GAIN (what the antenna
            % delivers) rather than directivity (what its shape alone
            % would give if lossless). Equals 1 for every analytic
            % element factor, leaving those unchanged.
            % Referenced to the pattern's OWN PEAK, not to sum|amp|.
            % db is 20*log10(|E|/ref) with ref = sum|amp|, which peaks at
            % 0 dB only when every element adds coherently AND the
            % element pattern is at its maximum there -- i.e. at
            % broadside. Steer a directive element and the peak of the
            % normalised curve drops by the element rolloff, so adding
            % 10log10(Dpk) on top put the whole absolute display that far
            % low. Measured on 8x8, 0.5 lambda, cos^1.5 steered to 60
            % deg: the cut peaked at 16.67 dBi against a title (and an
            % independent 4*pi*U/int-U calculation) of 20.77 / 20.35.
            % Subtracting the curve's own peak makes the displayed peak
            % equal Dpk*eff by construction, which is what the title
            % states, and changes nothing at broadside where that peak is
            % already 0 dB.
            switch ddShow.Value
                case 'Array factor only',   PradThis = Prad_af;
                case 'Element factor only', PradThis = Prad_ef;
                otherwise,                  PradThis = Prad_tot;
            end
            if cbAbs.Value
                gOff = safeOff(effThis, PradThis);
                if ~isfinite(gOff)
                    % Nothing is radiating, so there is no absolute level
                    % to shift to. Stay peak-relative and say why under
                    % the title rather than drawing a floor and calling it
                    % dBi.
                    gOff = 0;
                    caveats3{end+1} = 'NO RADIATED POWER: absolute level unavailable';
                end
            else
                gOff = 0;
            end
            % Unit-cell mode has to use its OWN absolute basis here too.
            % Dpk comes from the pattern integral, and for a periodic
            % unit cell that integral does not equal the radiated power
            % -- it is the exact quantity the unit-cell path exists to
            % avoid, and impEffLin() is derived from the same bad
            % integral. Left alone, ticking "Show absolute level" put a
            % peak on the 3D surface that disagreed with the array gain
            % printed in the title directly above it. Measured on a real
            % 0.5 lambda cell: the integral basis gave 7.54 dBi where the
            % file's own realized gain was 4.63.
            %
            % The array factor is deliberately NOT rebased: it is a
            % dimensionless interference pattern with no absolute level
            % to get wrong.
            if cbAbs.Value
                ucAbs = unitCellFigures();
                if ~isempty(ucAbs)
                    switch ddShow.Value
                        case 'Array factor only'
                            % left as-is, see above
                        % Rebased on the curve's own peak, as above.
                        % rg and gArr are absolute levels; db is
                        % peak-relative only at broadside. gArr already
                        % carries the element rolloff toward the steer
                        % angle, so adding it to a curve that carries the
                        % same rolloff counted it twice.
                        % rg is the embedded element's PEAK realized
                        % gain and efDb is unit-peak, so rg alone is the
                        % right offset there. For the array, the plotted
                        % field is sum(w_n E_n)/ref while the physical
                        % normalisation is per unit input power, so the
                        % offset is rg + 10log10(ref^2/sum|w|^2).
                        %
                        % NOT gArr - max(db): gArr is the gain toward the
                        % STEER direction, and pinning the pattern's
                        % global maximum to it drags the whole curve down
                        % whenever the beam is not at the steer angle --
                        % a 0/180 pair cancels at broadside, and that
                        % near-zero directional value was being applied
                        % to every direction.
                        case 'Element factor only'
                            gOff = ucAbs.rg;
                        otherwise
                            pInUC = sum(normalizedAmp.^2);
                            if pInUC > 0
                                gOff = ucAbs.rg + 10*log10(ref^2/pInUC);
                            else
                                gOff = ucAbs.rg;
                            end
                    end
                end
            end
            % Note: afDb/totDb/efDb were only needed to pick `db` above; they
            % are not read again anywhere after this point (only `db` feeds
            % the 3D surface), so only db needs the offset applied.
            db = db + gOff;

            % A dB radius makes sidelobes visible, but it nearly erases
            % shallow element-pattern waists: a crossed dipole is only
            % 3 dB down at the equator, so a 40 dB radius window makes
            % its equatorial radius 0.925 (almost a sphere). Linear power
            % gives its physical normalized power radius of 0.5. Keep
            % both older scales available for inspecting weak lobes.
            % Framed on the PLOTTED peak, not on gOff. gOff used to be
            % the peak level itself; it is now a power-integral offset,
            % and the two are equal only when the pattern happens to peak
            % at its own directivity. Using it as the window top clipped
            % the display -- a single dual-feed patch at q = 1.5 peaks
            % near 8.08 dBi against a gOff of about 7.07, so the top
            % decibel of the beam fell outside the axis.
            %
            % Computed BEFORE the branch and from db alone: db already
            % carries gOff (added just above), and clim below needs this
            % whichever radius scaling is selected. Defining it only in
            % the dB branch left it undefined for "Linear magnitude",
            % which is what the fuzz section caught.
            dbTop = max(db(:));
            switch ddScale.Value
                case 'Linear power'
                    R = 10.^((db-dbTop)/10);
                case 'Linear magnitude'
                    R = 10.^((db-dbTop)/20);
                otherwise
                    R = max(db - (dbTop-DR), 0)/DR;   % top -> 1, DR below -> 0
            end

            X = R.*sind(TH).*cosd(PH);
            Y = R.*sind(TH).*sind(PH);
            Z = R.*cosd(TH);
            prevView = get(ax3D,'View');       % keep the user's orientation
            cla(ax3D);
            surf(ax3D,X,Y,Z,db,'EdgeColor','none');
            try
                colormap(ax3D,turbo);
            catch
                colormap(ax3D,parula);   % turbo introduced in R2019b; parula is the older default
            end
            % A shallow pattern should use its actual colour span. With
            % the old fixed -40..0 dB colours, a crossed dipole's entire
            % 3 dB variation was rendered red/orange even after its
            % silhouette was fixed. Arrays with deep nulls still use the
            % selected DR floor. Isotropic needs a nonzero colour span.
            colourFloor = max(dbTop-DR,min(db(:)));
            if dbTop-colourFloor < 1e-6, colourFloor = dbTop-1; end
            clim(ax3D,[colourFloor dbTop]);
            % Same validity rule as the cut axis: an absolute unit only
            % goes on when there is radiated power to be absolute about.
            % Computed ONCE and shared by the colorbar unit and the title
            % tag below. They used to be decided separately -- the
            % colorbar tested the offset, the title tested only the
            % checkbox -- so a pattern that radiates nothing came out
            % with a plain "dB" colorbar under an "[absolute
            % directivity]" title. Every absolute-level annotation on
            % this axes now reads the same flag, which is the only way
            % they cannot drift apart again.
            absOK = cbAbs.Value && isfinite(safeOff(effThis, PradThis));
            cb3lbl = 'dB';
            if absOK
                cb3lbl = 'dBi';
            end
            axis(ax3D,'equal'); grid(ax3D,'on');
            xlabel(ax3D,'x'); ylabel(ax3D,'y'); zlabel(ax3D,'z');
            if isequal(prevView,[0 90]) || all(prevView==0)
                view(ax3D,135,25);             % first draw: sensible default
            else
                view(ax3D,prevView);           % otherwise keep what the user set
            end
            cb = colorbar(ax3D); cb.Label.String = cb3lbl;
            drawSteeringPointer(ax3D);
            if absOK
                % AF-only is always a dimensionless interference pattern;
                % its absolute dBi quantity is directivity even when the
                % selected ELEMENT came from an imported realized-gain file.
                % EF/Total for an import are realized gain because their
                % offsets include the imported element's efficiency / peak.
                if strcmp(ddShow.Value,'Array factor only')
                    absTag = ', absolute AF directivity';
                elseif strcmp(S.efType,'Imported (CST far-field)') || effThis < EFF_TOL
                    absTag = ', absolute gain (realized)';
                else
                    absTag = ', absolute directivity';
                end
            else
                % Blank rather than "normalized to own peak": that is the
                % default state, and a tag announcing the default is noise.
                % The absolute tags above stay, since those DO mark a
                % non-default reading. This branch is also where a TICKED
                % checkbox lands when there is no radiated power to be
                % absolute about; the caveat line says so instead.
                absTag = '';
            end
            % The title is ONE line: what is drawn and where the beam is
            % steered. The array's own figures (directivity or gain, HPBW,
            % SLL) are on the result cards and in Details -- the old second
            % line repeated them and, with the frequency line, pushed the
            % title into the card row at laptop size. Only the AF-only and
            % EF-only views keep a figure here, because no card shows the
            % directivity of a single factor.
            fig3 = '';
            if ~strcmp(ddShow.Value,'Total (EF x AF)')
                if isnan(Dpk) || Dpk <= 0
                    fig3 = 'D unavailable';
                elseif effThis < EFF_TOL
                    % Both numbers side by side: the imported file's
                    % efficiency is the whole reason they differ, and hiding
                    % either one is how a 2 dB discrepancy against CST turns
                    % into a mystery.
                    fig3 = sprintf('D = %.2f dBi, G = %.2f dBi', ...
                        10*log10(Dpk), 10*log10(Dpk*effThis));
                else
                    fig3 = sprintf('D = %.2f dBi', 10*log10(Dpk));
                end
                % In unit-cell mode the pattern-shape directivity of the
                % embedded element is not the directivity anyone means (the
                % lattice caps its gain at 4*pi*A_cell/lambda^2), so the
                % element view states the file's realized gain and the D
                % that follows from the entered efficiency. The AF stays
                % dimensionless, so its integrated directivity is kept.
                ucT = unitCellFigures();
                if ~isempty(ucT) && strcmp(ddShow.Value,'Element factor only')
                    fig3 = sprintf('unit cell RG = %.2f dBi, D = %.2f dBi', ucT.rg, ucT.d);
                end
            end
            % Polarization named in the title: the surface is polCombine's
            % output, so an RHCP/LHCP readout shows a genuinely different
            % pattern from Total on any CP array, not a rescaled one.
            pTag3 = polTag();
            if ~isempty(pTag3), ttl = [ttl ', ' pTag3 ' component']; end
            head3 = [ttl absTag];
            if ~isempty(fig3), head3 = [head3 '  ·  ' fig3]; end
            if size(S.el,1)==1
                head3 = [head3 '  ·  single element'];
                if abs(S.theta_s)>1e-9
                    caveats3{end+1} = 'ONE ELEMENT: steering command cannot move its pattern';
                end
            else
                head3 = [head3 '  ·  steered ' angleText(S.theta_s,S.phi_s)];
            end
            % Caveats, one per line under the title. The substituted
            % element and the missing import stay on the plot for as long
            % as they are what is drawn; the dialogs fire only once.
            if S.efFallback
                caveats3{end+1} = 'CUSTOM FORMULA FAILED - element substituted';
            end
            if S.impNeedsPattern
                caveats3{end+1} = 'LOADED DESIGN NEEDS AN IMPORTED PATTERN - showing Isotropic';
            end
            if magOnlyMixedRot()
                caveats3{end+1} = ['NOT VALID: magnitude-only import ' ...
                    'with mixed element rotations - polarisation is invented'];
            end
            % Frequency context only when the operating frequency differs
            % from the design one: at the design frequency, retuned and
            % frozen phases are the same thing. Under frozen phases the
            % steering-term prediction is reported, not a measured peak:
            % manual/rotation phases, the element factor, the taper or
            % competing lobes can move the actual maximum.
            fDiffers = abs(S.freqOpGHz - S.freqGHz) > ...
                1e-12*max([1 abs(S.freqGHz) abs(S.freqOpGHz)]);
            if S.squintMode && fDiffers
                fTxt = sprintf('At %.4g GHz (design %.4g GHz)', S.freqOpGHz, S.freqGHz);
                if ~S.retunePhase && size(S.el,1) > 1
                    sB = squintSinTheta(S.theta_s);
                    if abs(sB) <= 1
                        fTxt = sprintf('%s; phases fixed: steering term squints to %s', ...
                            fTxt, angleText(asind(sB),S.phi_s));
                    else
                        fTxt = sprintf(['%s; phases fixed: steering-term ' ...
                            'prediction OUTSIDE VISIBLE SPACE'], fTxt);
                    end
                elseif ~S.retunePhase
                    fTxt = sprintf('%s; phases fixed (single element: no AF squint)', fTxt);
                end
                caveats3{end+1} = fTxt;
            end
            setPlotTitle(ax3D, head3, caveats3);
        end

        % ---- principal cut: Theta sweep (fixed phi_s) or Phi sweep
        % (fixed theta), user-selected via ddCutMode ----
        isPhiCut = strcmp(S.cutMode, 'Phi cut (fixed theta)');
        if isPhiCut
            % azimuth cut: phi sweeps the full circle, theta held fixed.
            % No "negative side" flip needed (phi already wraps 0-360
            % naturally), unlike the theta-cut branch below.
            cutX = 0:0.5:360;
            thFix = S.cutFixedTheta;
            uc = sind(thFix)*cosd(cutX);
            vc = sind(thFix)*sind(cutX);
            thcForEF = thFix*ones(size(cutX));
            phcForEF = cutX;
        else
            % elevation cut: theta sweeps -90..90, phi_s fixed. Negative
            % theta automatically traverses the phi_s+180 side, so the
            % beam is shown whole even when it sits at broadside.
            % In full-sphere mode this extends to -180..180, closing the
            % great circle: the abs()/phi-flip mapping below stays exact
            % past 90 deg (cutX=150 -> theta=150 at phi_s; cutX=-150 ->
            % theta=150 at phi_s+180), and sind(cutX) still gives the
            % correct signed u,v there.
            if S.fullSphere, cutX = -180:0.1:180; else, cutX = -90:0.1:90; end
            % cutPhiVal(), not S.phi_s: the cut plane is independent of
            % where the beam is steered unless the user links them.
            cutPh = cutPhiVal();
            uc  = sind(cutX)*cosd(cutPh);
            vc  = sind(cutX)*sind(cutPh);
            thcForEF = abs(cutX);
            phcForEF = cutPh + 180*(cutX < 0);   % phi flips for negative theta
        end

        Ecut_th = zeros(size(cutX));  Ecut_ph = zeros(size(cutX));
        Acut = zeros(size(cutX));
        % Same rotation-keyed cache as the main grid loop above, for the
        % same reason: the element pattern depends only on rotation, so
        % evaluating it per ELEMENT repeats identical work for every
        % element sharing one. This loop was doing that -- on an imported
        % pattern each call is four scatteredInterpolant evaluations over
        % the whole cut, so a 1000-element array with four distinct
        % rotations did 250x more work than needed. Guarded on the same
        % two conditions (fewer distinct rotations than elements, and a
        % byte budget) so a design with a unique rotation per element
        % cannot turn this into an allocation problem; the cut is 1D, so
        % the budget is reached far later than on the 2D grid.
        bytesPerRotC = 2*numel(cutX)*16;
        useEFCacheC  = (nUR < numel(rotAll)) && (nUR*bytesPerRotC <= 256e6);
        if useEFCacheC
            EthCacheC = cell(nUR,1); EphCacheC = cell(nUR,1);
            for iu = 1:nUR
                [EthCacheC{iu}, EphCacheC{iu}] = elementFactor(thcForEF,phcForEF,uniqRot(iu));
            end
        end
        for n = 1:size(S.el,1)
            xn = S.el(n,1); yn = S.el(n,2);
            amp = normalizedAmp(n); ph0 = S.el(n,4); rot = S.el(n,5);
            phFeed = -S.k*phaseFreqRatio()*(xn*us + yn*vs) + deg2rad(ph0);
            if S.seqPhase, phFeed = phFeed + seqSgn*deg2rad(rot); end
            st = amp*exp(1j*phFeed)*exp(1j*S.k*fr*(xn*uc + yn*vc));
            Acut = Acut + st;
            if useEFCacheC
                Eth = EthCacheC{rotIdx(n)};  Eph = EphCacheC{rotIdx(n)};
            else
                [Eth, Eph] = elementFactor(thcForEF,phcForEF,rot);
            end
            Ecut_th = Ecut_th + st.*Eth;
            Ecut_ph = Ecut_ph + st.*Eph;
        end
        Ecut = polCombine(Ecut_th, Ecut_ph);
        [EFcut_th, EFcut_ph] = elementFactor(thcForEF,phcForEF,refRot);
        EFcut = polCombine(EFcut_th, EFcut_ph);

        % Per-field offsets, reusing Dpk_af/Dpk_ef/Dpk_tot computed above
        % from the full 2D pattern -- a field's directivity doesn't
        % depend on which view (3D surface vs. this 1D cut) you're
        % reading it through, so the same three values apply here too.
        % This is the fix for the bug where all three curves previously
        % shared whichever ONE offset the 3D dropdown happened to have
        % selected.
        % Declared before the branch so every path through the cut
        % drawing has it defined, including the axial-ratio branch which
        % skips the offset block entirely.
        noPower = false;
        if cbAbs.Value
            % Same efficiency treatment as the 3D surface above: EF and
            % Total become real gain for an imported pattern, AF stays
            % pure directivity (no physical antenna in an isotropic-
            % element array factor). eLin == 1 for analytic elements.
            % Each curve is rebased on the peak of its OWN full-sphere
            % grid, for the reason given at the 3D surface above: the
            % normalised curves are referenced to sum|amp|, which is
            % their peak only at broadside. The peaks come from the 3D
            % grids (afDb/efDb/totDb), not from the cut, because the cut
            % plane need not contain the global peak -- taking the cut's
            % own maximum would re-normalise every off-plane cut to
            % itself and quietly hide the loss.
            % Same power-integral basis as the 3D surface above, and
            % for the same reasons -- see the Prad note beside the
            % fieldDirectivity calls.
            eLin = impEffLin();
            gOff_af  = safeOff(1,    Prad_af);
            gOff_ef  = safeOff(eLin, Prad_ef);
            gOff_tot = safeOff(eLin, Prad_tot);
            % safeOff returns NaN when nothing is radiating. A NaN
            % offset would wipe the whole curve, so fall back to a
            % peak-relative plot; the title carries the explanation.
            noPower = ~isfinite(gOff_tot);
            if ~isfinite(gOff_af),  gOff_af  = 0; end
            if ~isfinite(gOff_ef),  gOff_ef  = 0; end
            if ~isfinite(gOff_tot), gOff_tot = 0; end
            % Unit-cell mode rebases these the same way the 3D surface
            % is rebased above, and for the same reason: Dpk_* and
            % impEffLin() both come from the pattern integral, which a
            % periodic unit cell invalidates. Fixing only the 3D view
            % left the cut plot 0.63 dB away from the array gain in the
            % title -- one number for the same array, read two ways.
            % AF keeps its own directivity: it is dimensionless.
            ucCut = unitCellFigures();
            if ~isempty(ucCut)
                gOff_ef = ucCut.rg;
                pInUC2  = sum(normalizedAmp.^2);
                if pInUC2 > 0
                    gOff_tot = ucCut.rg + 10*log10(ref^2/pInUC2);
                else
                    gOff_tot = ucCut.rg;
                end
            end
        else
            gOff_af = 0; gOff_ef = 0; gOff_tot = 0;
        end
        % The SAME three full-grid floors the 3D block computed -- not a
        % floor derived from this cut.
        %
        % A cut-local floor breaks on a cut that is entirely dead, which
        % says nothing about whether the antenna radiates elsewhere. The
        % fallback then handed back a flat 1e-4 (-80 dB) and the absolute
        % offset -- correctly computed from the whole pattern, and large
        % precisely because the field is small -- was added on top, so a
        % 1e-7-scaled cos(theta) element whose real peak is 7.78 dBi drew
        % its dead theta = 90 cut at about +67.78 dBi. Strong radiation
        % invented along a line that radiates nothing.
        %
        % Against the full-grid reference the same cut lands 80 dB below
        % the pattern's true peak, which is what "at or below the display
        % floor" should look like, and the title below says so outright.
        afCut  = 20*log10(max(abs(Acut)/ref, afFloor))  + gOff_af;
        totCut = 20*log10(max(abs(Ecut)/ref, totFloor)) + gOff_tot;
        efCut  = 20*log10(max(abs(EFcut),    efFloor))  + gOff_ef;
        % TWO different states, which the first version conflated:
        %
        %   cutAllFloored  every plotted sample sits at the display floor.
        %                  The curve is the floor, not data -- but the
        %                  field may be perfectly real and merely weak.
        %   cutTrulyDead   the field is identically zero. Nothing radiates.
        %
        % Calling the first one "NO RADIATION" was wrong and reachable:
        % the LHCP readout of a near-pure-RHCP element (E_theta = 1,
        % E_phi = -1j*(1-1e-6)) is the 1e-6 cross-pol, a genuine -77 dB
        % component, and it was labelled as no radiation at all. A weak
        % cross-pol is exactly what someone selects that readout to see.
        cutAllFloored = ~any(abs(Ecut(:))/ref > totFloor);
        cutTrulyDead  = ~any(abs(Ecut(:)) > 0);
        % The plotted ceiling, computed once as soon as the curve exists
        % and used for the axis limits, the peak line's baseline and the
        % grating-warning placement alike. All three previously used
        % gOff_tot, which stopped being a level when it became a
        % power-integral offset -- the baseline and the warning then sat
        % at an arbitrary height relative to the curve they annotate.
        % Defined BEFORE the axial-ratio branch, which reaches the axis
        % limits without passing through the gain-curve path.
        cutTop = max(totCut(isfinite(totCut)));
        if isempty(cutTop), cutTop = gOff_tot; end

        % What the cut below shows, for the pinned reference (refMismatch).
        % abs follows noPower: a request for dBi with nothing radiating
        % draws relative levels.
        if isPhiCut, cutPlaneKey = S.cutFixedTheta; else, cutPlaneKey = cutPhiVal(); end
        cutKind = 'gain'; if cbARCut.Value, cutKind = 'ar'; end
        S.cutView = struct('kind',cutKind,'isPhi',isPhiCut,'plane',cutPlaneKey, ...
            'conv',S.angleConvention,'abs',cbAbs.Value && ~noPower, ...
            'pol',ddPol.Value,'full',S.fullSphere, ...
            'gain',strcmp(S.efType,'Imported (CST far-field)'));
        % cla keeps the legend: without the explicit 'off' the axial-ratio
        % view below inherited an empty legend box from the gain curves.
        cla(axCut); legend(axCut,'off');
        P = pal();   % every coloured object below carries a role Tag (recolorPlots)
        if cbARCut.Value
            % ---- Axial ratio along this cut, instead of the usual
            % AF/EF/Total gain curves -- kept as a separate branch (not
            % overlaid on the same axes as the gain curves) because AR
            % (0 to ARCEIL_DB) and gain (peak-relative or absolute dBi)
            % don't share a sensible y-axis scale. HPBW/SLL/FNBW/grating-
            % lobe analysis is gain-pattern-specific and doesn't apply
            % here, so it's skipped entirely in this mode.
            E_R_cut = (Ecut_th + 1j*Ecut_ph)/sqrt(2);
            E_L_cut = (Ecut_th - 1j*Ecut_ph)/sqrt(2);
            magSumCut = abs(E_R_cut) + abs(E_L_cut);
            magDifCut = abs(abs(E_R_cut) - abs(E_L_cut));
            % Same scale-invariant ratio form as the 3D grid above.
            rCut = magDifCut ./ max(magSumCut, realmin);
            arCutDb = min(max(-20*log10(max(rCut, 10^(-ARCEIL_DB/20))), 0), ARCEIL_DB);
            % Same gain-relative null mask as the 3D grid, and the SAME
            % fix: uses a TRUE both-polarization power reference built
            % directly from Ecut_th/Ecut_ph, not totCut -- totCut goes
            % through polCombine (via Ecut), which depends on the
            % Polarization readout dropdown. Reusing totCut here would
            % falsely mask real, meaningful AR data anywhere the field is
            % purely the polarization sense OPPOSITE whatever RHCP/LHCP
            % component the dropdown happens to have selected.
            % See the 3D block: an identically zero field has no
            % polarisation state, and the relative null test cannot see
            % it because every sample sits on the floor together.
            pTrueCut = fieldPower(Ecut_th,Ecut_ph)/ref^2;
            % Scale-relative, exactly as the 3D mask -- and against the
            % FULL pattern's maximum (pRefMax, computed above), not this
            % cut's own. The cut that most needs this test is the one
            % that is dead along its whole length, and normalising that
            % by its own maximum would compare zero against zero.
            if pRefMax > 0
                noFieldCut = pTrueCut <= pRefMax * 1e-20;
            else
                noFieldCut = true(size(pTrueCut));
            end
            % Power ratio against this CUT's own maximum -- unfloored, for
            % the reason given at the 3D mask. Kept cut-local rather than
            % switched to pRefMax so the mask still means "far below the
            % strongest thing on this cut", which is what it has always
            % meant; only the floor is gone.
            cutRefMax = max(pTrueCut);
            if cutRefMax > 0
                nullMaskCut = pTrueCut < cutRefMax * 10^(-ARNULL_DB/10);
            else
                nullMaskCut = false(size(pTrueCut));   % noFieldCut owns this case
            end
            arCutDb(nullMaskCut) = ARCEIL_DB;
            arCutDb = min(arCutDb, ARCEIL_DB);
            arCutDb(noFieldCut) = ARCEIL_DB;   % see the 3D block

            plot(axCut,cutX,arCutDb,'LineWidth',2.0,'Color',P.traceAR, ...
                'Tag','cutAR','HitTest','off');
            hold(axCut,'on');
            yline(axCut,3,'--','Color',P.muted,'LineWidth',1,'Tag','cutARLimit', ...
                'Label','3 dB (common "good CP" threshold)','LabelHorizontalAlignment','left');
            hold(axCut,'off');
            grid(axCut,'on'); ylim(axCut,[0 ARCEIL_DB]);
            if isPhiCut
                xlim(axCut,[0 360]);
            else
                xlim(axCut,[min(cutX) max(cutX)]);
            end
            labelCutAngleAxis(axCut,isPhiCut,S.cutFixedTheta);
            ylabel(axCut,'Axial ratio (dB)');
            S.cutX = cutX; S.cutDb = arCutDb;   % click-to-inspect probe reads whatever's plotted
            [arMin,arMinIdx] = min(arCutDb);
            arCutTitle = ['Axial ratio, ' cutTitleText(isPhiCut,S.cutFixedTheta)];
            arCaveats = {};
            if all(noFieldCut)
                % Nothing radiates anywhere on this cut, so there is no
                % axial ratio to report -- saying "min = 0.00 dB" would
                % read as perfect circular polarisation.
                arCaveats{end+1} = 'NO FIELD ON THIS CUT';
            else
                arCutTitle = sprintf('%s  ·  min %.2f dB at %s', ...
                    arCutTitle, arMin, cutDirectionText(cutX(arMinIdx)));
            end
            if ~impBasisOK()
                arCaveats{end+1} = 'NOT VALID: imported pattern is not in a Theta/Phi basis';
            end
            refNote = refCaveat();
            if ~isempty(refNote), arCaveats{end+1} = refNote; end
            setPlotTitle(axCut, arCutTitle, arCaveats);
            % No gain curve, so no beamwidth or sidelobes: the cards say
            % why they are blank rather than looking out of date.
            S.cutMetrics.why = 'axial-ratio cut';
        else
        hAF = plot(axCut,cutX,afCut, 'LineWidth',1.8,'Color',P.traceAF, ...
            'LineStyle','--','DisplayName','Array factor','Tag','cutAF', ...
            'HitTest','off');
        hold(axCut,'on');
        hEF = plot(axCut,cutX,efCut, 'LineWidth',1.4,'Color',P.traceEF, ...
            'DisplayName','Element factor','Tag','cutEF','HitTest','off');
        hTot = plot(axCut,cutX,totCut,'LineWidth',2.0,'Color',P.traceTotal, ...
            'DisplayName','Total','Tag','cutTotal','HitTest','off');
        hRef = plotRefTrace();   % the pinned reference, when it applies here
        % Beside the plot, without a box: inside the axes any corner
        % covers sidelobes of some pattern, and above it the legend sat
        % over the title, taking height this short plot cannot spare.
        legend(axCut,[hAF hEF hTot hRef],'Location','eastoutside','Box','off', ...
            'AutoUpdate','off');
        % framed around the Total curve's OWN offset (gOff_tot), not the
        % 3D-surface-specific gOff, so the axis centers correctly on the
        % main curve regardless of what the 3D dropdown happens to show
        % Framed on the Total curve's own plotted peak, for the same
        % reason as the 3D window above: gOff_tot is an offset now, not a
        % level, so using it as the ceiling clipped the beam.
        % The headroom above the curves is a clear band for the peak and
        % strong-lobe labels: no trace reaches it, so they never sit on a
        % curve (and, inside the axes, never on the title). It starts
        % above the HIGHEST curve, not the Total peak: steer a directive
        % element and the normalised AF and EF curves stay at 0 dB while
        % Total drops by the element rolloff (or a pinned cut misses the
        % beam). Capped at half the range so an AF far above an absolute
        % Total does not squash the plot. A pinned reference trace counts
        % as a curve: pinned on a larger array, it stands above this one.
        curvesTop = max([cutTop; afCut(isfinite(afCut))'; efCut(isfinite(efCut))']);
        if ~isempty(hRef), curvesTop = max([curvesTop hRef.YData(isfinite(hRef.YData))]); end
        curvesTop = min(curvesTop, cutTop + DR/2);
        grid(axCut,'on'); ylim(axCut,[cutTop-DR curvesTop+max(2,0.14*DR)]);
        if isPhiCut
            xlim(axCut,[0 360]);
        else
            % Follows cutX rather than a hardcoded +-90: in full-sphere
            % mode the cut is computed out to +-180, and a fixed +-90
            % window would silently hide half the curve it just built.
            xlim(axCut,[min(cutX) max(cutX)]);
        end
        labelCutAngleAxis(axCut,isPhiCut,S.cutFixedTheta);
        % Names whichever quantity is actually plotted. gOff_ef/gOff_tot
        % above fold in the imported pattern's efficiency, so for such a
        % file these curves are GAIN; for every analytic element factor
        % (efficiency 1) they remain directivity. Calls impEffLin()
        % rather than reusing eLin, which only exists on the cbAbs
        % branch above.
        ylabel(axCut,'dB');
        % noPower, not cbAbs alone. Asking for absolute levels does not
        % make them exist: with nothing radiating the offset was replaced
        % by zero and the floored curve kept an absolute axis label, so
        % the plot still claimed dBi for a field that is not there. The
        % title, the axis label and the colorbar now all read the same
        % flag.
        if cbAbs.Value && ~noPower
            % Same correction as the 3D title: in unit-cell mode these
            % curves are offset by the file's realized GAIN, so they are
            % gain -- but impEffLin() returns 1 there by design, which
            % sent this to the "directivity" branch and mislabelled the
            % axis of a gain plot.
            if strcmp(S.efType,'Imported (CST far-field)')
                % Three unlike curves share this axes: AF keeps its own
                % directivity basis, while imported EF/Total use realized
                % gain. Calling the whole axis "realized gain" mislabeled AF.
                ylabel(axCut,'dBi (AF=directivity; imported EF/Total=realized gain)');
            else
                ylabel(axCut,'dBi (absolute directivity)');
            end
        end

        % Metrics are measured on the UNFLOORED curve. totCut carries the
        % display floor (1e-4 of the pattern's own peak), which flattens
        % every null to one level -- and the null structure is exactly
        % what FNBW and the sidelobe search walk. Beamwidth and SLL sit
        % far above the floor so they were already right, but measuring
        % them on clipped data was only accidentally correct: a design
        % whose features approach the floor would have been measured off
        % the floor itself. Clip for drawing, measure the real thing.
        %
        % -Inf where the field is exactly zero is fine here: every use is
        % a comparison or a local-extremum walk, and an exact zero IS the
        % deepest possible null.
        theta = cutX;
        cutDb = 20*log10(abs(Ecut)/ref) + gOff_tot;
        S.cutRaw = cutDb;
        S.cutTheta = thcForEF; S.cutPhiSamples = phcForEF;
        S.cutX = cutX; S.cutDb = totCut;   % the probe reads what is PLOTTED

        % Phi cuts and FULL-SPHERE theta cuts are periodic. The latter
        % join at -180/+180 (the same rearward direction). Only an upper-
        % hemisphere theta cut has two distinct, finite endpoints.
        periodicCut = isPhiCut || S.fullSphere;
        if ~any(isfinite(cutDb))
            hp = NaN; fnbw = NaN; sllRel = NaN;
            hold(axCut,'off');
        else
        [pk,ipk] = max(cutDb);
        nC = numel(cutDb);
        if periodicCut, nU = nC - 1; else, nU = nC; end
        if ipk > nU, ipk = 1; end

        % ---- main lobe: walk downhill from the peak to the first null
        % Each walk reports whether it stopped at a real turning point or
        % simply ran out of samples. Running off the end is NOT a null:
        % a single isotropic element falls monotonically to the horizon
        % and used to be credited with a first-null beamwidth despite
        % having no nulls at all.
        [rn, nStepR, whyR] = walkDown(cutDb, ipk, +1, nU, periodicCut);
        [ln, nStepL, whyL] = walkDown(cutDb, ipk, -1, nU, periodicCut);
        % A walk that closed the circle means the whole cut is one
        % downhill run from the peak -- there is no second lobe anywhere,
        % so there is nothing to call a sidelobe.
        fullCircle = strcmp(whyR,'circle') || strcmp(whyL,'circle') || ...
                     (nStepR + nStepL + 1 >= nU);

        % ---- HPBW, by interpolation on the -3 dB crossings
        hp = NaN;
        [raOK, ra] = crossOut(cutDb, theta, ipk, +1, nU, periodicCut, pk-3);
        [laOK, la] = crossOut(cutDb, theta, ipk, -1, nU, periodicCut, pk-3);
        if raOK && laOK
            if periodicCut
                hp = mod(ra - la, 360);
            else
                hp = ra - la;
            end
            if periodicCut && ra < la
                % The half-power span wraps the 0/360 seam. Drawn as one
                % segment it would run the LONG way round the plot,
                % crossing the whole pattern instead of marking the beam.
                plot(axCut,[la theta(end)],[pk-3 pk-3],'ro-','LineWidth',1.4,'HitTest','off');
                plot(axCut,[theta(1) ra], [pk-3 pk-3],'ro-','LineWidth',1.4,'HitTest','off');
            else
                plot(axCut,[la ra],[pk-3 pk-3],'ro-','LineWidth',1.4,'HitTest','off');
            end
        end
        % A dotted line through the peak, labelled with the exact angle
        % it occurs at -- answers "where is the main beam pointing on this
        % cut", directly and precisely, rather than having to eyeball it
        % off the curve. An xline keeps its label inside the axes, along
        % the floor (where the legend used to sit); the old label floated
        % above the peak and ran into the title once the beam was steered.
        % When the beam does not point where it was steered, say BOTH.
        % These are genuinely different quantities: the array factor
        % peaks exactly at theta_s, but multiplying by a directive
        % element pattern pulls the product's peak back toward broadside,
        % because past the beam centre the element is falling faster than
        % the array factor is rising. Measured on an 8x8 at 0.5 lambda
        % steered to 60 deg -- isotropic element 60.0 deg (no pull at
        % all), cos^1 57.0, cos^2 55.0, cos^4 52.2, cos^6 50.2 -- and the
        % same cos^2 element pulls 12.3 deg on a 4x4 but only 0.5 deg on
        % a 32x32, since a narrower array factor resists it.
        %
        % Correct physics, correctly computed, but previously silent: the
        % plot read 55 while the steering box read 60, with nothing to
        % say why, which looks exactly like a fault. Suppressed below
        % half a degree so the common case stays uncluttered, and skipped
        % for a phi cut, where the x axis is azimuth and theta_s is not
        % the quantity being plotted.
        cutSteer=signedSteerOnCut(S.theta_s,S.phi_s,cutPhiVal());
        pkLbl = cutPointText(theta(ipk));
        if size(S.el,1)>1 && ~isPhiCut && isfinite(cutSteer) && ...
                abs(theta(ipk)-cutSteer)>=0.5
            pkLbl = sprintf('%s (steered %s)', pkLbl, cutPointText(cutSteer));
        end
        pkMark = xline(axCut, theta(ipk), ':', ['Peak ' pkLbl], ...
            'Color',P.muted, 'LineWidth',1.1, 'FontSize',9, 'FontWeight','bold', ...
            'LabelVerticalAlignment','top', 'LabelOrientation','horizontal', ...
            'LabelHorizontalAlignment',xlineLabelSide(axCut,theta(ipk),NaN), ...
            'Tag','peakLine', 'HitTest','off');
        hold(axCut,'off');

        % peak sidelobe: exclude the main lobe out to its FIRST NULLS,
        % not merely out to the -3 dB points (those shoulders are still
        % part of the main lobe and would be reported as the sidelobe).
        % The span is the one already walked above, so the exclusion
        % wraps with the cut and cannot leave half a seam-straddling beam
        % behind to be found as a sidelobe.
        % Built from the STEP COUNTS. Walking from ln to rn by index
        % cannot tell a zero-width interval from one that wraps the whole
        % circle -- both have ln == rn -- and the index form excluded a
        % single sample in exactly the cases where it should have
        % excluded everything.
        mask = true(1,nC);
        if fullCircle
            mask(:) = false;
        else
            i = ipk;
            mask(i) = false;
            for k = 1:nStepR
                i = i + 1; if i > nU, i = 1; end
                mask(i) = false;
            end
            i = ipk;
            for k = 1:nStepL
                i = i - 1; if i < 1, i = nU; end
                mask(i) = false;
            end
        end
        if periodicCut, mask(nC) = mask(1); end
        if any(mask), sll = max(cutDb(mask)); else, sll = NaN; end
        % Displayed as sll-pk (relative to the main lobe), not sll alone:
        % SLL is, by universal antenna-engineering convention, always
        % expressed below the main lobe peak -- CST reports it the same
        % way even alongside an absolute "Main lobe magnitude (dBi)"
        % figure. Since gOff (the absolute-gain offset) is added
        % uniformly to every point in cutDb, it cancels exactly in this
        % subtraction, so sllRel is correct whether or not "Show
        % absolute gain" is checked.
        sllRel = sll - pk;

        % first-null beamwidth: the same two null indices used to exclude
        % the main lobe above, at zero extra cost -- but only when BOTH
        % walks actually found a null. A walk that ran to the end of the
        % cut found the edge of the plot, not a null, and reporting the
        % distance between two edges as a beamwidth is meaningless.
        % Summed as two arcs THROUGH the peak, not as the difference of
        % the two endpoint angles: on a periodic cut those endpoints
        % coincide when the lobe spans the circle, and the difference
        % then reads 0 deg for a beam that fills the whole plot.
        if ~strcmp(whyR,'null') || ~strcmp(whyL,'null') || fullCircle
            fnbw = NaN;
        elseif periodicCut
            fnbw = mod(theta(rn)-theta(ipk),360) + mod(theta(ipk)-theta(ln),360);
        else
            fnbw = theta(rn) - theta(ln);
        end

        % grating-lobe warning: a secondary lobe within 6 dB of the main
        % peak, well outside the main lobe's own null-to-null span, is a
        % strong grating-lobe candidate rather than an ordinary sidelobe
        % (an ordinary uniform sidelobe never gets closer than ~13 dB).
        gratingTxt = '';
        if any(mask) && ~fullCircle
            idxOutside = find(mask);
            [glLevel, relIdx] = max(cutDb(idxOutside));
            glIdx = idxOutside(relIdx);
            if glLevel > pk - 6
                % "+ 0" turns the -0 of a lobe level with the peak into 0.
                gratingTxt = sprintf('Strong lobe %.1f dB', round(glLevel-pk,1) + 0);
            end
        end
        delete(findall(axCut,'Tag','gratingWarn'));
        if ~isempty(gratingTxt)
            % Marked where it is, with its level below the peak, in the
            % "bad" colour. It used to be a banner at the top of the axes,
            % on the top gridline and into the title. The Grating lobes
            % card and Details say whether it is a true grating lobe.
            hold(axCut,'on');
            glMark = xline(axCut, theta(glIdx), '--', gratingTxt, ...
                'Color',P.bad, 'LineWidth',1.2, 'FontSize',9, 'FontWeight','bold', ...
                'LabelVerticalAlignment','top', 'LabelOrientation','horizontal', ...
                'LabelHorizontalAlignment',xlineLabelSide(axCut,theta(glIdx),theta(ipk)), ...
                'Tag','gratingWarn', 'HitTest','off');
            hold(axCut,'off');
            % The two labels point away from each other where the axes
            % leave room; where one has to point at the other's line, the
            % lobe's label moves down to the floor so they cannot meet.
            pkMark.LabelHorizontalAlignment = ...
                xlineLabelSide(axCut,theta(ipk),theta(glIdx));
            toward = @(h,xOther) xor(strcmp(h.LabelHorizontalAlignment,'right'), ...
                xOther < h.Value);
            if toward(glMark,theta(ipk)) || toward(pkMark,theta(glIdx))
                glMark.LabelVerticalAlignment = 'bottom';
            end
        end

        end % finite radiation required for beam metrics
        % Kept for the result cards, which refreshInfo fills after this.
        % The plane names the polarization too: like the title, these
        % figures describe whichever component the readout selects.
        cutPlane = cutPlaneText(isPhiCut, S.cutFixedTheta);
        if ~isempty(polTag()), cutPlane = [cutPlane ', ' polTag()]; end
        S.cutMetrics = struct('hpbw',hp,'fnbw',fnbw,'sll',sllRel, ...
            'measured',any(isfinite(cutDb)),'plane',cutPlane, ...
            'offBeam',~isPhiCut && ...
            ~isfinite(signedSteerOnCut(S.theta_s,S.phi_s,cutPhiVal())), ...
            'why','');
        if ~S.cutMetrics.measured, S.cutMetrics.why = 'no field on cut'; end

        % The title names the cut (and the polarization, which HPBW, FNBW
        % and SLL are measured on -- see polTag for the CP case where they
        % describe cross-pol leakage). The numbers are on the cards; what
        % used to follow them in the title is kept only where it changes
        % how the plot must be read, as caveat lines under it.
        cutTtl = cutTitleText(isPhiCut, S.cutFixedTheta);
        pTag = polTag();
        if ~isempty(pTag), cutTtl = [cutTtl ', ' pTag ' component']; end
        cutCaveats = {};
        if ~isPhiCut && ~isfinite(signedSteerOnCut(S.theta_s,S.phi_s,cutPhiVal()))
            % The cut plane is NOT the steering plane: the beam peak
            % generally does not lie in this cut, so HPBW/SLL describe
            % whatever the cut intersects, not the main beam. The signed
            % elevation cut contains both azimuth half-planes, so a cut
            % 180 deg from the steer azimuth still contains the beam, and
            % broadside lies on every such plane.
            if useAzEl()
                [steerAz,~] = azEl(S.theta_s,S.phi_s);
                cutCaveats{end+1} = sprintf(['Cut plane pinned off the beam ' ...
                    '(steering azimuth %.0f°): HPBW and SLL are not the main beam''s'], steerAz);
            else
                cutCaveats{end+1} = sprintf(['Cut plane pinned off the beam ' ...
                    '(steering φ = %.0f°): HPBW and SLL are not the main beam''s'], S.phi_s);
            end
        end
        if magOnlyMixedRot()
            % Stronger than the plain mixed-rotation note: that one says
            % the EF x AF factorisation no longer holds, which is a caveat
            % about how to read the curves. This says the numbers are
            % wrong.
            cutCaveats{end+1} = 'NOT VALID: magnitude-only import with mixed rotations';
        elseif mixedRot
            cutCaveats{end+1} = 'Mixed element rotations: Total ≠ EF × AF';
        end
        % Nothing radiating means no absolute level exists: said, rather
        % than letting the axis imply dBi over a floored curve.
        if noPower
            cutCaveats{end+1} = 'ABSOLUTE LEVEL UNAVAILABLE: nothing radiates';
        end
        if cutTrulyDead
            cutCaveats{end+1} = 'NO RADIATION ON THIS CUT (field is identically zero)';
        elseif cutAllFloored
            % The curve is the floor, so HPBW/SLL are being read off a
            % flat line -- but the field itself is real, just below what
            % the plot can show. Naming the two apart matters: one says
            % "there is nothing here", the other "there is something here
            % and it is too weak to draw".
            cutCaveats{end+1} = ['ENTIRELY BELOW DISPLAY FLOOR ' ...
                '(curve is the floor; field is nonzero)'];
        end
        refNote = refCaveat();   % change vs the reference, or why it is hidden
        if ~isempty(refNote), cutCaveats{end+1} = refNote; end
        setPlotTitle(axCut, cutTtl, cutCaveats);
        end

        % The info panel reports two things ONLY this function can
        % produce -- the grating/side-lobe search (S.glTheta/glPhi/
        % glRelDb) and the Metrics block's Calculated directivity /
        % Aperture efficiency (S.DpkTot) -- so it has to be refreshed
        % after they are written, here, not before.
        %
        % Two real staleness bugs this fixes, both confirmed:
        %  - refreshAll() calls refreshInfo() BEFORE maybeCompute(), so
        %    every readout driven through it was showing the PREVIOUS
        %    recompute's grating lobe and directivity.
        %  - assignSteer() never called refreshInfo() at all, so steering
        %    the beam left the panel frozen entirely. Measured on an 8x8
        %    at d=0.8 lambda (grating lobe expected from theta_s=14.5
        %    deg): steering to 20, 40 and 60 deg all kept reporting the
        %    broadside "highest sidelobe: theta=13 deg, 12.8 dB below
        %    peak" and never once flagged the grating lobe.
        % Placed at the end of the function rather than patched into each
        % caller because every path that recomputes reaches here, and
        % refreshInfo() only reads state -- it cannot re-enter this.
        refreshInfo();
        % An open phase map is a view of this design, so it is redrawn
        % with it. Guarded on the handle still being live: the window may
        % have been closed since, and redraw() reads only S -- it cannot
        % re-enter computePattern.
        refreshPhaseMaps();

        % ---- validity indicators, FINALISED after every evaluation ----
        % The 3D title's marker is applied where that title is built,
        % which is BEFORE the cut's own element-factor evaluation. A
        % formula that fails only on the finer cut grid therefore set the
        % flag after the marker had already been decided, and the surface
        % carried no warning while the cut beside it was drawn from a
        % substituted element. Re-stamping here -- after the cut, after
        % everything -- is the only point where the flag is final.
        %
        % addCaveat adds the line only where it is not already present,
        % so a compute that failed on the SURFACE grid does not get it
        % twice.
        if S.efFallback
            addCaveat(ax3D, 'CUSTOM FORMULA FAILED - element substituted');
            addCaveat(axCut, 'CUSTOM FORMULA FAILED - element substituted');
        end

        % A Custom formula that failed on the real grid has just been
        % replaced by an isotropic element. Announced ONCE per distinct
        % formula: elementFactor runs many times per compute and the
        % compute itself reruns on every control change, so alerting
        % unconditionally would make the app unusable -- but never
        % alerting is how an 80 dB substitution stayed invisible.
        if S.efFallback
            % Names the component that actually failed. Saying "replaced
            % by isotropic" was only true when E_theta was the casualty:
            % if E_phi alone fails it is zeroed while a perfectly good
            % E_theta survives, which is a LINEAR element, not isotropic.
            if all(S.efFallbackWhich)
                whatTxt = ['Both E_theta and E_phi failed. The element is now ' ...
                    'ISOTROPIC (E_theta = 1, E_phi = 0).'];
            elseif S.efFallbackWhich(1)
                whatTxt = ['E_theta failed and was replaced by 1. Your E_phi is ' ...
                    'still in use, so the element is neither what you typed nor ' ...
                    'a standard one.'];
            else
                whatTxt = ['E_phi failed and was replaced by 0. Your E_theta is ' ...
                    'still in use, so the element has become purely LINEAR -- not ' ...
                    'isotropic, and not what you typed.'];
            end
            key = [S.customFormula '|' S.customFormulaPh];
            if ~strcmp(key, S.efFallbackWarned)
                S.efFallbackWarned = key;
                uialert(fig, ['Your Custom element formula could not be evaluated ' ...
                    'over the full theta/phi grid -- it produced Inf, NaN, or the ' ...
                    'wrong shape somewhere on it.' newline newline whatTxt newline newline ...
                    'Everything shown -- directivity, gain, beamwidth, sidelobes -- ' ...
                    'describes the SUBSTITUTED element, not the formula you typed.' ...
                    newline newline ...
                    'A common cause is a division that is singular at a grid ' ...
                    'point: 1./(th-30) is finite on the small grid the field is ' ...
                    'checked against as you type, but theta = 30 deg IS on the ' ...
                    'plotting grid.'], 'Custom formula failed');
            end
        end
    end
end

% Stateless helpers stay in this file but outside the shared UI workspace.
% This keeps MATLAB below its nested-function parser size limit as the app
% grows; these functions depend only on their explicit arguments.
function [iEnd, nSteps, stopWhy] = walkDown(db, i0, dir, n, wrap)
    %WALKDOWN  Walk downhill from i0 to the first turning point.
    %   iEnd    where it stopped
    %   nSteps  how many samples it moved
    %   stopWhy 'null'   turned back up: a genuine local minimum
    %           'edge'   ran out of samples on a finite cut
    %           'circle' went the whole way round a periodic cut
    %                    without ever turning up
    %
    %   The step COUNT matters as much as the endpoint. On a periodic
    %   cut an index alone cannot say whether an interval has zero
    %   width or spans the entire circle -- both leave the left and
    %   right walks sitting on the same sample. A flat phi cut ended
    %   with ln == rn, only one sample got excluded as "main lobe",
    %   and every remaining sample was then classified as a sidelobe:
    %   SLL 0 dB and a grating-lobe warning on a pattern with no
    %   lobes at all.
    iEnd = i0; nSteps = 0; stopWhy = 'null';
    for k = 1:n
        nx = iEnd + dir;
        if nx > n
            if wrap, nx = 1; else, stopWhy = 'edge'; return; end
        elseif nx < 1
            if wrap, nx = n; else, stopWhy = 'edge'; return; end
        end
        if db(nx) > db(iEnd), return; end     % turned back up: a null
        iEnd = nx; nSteps = nSteps + 1;
    end
    % Never turned up anywhere: a constant or monotone-round pattern.
    if wrap, stopWhy = 'circle'; else, stopWhy = 'edge'; end
end

function [ok, ang] = crossOut(db, ang0, i0, dir, n, wrap, lvl)
    %CROSSOUT  Angle at which the pattern first falls through lvl,
    %   walking outward from the peak. Interpolated between the two
    %   samples that straddle it. ok is false when the walk leaves
    %   the cut without ever crossing -- there is no half-power point
    %   in that direction, and inventing one from the end sample
    %   would report a beamwidth for a beam that is still open.
    ok = false; ang = NaN;
    i = i0;
    for k = 1:n
        nx = i + dir;
        if nx > n
            if wrap, nx = 1; else, return; end
        elseif nx < 1
            if wrap, nx = n; else, return; end
        end
        if db(nx) <= lvl
            if db(i) == db(nx), ang = ang0(nx); else
                if ~isfinite(db(i)) || ~isfinite(lvl), return; end
                % Referencing the upper sample avoids overflow; -Inf
                % is an exact zero power endpoint, not an infinite slope.
                pNext = 10.^((db(nx)-db(i))/10);
                pLevel = 10.^((lvl-db(i))/10);
                t = (pLevel-1)/(pNext-1);
                a0 = ang0(i); a1 = ang0(nx);
                if wrap
                    % interpolate the SHORT way round the circle
                    d = mod(a1 - a0 + 180, 360) - 180;
                    ang = ang0(1) + mod(a0 + t*d - ang0(1), 360);
                else
                    ang = a0 + t*(a1 - a0);
                end
            end
            ok = true; return;
        end
        i = nx;
    end
end

function [strip,buttons,closers,helpIcon] = buildTaskStrip(parent,designNames,resultNames,onSelect,onClose)
    nDesign = numel(designNames); nResult = numel(resultNames);
    names = [designNames resultNames];
    strip = uigridlayout(parent,[1 nDesign+nResult+1]);
    designWidth = 104;
    strip.ColumnWidth = [repmat({designWidth},1,nDesign) repmat({0},1,nResult) {'1x'}];
    strip.Padding = [0 0 0 0];
    strip.ColumnSpacing = 0;
    strip.BackgroundColor = [0.025 0.255 0.43];
    buttons = gobjects(numel(names),1);
    closers = gobjects(nResult,1);
    for i = 1:numel(names)
        % Every tab has the same full-height, square-cornered artwork.
        % The underlying button preserves the programmatic tab callback;
        % the image receives pointer clicks without the platform's rounded
        % uibutton border. Result tabs add a close glyph inside this shell.
        shell = uipanel(strip,'BorderType','none', ...
            'BackgroundColor',strip.BackgroundColor, ...
            'AutoResizeChildren','off');
        shell.Layout.Row = 1;
        shell.Layout.Column = i;
        if i <= nDesign
            buttons(i) = uibutton(shell,'Text',names{i},'FontSize',12, ...
                'Tag',['tab' names{i}],'ButtonPushedFcn',onSelect, ...
                'Tooltip',['Show the ' names{i} ' task settings and actions.']);
            width = designWidth;
        else
            j = i-nDesign;
            buttons(i) = uibutton(shell,'Text',names{i},'FontSize',11, ...
                'HorizontalAlignment','center', ...
                'Tag',['tab' names{i}],'ButtonPushedFcn',onSelect, ...
                'Tooltip',['Show the ' names{i} ' result.']);
            buttons(i).Visible = 'off';
            width = 100;
        end
        btn = buttons(i);
        visual = uiimage(shell, ...
            'ImageSource',tabArtworkFile(names{i},width,false,i>nDesign), ...
            'Tag',['visual' btn.Tag], ...
            'Tooltip',btn.Tooltip, ...
            'ImageClickedFcn',@(~,e)onSelect(btn,e));
        if i > nDesign
            visual.Visible = 'off';
            closers(j) = uiimage(shell,'ImageSource',closeTabGlyphFile(), ...
                'Tag',['close' names{i}], ...
                'Tooltip',['Close the ' names{i} ' result tab'], ...
                'ImageClickedFcn',onClose);
            closers(j).Visible = 'off';
            layoutClosableTab(btn,visual,closers(j),width);
        else
            btn.Position = [0 0 width 32];
            visual.Position = [0 0 width 32];
        end
    end
    filler = uipanel(strip,'BorderType','none');
    filler.Layout.Column = nDesign+nResult+1;
    filler.BackgroundColor = strip.BackgroundColor;
    helpRow = uigridlayout(filler,[1 2]);
    helpRow.ColumnWidth = {'1x',86};
    helpRow.Padding = [0 0 8 0];
    helpRow.ColumnSpacing = 0;
    helpRow.BackgroundColor = strip.BackgroundColor;
    helpIcon = uiimage(helpRow, ...
        'ImageSource',ribbonHelpArtworkFile(), ...
        'Tag','ribbonHelpIcon', ...
        'Tooltip','Open Quick start help; User guide and About are in the Help menu');
    helpIcon.Layout.Column = 2;
end

function path = ribbonHelpArtworkFile()
    % The blue chevron and circular question mark echo the reference
    % toolstrip artwork while remaining sharp on high-DPI screens.
    iconDir = fullfile(tempdir,'padRibbonHelp_v1');
    if ~isfolder(iconDir), mkdir(iconDir); end
    path = fullfile(iconDir,'help.svg');
    if isfile(path), return; end
    fid = fopen(path,'w');
    assert(fid>0,'phasedArrayDesigner:helpArtwork', ...
        'Could not create the Help artwork.');
    cleanup = onCleanup(@()fclose(fid));
    fprintf(fid,['<svg xmlns="http://www.w3.org/2000/svg" ' ...
        'width="86" height="32" viewBox="0 0 86 32">' ...
        '<defs><radialGradient id="orb" cx="42%%" cy="32%%" r="65%%">' ...
        '<stop offset="0" stop-color="#ffffff"/>' ...
        '<stop offset="0.72" stop-color="#d5e2ea"/>' ...
        '<stop offset="1" stop-color="#8aa8bb"/>' ...
        '</radialGradient></defs>' ...
        '<path d="M14 0H86V32H14L0 16Z" fill="#5181aa"/>' ...
        '<circle cx="50" cy="17" r="12.3" fill="#254760" opacity=".4"/>' ...
        '<circle cx="50" cy="15" r="12.3" fill="url(#orb)" ' ...
        'stroke="#6f91a6" stroke-width="1.2"/>' ...
        '<text x="50" y="22" text-anchor="middle" ' ...
        'font-family="Arial,Helvetica,sans-serif" font-weight="bold" ' ...
        'font-size="22" fill="#355a76">?</text></svg>']);
end

function pane = makeResultPane(host,strip,button,closer,nDesign,index,width,tag)
    pane = uipanel(host,'BorderType','none','Tag',tag,'Visible','off');
    pane.Layout.Row = 1;
    pane.Layout.Column = 1;
    column = nDesign+index;
    widths = strip.ColumnWidth;
    widths(column) = {width};
    strip.ColumnWidth = widths;
    button.Visible = 'on';
    visual = findall(button.Parent,'Tag',['visual' button.Tag]);
    visual.Visible = 'on';
    closer.Visible = 'on';
    layoutClosableTab(button,visual,closer,width);
end

function removeResultPane(pane,strip,button,closer,nDesign,index)
    if isgraphics(pane), delete(pane); end
    column = nDesign+index;
    widths = strip.ColumnWidth;
    widths(column) = {0};
    strip.ColumnWidth = widths;
    button.Visible = 'off';
    visual = findall(button.Parent,'Tag',['visual' button.Tag]);
    visual.Visible = 'off';
    closer.Visible = 'off';
end

function layoutClosableTab(button,visual,closer,width)
    % All three layers fill one grid-managed shell; the close icon is the
    % only distinct hit target and sits inside the common flat tab face.
    button.Position = [0 0 width 32];
    visual.Position = [0 0 width 32];
    closer.Position = [width-23 8 16 16];
end

function path = closeTabGlyphFile()
    persistent glyphPath
    if isempty(glyphPath) || ~isfile(glyphPath)
        iconDir = fullfile(tempdir,'padTabIcons_v1');
        if ~isfolder(iconDir), mkdir(iconDir); end
        glyphPath = fullfile(iconDir,'close.png');
        [x,y] = meshgrid(1:40);
        distance = min(abs(x-y),abs(x+y-41))/sqrt(2);
        inside = x >= 9 & x <= 32 & y >= 9 & y <= 32;
        opacity = uint8(255*max(0,min(1,(2.2-distance)/1.2)).*inside);
        rgb = repmat(reshape(uint8([203 216 230]),1,1,3),40,40);
        imwrite(rgb,glyphPath,'Alpha',opacity);
    end
    path = glyphPath;
end

function path = tabArtworkFile(name,width,selected,closable)
    % Flat tab face rendered at its display size, avoiding MATLAB's
    % rounded uibutton frame and the extra outline on closable tabs.
    iconDir = fullfile(tempdir,'padFlatTabs_v2');
    if ~isfolder(iconDir), mkdir(iconDir); end
    key = regexprep(name,'[^A-Za-z0-9]','');
    path = fullfile(iconDir,sprintf('%s_%d_%d.svg',key,width,selected));
    if isfile(path), return; end
    if selected
        fill = '#202020'; weight = '600';
    else
        fill = '#084976'; weight = '400';
    end
    labelX = width/2;
    if closable, labelX = (width-20)/2; end
    fid = fopen(path,'w');
    assert(fid>0,'phasedArrayDesigner:tabArtwork', ...
        'Could not create the ribbon tab artwork.');
    cleanup = onCleanup(@()fclose(fid));
    fprintf(fid,['<svg xmlns="http://www.w3.org/2000/svg" ' ...
        'width="%d" height="32" viewBox="0 0 %d 32">' ...
        '<rect width="%d" height="32" fill="%s"/>' ...
        '<rect x="%d" width="1" height="32" fill="#337da9"/>' ...
        '<text x="%.1f" y="21" text-anchor="middle" ' ...
        'font-family="Arial,Helvetica,sans-serif" font-size="13" ' ...
        'font-weight="%s" fill="#ffffff">%s</text></svg>'], ...
        width,width,width,fill,width-1,labelX,weight,name);
end

function paintResultTabs(strip,buttons,closeButtons,selected)
    strip.BackgroundColor = [0.025 0.255 0.43];
    for i = 1:numel(buttons)
        if i == selected
            buttons(i).BackgroundColor = [0.125 0.125 0.125];
            buttons(i).FontWeight = 'bold';
        else
            buttons(i).BackgroundColor = [0.025 0.255 0.43];
            buttons(i).FontWeight = 'normal';
        end
        buttons(i).FontColor = [1 1 1];
        visual = findall(buttons(i).Parent,'Tag',['visual' buttons(i).Tag]);
        if ~isempty(visual)
            width = round(buttons(i).Position(3));
            visual.ImageSource = tabArtworkFile(buttons(i).Text,width, ...
                i==selected,i>4);
        end
    end
    for i = 1:numel(closeButtons)
        if i+4 == selected
            closeButtons(i).BackgroundColor = [0.125 0.125 0.125];
        else
            closeButtons(i).BackgroundColor = [0.025 0.255 0.43];
        end
    end
end

function [popup,btns,models] = buildElementGalleryPopup(owner,iconFor,choose)
    popup = uipanel(owner,'Tag','elementGalleryPopup', ...
        'Position',[10 10 475 475],'Visible','off');
    whole = uigridlayout(popup,[4 1]);
    whole.RowHeight = {18,290,18,105};
    whole.Padding = [10 8 10 8];
    whole.RowSpacing = 5;
    hA = uilabel(whole,'Text','ANTENNAS','FontWeight','bold');
    hA.Layout.Row = 1;
    gA = uigridlayout(whole,[3 3]);
    gA.Layout.Row = 2;
    gA.RowHeight = {'1x','1x','1x'};
    gA.ColumnWidth = {'1x','1x','1x'};
    gA.Padding = [0 0 0 0];
    gA.RowSpacing = 4; gA.ColumnSpacing = 4;
    hP = uilabel(whole,'Text','POLARIZED ANTENNAS','FontWeight','bold');
    hP.Layout.Row = 3;
    gP = uigridlayout(whole,[1 3]);
    gP.Layout.Row = 4;
    gP.ColumnWidth = {'1x','1x','1x'};
    gP.Padding = [0 0 0 0];
    gP.ColumnSpacing = 4;
    entries = { ...
        gA,'Isotropic','Isotropic','efIsotropic'; ...
        gA,'Cosine','cos^q(theta)','efCosine'; ...
        gA,'Cardioid','Cardioid','efCardioid'; ...
        gA,'Gaussian','Gaussian','efGaussian'; ...
        gA,'Sinc','Sinc','efSinc'; ...
        gA,'Patch','Patch (cos^q x lin pol)','efPatch'; ...
        gA,'Custom formula','Custom (formula)','efCustom'; ...
        gA,'CST far-field','Imported (CST far-field)','download'; ...
        gP,'Short dipole','Short dipole (z-axis)','efDipole'; ...
        gP,'Crossed dipole','Crossed dipole (RHCP)','efCrossed'; ...
        gP,'3GPP shape','3GPP TR 38.901 shape','ef3GPP'};
    btns = gobjects(size(entries,1),1);
    models = entries(:,3);
    for k = 1:size(entries,1)
        btns(k) = uibutton(entries{k,1},'Text',entries{k,2}, ...
            'Icon',iconFor(entries{k,4}), ...
            'IconAlignment','top','FontSize',10, ...
            'Tag',['gallery' regexprep(entries{k,4},'[^a-zA-Z0-9]','')], ...
            'Tooltip',entries{k,3}, ...
            'ButtonPushedFcn',@(s,e)choose(s));
    end
end

function [label,icon] = elementPreviewSpec(model)
    switch model
        case 'Isotropic',                   label = 'Isotropic';      icon = 'efIsotropic';
        case 'cos^q(theta)',                label = 'Cosine';         icon = 'efCosine';
        case 'Cardioid',                    label = 'Cardioid';       icon = 'efCardioid';
        case 'Gaussian',                    label = 'Gaussian';       icon = 'efGaussian';
        case 'Sinc',                        label = 'Sinc';           icon = 'efSinc';
        case 'Patch (cos^q x lin pol)',     label = 'Patch';          icon = 'efPatch';
        case 'Custom (formula)',            label = 'Custom';         icon = 'efCustom';
        case 'Imported (CST far-field)',    label = 'CST far-field';  icon = 'download';
        case 'Short dipole (z-axis)',       label = 'Short dipole';   icon = 'efDipole';
        case 'Dipole (linear pol)',         label = 'Plane dipole';   icon = 'efDipole';
        case 'Crossed dipole (RHCP)',       label = 'Crossed dipole'; icon = 'efCrossed';
        case '3GPP TR 38.901 shape',       label = '3GPP shape';     icon = 'ef3GPP';
        otherwise,                         label = model;            icon = 'efGallery';
    end
end

function setElementPreviewButton(button,model,iconFor)
    [label,icon] = elementPreviewSpec(model);
    button.Text = label;
    button.Icon = iconFor(icon);
    button.Tooltip = model;
    button.UserData = model;
end

function positionElementGalleryPopup(popup,owner,anchorGrid)
    % A shown UIFigure reports live grid coordinates. Headless MATLAB can
    % leave layout positions at placeholders, so fall back to the ribbon.
    figSize = owner.Position(3:4);
    popupW = min(475,max(300,figSize(1)-16));
    popupH = min(475,max(300,figSize(2)-190));
    anchor = getpixelposition(anchorGrid,true);
    if anchor(2) < figSize(2)/2
        anchor(1) = 165;
        anchor(2) = figSize(2)-280;
    end
    x = max(8,min(anchor(1),figSize(1)-popupW-8));
    y = max(8,min(anchor(2)-popupH,figSize(2)-popupH-8));
    popup.Position = [x y popupW popupH];
end

function [older,previous,current,openButton] = buildShapeRibbon(group,iconFor,choose)
    group.ColumnWidth = {'1x','1x','1x',24};
    older = uibutton(group,'Text','Diamond', ...
        'Icon',iconFor('shpDiamond'),'IconAlignment','top', ...
        'FontSize',10,'Tag','shapePreviewOlder', ...
        'Tooltip','Select an earlier array shape', ...
        'ButtonPushedFcn',@(s,e)choose(s.UserData));
    previous = uibutton(group,'Text','Circle', ...
        'Icon',iconFor('shpCircle'),'IconAlignment','top', ...
        'FontSize',10,'Tag','shapePreviewPrevious', ...
        'Tooltip','Select the previous array shape', ...
        'ButtonPushedFcn',@(s,e)choose(s.UserData));
    current = uibutton(group,'Text','Custom', ...
        'Icon',iconFor('shpCustomActive'),'IconAlignment','top', ...
        'FontSize',10,'Tag','shapePreviewCurrent', ...
        'Tooltip','Current array shape; click to open the gallery', ...
        'ButtonPushedFcn',@(s,e)choose(s.UserData));
    openButton = uibutton(group,'Text','▼','FontSize',12, ...
        'Tag','btnShapeGallery', ...
        'Tooltip','Open the array shape gallery below this ribbon', ...
        'ButtonPushedFcn',@(~,~)[]);
end

function [index,taper,axesBox,annotation] = ...
        buildGeometryDisplayControls(group,settings,change)
    group.RowHeight = {'1x','1x'};
    group.ColumnWidth = {110,'1x'};
    index = uicheckbox(group,'Text','Show Index', ...
        'Value',settings.indices,'Tag','cbLayoutIndex', ...
        'Tooltip','Show element numbers when there is room to read them.', ...
        'ValueChangedFcn',@(s,e)change('indices',s.Value));
    taper = uicheckbox(group,'Text','Show Taper', ...
        'Value',settings.taper,'Tag','cbLayoutTaper', ...
        'Tooltip','Scale marker sizes by excitation amplitude.', ...
        'ValueChangedFcn',@(s,e)change('taper',s.Value));
    axesBox = uicheckbox(group,'Text','Local Coordinates', ...
        'Value',settings.localAxes,'Tag','cbLayoutAxes', ...
        'Tooltip','Draw the array-plane +x and +y directions.', ...
        'ValueChangedFcn',@(s,e)change('localAxes',s.Value));
    annotation = uicheckbox(group,'Text','Show Annotation', ...
        'Value',settings.annotation,'Tag','cbLayoutAnnotation', ...
        'Tooltip','Show the selected element amplitude, phase and rotation.', ...
        'ValueChangedFcn',@(s,e)change('annotation',s.Value));
    index.Layout.Row = 1; index.Layout.Column = 1;
    taper.Layout.Row = 2; taper.Layout.Column = 1;
    axesBox.Layout.Row = 1; axesBox.Layout.Column = 2;
    annotation.Layout.Row = 2; annotation.Layout.Column = 2;
end

function [popup,buttons] = buildShapeGalleryPopup(owner,names,iconFor,choose)
    popup = uipanel(owner,'Tag','shapeGalleryPopup', ...
        'Position',[10 10 420 250],'Visible','off');
    whole = uigridlayout(popup,[2 1]);
    whole.RowHeight = {20,'1x'};
    whole.Padding = [8 8 8 8];
    whole.RowSpacing = 5;
    uilabel(whole,'Text','CHOOSE ARRAY SHAPE','FontWeight','bold');
    tiles = uigridlayout(whole,[2 3]);
    tiles.Layout.Row = 2;
    tiles.ColumnWidth = {'1x','1x','1x'};
    tiles.RowHeight = {'1x','1x'};
    tiles.Padding = [0 0 0 0];
    tiles.ColumnSpacing = 5; tiles.RowSpacing = 5;
    buttons = gobjects(numel(names),1);
    for i = 1:numel(names)
        buttons(i) = uibutton(tiles,'Text',names{i}, ...
            'Icon',iconFor(['shp' names{i}]), ...
            'IconAlignment','top','FontSize',11, ...
            'Tag',['btnShape' names{i}], ...
            'Tooltip',['Populate the uniform grid inside the ' names{i} ' boundary.'], ...
            'ButtonPushedFcn',@(s,e)choose(s.Text));
    end
end

function positionShapeGalleryPopup(popup,owner,anchorGrid)
    figSize = owner.Position(3:4);
    popupW = min(420,figSize(1)-16);
    popupH = min(250,figSize(2)-190);
    anchor = getpixelposition(anchorGrid,true);
    if anchor(2) < figSize(2)/2
        anchor(1) = 16;
        anchor(2) = figSize(2)-180;
    end
    x = max(8,min(anchor(1),figSize(1)-popupW-8));
    y = max(8,min(anchor(2)-popupH,figSize(2)-popupH-8));
    popup.Position = [x y popupW popupH];
end

function toggleShapeGallery(popup,owner,anchorGrid)
    if strcmp(popup.Visible,'on')
        popup.Visible = 'off';
        return;
    end
    elementPopup = findall(owner,'Tag','elementGalleryPopup');
    if ~isempty(elementPopup), elementPopup.Visible = 'off'; end
    positionShapeGalleryPopup(popup,owner,anchorGrid);
    popup.Visible = 'on';
end

function P = paletteForFigure(fig)
        %PAL  Named colours, by role, for the window's CURRENT theme.
        %   Components left on their defaults follow the theme by
        %   themselves. Anything given an explicit colour does not, and a
        %   literal RGB is right for one theme only: the dark-grey peak
        %   label and element numbers vanished on the dark axes, and the
        %   light-grey dividers became glaring bars. So every explicit
        %   colour in the app is looked up here by what it is FOR, and
        %   applyTheme re-applies them when the theme changes.
        %
        %   Popups are brought onto the main window's theme (matchTheme)
        %   so these answers hold for them too. Values were picked for
        %   contrast against that theme's backgrounds: text roles (ink,
        %   muted, good, warn, bad, accentText on accent) at least 4.5:1,
        %   plot traces at least 3:1 -- testPhasedArrayDesignerTheme
        %   measures this on the live plots, so a new role should keep
        %   to it.
        isDark = false;
        try
            isDark = strcmpi(fig.Theme.BaseColorStyle,'dark');
        catch
            % No figure themes on this release: always the light look.
        end
        if isDark
            P.accent     = [0.16 0.44 0.78];  % selection highlight
            P.accentText = [1 1 1];           % text on the accent
            P.ink        = [0.85 0.85 0.85];  % primary text on plots
            P.muted      = [0.62 0.62 0.62];  % captions, secondary text
            P.rule       = [0.32 0.32 0.32];  % dividers
            P.panelBg    = [0.13 0.13 0.13];  % theme window background
            P.axesBg     = [0.07 0.07 0.07];  % theme axes background
            P.btnBg      = [0.13 0.13 0.13];
            P.btnFg      = [0.85 0.85 0.85];
            P.good       = [0.36 0.80 0.46];
            P.warn       = [0.98 0.72 0.26];
            P.bad        = [1.00 0.45 0.42];
            P.refTrace   = [0.66 0.66 0.66];  % pinned reference curve
            P.overlayTrace = [1.00 0.60 0.22];% external (CST) overlay
            P.gridLines  = [0.30 0.30 0.30];
            P.traceTotal = [0.30 0.64 0.98];
            P.traceAF    = [0.78 0.52 0.94];
            P.traceEF    = [0.62 0.62 0.62];
            P.traceAR    = [0.36 0.80 0.46];
            P.elemFace   = [0.24 0.55 0.90];  % layout element body
            P.elemEdge   = [0.05 0.18 0.32];
            P.selEdge    = [1.00 0.80 0.20];  % selected element outline
            P.latticeDot = [0.42 0.42 0.42];  % empty lattice positions
            P.pick       = [1.00 0.55 0.20];  % 9-point candidates
            P.pickGuide  = [0.72 0.46 0.30];
            P.rotArrow   = [1.00 0.38 0.38];
            P.cellNeutral= [0.36 0.36 0.36];  % map cell holding several
            P.readOnlyBg = [0.19 0.19 0.19];  % computed table columns
            wash = 0.22;                      % phase map: toward the bg
        else
            P.accent     = [0.13 0.45 0.78];
            P.accentText = [1 1 1];
            P.ink        = [0.13 0.13 0.13];
            P.muted      = [0.42 0.42 0.42];
            P.rule       = [0.80 0.80 0.80];
            P.panelBg    = [0.96 0.96 0.96];
            P.axesBg     = [1 1 1];
            P.btnBg      = [0.96 0.96 0.96];
            P.btnFg      = [0.13 0.13 0.13];
            P.good       = [0.08 0.50 0.16];
            P.warn       = [0.60 0.34 0.00];  % 4.5:1 on the window, too
            P.bad        = [0.80 0.12 0.12];
            P.refTrace   = [0.40 0.40 0.40];
            P.overlayTrace = [0.85 0.38 0.02];
            P.gridLines  = [0.85 0.85 0.85];
            P.traceTotal = [0.00 0.45 0.74];
            P.traceAF    = [0.49 0.18 0.56];
            P.traceEF    = [0.45 0.45 0.45];
            P.traceAR    = [0.10 0.62 0.28];
            P.elemFace   = [0.10 0.45 0.75];
            P.elemEdge   = [0 0 0];
            P.selEdge    = [0.95 0.70 0.05];
            P.latticeDot = [0.72 0.72 0.72];
            P.pick       = [0.90 0.40 0.10];
            P.pickGuide  = [0.92 0.62 0.42];
            P.rotArrow   = [0.88 0.10 0.10];
            P.cellNeutral= [0.84 0.84 0.84];
            P.readOnlyBg = [0.93 0.93 0.93];  % muted text on it: >= 4.5:1
            wash = 0.30;
        end
        % Text drawn ON a coloured fill picks between these by the fill's
        % brightness (inkOn), not by the theme. Pure black and white:
        % with both extremes available every fill gets at least 4.58:1,
        % which a softened near-black cannot promise for mid-bright hues.
        P.onLight = [0 0 0];
        P.onDark  = [1 1 1];
        % Phase is cyclic -- 359 deg and 1 deg are nearly the same
        % excitation -- so its colour scale must wrap too; a linear
        % light-to-dark ramp painted them as opposites. hsv wraps, but
        % at full saturation it is harsh behind text, so it is washed
        % a little toward the axes background.
        P.phaseMap = (1-wash)*hsv(360) + wash*P.axesBg;
    end

function inkC = inkOn(P, bg)
        %INKON  Readable text colour on a fill of colour bg.
        %   Whichever of P.onLight / P.onDark has the higher WCAG contrast
        %   against the fill. A plain luma threshold on the gamma-encoded
        %   values put white on the washed reds and purples, where it
        %   falls to ~3:1; relative luminance is linear light, so the
        %   sRGB curve is undone first.
        lin = bg(:);
        lo = lin <= 0.04045;
        lin(lo) = lin(lo)/12.92;
        lin(~lo) = ((lin(~lo)+0.055)/1.055).^2.4;
        lum = [0.2126 0.7152 0.0722]*lin;
        % Against white (L = 1) the ratio is 1.05/(L+0.05); against
        % black (L = 0) it is (L+0.05)/0.05. They cross at L = 0.179.
        if (lum+0.05)/0.05 >= 1.05/(lum+0.05)
            inkC = P.onLight;
        else
            inkC = P.onDark;
        end
    end

function drawElementLayout(ax,S,view,P,amp,latticePos)
        cla(ax);
        hold(ax,'on');
        % lattice -- background dots always show the FULL (un-shape-
        % masked) MxN lattice, so the boundary that Array shape trims
        % away is visible as empty dots around the placed elements.
        % ndgrid (not meshgrid) so GR(i,j)=i, GC(i,j)=j directly -- no
        % row/col transpose to get wrong.
        [GR,GC] = ndgrid(1:S.M,1:S.N);
        if strcmp(S.mode,'Uniform grid')
            [GX,GY] = latticePos(GR,GC,S.dx,S.dy,S.gridAngle,S.stagger);
        else
            % Sparse / Sub-position modes place elements on the plain
            % orthogonal dx,dy grid regardless of Grid angle/Row stagger
            % (see onGridClick / onGeom's else-branch) -- draw the SAME
            % plain grid here. Otherwise, with a non-90 Grid angle or
            % nonzero stagger set, the reference dots would visibly skew
            % away from where clicks and placed elements actually land.
            [GX,GY] = latticePos(GR,GC,S.dx,S.dy,90,0);
        end
        plot(ax,GX(:),GY(:),'.','Color',P.latticeDot, ...
            'MarkerSize',9,'HitTest','off');
        % pending 9-point cluster
        if ~isempty(S.pending)
            cx = (S.pending(2)-1)*S.dx; cy = (S.pending(1)-1)*S.dy;
            ox = S.subOff*[-1 0 1]; oy = S.subOff*[-1 0 1];
            [PX,PY] = meshgrid(cx+ox, cy+oy);
            % faint guides from the parent cell out to each candidate
            for q = 1:numel(PX)
                plot(ax,[cx PX(q)],[cy PY(q)],'-', ...
                    'Color',P.pickGuide,'LineWidth',0.8,'HitTest','off');
            end
            plot(ax,PX(:),PY(:),'o','Color',P.pick, ...
                'MarkerSize',11,'LineWidth',1.6,'HitTest','off');
            plot(ax,cx,cy,'s','Color',P.pick, ...
                'MarkerSize',16,'LineWidth',2,'HitTest','off');
        end
        % Element body and rotation markers; feed details are in the phase map.
        isPatchShape = strcmp(S.efType,'Patch (cos^q x lin pol)') || ...
            strcmp(S.efType,'Imported (CST far-field)');
        if ~isempty(S.el)
            numOff = repmat(0.2*min(S.dx,S.dy), size(S.el,1), 1);
            for n = 1:size(S.el,1)
                if view.taper, ms = 5 + 8*amp(n); else, ms = 13; end
                isSel = ismember(n,S.sel);
                if isSel
                    ec = P.selEdge;  lw = 2.5;
                else
                    ec = P.elemEdge; lw = 0.5;
                end
                cx = S.el(n,1); cy = S.el(n,2);
                if isPatchShape
                    % Rectangular patch BODY, replacing the circle marker
                    % entirely (not drawn alongside/around it) -- this is
                    % what makes the layout actually look like a patch
                    % array instead of a generic point-source array.
                    if view.taper
                        hs = (0.16 + 0.10*amp(n)) * min(S.dx,S.dy);
                    else
                        hs = 0.26 * min(S.dx,S.dy);
                    end
                    fill(ax, cx+hs*[-1 1 1 -1], cy+hs*[-1 -1 1 1], ...
                        P.elemFace, 'EdgeColor',ec, 'LineWidth',lw, 'HitTest','off');
                    numOff(n) = max(numOff(n), hs + 0.03*min(S.dx,S.dy));
                else
                    plot(ax,cx,cy,'o','MarkerSize',ms, ...
                        'MarkerFaceColor',P.elemFace,'MarkerEdgeColor',ec, ...
                        'LineWidth',lw,'HitTest','off');
                end
                if S.el(n,5) ~= 0
                    % rotation arrow: 0 deg = up (+y), positive = clockwise.
                    % Only drawn when the element is actually rotated.
                    L    = 0.50*min(S.dx,S.dy);
                    az   = 90 - S.el(n,5);
                    xt   = S.el(n,1) + L*cosd(az);
                    yt   = S.el(n,2) + L*sind(az);
                    head = 0.35*L;
                    plot(ax,[S.el(n,1) xt],[S.el(n,2) yt], ...
                        '-','Color',P.rotArrow,'LineWidth',2,'HitTest','off');
                    plot(ax,[xt xt+head*cosd(az+150)],[yt yt+head*sind(az+150)], ...
                        '-','Color',P.rotArrow,'LineWidth',2,'HitTest','off');
                    plot(ax,[xt xt+head*cosd(az-150)],[yt yt+head*sind(az-150)], ...
                        '-','Color',P.rotArrow,'LineWidth',2,'HitTest','off');
                end
            end
            % Element numbers, in the theme's text colour, sit at each
            % body's upper-LEFT corner. Centred on the marker's top edge,
            % as they were, half of each was dark-on-blue and unreadable.
            % Straight above is no better once elements are rotated: a
            % sequential-rotation arrow pointing down from the element
            % above lands exactly there, and the up/down/left/right
            % arrows of a sequentially rotated array all miss the
            % upper-left corner. Drawn last, in one call, on an
            % axes-coloured backing that hides the grid line otherwise
            % running through the digits. Only while there is room: past
            % 10 cells a side (or 100 elements) the backings would cover
            % the neighbouring bodies, and the click readout and the
            % table still name every element.
            if view.indices && numel(numOff) <= 100 && max(S.M,S.N) <= 10
                numTxt = arrayfun(@(v) sprintf('%d',v), (1:numel(numOff))', ...
                    'UniformOutput',false);
                text(ax, S.el(:,1)-numOff, S.el(:,2)+numOff, numTxt, ...
                    'HorizontalAlignment','right','VerticalAlignment','bottom', ...
                    'FontSize',8,'Color',P.ink,'BackgroundColor',ax.Color, ...
                    'Margin',0.5,'HitTest','off');
            end
        end
        pad = 0.6*max(S.dx,S.dy) + S.subOff;
        % bounding box from the (possibly skewed) background lattice
        % dots, not the plain (N-1)*dx,(M-1)*dy formula -- a non-90 grid
        % angle shifts the true extents away from that orthogonal box
        xhi = max(GX(:)); xlo = min(GX(:));
        yhi = max(GY(:)); ylo = min(GY(:));
        if ~isempty(S.el)
            xlo = min(xlo,min(S.el(:,1))); xhi = max(xhi,max(S.el(:,1)));
            ylo = min(ylo,min(S.el(:,2))); yhi = max(yhi,max(S.el(:,2)));
        end
        xlim(ax,[xlo-pad xhi+pad]);
        ylim(ax,[ylo-pad yhi+pad]);
        drawLayoutDecorations(ax,view,P,xlo,xhi,ylo,yhi,pad);
        hold(ax,'off');
        titleTxt = sprintf('%s   |   %d elements', S.mode, size(S.el,1));
        if strcmp(S.efType,'Imported (CST far-field)')
            titleTxt=sprintf('%s | Imported element',titleTxt);
        end
        title(ax,titleTxt);
    end

function drawLayoutDecorations(ax,view,P,xlo,xhi,ylo,yhi,pad)
    if view.localAxes
        xs = (xhi-xlo)+2*pad; ys = (yhi-ylo)+2*pad;
        x0 = xlo-pad+0.12*xs; y0 = ylo-pad+0.12*ys;
        L = 0.08*min(xs,ys);
        plot(ax,[x0 x0+L],[y0 y0],'-','Color',P.ink, ...
            'LineWidth',1.6,'HitTest','off');
        plot(ax,[x0 x0],[y0 y0+L],'-','Color',P.ink, ...
            'LineWidth',1.6,'HitTest','off');
        text(ax,x0+1.12*L,y0,'x','Color',P.ink, ...
            'FontSize',9,'HitTest','off');
        text(ax,x0,y0+1.12*L,'y','Color',P.ink, ...
            'FontSize',9,'HitTest','off');
    end
end

function updateLayoutPanelTitle(panel,S,view)
    if view.annotation && ~isempty(S.el) && ~isempty(S.sel)
        n = S.sel(1);
        if n >= 1 && n <= size(S.el,1)
            panel.Title = sprintf('Element layout  |  #%d  Amp %.3g  Phase %.3g°  Rot %.3g°', ...
                n,S.el(n,3),S.el(n,4),S.el(n,5));
            return;
        end
    end
    panel.Title = 'Element layout';
end

    function [RGB,A] = ribbonGlyph(name)
        %RIBBONGLYPH  Draw one 128px glyph as an RGB image plus alpha mask.
        %   Composed from algebraic masks (rectangles, rounded
        %   rectangles, discs, rings, triangles) supersampled 4x and
        %   box-filtered down, rather than rendered from a figure: this
        %   runs identically headless, needs no graphics context, and is
        %   deterministic.
        %
        %   Each glyph is a stack of {mask, colour} LAYERS painted in
        %   order, so a glyph can be two-tone the way toolstrip icons
        %   normally are -- a steel floppy with a light shutter, an amber
        %   folder with a lighter front. A single flat colour per glyph
        %   read as a set of monochrome silhouettes, which is what made
        %   the first version look nothing like a real ribbon.
        %
        %   Colours are picked to carry meaning rather than to decorate:
        %   green for the action that produces a result, amber for the
        %   one that reads from disk, steel for the one that writes,
        %   blue for the analysis pair.
        % 128 px source, not 32. A button's icon area is around 40 px and
        % more than that on a HiDPI display, so a 32 px glyph was being
        % scaled UP -- which is what made these look blocky. Downscaling
        % a larger source is what stays crisp at any button size. 4x
        % supersampling on top of that still gives 16 samples per output
        % pixel for the antialiasing.
        N = 128; SS = 4; n = N*SS;
        px = (0.5:n-0.5)/n;
        [X,Y] = meshgrid(px,1-px);          % Y increases upward
        steel  = [0.52 0.58 0.65];
        light  = [0.88 0.91 0.94];
        amber  = [0.87 0.63 0.16];
        amberL = [0.98 0.78 0.34];
        green  = [0.30 0.75 0.38];
        blue   = [0.25 0.65 0.92];
        switch name
            case 'new'       % yellow plus, as in the reference ribbon
                plus = rr(X,Y,.5,.5,.22,.88,.03) | ...
                    rr(X,Y,.5,.5,.88,.22,.03);
                L = { plus, [0.24 0.19 0.07]; ...
                      rr(X,Y,.5,.5,.17,.83,.02) | ...
                      rr(X,Y,.5,.5,.83,.17,.02), [1.00 0.81 0.31] };
            case 'save'      % pale floppy with dark shutter and label slot
                L = { rr(X,Y,.5,.5,.83,.87,.025), [0.09 0.10 0.11]; ...
                      rr(X,Y,.5,.5,.77,.81,.015), light; ...
                      rct(X,Y,.28,.72,.61,.88), [0.09 0.10 0.11]; ...
                      rct(X,Y,.32,.61,.65,.85), light; ...
                      rct(X,Y,.63,.69,.65,.85), steel; ...
                      rr(X,Y,.5,.31,.59,.33,.02), [0.09 0.10 0.11]; ...
                      rr(X,Y,.5,.29,.49,.23,.01), light };
            case 'import'    % narrow curved mint arrow in the reference
                arrow = stroke(X,Y,[.42 .81],[.50 .36],.13) | ...
                    stroke(X,Y,[.24 .46],[.50 .19],.13) | ...
                    stroke(X,Y,[.76 .46],[.50 .19],.13) | ...
                    ell(X,Y,.40,.82,.09,.06,-.3);
                L = { arrow, [0.71 0.96 0.70] };
            case 'open'      % folder: darker back with tab, lighter front
                L = { rct(X,Y,.08,.50,.68,.82) | rr(X,Y,.5,.50,.84,.44,.05), amber; ...
                      rr(X,Y,.5,.36,.84,.36,.05),            amberL };
            case 'play'      % solid right-pointing triangle
                L = { tri(X,Y,[.28 .12; .28 .88; .86 .50]),  green };
            case 'polar2d'   % white polar dial with a blue bulb-shaped cut
                bearing = atan2(Y-.52,X-.5);
                ticks = rng(X,Y,.5,.52,.47,.016) & ...
                    abs(sin(8*bearing)) < .13;
                dial = stroke(X,Y,[.20 .52],[.80 .52],.012) | ...
                    stroke(X,Y,[.50 .21],[.50 .83],.012) | ...
                    rng(X,Y,.5,.52,.31,.012);
                bulb = (rng(X,Y,.5,.61,.205,.026) & Y>=.59) | ...
                    stroke(X,Y,[.30 .61],[.42 .40],.029) | ...
                    stroke(X,Y,[.70 .61],[.58 .40],.029) | ...
                    stroke(X,Y,[.42 .40],[.58 .40],.029) | ...
                    stroke(X,Y,[.50 .40],[.50 .23],.031);
                L = { ticks, [0.39 0.46 0.53]; ...
                      disc(X,Y,.5,.52,.41), [0.97 0.98 0.99]; ...
                      rng(X,Y,.5,.52,.41,.018) | dial, [0.58 0.66 0.73]; ...
                      bulb, [0.24 0.66 0.92]; ...
                      disc(X,Y,.5,.23,.027), [0.24 0.66 0.92] };
            case 'wave'      % scan loss: smooth total curve
                axesMask = stroke(X,Y,[.20 .23],[.20 .79],.055) | ...
                    stroke(X,Y,[.20 .23],[.83 .23],.055);
                t = max(0,min(1,(X-.23)/.58));
                inPlot = X>=.23 & X<=.81;
                total = inPlot & abs(Y-(.76-.48*t.^3))<=.027;
                L = { disc(X,Y,.5,.5,.48),                    light; ...
                      axesMask,                              [0.12 0.17 0.21]; ...
                      total,                                 [0.09 0.27 0.94] };
            case 'band'      % two scan responses across frequency
                axesMask = stroke(X,Y,[.20 .22],[.20 .79],.05) | ...
                    stroke(X,Y,[.20 .22],[.83 .22],.05);
                t = max(0,min(1,(X-.24)/.56));
                inPlot = X>=.24 & X<=.80;
                phaseShifter = inPlot & abs(Y-(.68-.32*t.^2))<=.026;
                trueDelay = inPlot & abs(Y-(.53+.035*t))<=.025;
                L = { disc(X,Y,.5,.5,.48),                    light; ...
                      axesMask,                              [0.12 0.17 0.21]; ...
                      phaseShifter,                          [0.09 0.27 0.94]; ...
                      trueDelay,                             green };
            case 'coverage'  % polar scan cone with a highlighted sector
                rings = rng(X,Y,.5,.5,.33,.027) | ...
                    rng(X,Y,.5,.5,.20,.027);
                radial = stroke(X,Y,[.50 .50],[.82 .50],.024) | ...
                    stroke(X,Y,[.50 .50],[.50 .82],.024) | ...
                    stroke(X,Y,[.50 .50],[.18 .50],.024) | ...
                    stroke(X,Y,[.50 .50],[.50 .18],.024);
                sector = disc(X,Y,.5,.5,.31) & Y>=.5 & X>=.5;
                marker = disc(X,Y,.72,.72,.043);
                L = { disc(X,Y,.5,.5,.48),                    light; ...
                      sector,                                [0.67 0.86 0.96]; ...
                      radial,                                [0.12 0.17 0.21]; ...
                      rings,                                 [0.09 0.27 0.94]; ...
                      marker,                                green };
            case 'peak'      % max scan: tiled phased array and beam fan
                tileEdges = false(n,n); tileFaces = false(n,n);
                origin = [.13 .28]; dc = [.18 -.035]; dr = [.10 .13];
                for row = 0:1
                    for col = 0:2
                        p0 = origin + col*dc + row*dr;
                        P = [p0; p0+dc; p0+dc+dr; p0+dr];
                        Pi = mean(P,1) + .86*(P-mean(P,1));
                        tileEdges = tileEdges | tri(X,Y,P(1:3,:)) | ...
                            tri(X,Y,P([1 3 4],:));
                        tileFaces = tileFaces | tri(X,Y,Pi(1:3,:)) | ...
                            tri(X,Y,Pi([1 3 4],:));
                    end
                end
                % A few broad, rounded lobes carry the reference image's
                % yellow/green fan without turning into noise at 40 px.
                backPetals = ell(X,Y,.28,.65,.22,.075,2.08) | ...
                    ell(X,Y,.40,.72,.23,.077,1.82) | ...
                    ell(X,Y,.52,.75,.24,.080,1.57) | ...
                    ell(X,Y,.64,.70,.22,.077,1.30) | ...
                    ell(X,Y,.75,.62,.20,.073,1.02);
                frontPetals = ell(X,Y,.36,.55,.15,.068,1.90) | ...
                    ell(X,Y,.50,.58,.17,.072,1.56) | ...
                    ell(X,Y,.63,.53,.14,.065,1.22);
                L = { disc(X,Y,.5,.5,.48),                    light; ...
                      tileEdges,                             [0.16 0.21 0.23]; ...
                      tileFaces,                             [0.78 0.64 0.35]; ...
                      backPetals,                            [1.00 0.79 0.07]; ...
                      frontPetals,                           [0.16 0.72 0.49]; ...
                      ell(X,Y,.50,.60,.15,.055,1.57),        [1.00 0.89 0.08] };
            case 'doc'       % page with text lines
                L = { rr(X,Y,.5,.5,.66,.88,.06),             light; ...
                      rct(X,Y,.28,.72,.66,.73),              blue; ...
                      rct(X,Y,.28,.72,.47,.54),              blue; ...
                      rct(X,Y,.28,.60,.28,.35),              blue };
            case 'grid'      % lattice, echoing the app's own icon
                M = false(n,n);
                for i2 = 1:3
                    for j2 = 1:3
                        M = M | rr(X,Y,.20+(j2-1)*.30,.20+(i2-1)*.30,.20,.20,.05);
                    end
                end
                L = { M, blue };
            case 'trash'
                L = { rct(X,Y,.24,.76,.72,.80) | rct(X,Y,.42,.58,.82,.90), steel; ...
                      rr(X,Y,.5,.40,.56,.62,.06),            blue; ...
                      rct(X,Y,.44,.50,.16,.60),              light; ...
                      rct(X,Y,.56,.62,.16,.60),              light };
            case 'cross'
                L = { (abs((X-.5)-(Y-.5))<.075 | abs((X-.5)+(Y-.5))<.075) ...
                      & hypot(X-.5,Y-.5)<.40,                steel };
            case 'pin'       % push pin (pin reference): knob, flange, needle
                L = { stroke(X,Y,[.5 .52],[.5 .08],.06),     steel; ...
                      rr(X,Y,.5,.84,.36,.12,.05) | rct(X,Y,.41,.59,.60,.82) | ...
                      rr(X,Y,.5,.57,.62,.09,.04),            blue };
            case 'download'
                % Green download arrow entering a dark, open tray. The
                % pale disc preserves the reference's contrast on the
                % dark ribbon and makes the tray legible at button size.
                tray = rr(X,Y,.5,.24,.68,.12,.035) | ...
                    rr(X,Y,.22,.36,.12,.32,.035) | ...
                    rr(X,Y,.78,.36,.12,.32,.035);
                arrow = rr(X,Y,.5,.65,.12,.38,.045) | ...
                    tri(X,Y,[.25 .59; .75 .59; .50 .34]);
                L = { disc(X,Y,.5,.5,.48),                    light; ...
                      tray,                                  [0.08 0.09 0.10]; ...
                      arrow,                                 [0.30 0.64 0.31] };
            case 'macro'     % script: the CSV page with code brackets
                L = { rr(X,Y,.5,.5,.66,.88,.06),             light; ...
                      stroke(X,Y,[.43 .67],[.28 .50],.075) | ...
                      stroke(X,Y,[.28 .50],[.43 .33],.075) | ...
                      stroke(X,Y,[.57 .67],[.72 .50],.075) | ...
                      stroke(X,Y,[.72 .50],[.57 .33],.075),  blue };
            % ---- element gallery: original schematic pattern glyphs ---
            case 'efGallery'
                L = { disc(X,Y,.5,.5,.47), light; ...
                      ell(X,Y,.32,.65,.17,.21,0) | ...
                      ell(X,Y,.68,.65,.17,.21,0) | ...
                      ell(X,Y,.32,.35,.17,.21,0) | ...
                      ell(X,Y,.68,.35,.17,.21,0), blue; ...
                      disc(X,Y,.5,.5,.105), green };
            case 'efIsotropic'
                L = { disc(X,Y,.5,.5,.42), [0.03 .48 .53]; ...
                      ell(X,Y,.40,.64,.16,.13,-.3), [0.27 .76 .76]; ...
                      rng(X,Y,.5,.5,.42,.026), steel };
            case 'efCosine'
                fan = ell(X,Y,.58,.5,.30,.36,0) & X>=.29;
                L = { stroke(X,Y,[.24 .17],[.24 .83],.055), steel; ...
                      fan, blue; ...
                      ell(X,Y,.52,.59,.17,.24,-.25), [0.40 .83 .96] };
            case 'efCardioid'
                a = atan2(Y-.5,X-.33); rad = hypot(X-.33,Y-.5);
                card = rad <= .25*(1+cos(a));
                L = { card, amber; ...
                      ell(X,Y,.58,.62,.19,.12,-.25) & card, amberL; ...
                      disc(X,Y,.33,.5,.055), steel };
            case 'efGaussian'
                L = { ell(X,Y,.55,.5,.37,.27,-.30), blue; ...
                      ell(X,Y,.55,.5,.28,.20,-.30), green; ...
                      ell(X,Y,.55,.5,.17,.13,-.30), amberL };
            case 'efSinc'
                L = { ell(X,Y,.50,.50,.28,.30,0), blue; ...
                      ell(X,Y,.14,.50,.10,.18,0) | ...
                      ell(X,Y,.86,.50,.10,.18,0), green; ...
                      ell(X,Y,.50,.57,.15,.16,0), amberL };
            case 'efPatch'
                % A broadside radiation lobe over a gold patch. The
                % supplied gallery has no patch tile, so use its white
                % card and yellow/green/cyan field palette here.
                boardPoly = [.18 .27; .82 .27; .72 .41; .28 .41];
                board = tri(X,Y,boardPoly(1:3,:)) | ...
                    tri(X,Y,boardPoly([1 3 4],:));
                gridLines = stroke(X,Y,[.29 .31],[.72 .31],.012) | ...
                    stroke(X,Y,[.27 .36],[.73 .36],.012) | ...
                    stroke(X,Y,[.38 .27],[.42 .41],.012) | ...
                    stroke(X,Y,[.50 .27],[.50 .41],.012) | ...
                    stroke(X,Y,[.62 .27],[.58 .41],.012);
                L = { true(n,n), [1 1 1]; ...
                      board, [0.49 0.34 0.12]; ...
                      board & Y>=.30, [0.93 0.71 0.22]; ...
                      gridLines & board, [0.62 0.43 0.13]; ...
                      ell(X,Y,.50,.46,.22,.17,0), [0.10 0.23 0.82]; ...
                      ell(X,Y,.50,.51,.24,.19,0), [0.11 0.66 0.86]; ...
                      ell(X,Y,.50,.57,.23,.23,0), [0.16 0.79 0.50]; ...
                      ell(X,Y,.50,.64,.21,.24,0), [0.99 0.75 0.02]; ...
                      ell(X,Y,.49,.69,.17,.18,0), [1.00 0.91 0.03]; ...
                      ell(X,Y,.44,.77,.075,.07,0), [1.00 0.98 0.54] };
            case 'efCustom'
                axesMask = stroke(X,Y,[.22 .24],[.22 .79],.05) | ...
                    stroke(X,Y,[.22 .24],[.83 .24],.05);
                curve = abs(Y-(.44+.23*sin(10*(X-.23))))<.026 & ...
                    X>=.23 & X<=.82;
                L = { rr(X,Y,.5,.5,.84,.84,.08), light; ...
                      axesMask, steel; ...
                      curve, blue };
            case 'efDipole'
                wire = stroke(X,Y,[.5 .24],[.5 .76],.06);
                L = { ell(X,Y,.26,.5,.20,.34,0) | ...
                      ell(X,Y,.74,.5,.20,.34,0), blue; ...
                      wire, steel; ...
                      disc(X,Y,.5,.5,.055), amberL };
            case 'efCrossed'
                crossWire = stroke(X,Y,[.23 .5],[.77 .5],.05) | ...
                    stroke(X,Y,[.5 .23],[.5 .77],.05);
                petals = ell(X,Y,.5,.76,.15,.18,0) | ...
                    ell(X,Y,.5,.24,.15,.18,0) | ...
                    ell(X,Y,.76,.5,.18,.15,0) | ...
                    ell(X,Y,.24,.5,.18,.15,0);
                L = { petals, blue; ...
                      crossWire, steel; ...
                      disc(X,Y,.5,.5,.07), green };
            case 'ef3GPP'
                fan = ell(X,Y,.50,.64,.35,.25,0) & Y>=.43;
                L = { rr(X,Y,.5,.26,.56,.13,.03), steel; ...
                      fan, blue; ...
                      ell(X,Y,.5,.64,.25,.18,0) & Y>=.43, green; ...
                      ell(X,Y,.5,.64,.13,.11,0) & Y>=.43, amberL };
            % ---- array-shape gallery -----------------------------------
            % Depict the populated ELEMENTS, as in an array gallery,
            % rather than solid silhouettes of the trimming boundaries.
            % A dark thumbnail tile preserves the current dark ribbon;
            % selected shapes use brighter dots on their blue button.
            case {'shpCustom','shpDiamond','shpHexagon', ...
                  'shpOctagon','shpCircle','shpEllipse', ...
                  'shpCustomActive','shpDiamondActive', ...
                  'shpHexagonActive','shpOctagonActive', ...
                  'shpCircleActive','shpEllipseActive'}
                activeIcon = endsWith(name,'Active');
                shapeName = name(4:end);
                if activeIcon, shapeName = shapeName(1:end-6); end
                if activeIcon, dotColor = [0.48 0.84 1.00];
                else, dotColor = [0.18 0.67 0.98]; end
                smoothBody = false(n,n); smoothEdge = false(n,n);
                if strcmp(shapeName,'Circle')
                    smoothBody = disc(X,Y,.5,.5,.38);
                    smoothEdge = rng(X,Y,.5,.5,.38,.025);
                elseif strcmp(shapeName,'Ellipse')
                    smoothBody = ell(X,Y,.5,.5,.43,.28,0);
                    smoothEdge = smoothBody & ...
                        ~ell(X,Y,.5,.5,.405,.255,0);
                end
                if activeIcon, edgeColor = [0.31 0.69 0.88];
                else, edgeColor = [0.16 0.46 0.65]; end
                L = { rr(X,Y,.5,.5,.88,.88,.045), [0.29 0.36 0.42]; ...
                      rr(X,Y,.5,.5,.85,.85,.035), [0.08 0.12 0.16]; ...
                      smoothBody, [0.09 0.19 0.26]; ...
                      smoothEdge, edgeColor; ...
                      shapeLatticeDots(X,Y,shapeName), dotColor };
            otherwise
                L = { disc(X,Y,.5,.5,.40),                   steel };
        end
        % Paint in order: a later layer overwrites the colour AND claims
        % the alpha wherever it covers an earlier one.
        Am = zeros(n,n); C = zeros(n,n,3);
        for k = 1:size(L,1)
            m = L{k,1}; col = L{k,2};
            Am(m) = 1;
            for ch = 1:3
                Ck = C(:,:,ch); Ck(m) = col(ch); C(:,:,ch) = Ck;
            end
        end
        A = boxdown(Am,SS);
        RGB = zeros(N,N,3);
        for ch = 1:3, RGB(:,:,ch) = boxdown(C(:,:,ch),SS); end
        % Un-premultiply: boxdown averaged colour against black in the
        % edge pixels, which would fringe every glyph with a dark halo
        % once the alpha channel faded them out.
        safe = max(A,1e-6);
        for ch = 1:3, RGB(:,:,ch) = min(RGB(:,:,ch)./safe,1); end
    end

    function dots = shapeLatticeDots(X,Y,shapeName)
        % Quantize the raster to its nearest lattice site rather than
        % looping over every supersampled circular dot mask. The octagon
        % uses deliberate row widths: analytic clipping of only nine rows
        % made it nearly square. The circular options use fewer dots
        % inside their smooth boundary, so the dot silhouette does not
        % compete with the circle or ellipse at ribbon-icon size.
        if strcmp(shapeName,'Circle') || strcmp(shapeName,'Ellipse')
            if strcmp(shapeName,'Circle')
                firstX = .205; pitchX = .0983; nCol = 7;
                firstY = .205; pitchY = .0983; nRow = 7;
                radiusX = .32; radiusY = .32; dotRadius = .026;
            else
                firstX = .17; pitchX = .0825; nCol = 9;
                firstY = .32; pitchY = .09; nRow = 5;
                radiusX = .37; radiusY = .23; dotRadius = .026;
            end
            row = round((Y-firstY)/pitchY);
            col = round((X-firstX)/pitchX);
            cx = firstX + col*pitchX;
            cy = firstY + row*pitchY;
            valid = col>=0 & col<nCol & row>=0 & row<nRow;
            inside = ((cx-.5)/radiusX).^2 + ...
                ((cy-.5)/radiusY).^2 <= 1;
            dots = valid & inside & hypot(X-cx,Y-cy)<=dotRadius;
            return;
        end
        first = .17; pitch = .0825;
        col = round((X-first)/pitch);
        row = round((Y-first)/pitch);
        cx = first + col*pitch;
        cy = first + row*pitch;
        valid = col>=0 & col<=8 & row>=0 & row<=8;
        dx = cx-.5; dy = cy-.5;
        switch shapeName
            case 'Custom'
                inside = true(size(X));
            case 'Diamond'
                inside = abs(dx)+abs(dy) <= .36;
            case 'Hexagon'
                inside = poly(cx,cy,.5,.5,.43,6,pi/6);
            case 'Octagon'
                widths = [5 7 9 9 9 9 9 7 5];
                rowIndex = min(9,max(1,row+1));
                inside = abs(col-4) <= (widths(rowIndex)-1)/2;
            otherwise
                inside = false(size(X));
        end
        dots = valid & inside & hypot(X-cx,Y-cy)<=.0245;
    end

    function m = rct(X,Y,x0,x1,y0,y1), m = X>=x0 & X<=x1 & Y>=y0 & Y<=y1; end
    function m = stroke(X,Y,a,b,width)
        % Rounded line segment, independent of image resolution.
        d = b-a;
        t = max(0,min(1,((X-a(1))*d(1)+(Y-a(2))*d(2))/sum(d.^2)));
        m = hypot(X-(a(1)+t*d(1)),Y-(a(2)+t*d(2))) <= width/2;
    end
    function m = ell(X,Y,cx,cy,a,b,rot)
        xr =  (X-cx)*cos(rot) + (Y-cy)*sin(rot);
        yr = -(X-cx)*sin(rot) + (Y-cy)*cos(rot);
        m = (xr/a).^2 + (yr/b).^2 <= 1;
    end
    function m = poly(X,Y,cx,cy,R,nSides,rot)
        % Regular n-gon by the closed form "distance to the nearest edge
        % plane": fold the polar angle into one wedge, then compare the
        % radius against the apothem. Vectorises over the whole
        % supersampled grid, which inpolygon would not.
        th = atan2(Y-cy,X-cx) - rot;
        a  = 2*pi/nSides;
        m  = hypot(X-cx,Y-cy) .* cos(mod(th+a/2,a)-a/2) <= R*cos(pi/nSides);
    end
    function m = disc(X,Y,cx,cy,rad), m = hypot(X-cx,Y-cy) <= rad; end
    function m = rng(X,Y,cx,cy,rad,t)
        d = hypot(X-cx,Y-cy); m = d <= rad+t/2 & d >= rad-t/2;
    end
    function m = rr(X,Y,cx,cy,w,h,rad)
        dx = max(abs(X-cx)-(w/2-rad),0);
        dy = max(abs(Y-cy)-(h/2-rad),0);
        m  = hypot(dx,dy) <= rad;
    end
    function m = tri(X,Y,P)
        % Inside-triangle test by consistent sign of the three edge
        % cross-products. Written out rather than using inpolygon so it
        % stays vectorised over the whole supersampled grid.
        s = @(a,b) (b(1)-a(1)).*(Y-a(2)) - (b(2)-a(2)).*(X-a(1));
        d1 = s(P(1,:),P(2,:)); d2 = s(P(2,:),P(3,:)); d3 = s(P(3,:),P(1,:));
        m = ~(((d1<0)|(d2<0)|(d3<0)) & ((d1>0)|(d2>0)|(d3>0)));
    end
    function out = boxdown(in,f)
        nn = size(in,1)/f;
        out = squeeze(mean(reshape(mean(reshape(in,f,[]),1),nn,f,nn),2));
    end


    function ffS = loadCstFarfieldASCII(filename)
        % Supported contract: CST realized gain in dBi, complex Theta/Phi,
        % a complete regular spherical grid. Never synthesize missing angles.
        rawText = fileread(filename);
        inLines = regexp(rawText,'\r\n|\n|\r','split');
        hdr = ''; hdrIdx = 0;
        for pr = 1:numel(inLines)
            if contains(inLines{pr},'Theta','IgnoreCase',true) && ...
                    contains(inLines{pr},'Phi','IgnoreCase',true) && ...
                    contains(inLines{pr},'Abs(','IgnoreCase',true)
                hdr = inLines{pr}; hdrIdx = pr; break;
            end
        end
        % The closing bracket MUST be escaped: [^\]] not [^]].
        %
        % In PCRE, "[^]]" reads as "any character except ]". MATLAB's
        % engine does not agree -- the ] immediately after [^ closes an
        % empty set, so the pattern never matches anything. numel(desc)
        % came back 0 for every file, the "~= 5" test below fired every
        % time, and the import refused EVERY export before any other
        % check ran. Measured on a real CST header: 0 descriptors with
        % [^]], 5 with [^\]].
        desc = regexpi(hdr,'(Abs|Phase)\(\s*([^)]*?)\s*\)\s*\[\s*([^\]]*?)\s*\]','tokens');
        if numel(desc) ~= 5
            error('PAD:ImportSchema',['Export realized gain with Abs(total), Abs(Theta), Phase(Theta), ' ...
                'Abs(Phi), Phase(Phi), and optional axial ratio. Unrecognized columns are not guessed.']);
        end
        names = cellfun(@(d)lower(strtrim(d{2})),desc,'UniformOutput',false);
        kinds = cellfun(@(d)lower(d{1}),desc,'UniformOutput',false);
        units = cellfun(@(d)lower(strtrim(d{3})),desc,'UniformOutput',false);
        if ~strcmp(kinds{1},'abs') || ...
                ~any(strcmp(names{1},{'grlz','realized gain','realised gain','realizedgain'})) || ...
                ~isequal(kinds,{'abs','abs','phase','abs','phase'}) || ...
                ~all(strcmp(units([1 2 4]),'dbi')) || ...
                ~all(ismember(units([3 5]),{'deg','deg.','degrees'})) || ...
                ~isequal(sort(names([2 4])),{'phi','theta'}) || ...
                ~strcmp(names{2},names{3}) || ~strcmp(names{4},names{5})
            error('PAD:ImportQuantity',['Only realized-gain dBi exports with complex Theta/Phi components are supported. ' ...
                'Directivity, ordinary gain, linear fields and other component bases must be re-exported.']);
        end
        % No other descriptors may shift these numeric column indices.
        if numel(regexpi(hdr,'(?:Abs|Phase)\(')) ~= 5
            error('PAD:ImportSchema','Ambiguous column descriptors.');
        end
        if isempty(regexpi(hdr,'Theta\s*\[\s*deg\.?\s*\]\s*Phi\s*\[\s*deg\.?\s*\]','once'))
            error('PAD:ImportAngles','Expected Theta and Phi as the first two columns, in degrees.');
        end
        nExpected = 7 + double(contains(lower(hdr),'ax.ratio'));
        rows = cell(0,1);
        for pr = hdrIdx+1:numel(inLines)
            ln = strtrim(inLines{pr});
            if isempty(ln) || ~isempty(regexp(ln,'^[-=\s]+$','once')), continue; end
            vals = str2double(strsplit(ln));
            if numel(vals) ~= nExpected || ~isreal(vals) || any(isnan(vals))
                error('PAD:ImportRow','Malformed numeric data at line %d. No rows were discarded.',pr);
            end
            if any(~isfinite(vals([1 2 5 7]))) || ...
                    any(vals([3 4 6]) == Inf)
                error('PAD:ImportRow','Invalid angle, phase or gain at line %d.',pr);
            end
            rows{end+1,1} = vals; %#ok<AGROW>
        end
        if isempty(rows), error('PAD:ImportEmpty','No far-field samples found.'); end
        data = vertcat(rows{:});
        th = data(:,1); ph = data(:,2);
        if any(abs(th)>180) || any(abs(ph)>360)
            error('PAD:ImportAngles','Theta must be within -180..180 and phi within -360..360 degrees.');
        end
        neg = th<0; th(neg) = -th(neg); ph(neg) = ph(neg)+180; ph = mod(ph,360);
        if strcmp(names{2},'theta'), ac=[4 6]; pc=[5 7]; else, ac=[6 4]; pc=[7 5]; end
        et = 10.^(data(:,ac(1))/20).*exp(1j*deg2rad(data(:,pc(1))));
        ep = 10.^(data(:,ac(2))/20).*exp(1j*deg2rad(data(:,pc(2))));
        et(neg)=-et(neg); ep(neg)=-ep(neg);
        total = 10.^(data(:,3)/10);
        pw = abs(et).^2+abs(ep).^2;
        if any(~isfinite(pw)) || any(~isfinite(total)) || ~any(pw>0) || ~any(total>0)
            error('PAD:ImportPower','The export has no usable finite positive power.');
        end
        good = max(total,pw)>max([total;pw])*1e-7;
        if any(abs(10*log10(pw(good)./total(good)))>0.15)
            error('PAD:ImportComponents',['Theta/Phi powers disagree with total realized gain. ' ...
                'Export actual complex components, not placeholder copies of the total.']);
        end
        [keys,~,group] = unique(round([th ph],6),'rows');
        cnt = accumarray(group,1);
        etMean = accumarray(group,real(et))./cnt + 1j*accumarray(group,imag(et))./cnt;
        epMean = accumarray(group,real(ep))./cnt + 1j*accumarray(group,imag(ep))./cnt;
        sc = accumarray(group,hypot(abs(et),abs(ep)),[],@max);
        dev = accumarray(group,hypot(abs(et-etMean(group)),abs(ep-epMean(group))),[],@max);
        if any(cnt>1 & dev>max(0.05*sc,1e-6*max(sc)))
            error('PAD:ImportDuplicate','Duplicate directions have conflicting complex fields.');
        end
        th=keys(:,1); ph=keys(:,2); et=etMean; ep=epMean;
        tu=unique(th); pu=unique(ph);
        if numel(tu)<3 || abs(tu(1))>1e-6 || abs(tu(end)-180)>1e-6
            error('PAD:ImportTheta','A complete theta range 0..180 degrees is required. Missing radiation is not assumed zero.');
        end
        gaps=diff([pu;pu(1)+360]);
        if numel(pu)<8 || max(gaps)>45+1e-6 || ...
                max(abs(gaps-median(gaps)))>1e-5 || ...
                max(abs(diff(tu)-median(diff(tu))))>1e-5 || ...
                numel(th)~=numel(tu)*numel(pu)
            error('PAD:ImportCoverage',['A complete regular Theta/Phi grid is required. ' ...
                'Re-export theta 0..180 and a full phi circle, with no missing rows or angular gaps.']);
        end
        % Periodic copies support a shifted azimuth origin as well as phi=0.
        tx=[th;th;th]; px=[ph-360;ph;ph+360]; ex=[et;et;et]; ey=[ep;ep;ep];
        ffS=struct();
        ffS.FReEth=scatteredInterpolant(tx,px,real(ex),'linear','none');
        ffS.FImEth=scatteredInterpolant(tx,px,imag(ex),'linear','none');
        ffS.FReEph=scatteredInterpolant(tx,px,real(ey),'linear','none');
        ffS.FImEph=scatteredInterpolant(tx,px,imag(ey),'linear','none');
        pw=abs(et).^2+abs(ep).^2;
        [pk,ix]=max(pw);
        ffS.peak=[10*log10(pk),th(ix),ph(ix)];
        if th(ix)==0, ffS.peak(3)=0; end
        ffS.peakMag=sqrt(pk);
        ffS.thMin=0; ffS.thMax=180; ffS.phiFull=true; ffS.thetaFull=true;
        ffS.compThetaPhi=true; ffS.compNames={'Theta','Phi'}; ffS.noComponents=false;
        ffS.quantity='realized_gain'; ffS.schemaVersion=2; ffS.dupWarn='';
        ffS=refreshImportedPowerMetadata(ffS);
    end

    function ffS = refreshImportedPowerMetadata(ffS)
        % Derived cache: recompute on import and config restore from the
        % validated fields. This is not the entered unit-cell efficiency.
        ti=linspace(0,180,181); pi_=linspace(0,360,361); [tt,pp]=meshgrid(ti,pi_);
        ei=ffS.FReEth(tt,pp)+1j*ffS.FImEth(tt,pp);
        ej=ffS.FReEph(tt,pp)+1j*ffS.FImEph(tt,pp);
        sampledPower=abs(ei).^2+abs(ej).^2;
        integral=trapz(deg2rad(pi_(:)),trapz(deg2rad(ti), ...
            sampledPower.*sind(tt),2))/(4*pi);
        % The app's 1-degree theta / 2-degree phi grid cannot resolve
        % every feature in a finer CST export. Reject a field when that
        % grid misses appreciable source power or the source peak.
        pts=ffS.FReEth.Points;
        base=find(pts(:,2)>=0 & pts(:,2)<360);
        [sorted,order]=sortrows(pts(base,:));
        thNative=unique(sorted(:,1)); phNative=unique(sorted(:,2));
        nTheta=numel(thNative); nPhi=numel(phNative);
        nativeIdx=base(order);
        etNative=reshape(ffS.FReEth.Values(nativeIdx)+ ...
            1j*ffS.FImEth.Values(nativeIdx),nPhi,nTheta);
        epNative=reshape(ffS.FReEph.Values(nativeIdx)+ ...
            1j*ffS.FImEph.Values(nativeIdx),nPhi,nTheta);
        nativePower=abs(etNative).^2+abs(epNative).^2;
        nativeThetaIntegral=trapz(deg2rad(thNative).', ...
            nativePower.*sind(thNative).',2);
        nativeIntegral=trapz(deg2rad([phNative;phNative(1)+360]), ...
            [nativeThetaIntegral;nativeThetaIntegral(1)])/(4*pi);
        sourceFinerThanGrid=median(diff(thNative))<1-1e-9 || ...
            median(diff([phNative;phNative(1)+360]))<2-1e-9;
        if sourceFinerThanGrid && nativeIntegral>0 && integral>0
            integralErrorDb=abs(10*log10(integral/nativeIntegral));
        elseif sourceFinerThanGrid
            integralErrorDb=Inf;
        else
            integralErrorDb=0; % sparse source nodes do not integrate sin(theta) reliably
        end
        sampledPeak=max(sampledPower(:));
        peakErrorDb=10*log10(ffS.peakMag^2/max(sampledPeak,realmin));
        if integralErrorDb>1 || peakErrorDb>1
            error('PAD:ImportResolution',[ ...
                'This imported field has angular detail that the 1-degree theta / 2-degree phi ' ...
                'calculation grid cannot resolve (power integral differs by %.1f dB; ' ...
                'sampled peak is %.1f dB low). No pattern was loaded. ' ...
                'Use a solver that evaluates the CST data on a finer angular grid.'], ...
                integralErrorDb,peakErrorDb);
        end
        ffS.integralEfficiency=integral;
        ffS.eff=integral;
        if ~isfinite(integral) || integral<=0 || integral>1
            ffS.eff=NaN; % no clipping and no invented lossless efficiency
        end
    end

    function [Eth,Eph] = evalImportedFF(ffS, TH, PH)
        % The strict import covers the full sphere and includes periodic
        % azimuth copies. Out-of-domain data is an error, never zero-filled.
        if any(TH(:)<0 | TH(:)>180) || any(~isfinite(TH(:))) || any(~isfinite(PH(:)))
            error('PAD:FieldAngles','Field query angles are invalid.');
        end
        PH=mod(PH,360);
        Eth=(ffS.FReEth(TH,PH)+1j*ffS.FImEth(TH,PH))/ffS.peakMag;
        Eph=(ffS.FReEph(TH,PH)+1j*ffS.FImEph(TH,PH))/ffS.peakMag;
        if any(~isfinite(Eth(:)))||any(~isfinite(Eph(:)))
            error('PAD:FieldCoverage','Imported field has missing or invalid samples at the requested angles.');
        end
    end

    function [Eth,Eph] = patchEF(TH,PH,azArg,q)
        % Single-feed Patch (cos^q x lin pol) element pattern at an
        % arbitrary orientation azArg, factored out of elementFactor so
        % the dual-feed CP case can evaluate it twice (once per port,
        % 90 deg apart) without duplicating the formula.
        ct = max(cosd(TH),0);
        % Hemisphere mask applied EXPLICITLY rather than relying on
        % ct.^(q/2) to vanish behind the array. It does vanish for every
        % q > 0, but MATLAB evaluates 0^0 as 1, so at exactly q = 0 the
        % clamp is defeated and the element radiates backwards at full
        % strength. The discontinuity is the tell: q = 0 read 0.00 dBi
        % (omnidirectional) while q = 0.001 read 3.05 dBi (hemispherical).
        front = ct > 0;
        Eth = front.*(ct.^(q/2)).*cosd(PH-azArg);
        Eph = front.*(ct.^(q/2)).*ct.*(-sind(PH-azArg));
    end

    function fieldValue = formulaFieldValue(fieldValue,fieldSize)
        if ~(isnumeric(fieldValue)||islogical(fieldValue))
            error('PAD:FormulaType','A field formula must return numeric or logical values, not text or objects.');
        end
        fieldValue=double(fieldValue);
        if isscalar(fieldValue)
            fieldValue=repmat(fieldValue,fieldSize);
        elseif ~isequal(size(fieldValue),fieldSize)
            error('PAD:FormulaShape','A field formula must return a scalar or one value per angle.');
        end
        if any(~isfinite(fieldValue(:))) || any(abs(fieldValue(:))>1e100)
            error('PAD:FormulaValue','A field formula returned nonfinite or excessively large values.');
        end
    end
