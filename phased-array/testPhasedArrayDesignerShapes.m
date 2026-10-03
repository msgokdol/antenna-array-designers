function testPhasedArrayDesignerShapes
% Circle and Ellipse must make different arrays on the default square grid.
addpath(fileparts(mfilename('fullpath')),'-begin');
close(findall(0,'Type','figure'),'force');
cleanup = onCleanup(@()close(findall(0,'Type','figure'),'force')); %#ok<NASGU>
phasedArrayDesigner;
fig = findall(0,'Type','figure','Name','Planar Phased Array Designer');
table = findall(fig,'Tag','tblElements');
circleButton = findall(fig,'Tag','btnShapeCircle');
ellipseButton = findall(fig,'Tag','btnShapeEllipse');
circleButton.ButtonPushedFcn(circleButton,[]);
circleXY = table.Data(:,1:2);
ellipseButton.ButtonPushedFcn(ellipseButton,[]);
ellipseXY = table.Data(:,1:2);
assert(~isequal(circleXY,ellipseXY), ...
    'Circle and Ellipse populated identical default-grid cells.');
assert(max(ellipseXY(:,1))-min(ellipseXY(:,1)) > ...
    max(circleXY(:,1))-min(circleXY(:,1)), ...
    'Ellipse should visibly span more columns than Circle.');

% Changing aspect ratio must remask the oval even when M and N stay fixed.
dx = findall(fig,'Tag','spDx');
dx.Value = .75;
dx.ValueChangedFcn(dx,[]);
assert(size(table.Data,1) ~= size(ellipseXY,1), ...
    'Ellipse did not update when dx/dy changed.');

% A one-column grid has no 2D oval boundary and keeps its full line.
n = findall(fig,'Tag','spN');
n.Value = 1;
n.ValueChangedFcn(n,[]);
shape = findall(fig,'Tag','ddShape');
assert(size(table.Data,1)==8 && strcmp(shape.Value,'Ellipse'), ...
    'Ellipse should keep all elements of a one-column array.');
disp('Circle/Ellipse geometry and aspect-ratio checks passed.');
end
