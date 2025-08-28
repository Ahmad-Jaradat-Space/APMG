function printjpg(h,outfilename,renderer)

if nargin == 2
    renderer = '-zbuffer';
end
pause(0.01)
set(h, 'PaperUnits','centimeters');
set(h, 'Units','centimeters');
pos=get(h,'Position');
set(h, 'PaperSize', [pos(3) pos(4)]);
set(h, 'PaperPositionMode', 'manual');
set(h, 'PaperPosition',[0 0 pos(3) pos(4)]);
pause(0.01)
print(h,'-djpeg','-r300',renderer,outfilename);
