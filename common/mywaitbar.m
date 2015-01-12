function mywaitbar(x,h,title)
if nargin < 1
    error('Input arguments not valid');
end
if nargin > 1
    hTitle = get(h,'title');
    set(hTitle,'String',title);
    set(h,'visible','on');
end
fractioninput = x;
x = max(0,min(100*x,100));
xline = [100 0 0 100 100];
yline = [0 0 1 1 0];
if fractioninput == 0    
    cla(h)
    pause(eps) 
    xpatch = [0 x x 0];
    ypatch = [0 0 1 1];
    patch(xpatch,ypatch,'b','EdgeColor','b','EraseMode','none','parent',h);
else
    p = findobj(gcf,'Type','patch');
    xpatch = [0 x x 0];
    set(p,'XData',xpatch,'parent',h);
    line(xline,yline,'EraseMode','none','color','k','parent',h);
end
if fractioninput==1
    p = findobj(gcf,'Type','patch');
    set(p,'XData',[]);
    set(h,'visible','off')
    line(xline,yline,'color',[1 1 1],'parent',h);
end
drawnow;