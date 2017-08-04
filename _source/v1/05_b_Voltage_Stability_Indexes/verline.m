function verline(x,linetype,label);

    g=ishold(gca);
    hold on
    y=get(gca,'ylim');
    h=plot([x x],y,linetype);
    
    if exist('label')
        xx=get(gca,'xlim');
        xrange=xx(2)-xx(1);
        xunit=(x-xx(1))/xrange;
        if xunit<0.8
            text(x+0.01*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        else
            text(x-.05*xrange,y(1)+0.1*(y(2)-y(1)),label,'color',get(h,'color'))
        end
    end     