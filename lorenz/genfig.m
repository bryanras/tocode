
% This is a script for generating a figure for use in the paper.
% I'm just sick of typing this shit over and over again.

% First, make the plot. Make sure the hold is off.
semilogy(kjjs(:,1),kjjs(:,2),'.-');

% Now change the properties.
ax=gca;
ylims=[4e3,3e5];
set(ax,'XLim',[13,25])
set(ax,'XDir','reverse')
set(ax,'XTick',[14.1,16.0,24.2]);
set(ax,'XGrid','on')
set(ax,'YLim',ylims)
set(ax,'YTickMode','auto')
set(ax,'YTickLabelMode','auto')
set(ax,'Box','off')

% Add some labels.
xlabel('\lambda','FontSize',12)
ylabel('Condition # of Jacobian, \kappa(J)','FontSize',12)

% Now add the stuff on the right.
%axt=axes;
%set(axt,'Color','none')

% Plot the vector of Delta-y's.
%hold on
%plot( [24.5,24.5],ylims, [24.2,24.2],ylims, [24.1,24.1],ylims,...
%	 	[16.0,16.0],ylims, [15.95,15.95],ylims, [14.1,14.1],ylims), ...
%		[14.08,14.08],ylims, [13.93,13.93],ylims, ...
%		[13.929,13.929],ylims, [13.923,13.923],ylims);

%hold off

% This is much the same as before.
%set(axt,'YScale','log')
%set(axt,'Xdir','reverse')
%set(axt,'XLim',[13,25])
%set(axt,'YLim',[0,1])
%set(axt,'YTick',[])


