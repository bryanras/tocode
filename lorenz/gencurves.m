
% Yet another figure generator.
colors = ['k', 'c', 'y', 'm', 'r', 'g', 'b'];
colors = [colors,colors,colors,colors];			% Kludge.  I don't care.

if 1==1
	plot3(huges(:,1),huges(:,2),huges(:,3),colors(1));
	hold on
	for ii=2:size(huges,2)/3
		plot3(huges(:,3*ii-2),huges(:,3*ii-1),huges(:,3*ii), ...
			colors(mod(ii-1,length(colors))+1))
	end

	hold off
	xlabel('x_1','FontSize',12);
	ylabel('x_2','FontSize',12);
	zlabel('x_3','FontSize',12,'Rotation',0);

else 
	plot(huges(:,2),huges(:,3),colors(1));
	hold on
	for ii=2:size(huges,2)/3
		plot(huges(:,3*ii-1),huges(:,3*ii),colors(ii))
	end

	hold off
	xlabel('x_2','FontSize',12);
	ylabel('x_3','FontSize',12,'Rotation',0);
	set(gca,'XTick',[-15:10:5])
	set(gca,'YTick',[-5:10:30])

end
