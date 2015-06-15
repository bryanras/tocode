

% Script file for generating some figures.  Set slc and scn your own damn self.

slc=[22*45+1:45^2,1:12*45];  % Or maybe not.

torit(huge(slc,scn),35,ptsj,1,1);
cl=centerline(huge(:,scn),ptsj);
hold on
plot3(cl(:,1),cl(:,2),cl(:,3));
hold off



