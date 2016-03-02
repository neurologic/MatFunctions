function [colvec colsize rowvec rowsize] = subplotinds(ncols,nitems,border)

% subplot('Position',positionVector) creates a new axes at the position ...
% specified by positionVector. The positionVector is a four-element ...
% vector of the form [left,bottom,width,height], such that the entries ...
% are normalized values between 0.0 to 1.0. If the position vector ...
% specifies an axes that overlaps any previous axes, then the new axes ...
% replaces the existing ones.

% border=0.03;

nrows=ceil(nitems/ncols);
colsize=(1-(border*(ncols+1)))/ncols;  
colind=[(border):(colsize+border):1-(colsize+border)];
colvec=reshape(repmat(colind,nrows,1),1,(nrows*ncols));
colvec=reshape(colvec,nrows,ncols);
rowsize=(1-(border*(nrows+1)))/nrows;
rowind=[(border):(rowsize+border):1-(rowsize+border)];
rowvec=repmat(rowind,1,ncols);
rowvec=flipud(reshape(rowvec,nrows,ncols));