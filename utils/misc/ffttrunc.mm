function [shift,ph] = gridcorr(nav,k,gridsize)

	npts = size(nav,2);

	fprintf('\tgridcorr: gridding navigators\n');
	% grid navigators
	gridnav = grid_navs(nav,k,gridsize);

	fprintf('\tgridcorr: finding shift (center-of-mass)\n');
	% take centroid of gridded navigator to be shifted center of k-space
	absnav = abs(gridnav.');
	s = sum(absnav,2);

	t = repmat(1:gridsize,1,gridsize)';
	cx = absnav*t; 	cx = round(cx./s);
	t = repmat(1:gridsize,gridsize,1); t = t(:);
	cy = absnav*t;	cy = round(cy./s);

	shift = (cy - gridsize/2) + i*(cx - gridsize/2);
	% shift = cy + i*cx;

	fprintf('gridcorr: finding phase\n');
	in = 0:npts-1;
	in = in'*gridsize*gridsize + cy*gridsize + cx;
	ph = mod(phase(gridnav(in).'),2*pi);
