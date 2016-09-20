function sp = compute_bspline(knots,k,x,y,w)
%SPAP2 Least squares spline approximation.
%
%   SPAP2(KNOTS,K,X,Y)  returns the B-form of the least-squares approximation
%   f  to the data X, Y by splines of order K with knot sequence KNOTS.  
%   The spline approximates, at the data site X(j), the given data value
%   Y(:,j), j=1:length(X).
%   The data values may be scalars, vectors, matrices, or even ND-arrays.
%   Data points with the same site are replaced by their (weighted) average.
%   
%    f  is the spline of order K with knot sequence KNOTS for which
%
%   (*)    Y = f(X)
%
%   in the mean-square sense, i.e., the one that minimizes
%
%   (**)   sum_j W(j) |Y(:,j) - f(X(j))|^2 ,
%
%   with  W = ones(size(X)). Other weights can be specified by an optional
%   additional argument, i.e., by using
%
%   SPAP2(KNOTS,K,X,Y,W) which returns the spline  f  of order K with knot
%   sequence KNOTS that minimizes (**). A better choice than the default
%   W = ones(size(X))  would be the composite trapezoid rule weights 
%   W = ([dx;0]+[0;dx]).'/2 , with dx = diff(X(:))  (assuming that X is
%   strictly increasing).
%
%   If the data sites satisfy the Schoenberg-Whitney conditions
%
%   (***)   KNOTS(j) < X(j) < KNOTS(j+K) ,
%                               j=1:length(X)=length(KNOTS)-K ,
%
%   (with equality permitted at knots of multiplicity K), then  f  is
%   the unique spline of that order satisfying  (***)  exactly.  No
%   spline is returned unless (***) is satisfied for some subsequence
%   of X.
%   
%   Since it might be difficult to supply such a knot sequence for the given
%   data sites, it is also possible to specify KNOTS as a positive integer,
%   in which case, if possible, a knot sequence will be supplied that satisfies
%   (***) for some subsequence of X and results in a spline consisting of KNOTS 
%   polynomial pieces.
%
%   If Y is a matrix or, more generally, an ND array, of size [d1,...,ds,n] say,
%   then Y(:,...,:,j) is the value being approximated at X(j), and the
%   resulting spline is correspondingly [d1,...,ds]-valued. In that case, the
%   expression  |Y(:,j) - f(X(j))|^2  in the error measure (**) is meant as
%   the sum of squares of all the d1*...*ds entries of  Y(:,j)-f(X(j)) .
%   
%   It is also possible to fit to gridded data:
%
%   SPAP2( {KNOTS1,...,KNOTSm}, [K1,...,Km], {X1,...,Xm}, Y ) returns
%   the m-variate tensor-product spline of coordinate order Ki and with knot 
%   sequence KNOTSi in the i-th variable, i=1,...,m, for which
%
%   Y(:,...,:,i1,...,im) = f(X1(i1),...,Xm(im)),  all i := (i1,...,im) 
%
%   in the (possibly weighted) mean-square sense.
%   Note the possibility of fitting to vector-valued and even ND-valued data.
%   However, in contrast to the univariate case, if the data to be fitted are
%   scalar-valued, then the input array Y is permitted to be m-dimensional,
%   in which case
%   Y(i1,...,im) = f(X1(i1),...,Xm(im)),  all i := (i1,...,im) 
%   in the (possibly weighted) mean-square sense.
%
%   Example 1:
%
%      spap2(augknt(x([1 end]),2),2,x,y);
%
%   provides the least-squares straight-line fit to data x,y, assuming that
%   all the sites x(j) lie in the interval [x(1) .. x(end)], while
%
%      spap2(1,2,x,y);
%
%   accomplishes this without that assumption, and, with that assumption,
%
%      w = ones(size(x)); w([1 end]) = 100;
%      spap2(1,2,x,y,w);
%
%   forces that fit to come very close to the first and last data point.
%
%   Example 2: The statements
%
%      x = -2:.2:2; y=-1:.25:1; [xx, yy] = ndgrid(x,y); 
%      z = exp(-(xx.^2+yy.^2)); 
%      sp = spap2({augknt([-2:2],3),2},[3 4],{x,y},z);
%      fnplt(sp)
%
%   produce the picture of an approximant to a bivariate function. 
%   Use of MESHGRID instead of NDGRID here would produce an error.
%
%   See also SPAPI, SPAPS.

%   Copyright 1987-2008 The MathWorks, Inc.
%   $Revision: 1.21.4.3 $

if nargin<5, w = []; end

if iscell(knots) % gridded data are to be fitted by tensor product splines

   if ~iscell(x)
      error('SPLINES:SPAP2:Xnotcell', ...
            'If KNOTS is a cell-array, then also X must be one.')
   end
   m = length(knots);
   if m~=length(x)
      error('SPLINES:SPAP2:wrongsizeX', ...
            'If KNOTS is a cell-array, then X must be one of the same length.')
   end
   sizey = size(y);
   if length(sizey)<m
     error('SPLINES:SPAP2:wrongsizeY', ...
          ['If KNOTS is a cell-array of length m, then Y must have', ...
            ' at least m dimensions.'])
   end

   if length(sizey)==m,  % grid values of a scalar-valued function
     if issparse(y), y = full(y); end 
     sizey = [1 sizey]; 
   end

   sizeval = sizey(1:end-m); sizey = [prod(sizeval), sizey(end-m+(1:m))];
   y = reshape(y, sizey); 

   if iscell(k), k = cat(2,k{:}); end
   if length(k)==1, k = repmat(k,1,m); end
   if isempty(w), w = cell(1,m); end
   
   v = y; sizev = sizey;
   for i=m:-1:1   % carry out coordinatewise least-squares fitting
      [knots{i},v,sizev(m+1),k(i)] = spbrk(spap21(knots{i}, k(i), x{i}, ...
                      reshape(v,prod(sizev(1:m)),sizev(m+1)),w{i}));
      v = reshape(v,sizev);
      if m>1
         v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
      end
   end
   % At this point, V contains the tensor-product B-spline coefficients.
   % It remains to put together the spline:
   sp = spmak(knots, v, sizev);
   if length(sizeval)>1, sp = fnchg(sp,'dz',sizeval); end

else             % univariate spline interpolation
   sp = spap21(knots,k,x,y,w);
end

function sp = spap21(knots,k,x,y,w)
%SPAP21 Univariate least squares spline approximation.

if isempty(w), [x,y,sizeval]   = chckxywp(x,y,1);
else           [x,y,sizeval,w] = chckxywp(x,y,1,w);
end
nx = length(x);

if length(knots)==1 % we are to use a spline with KNOTS pieces
   k = min(k,nx); maxpieces = nx-k+1;
   if knots<1||knots>maxpieces
      warning('SPLINES:SPAP2:wrongknotnumber', ...
      ['The number of polynomial pieces to be used must be\n',...
       'positive but not larger than length(x)-k+1 = %g.\n'], maxpieces)
      knots = max(1,min(maxpieces,knots));
   end
   if knots==1&&k==1
      if nx==1, knots = [x(1) x(1)+1];
      else      knots = x([1 end]).';
      end
   else
      knots = create_knots_acceptable(x(round(linspace(1,nx,knots-1+k))),k);
   end
end

%  Generate the collocation matrix and divide it into the possibly reordered
%  sequence of given values to generate the B-spline coefficients of the
%  interpolant, then put it all together into SP. But trap any error from
%  SLVBLK, in order to provide a more helpful error message.

try
   if isempty(w)
      sp = spmak(knots,slvblk(spcol(knots,k,x,'slvblk','noderiv'),y).');
   else
      sp = spmak(knots,slvblk(spcol(knots,k,x,'slvblk','noderiv'),y,w).');
   end
catch laster
  if ~isempty(findstr('SLVBLK',laster.identifier))
     error('SPLINES:SPAP2:noSWconds',...
     ['\nThe given knots and order are incompatible with the given', ...
             ' data sites.\n', ...
       'No subsequence of the data sites satisfies the',...
       ' Schoenberg-Whitney\n',...
       '   conditions of the given order wrto the given knots.\n', ...
       'If you have trouble coming up with a good knot sequence for the', ...
       ' data\n', ...
       '   sites, try using spap2(l,k,x,y), with l the desired number of\n',...
       '   polynomial pieces in the spline fit.'],0)
  else
     error('SPLINES:SPAP2:lasterr',laster.message)
  end
end
if length(sizeval)>1, sp = fnchg(sp,'dz',sizeval); end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISSING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function 1
function [x,y,sizeval,w,origint,p,tolred] = chckxywp(x,y,nmin,w,p,adjtol)
% make sure X is a vector:
if iscell(x)||length(find(size(x)>1))>1
   error('SPLINES:CHCKXYWP:Xnotvec','X must be a vector.'), end

% make sure X is real:
if ~all(isreal(x))
   x = real(x);
   warning('SPLINES:CHCKXYWP:Xnotreal', ...
           'Imaginary part of complex data sites ignored.')
end

% deal with NaN's and Inf's among the sites:
nanx = find(~isfinite(x));
if ~isempty(nanx)
   x(nanx) = [];
   warning('SPLINES:CHCKXYWP:NaNs', ...
           'All data points with NaN or Inf as their site will be ignored.')
end

n = length(x);
if nargin>2&&nmin>0, minn = nmin; else minn = 2; end
if n<minn
   error('SPLINES:CHCKXYWP:toofewpoints', ...
   'There should be at least %g data sites.',minn), end

% re-sort, if needed, to ensure nondecreasing site sequence:
tosort = false;
if any(diff(x)<0), tosort = true; [x,ind] = sort(x); end

nstart = n+length(nanx);
% if Y is ND, reshape it to a matrix by combining all dimensions but the last:
sizeval = size(y);
yn = sizeval(end); sizeval(end) = []; yd = prod(sizeval);
if length(sizeval)>1
   y = reshape(y,yd,yn);
else
   % if Y happens to be a column matrix, of the same length as the original X,
   % then change Y to a row matrix
   if yn==1&&yd==nstart
      yn = yd; y = reshape(y,1,yn); yd = 1; sizeval = yd;
   end
end
y = y.'; x = reshape(x,n,1);

% make sure that sites, values and weights match in number:

if nargin>2&&~nmin % in this case we accept two more data values than
                   % sites, stripping off the first and last, and returning
		   % them separately, in W, for use in CSAPE1.
   switch yn
   case nstart+2, w = y([1 end],:); y([1 end],:) = [];
      if ~all(isfinite(w)),
         error('SPLINES:CHCKXYWP:InfY', ...
	 'Some of the end condition values fail to be finite.')
      end
   case nstart, w = [];
   otherwise
      error('SPLINES:CHCKXYWP:XdontmatchY', ...
           ['The number of sites, %g, does not match the number of', ...
          ' values, %g.'], nstart, yn)
   end
else
   if yn~=nstart
      error('SPLINES:CHCKXYWP:XdontmatchY', ...
           ['The number of sites, %g, does not match the number of', ...
          ' values, %g.'], nstart, yn)
   end
end

nonemptyw = nargin>3&&~isempty(w);
if nonemptyw
   if length(w)~=nstart
      error('SPLINES:CHCKXYWP:weightsdontmatchX', ...
       ['The number of weights, %g, does not match the number of', ...
       ' sites, %g.'], length(w), nstart)
   else
      w = reshape(w,1,nstart);
   end
end

roughnessw = exist('p','var')&&length(p)>1;
if roughnessw
   if tosort
      warning('SPLINES:CHCKXYWP:cantreorderrough', ...
           'Since data sites are not ordered, roughness weights are ignored.')
      p = p(1);
   else
      if length(p)~=nstart
         error('SPLINES:CHCKXYWP:rweightsdontmatchX', ...
	 ['The number of roughness weights is incompatible with the', ...
	        ' number of sites, %g.'], nstart)
      end
   end
end

%%% remove values and error weights corresponding to nonfinite sites:
if ~isempty(nanx), y(nanx,:) = []; if nonemptyw, w(nanx) = []; end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(nanx,2)) = [];
   end
end
if tosort, y = y(ind,:); if nonemptyw, w = w(ind); end, end

% deal with nonfinites among the values:
nany = find(sum(~isfinite(y),2));
if ~isempty(nany)
   y(nany,:) = []; x(nany) = []; if nonemptyw, w(nany) = []; end
   warning('SPLINES:CHCKXYWP:NaNs', ...
           'All data points with NaNs or Infs in their value will be ignored.')
   n = length(x);
   if n<minn
      error('SPLINES:CHCKXYWP:toofewX', ...
      'There should be at least %g data sites.',minn), end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(nany,2)) = [];
   end
end

if nargin==3&&nmin, return, end % for SPAPI, skip the averaging

if nargin>3&&isempty(w) %  use the trapezoidal rule weights:
   dx = diff(x);
   if any(dx), w = ([dx;0]+[0;dx]).'/2;
   else,       w = ones(1,n);
   end
   nonemptyw = ~nonemptyw;
end

tolred = 0;
if ~all(diff(x)) % conflate repeat sites, averaging the corresponding values
                 % and summing the corresponding weights
   mults = knt2mlt(x);
   for j=find(diff([mults;0])<0).'
      if nonemptyw
         temp = sum(w(j-mults(j):j));
	 if nargin>5
	    tolred = tolred + w(j-mults(j):j)*sum(y(j-mults(j):j,:).^2,2); 
	 end
         y(j-mults(j),:) = (w(j-mults(j):j)*y(j-mults(j):j,:))/temp;
         w(j-mults(j)) = temp;
         if nargin>5
	    tolred = tolred - temp*sum(y(j-mults(j),:).^2);
	 end
      else
         y(j-mults(j),:) = mean(y(j-mults(j):j,:),1);
      end
   end
      
   repeats = find(mults);
   x(repeats) = []; y(repeats,:) = []; if nonemptyw, w(repeats) = []; end
   if roughnessw  % as a first approximation, simply ignore the
                  % specified weight to the left of any ignored point.
      p(max(repeats,2)) = [];
   end
   n = length(x);
   if n<minn, error('SPLINES:CHCKXYWP:toofewX', ...
      'There should be at least %g data sites.',minn), end
end

if nargin<4, return, end


% remove all points corresponding to relatively small weights (since a
% (near-)zero weight in effect asks for the corresponding datum to be dis-
% regarded while, at the same time, leading to bad condition and even
% division by zero).
origint = []; % this will be set to x([1 end]).' in case the weight for an end
             % data point is near zero, hence the approximation is computed
             % without that endpoint.
if nonemptyw
   ignorep = find( w <= (1e-13)*max(abs(w)) );
   if ~isempty(ignorep)
      if ignorep(1)==1||ignorep(end)==n, origint = x([1 end]).'; end
      x(ignorep) = []; y(ignorep,:) = []; w(ignorep) = []; 
      if roughnessw
                     % as a first approximation, simply ignore the
                     % specified weight to the left of any ignored point.
         p(max(ignorep,2)) = [];
      end
      n = length(x);
      if n<minn
        error('SPLINES:CHCKXYWP:toofewX', ...
	     ['There should be at least %g data points with positive',...
	       ' weights.'],minn)
      end
   end
end


%% Function 2
function colloc = spcol(knots,k,tau,varargin)
if ~isempty(find(diff(knots)<0))
   error('SPLINES:SPCOL:knotsdecreasing',...
   'The knot sequence KNOTS should be nondecreasing.')
end
if ~isempty(find(diff(tau)<0))
   error('SPLINES:SPCOL:TAUdecreasing', ...
   'The point sequence TAU should be nondecreasing.')
end

%  Compute the number  n  of B-splines of order K supported by the given
%  knot sequence and return an empty matrix in case there aren't any.

npk=length(knots); n=npk-k;
if n<1, warning('SPLINES:SPCOL:noBsplines', ...
                'There are no B-splines for the given input.')
   colloc = zeros(length(tau),0); return
end

% Settle the options:
slvblk=0; noderiv=0;
for j=4:nargin
   argj = varargin{j-3};
   if ~isempty(argj)
      if ischar(argj)
         if     argj(1)=='s', slvblk=1;
            if length(argj)>1, if argj(2)=='p', slvblk=2; end, end
         elseif argj(1)=='n', noderiv=1;
         else error('SPLINES:SPCOL:wronginarg2',...
	 ['The second optional argument should be ''sl'', ',...
         '''sp'', ''n'', or a number.'])
         end
      else 
         switch j  % for backward compatibility
         case 4,  slvblk=1; 
         case 5,  noderiv=1;
         end
      end
   end
end

% If  NODERIV==0, remove all multiplicities from TAU and generate repetitions
% of rows instead of rows containing values of successive derivatives.
nrows = length(tau); tau = reshape(tau,1,nrows);
if noderiv
   index = 1:nrows; m = ones(1,nrows); nd = 1; pts = tau;
else
   index = [1 find(diff(tau)>0)+1];
   m = diff([index nrows+1]); nd = max(m);
   if nd>k
      error('SPLINES:SPCOL:multtoohigh',...
      'Point multiplicity should not exceed the given order %g.',k);
   end
   pts = tau(index);
end

%set some abbreviations
km1 = k-1;

%  augment knot sequence to provide a K-fold knot at each end, in order to avoid
% struggles near ends of basic interval,  [KNOTS(1) .. KNOTS(npk)] .
% The resulting additional B-splines, if any, will NOT appear in the output.

[augknot,addl] = create_knots_augment(knots,k); naug = length(augknot)-k;
pts = pts(:); augknot = augknot(:);

%  For each  i , determine  savl(i)  so that  K <= savl(i) < naug+1 , and,
% within that restriction,
%        augknot(savl(i)) <= pts(i) < augknot(savl(i)+1) .

savl = max(sorted(augknot(1:naug),pts), k);

b = zeros(nrows,k);

% first do those without derivatives
index1 = find(m==1);
if ~isempty(index1)
   pt1s = pts(index1); savls = savl(index1); lpt1 = length(index1);
   % initialize the  b  array.
   lpt1s = index(index1); b(lpt1s,1) = ones(lpt1,1);

   % run the recurrence simultaneously for all  pt1(i) .
   for j=1:km1
      saved = zeros(lpt1,1);
      for r=1:j
         tr = augknot(savls+r)-pt1s;
         tl = pt1s-augknot(savls+r-j);
         term = b(lpt1s,r)./(tr+tl);
         b(lpt1s,r) = saved+tr.*term;
         saved = tl.*term;
      end
      b(lpt1s,j+1) = saved;
   end
end

% then do those with derivatives, if any:
if nd>1
   indexm=find(m>1);ptss=pts(indexm);savls=savl(indexm);lpts=length(indexm);
   % initialize the  bb  array.
   %temp = [1 zeros(1,km1)]; bb = temp(ones(nd*lpts,1),:);
   bb = repmat([1 zeros(1,km1)],nd*lpts,1);
   lptss = nd*[1:lpts];

   % run the recurrence simultaneously for all  pts(i) .
   % First, bring it up to the intended level:
   for j=1:k-nd
      saved = zeros(lpts,1);
      for r=1:j
         tr = augknot(savls+r)-ptss;
         tl = ptss-augknot(savls+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptss,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptss,j+1) = saved;
   end

   % save the B-spline values in successive blocks in  bb .

   for jj=1:nd-1
      j = k-nd+jj; saved = zeros(lpts,1); lptsn = lptss-1;
      for r=1:j
         tr = augknot(savls+r)-ptss;
         tl = ptss-augknot(savls+r-j);
         term = bb(lptss,r)./(tr+tl);
         bb(lptsn,r) = saved+tr.*term;
         saved = tl.*term;
      end
      bb(lptsn,j+1) = saved; lptss = lptsn;
   end

   % now use the fact that derivative values can be obtained by differencing:

   for jj=nd-1:-1:1
      j = k-jj;
      temp = repmat([jj:nd-1].',1,lpts)+repmat(lptsn,nd-jj,1); lptss=temp(:);
      for r=j:-1:1
         temp = repmat((augknot(savls+r)-augknot(savls+r-j)).'/j,nd-jj,1);
         bb(lptss,r) = -bb(lptss,r)./temp(:);
         bb(lptss,r+1) = bb(lptss,r+1) - bb(lptss,r);
      end
   end

   % finally, combine appropriately with  b  by interspersing the multiple
   % point conditions appropriately:
   dtau = diff([tau(1)-1 tau(:).' tau(nrows)+1]);
   index=find(min(dtau(2:nrows+1),dtau(1:nrows))==0); % Determines all rows
                                                    % involving multiple tau.
   dtau=diff(tau(index));index2=find(dtau>0)+1;     % We need to make sure to
   index3=[1 (dtau==0)];                            % skip unwanted derivs:
   if ~isempty(index2)
             index3(index2)=1+nd-m(indexm(1:length(indexm)-1));end
   b(index,:)=bb(cumsum(index3),:);

   % ... and appropriately enlarge  savl
   index = cumsum([1 (diff(tau)>0)]);
   savl = savl(index);
end

% Finally, zero out all rows of  b  corresponding to TAU outside the basic
% interval,  [knots(1) .. knots(npk)] .

index = find(tau<knots(1)|tau>knots(npk));
if ~isempty(index)
   b(index,:) = 0;
end

% The first B-spline of interest begins at KNOTS(1), i.e., at  augknot(1+addl)
% (since  augknot's  first knot has exact multiplicity K). If  addl<0 ,
% this refers to a nonexistent index and means that the first  -addl  columns
% of the collocation matrix should be trivial.  This we manage by setting
savl = savl+max(0,-addl);

if slvblk     % return the collocation matrix in almost block diagonal form.
              % For this, make the blocks out of the entries with the same
              %  SAVL(i) , with  LAST  computed from the differences.
   % There are two issues, the change in the B-splines considered because of
   % the use of  AUGKNOT  instead of  KNOTS , and the possible drop of B-splines
   % because the extreme  TAU  fail to involve the extreme knot intervals.

   % SAVL(j) is the index in  AUGKNOT  of the left knot for  TAU(j) , hence the
   % corresponding row involves  B-splines to index  savl(j) wrto augknot, i.e.,
   % B-splines to index  savl(j)-addl  wrto  KNOTS.
   % Those with negative index are removed by cutting out their columns (i.e.,
   % shifting appropriately the blocks in which they lie). Those with index
   % greater than  n  will be ignored because of  last .

   last0 = max(0,savl(1)-max(0,addl)-k); % number of cols in trivial first block
   if addl>0   % if B-splines were added on the left, remove them now:
      width = km1+k;cc = zeros(nrows*width,1);
      index = min(k,savl-addl); 
      temp = +repmat(nrows*[0:km1],nrows,1);
    cc(repmat(([1-nrows:0]+nrows*index).',1,k)+repmat(nrows*[0:km1],nrows,1))=b;
      b(:)=cc(repmat([1-nrows:0].',1,k)+repmat(nrows*(k+[0:km1]),nrows,1));
      savl=savl+k-index;
   end
   ds=(diff(savl));
   index=[0 find(ds>0) nrows];
   rows=diff(index);
   nb=length(index)-1;
   last=ds(index(2:nb));
   if addl<0  nb=nb+1; rows=[0 rows]; last=[last0 last]; end
   if slvblk==1
      colloc=[41 nb rows k last n-sum(last) b(:).'];
   else   % return the equivalent sparse matrix (cf BKBRK)
      nr = (1:nrows).'; nc = 1:k; nrnc = nrows*k;
      ncc = zeros(1,nrows); ncc(1+cumsum(rows(1:(nb-1)))) = last;
      ncc(1) = last0; ncc = reshape(cumsum(ncc),nrows,1);
      ijs = [reshape(repmat(nr,1,k),nrnc,1), ...
           reshape(repmat(ncc,1,k)+repmat(nc,nrows,1), nrnc,1), ...
           reshape(b,nrnc,1)];
      index = find(ijs(:,2)>n);
      if ~isempty(index), ijs(index,:) = []; end
      colloc = sparse(ijs(:,1),ijs(:,2),ijs(:,3),nrows,n);
   end
else          % return the collocation matrix in standard matrix form
   width = max([n,naug])+km1+km1;
   cc = zeros(nrows*width,1);
   cc(repmat([1-nrows:0].',1,k)+ ...
              repmat(nrows*savl.',1,k)+repmat(nrows*[-km1:0],nrows,1))=b;
   % (This uses the fact that, for a column vector  v  and a matrix  A ,
   %  v(A)(i,j)=v(A(i,j)), all i,j.)
   colloc = reshape(cc(repmat([1-nrows:0].',1,n) + ...
                    repmat(nrows*(max(0,addl)+[1:n]),nrows,1)), nrows,n);
end


%% Function 3
function pointer = sorted(meshsites, sites)
[ignored,index] = sort([meshsites(:).' sites(:).']);
pointer = find(index>length(meshsites))-(1:length(sites));

%% Function 4
function x = slvblk(blokmat,b,w)
% If BLOKMAT is sparse, handle the problem sparsely:
if issparse(blokmat)
   if nargin>2&&~isempty(w)
      n = length(w); spw = sparse(1:n,1:n,sqrt(w));
      x = (spw*blokmat)\(spw*b);
   else
      x = blokmat\b;
   end
   return
end

% get the basic information
[nb,rows,ncols,last,blocks] = bkbrk(blokmat);

ne = sum(rows);nu = sum(last);
if any(cumsum(rows)<cumsum(last))||any(last>ncols)
   error('SPLINES:SLVBLK:matrixnot11', ...
   'The coefficient matrix has a nontrivial nullspace.')
end

[brow,bcol] = size(b);
if(ne~=brow)
   error('SPLINES:SLVBLK:wrongrightside',...
   'Matrix and right side are incompatible.')
end

blocks = [blocks b];
ccols = ncols+bcol;
if nargin>2, w = sqrt(w); blocks = repmat(w(:),1,ccols).*blocks; end

f = 1; l = 0; elim = 0;
for j=1:nb
   if (f<=l) % shift the rows still remaining from previous block
      blocks(f:l,:) = ...
         [blocks(f:l,elim+1:ncols) zeros(l+1-f,elim),blocks(f:l,ncols+1:ccols)];
   end
   l = l+rows(j);

   elim = last(j);
   % ideally, one would now use
   %   [q,r] = qr(blocks(f:l,1:elim));
   % followed up by
   %   blocks(f:l,:) = q'*blocks(f:l,:);
   %   f = f+elim;
   % but, unfortunately, this generates the possibly very large square matrix q
   % The unhappy alternative is to do the elimination explicitly here, using
   % Householder reflections (and an additional inner loop):
   for k=1:elim
      a = norm(blocks(f:l,k));
      vv = abs(blocks(f,k))+a;
      c = vv*a;
      if blocks(f,k)<0, vv = -vv; end
      q = [vv;blocks(f+1:l,k)];
      blocks(f:l,:) = ...
       blocks(f:l,:)-repmat(q/c,1,ccols).*repmat(q'*blocks(f:l,:),l+1-f,1);
       %blocks(f:l,:)-((q/c)*ones(1,ccols)).*(ones(l+1-f,1)*(q'*blocks(f:l,:)));
      f = f+1;
   end
end

% now we are ready for back-substitution
x = zeros(f-elim-1+ncols,bcol);

for j=nb:-1:1
   elim = last(j); l = f-1; f = f-elim;
   % here is another occasion where empty matrices of various sizes would help;
   % instead, use an if statement:
   if elim<ncols, blocks(f:l,ncols+1:ccols) = blocks(f:l,ncols+1:ccols) ...
                    - blocks(f:l,elim+1:ncols)*x(f-1+[elim+1:ncols],:); end
   x(f:l,:) = blocks(f:l,1:elim) \ blocks(f:l,ncols+1:ccols);
end
x = x(1:nu,:);


%% Function 5
function [nbo,rows,ncols,last,blocks] = bkbrk(blokmat)
if blokmat(1)==41 % data type number for the spline block format is 41
   % Here are the details of this particular sparse format:
   % The matrix is sum(ROWS)-by-sum(LAST).
   % There are NB blocks. The i-th block has ROWS(i) rows and NCOLS columns.
   % The first column of the (i+1)st block is exactly LAST(i) columns to the
   % right of the first column of the i-th block.
   nb = blokmat(2);
   rows = blokmat(2+[1:nb]);
   ncols = blokmat(2+nb+1);
   last = blokmat(3+nb+[1:nb]);
   blocks = reshape(blokmat(3+2*nb+[1:sum(rows)*ncols]),sum(rows),ncols);

elseif blokmat(1)==40 % data type number for general almost block diagonal
                      % format is 40;
   nb = blokmat(2);
   rows = blokmat(2+[1:nb]);
   cols = blokmat(2+nb+[1:nb]);
   last = blokmat(2+2*nb+[1:nb]);
   row = cumsum([0,rows]);
   ne = sum(rows);ncols = max(cols);
   len = rows.*cols;
   index = cumsum([2+3*nb len]);
   blocks = zeros(ne,ncols);
   for j=1:nb
      block = reshape(blokmat(index(j)+[1:len(j)]),rows(j),cols(j));
      blocks(row(j)+[1:row(j+1)],[1:cols(j)]) = block;
   end
else
   error('SPLINES:BKBRK:unknownarg', ...
         'The argument does not appear to be an almost block diagonal matrix.')
end

if nargout==0 % print out the blocks
   if blokmat(1)==41 % generate COLS
      temp = cumsum([0 last]); temp = temp(nb+1)-temp;
      cols = min(temp(1:nb),ncols);
   end
   rowsum = cumsum([0 rows]);
   for j=1:nb
      fprintf(['block ',int2str(j),' has ',int2str(rows(j)),' row(s)\n'])
      disp(blocks(rowsum(j)+[1:rows(j)],1:cols(j)))
    fprintf([' next block is shifted over ',int2str(last(j)),' column(s)\n\n'])
   end
else
   nbo = nb;
end


%% Function 6
function spline = spmak(knots,coefs,sizec)
if nargin==0;
   knots = input('Give the vector of knots  >');
   coefs = input('Give the array of B-spline coefficients  >');
end

if nargin>2
   if numel(coefs)~=prod(sizec)
     error('SPLINES:SPMAK:coefsdontmatchsize', ...
           'The coefficient array is not of the explicitly specified size.')
   end
else
   if isempty(coefs)
      error('SPLINES:SPMAK:emptycoefs','The coefficient array is empty.')
   end
   sizec = size(coefs);
end

m = 1; if iscell(knots), m = length(knots); end
if length(sizec)<m
   error('SPLINES:SPMAK:coefsdontmatchknots', ...
        ['According to KNOTS, the function is %g-dimensional;\n',...
          'hence COEFS must be at least %g-dimensional.'],m,m)
end
if length(sizec)==m,  % coefficients of a scalar-valued function
   sizec = [1 sizec];
end

% convert ND-valued coefficients into vector-valued ones, retaining the
% original size in SIZEVAL, to be stored eventually in SP.DIM .
sizeval = sizec(1:end-m); sizec = [prod(sizeval), sizec(end-m+(1:m))];
coefs = reshape(coefs, sizec);

if iscell(knots), % we are putting together a tensor-product spline
   [knots,coefs,k,sizec] = chckknt(knots,coefs,sizec);
else            % we are putting together a univariate spline
   [knots,coefs,k,sizec] = chckknt({knots},coefs,sizec); knots = knots{1};
end

spline.form = 'B-';
spline.knots = knots;
spline.coefs = coefs;
spline.number = sizec(2:end);
spline.order = k;
spline.dim = sizeval;
% spline = [11 d n coefs(:).' k knots(:).'];

function [knots,coefs,k,sizec] = chckknt(knots,coefs,sizec)
%CHCKKNT check knots, omit trivial B-splines

for j=1:length(sizec)-1
   n = sizec(j+1); k(j) = length(knots{j})-n;
   if k(j)<=0, error('SPLINES:SPMAK:knotsdontmatchcoefs', ...
                     'There should be more knots than coefficients.'), end
   if any(diff(knots{j})<0)
      error('SPLINES:SPMAK:knotdecreasing',...
      'The knot sequence should be nondecreasing.')
   end
   if knots{j}(1)==knots{j}(end)
      error('SPLINES:SPMAK:extremeknotssame',...
      'The extreme knots should be different.')
   end

   % make sure knot sequence is a row matrix:
   knots{j} = reshape(knots{j},1,n+k(j));
   % throw out trivial B-splines:
   index = find(knots{j}(k(j)+(1:n))-knots{j}(1:n)>0);
   if length(index)<n
      oldn = n; n = length(index);
      knots{j} = reshape(knots{j}([index oldn+(1:k(j))]),1,n+k(j));
      coefs = ...
          reshape(coefs, [prod(sizec(1:j)),sizec(j+1),prod(sizec(j+2:end))]);
      sizec(j+1) = n; coefs = reshape(coefs(:,index,:),sizec);
   end
end


%% Function 7
function varargout = spbrk(sp,varargin)
if ~isstruct(sp)
  if sp(1)~=11&&sp(1)~=12
     error('SPLINES:SPBRK:fnotBform', ...
     'The input array does not seem to describe a function in B-form.')
  else
     di=sp(2);ni=sp(3);
     ci=reshape(sp(3+(1:di*ni)),di,ni);
     kk=sp(4+di*ni);ki=sp(4+di*ni+(1:kk+ni));
     sp = spmak(ki,ci);
  end
end

if length(sp.form)~=2||sp.form(1)~='B'
   error('SPLINES:SPBRK:snotBform',...
   'The input does not seem to describe a spline in B-form.')
end
if nargin>1 % we have to hand back one or more parts
   lp = max(1,nargout); % SPBRK(SP,PART) may be part of an expression
   if lp>length(varargin)
      error('SPLINES:SPBRK:moreoutthanin',...
            'Too many output arguments for the given input.')
   end
   varargout = cell(1,lp);
   for jp=1:lp
      part = varargin{jp};
      if ischar(part)
         if isempty(part)
	    error('SPLINES:SPBRK:partemptystr',...
	    'Part specification should not be an empty string.')
	 end
         switch part(1)
         case 'f',       out1 = [sp.form,'form'];
         case 'd',       out1 = sp.dim;
         case 'n',       out1 = sp.number;
         case {'k','t'}, out1 = sp.knots;
         case 'o',       out1 = sp.order;
         case 'c',       out1 = sp.coefs;
	 case 'v',       out1 = length(sp.order);
         case 'i', % this must be treated differently in multivariate case
            if iscell(sp.knots)
               for i=length(sp.knots):-1:1  % loop backward to avoid redef.
                  out1{i} = sp.knots{i}([1 end]);
               end
            else
               out1 = sp.knots([1 end]);
            end
	 case 'b', % this must be treated differently in multivariate case
	    if iscell(sp.knots)
               for i=length(sp.knots):-1:1  % loop backward to avoid redef.
                  out1{i} = knt2brk(sp.knots{i});
               end
            else
               out1 = knt2brk(sp.knots);
            end
         otherwise
            error('SPLINES:SPBRK:wrongpart',...
	    ['''',part,''' is not part of a B-form.'])
         end
      elseif isempty(part)
	 out1 = sp;
      else
         if iscell(part)  % we must be dealing with a tensor-product spline
            c = sp.coefs; knots = sp.knots; m = length(knots);
            sizec = size(c);
            if length(sizec)~=m+1 % trouble because of trailing singleton dims
               sizec = [sp.dim,sp.number]; c = reshape(c,sizec);
            end
            for i=m:-1:1
               dd = prod(sizec(1:m));
               spi = spcut(spmak(knots{i},reshape(c,dd,sp.number(i))), part{i});
               knots{i} = spi.knots; sizec(m+1) = spi.number;
               c = reshape(spi.coefs,sizec);
               if m>1
                  c = permute(c,[1,m+1,2:m]);
                  sizec(2:m+1) = sizec([m+1,2:m]);
               end
            end
            out1 = spmak(knots,c,sizec);

         else             % we must be dealing with a univariate spline
            out1 = spcut(sp,part);
         end
      end
      varargout{jp} = out1;
   end
else
   if nargout==0
     if iscell(sp.knots) % we have a multivariate spline and, at present,
                         % I can't think of anything clever to do; so...
       disp(sp)
     else
       disp('knots(1:n+k)'),disp(sp.knots),
       disp('coefficients(d,n)'),disp(sp.coefs),
       disp('number n of coefficients'),disp(sp.number),
       disp('order k'),disp(sp.order),
       disp('dimension d of target'),disp(sp.dim),
     end
   else
    varargout = {sp.knots,sp.coefs, sp.number, sp.order, sp.dim};
   end
end
function out1 = spcut(sp,interv)
%SPCUT change the basic interval

if isempty(interv)||ischar(interv), out1 = sp; return, end

sizei = size(interv);
if sizei(2)>1 % we are to change the basic interval
   tl = interv(1,1); tr = interv(1,2);
   if tl==tr
      warning('SPLINES:SPBRK:SPCUT:trivialinterval', ...
	         'No changes made since the given end points are equal.')
      out1 = sp; return
   end
   if tl>tr, tl = tr; tr = interv(1); end

   index = sorted(sp.knots,[tl,tr]); mults = knt2mlt(sp.knots);
   if tl<sp.knots(1),      m1 = 1;
   elseif tl==sp.knots(1), m1 = 0;
   else                    m1 = sp.order;
      if tl==sp.knots(index(1)), m1 = m1-mults(index(1))-1; end
   end
   if tr>sp.knots(end),      m2 = 1;
   elseif tr==sp.knots(end), m2 = 0;
   else                      m2 = sp.order;
      if tr==sp.knots(index(2)), m2 = m2-mults(index(2))-1; end
   end
   sp = fnrfn(sp, [repmat(tl,1,m1),repmat(tr,1,m2)]);
   index = sorted(sp.knots,[tl tr]);
   if sp.knots(end)>tr
      sp = spmak(sp.knots(1:index(2)),sp.coefs(:,1:(index(2)-sp.order)));
   end
   if sp.knots(1)<tl
      sp = spmak(sp.knots(index(1)-sp.order+1:end), ...
                 sp.coefs(:,index(1)-sp.order+1:end));
   end
   out1 = sp;
else
   error('SPLINES:SPBRK:partnotinterv',...
   'The given part, %g, does not specify an interval.',interv)
end