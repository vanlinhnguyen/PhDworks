function v = spline_getvals(f,varargin)
%FNVAL Evaluate a function.
%
%   V = FNVAL(F,X)  or  FNVAL(X,F)  provides the value at the points
%   in  X  of the function described by  F .
%
%   Roughly speaking, V is obtained by replacing each entry of X by the 
%   value of  f  there. This is exactly true in case  f  is scalar-valued
%   and univariate, and is the intent in all other cases, except that, for a 
%   d-valued m-variate function, this means replacing m-vectors by d-vectors.
%   The full details follow.
%
%   For a univariate  f :
%   If f is scalar-valued, then V is of the same size as X. 
%   If f is [d1,...,dr]-valued, and X has size [n1,...,ns], then V has size
%   [d1,...,dr, n1,...,ns], with V(:,...,:, j1,...,js) the value of  f  at 
%   X(j1,...,js), -- except that 
%   (1) n1 is ignored if it is 1 and s is 2, i.e., if X is a row vector; and
%   (2) MATLAB ignores any trailing singleton dimensions of X.
%
%   For an m-variate  f  with  m>1 ,  with  f  [d1,...,dr]-valued, X may be
%   either an array, or else a cell array {X1,...,Xm}.
%   If X is an array, of size [n1,...,ns] say, then n1 must equal m, and V has
%   size [d1,...,dr, n2,...,ns], with V(:,...,:, j2,...,js) the value of  f
%   at X(:,j2,...,js), -- except that
%   (1) d1, ..., dr is ignored in case  f  is scalar-valued, i.e., r==1==d1;
%   (2) MATLAB ignores any trailing singleton dimensions of X.
%   If X is a cell array, then it must be of the form {X1,...,Xm}, with Xj
%   a vector, of length nj, and, in that case, V has size
%   [d1,...,dr, n1,...,nm], with V(:,...,:, j1,...,jm) the value of  f
%   at (X1(j1), ..., Xm(jm)), -- except that
%   d1, ..., dr is ignored in case  f  is scalar-valued, i.e., r==1==d1.
%
%   By agreement, all piecewise polynomial functions in this toolbox are
%   continuous from the right. But FNVAL can be made to treat them as
%   continuous from the left by calling it with an optional third argument,
%   as follows.
%
%   FNVAL(F,X,LEFT)  or  FNVAL(X,F,LEFT)  takes the function to be
%   left-continuous if LEFT is a string that begins with 'l'.
%   If the function is m-variate and LEFT is an m-cell, then continuity
%   from the left is enforced in the i-th  variable if  LEFT{i}(1) is 'l'.
%
%   See also  PPUAL, RSVAL, SPVAL, STVAL, PPVAL.

%   Copyright 1987-2008 The MathWorks, Inc.
%   $Revision: 1.20.4.3 $

if ~isstruct(f)
   if isstruct(varargin{1})
      temp = f; f = varargin{1}; varargin{1} = temp;
   else
      f = fn2fm(f);
   end
end

try
   [m, sizeval] = fnbrk(f,'var','dim');
catch
   error('SPLINES:FNVAL:unknownfn',...
   ['Cannot handle the given form, ',f.form,'.'])
end
    % record, then adjust, size of site array and of function values.
sizex = size(varargin{1});
if ~iscell(varargin{1})
   if m>1
      if sizex(1)~=m
         error('SPLINES:FNVAL:wrongsizex', ...
	      ['Each X(:,j) must be a ',num2str(m),'-vector.'])
      end
   sizex(1) = [];
   end
   if length(sizex)>2
      varargin{1} = reshape(varargin{1},sizex(1),prod(sizex(2:end)));
   elseif length(sizex)==2&&sizex(1)==1, sizex = sizex(2); end
else
   if sizex(2)~=m
      if sizex(2)==1&&sizex(1)==m, varargin{1}=varargin{1}.';
      else
         error('SPLINES:FNVAL:wrongsizex', ...
	      ['X must be a cell array of size (1,',num2str(m),').'])
      end
   end
   sizex = cellfun('length',varargin{1});
end
if length(sizeval)>1, f = fnchg(f,'dz',prod(sizeval));
else if sizeval==1&&length(sizex)>1; sizeval = []; end
end

switch f.form(1)
case 'B',  ff = @spval;
case 'p',  ff = @ppual;
case 'r',  ff = @rsval;
case 's',  ff = @stval;
otherwise
   error('SPLINES:FNVAL:unknownfn','Unknown function type encountered.')
end

v = reshape(feval(ff,f,varargin{:}),[sizeval,sizex]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MISSING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%
%% Function 1
function varargout = fnbrk(fn,varargin)

if nargin>1
   np = max(1,nargout); % FNBRK(FN,PART) may be part of an expression
   if np <= length(varargin)
      varargout = cell(1,np);
   else
      error('SPLINES:FNBRK:moreoutthanin', ...
            'Too many output arguments for the given input.')
   end
end 

if ~isstruct(fn)    % this branch should eventually be abandoned
   switch fn(1)
   %
   % curves:
   %
      case 10, fnform = 'ppform, univariate, array format';
      case 11, fnform = 'B-form, univariate, array format';
      case 12, fnform = 'BBform, univariate, array format';
      case 15, fnform = 'polynomial in Newton form';
   %
   % surfaces:
   %
      case 20, fnform = 'ppform, bivariate tensor product, array format';
      case 21, fnform = 'B-form, bivariate tensor product, array format';
      case 22, fnform = 'BBform, bivariate, array format';
      case 24, fnform = 'polynomial in shifted power form, bivariate';
      case 25, fnform = 'thin-plate spline, bivariate';
   %
   % matrices:
   %
      case 40, fnform = 'almost block diagonal form';
      case 41, fnform = 'spline version of almost block diagonal form';
   % 42 = 'factorization of spline version of almost block diagonal form'
   %      (not yet implemented)
   
   %
   % multivariate:
   %
      case 94, fnform = ...
                  'polynomial in shifted normalized power form, multivariate';
      otherwise
         error('SPLINES:FNBRK:unknownform','Input is of unknown (function) form.')
   end
   
   if nargin>1 %  return some parts if possible
      switch fn(1)
      case 10, [varargout{:}] = ppbrk(fn,varargin{:});
      case {11,12}, [varargout{:}] = spbrk(fn,varargin{:});
      otherwise
         error('SPLINES:FNBRK:unknownpart',...
              ['Parts for ',fnform,' are not (yet) available.'])
      end
   else        % print available information
      if nargout
         error('SPLINES:FNBRK:partneeded','You need to specify a part to be returned.')
      else
         fprintf(['The input describes a ',fnform,'\n\n'])
         switch fn(1)
         case 10, ppbrk(fn);
         case {11,12}, spbrk(fn);
         otherwise
            fprintf('Its parts are not (yet) available.\n')
         end
      end
   end
   return
end   
 
     % we reach this point only if FN is a structure.

switch fn.form(1:2) 
case 'pp',        ffbrk = @ppbrk;
case 'rp',        ffbrk = @rpbrk;
case 'st',        ffbrk = @stbrk;
case {'B-','BB'}, ffbrk = @spbrk;
case 'rB',        ffbrk = @rsbrk;
otherwise
   error('SPLINES:FNBRK:unknownform','Input is of unknown (function) form.')
end

if nargin>1
   [varargout{:}] = ffbrk(fn,varargin{:});
else
   if nargout
      error('SPLINES:FNBRK:partneeded','You need to specify a part to be returned.')
   else
      fprintf(['The input describes a ',fn.form(1:2),'form\n\n'])
      ffbrk(fn)
   end
end


%% Function 2
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


%% Function 3
function v = spval(sp,x,left)

if ~isstruct(sp)
   error('SPLINES:SPVAL:fnnotstruct','SP must be a structure.'), end
if iscell(sp.knots)  % we are dealing with a tensor product spline

   [t,a,n,k,d] = spbrk(sp); m = length(t);

   if nargin>2 % set up left appropriately
      if ~iscell(left)
         temp = left; left = cell(1,m); [left{:}] = deal(temp);
      end
   else
      left = cell(1,m);
   end

   if iscell(x)  % evaluation on a mesh

      v = a; sizev = [d,n]; nsizev = zeros(1,m);

      for i=m:-1:1
         nsizev(i) = length(x{i}(:));
         v = reshape(...
         spval1(spmak(t{i},reshape(v,prod(sizev(1:m)),sizev(m+1))), ...
                 x{i},left{i}),   [sizev(1:m),nsizev(i)]);
         sizev(m+1) = nsizev(i);
         if m>1
            v = permute(v,[1,m+1,2:m]); sizev(2:m+1) = sizev([m+1,2:m]);
         end
      end
      if d>1
         v = reshape(v,[d,nsizev]);
      else
         v = reshape(v,nsizev);
      end
 
   else          % evaluation at scattered points;
                 % this will eventually be done directly here.
      v = ppual(sp2pp(sp),x);
      temp = cell2mat(fnbrk(sp,'interv')).'; mm = 2:2:(2*m); nx = size(x,2);
      v(min([x-repmat(temp(mm-1),1,nx); repmat(temp(mm),1,nx)-x])<0) = 0;
   end

else                 % we are dealing with a univariate spline
   if nargin<3, left = []; end
   v = spval1(sp,x,left);
end

function v = spval1(sp,x,left)
%SPVAL1 Evaluate univariate function in B-form.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);

%  Take apart spline:
[t,a,n,k,d] = spbrk(sp);
%  If there are no points to evaluate at, return empty matrix of appropriate
%  size:
if lx==0, v = zeros(d,0); return, end

%  Otherwise, augment the knot sequence so that first and last knot each
%  has multiplicity  >= K . (AUGKNT would not be suitable for this since
%  any change in T must be accompanied by a corresponding change in A.)

index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if ( addl>0 || addr>0 )
   npk = n+k; t = t([ones(1,addl) 1:npk npk(ones(1,addr))]);
   a = [zeros(d,addl) a zeros(d,addr)];
   n = n+addl+addr;
end

% For each data point, compute its knot interval:
if isempty(left)||left(1)~='l'
   [ignored, index] = histc(xs,[-inf,t(k+1:n),inf]);
   NaNx = find(index==0); index = min(index+(k-1),n);
else
   [ignored, index] = histc(-xs,[-inf,-fliplr(t(k+1:n)),inf]);
   NaNx = find(index==0); index = max(n+1-index,k);
end
if ~isempty(NaNx), index(NaNx) = k; end

% Now, all is ready for the evaluation.
if  k>1  % carry out in lockstep the first spline evaluation algorithm
         % (this requires the following initialization):
   dindex = reshape(repmat(index,d,1),d*lx,1);
   tx =reshape(t(repmat(2-k:k-1,d*lx,1)+repmat(dindex,1,2*(k-1))),d*lx,2*(k-1));
   tx = tx - repmat(reshape(repmat(xs,d,1),d*lx,1),1,2*(k-1));
   dindex = reshape(repmat(d*index,d,1)+repmat((1-d:0).',1,lx),d*lx,1);
   b = repmat(d*(1-k):d:0,d*lx,1)+repmat(dindex,1,k);
   a = a(:); b(:) = a(b);

   % (the following loop is taken from SPRPP)

   for r = 1:k-1
      for i = 1:k-r
         b(:,i) = (tx(:,i+k-1).*b(:,i)-tx(:,i+r-1).*b(:,i+1)) ./ ...
                  (tx(:,i+k-1)    -    tx(:,i+r-1));
      end
   end

   v = reshape(b(:,1),d,lx);
else     % the spline is piecewise constant, hence ...
   v = a(:,index);
   if ~isempty(NaNx), v(:,NaNx) = NaN; end
end

% Finally, zero out all values for points outside the basic interval:
index = find(x<t(1)|x>t(n+k));
if ~isempty(index)
   v(:,index) = zeros(d,length(index));
end
v = reshape(v,d*mx,nx);


