function [points,t] = spline_plot(f,varargin)
%FNPLT Plot a function.
%
%   FNPLT(F)  plots the function in F on its basic interval.
%
%   FNPLT(F,SYMBOL,INTERV,LINEWIDTH,JUMPS) plots the function F
%   on the specified INTERV = [a,b] (default is the basic interval), 
%   using the specified plotting SYMBOL (default is '-'), 
%   and the specified LINEWIDTH (default is 1), 
%   and using NaNs in order to show any jumps as actual jumps only 
%   in case JUMPS is a string beginning with 'j'.
%
%   The four optional arguments may appear in any order, with INTERV
%   the one of size [1 2], SYMBOL and JUMPS strings, and LINEWIDTH the
%   scalar. Any empty optional argument is ignored.
%
%   If the function in F is 2-vector-valued, the planar curve is
%   plotted.  If the function in F is d-vector-valued with d>2, the
%   space curve given by the first three components of F is plotted.
%
%   If the function is multivariate, it is plotted as a bivariate function,
%   at the midpoint of its basic intervals in additional variables, if any.
%
%   POINTS = FNPLT(F,...)   does not plot, but returns instead the sequence 
%   of 2D-points or 3D-points it would have plotted.
%
%   [POINTS,T] = FNPLT(F,...)  also returns, for a vector-valued F, the 
%   corresponding vector T of parameter values.
%
%   Example:
%      x=linspace(0,2*pi,21); f = spapi(4,x,sin(x));
%      fnplt(f,'r',3,[1 3])
%
%   plots the graph of the function in f, restricted to the interval [1 .. 3],
%   in red, with linewidth 3 .

%   Copyright 1987-2008 The MathWorks, Inc.
%   $Revision: 1.22.4.3 $

% interpret the input:
symbol=''; interv=[]; linewidth=[]; jumps=0;
for j=2:nargin
   arg = varargin{j-1};
   if ~isempty(arg)
      if ischar(arg)
         if arg(1)=='j', jumps = 1;
         else symbol = arg;
         end
      else
         [ignore,d] = size(arg);
         if ignore~=1
	    error('SPLINES:FNPLT:wrongarg',['arg',num2str(j),' is incorrect.']), end
         if d==1
            linewidth = arg;
         else
            interv = arg;
         end
      end
   end
end

% generate the plotting info:
if ~isstruct(f), f = fn2fm(f); end

% convert ND-valued to equivalent vector-valued:
d = fnbrk(f,'dz'); if length(d)>1, f = fnchg(f,'dim',prod(d)); end

switch f.form(1:2)
case 'st'
   if ~isempty(interv), f = stbrk(f,interv);
   else
      interv = stbrk(f,'interv');
   end
   npoints = 51; d = stbrk(f,'dim');
   switch fnbrk(f,'var')
   case 1
      x = linspace(interv{1}(1),interv{1}(2),npoints);
      v = stval(f,x);
   case 2
      x = {linspace(interv{1}(1),interv{1}(2),npoints), ...
                    linspace(interv{2}(1),interv{2}(2),npoints)};
      [xx,yy] = ndgrid(x{1},x{2});
      v = reshape(stval(f,[xx(:),yy(:)].'),[d,size(xx)]);
   otherwise
      error('SPLINES:FNPLT:atmostbivar', ...
            'Cannot handle st functions with more than 2 variables.')
   end
otherwise
   if ~strcmp(f.form([1 2]),'pp')
      givenform = f.form; f = fn2fm(f,'pp'); basicint = ppbrk(f,'interval');
   end
   
   if ~isempty(interv), f = ppbrk(f,interv); end
      
   [breaks,l,d] = ppbrk(f,'b','l','d');
   if iscell(breaks)
      m = length(breaks);
      for i=m:-1:3
         x{i} = (breaks{i}(1)+breaks{i}(end))/2;
      end
      npoints = 51;
      ii = 1; if m>1, ii = [2 1]; end
      for i=ii
         x{i}= linspace(breaks{i}(1),breaks{i}(end),npoints);
      end
      v = ppual(f,x);
      if exist('basicint','var') 
                         % we converted from B-form to ppform, hence must now
                         % enforce the basic interval for the underlying spline.
         for i=ii
            temp = find(x{i}<basicint{i}(1)|x{i}>basicint{i}(2));
            if d==1
               if ~isempty(temp), v(:,temp,:) = 0; end
               v = permute(v,[2,1]);
            else
               if ~isempty(temp), v(:,:,temp,:) = 0; end
               v = permute(v,[1,3,2]);
            end
         end
      end
   else     % we are dealing with a univariate spline
      npoints = 101;
      x = [breaks(2:l) linspace(breaks(1),breaks(l+1),npoints)];
      v = ppual(f,x); 
      if l>1 % make sure of proper treatment at jumps if so required
         if jumps
            tx = breaks(2:l); temp = repmat(NaN, d,l-1);
         else
            tx = []; temp = zeros(d,0);
         end
         x = [breaks(2:l) tx x]; 
         v = [ppual(f,breaks(2:l),'left') temp v];
      end
      [x,inx] = sort(x); v = v(:,inx);
   
      if exist('basicint','var') 
                         % we converted from B-form to ppform, hence must now
                         % enforce the basic interval for the underlying spline.
                         % Note that only the first d components are set to zero
                         % outside the basic interval, i.e., the (d+1)st 
                         % component of a rational spline is left unaltered :-)
         if jumps, extrap = repmat(NaN,d,1); else extrap = zeros(d,1); end
         temp = find(x<basicint(1)); ltp = length(temp);
         if ltp
            x = [x(temp),basicint([1 1]), x(ltp+1:end)];
            v = [zeros(d,ltp+1),extrap,v(:,ltp+1:end)];
         end
         temp = find(x>basicint(2)); ltp = length(temp);
         if ltp
            x = [x(1:temp(1)-1),basicint([2 2]),x(temp)];
            v = [v(:,1:temp(1)-1),extrap,zeros(d,ltp+1)];
         end
         %   temp = find(x<basicint(1)|x>basicint(2));
         %   if ~isempty(temp), v(temp) = zeros(d,length(temp)); end
      end
   end
   
   if exist('givenform','var')&&givenform(1)=='r' 
                                           % we are dealing with a rational fn:
                                           % need to divide by last component
      d = d-1;
      sizev = size(v); sizev(1) = d;
      % since fnval will replace any zero value of the denominator by 1,
      % so must we here, for consistency:
      v(d+1,v(d+1,:)==0) = 1;
      v = reshape(v(1:d,:)./repmat(v(d+1,:),d,1),sizev);
   end 
end

%  use the plotting info, to plot or else to output:
if nargout==0
   if iscell(x)
      switch d
      case 1
         [yy,xx] = meshgrid(x{2},x{1});
         surf(xx,yy,reshape(v,length(x{1}),length(x{2})))
      case 2
         v = squeeze(v); roughp = 1+(npoints-1)/5;
         vv = reshape(cat(1,...
              permute(v(:,1:5:npoints,:),[3,2,1]),...
              repmat(NaN,[1,roughp,2]),...
              permute(v(:,:,1:5:npoints),[2,3,1]),...
              repmat(NaN,[1,roughp,2])), ...
                    [2*roughp*(npoints+1),2]);
         plot(vv(:,1),vv(:,2))
      case 3
         v = permute(reshape(v,[3,length(x{1}),length(x{2})]),[2 3 1]);
         surf(v(:,:,1),v(:,:,2),v(:,:,3))
      otherwise
      end
   else
      if isempty(symbol), symbol = '-'; end
      if isempty(linewidth), linewidth = 2; end
      switch d
      case 1, plot(x,v,symbol,'linew',linewidth)
      case 2, plot(v(1,:),v(2,:),symbol,'linew',linewidth)
      otherwise
         plot3(v(1,:),v(2,:),v(3,:),symbol,'linew',linewidth)
      end
   end
else
   if iscell(x)
      switch d
      case 1
         [yy,xx] = meshgrid(x{2},x{1});
         points = {xx,yy,reshape(v,length(x{1}),length(x{2}))};
      case 2
         [yy,xx] = meshgrid(x{2},x{1});
         points = {xx,yy,reshape(v,[2,length(x{1}),length(x{2})])};
      case 3
         points = {squeeze(v(1,:)),squeeze(v(2,:)),squeeze(v(3,:))};
         t = {x{1:2}};
      otherwise
      end
   else
      if d==1, points = [x;v];
      else t = x; points = v(1:min([d,3]),:); end
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%% MISSING FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
function g = fn2fm(f,form,sconds)
   
if ~isstruct(f) && numel(f)<2 % reject empty F or F that's just a scalar
   error('SPLINES:FN2FM:notafunction','Input is not a function')
end

if nargin==1 % set FORM to convert the form in F into its current version
   if isstruct(f)
      form = f.form;
   else
      switch f(1)
      case 10,  form = 'pp';
      case 11,  form = 'B-';
      case 12,  form = 'BB';
      case 15,  form = 'NP';
      case 25,  form = 'st-tp00';
      otherwise
         error('SPLINES:FN2FM:unknownform','Input is of unknown (function) form.')
      end
   end
else
   if length(form)<2
      error('SPLINES:FN2FM:unknownform',...
           ['The specified form, ''',form,''', is not (yet) recognized.']),
   end 
end

switch form(1:2)
case 'pp'
   if isstruct(f)
      switch f.form(1:2)
      case 'pp',          g = f;
      case 'rp',          g = f; g.form = 'pp'; g.dim = g.dim+1;
      case 'rB',          f.form = 'B-'; f.dim = f.dim+1; g = sp2pp(f);
      case {'B-','BB'},   g = sp2pp(f);
      end
   else
      switch f(1)
      case 10,      [b,c,l,k,d] = ppbrk(f); g = ppmak(b,c,d);
      case {11,12}, g = sp2pp(f);
      case 15
       error('SPLINES:FN2FM:nonnewtonhere',...
             'Use NP2PP directly to convert from Newton form to ppform.')
      end
   end
case {'sp','B-'}
   form = 'B-';
   if isstruct(f)
      switch f.form(1:2)
      case 'pp'
         if nargin<3, g = pp2sp(f);
         else         g = pp2sp(f,sconds);
         end
      case 'rp',      f.form = 'pp'; f.dim = f.dim+1; g = pp2sp(f);
      case 'rB',      g = f; g.form = 'B-'; g.dim = g.dim+1;
      case 'B-',      g = f;
      case 'BB',      g = pp2sp(sp2pp(f));
      end
   else
      switch f(1)
      case 10, g = pp2sp(f);
      case 11, [knots,coefs] = spbrk(f); g = spmak(knots,coefs);
      case 12, g = pp2sp(sp2pp(f));
      end
   end
case {'bb','BB'}
   form = 'BB';
   if isstruct(f)
      switch f.form(1:2)
      case 'pp', g = sp2bb(pp2sp(f));          
      case 'rp', f.form = 'pp'; f.dim = f.dim+1; g = sp2bb(pp2sp(f));
      case 'rB', f.form = 'B-'; f.dim = f.dim+1; g = sp2bb(f);
      case 'B-', g = sp2bb(f);
      case 'BB',  g = f;
      end
   else
      switch f(1)
      case 10, g = sp2bb(pp2sp(f));
      case 11, g = sp2bb(f);
      case 12, [knots,coefs] = spbrk(f); g = spmak(knots,coefs); g.form = 'BB';
      end
   end
case 'st' 
   if isstruct(f)
      switch f.form(1:2)
      case 'st', g = f;
      end
   else
      switch f(1)
      case 25, [ce,co] = stbrk(f); g = stmak(ce,co,'tp00');
      end
   end
case 'rp'  % switch first into ppform if need be
   if isstruct(f)
      switch f.form(1:2)
      case 'pp',      g = f;
      case 'rp',      g = f; g.dim = g.dim+1;
      case 'rB',      g = sp2pp(fn2fm(f,'B-'));
      case {'B-','BB','sp','bb'}, g = sp2pp(f);
      end
   else
      switch f(1)
      case 10,      [b,c,l,k,d] = ppbrk(f); g = ppmak(b,c,d);
      case {11,12}, g = sp2pp(f);
      end
   end
   if exist('g','var')
      g.form = 'rp'; 
      if length(g.dim)>1
         warning('SPLINES:FN2FM:noNDrats', ...
	 ['While the given function has values of size [', num2str(g.dim), ...
	  '], the function returned is ', num2str(prod(g.dim)-1), ...
	  '-vector valued.'])
	 g.dim = prod(g.dim);
      end
      g.dim = g.dim-1;
      if g.dim<1
         error('SPLINES:FN2FM:ratneedsmoredim',...
               'A rational spline must have more than one component.')
      end
   end
case 'rB'  % switch first into B-form if need be
   if isstruct(f)
      switch f.form(1:2)
      case 'pp'
         if nargin<3, g = pp2sp(f);
         else         g = pp2sp(f,sconds);
         end
      case 'rp',      g = fn2fm(fn2fm(f,'pp'),'B-');
      case 'rB',      g = f; g.dim = g.dim+1;
      case {'B-','BB','sp','bb'}, g = f;
      end   
   else 
      switch f(1)
      case 10, g = pp2sp(f);
      case 11, g = spmak(fnbrk(f,'knots'),fnbrk(f,'coefs'));
      case 12, g = pp2sp(sp2pp(f));
      end
   end
   if exist('g','var')
      g.form = 'rB';
      if length(g.dim)>1
         warning('SPLINES:FN2FM:noNDrats', ...
	 ['While the given function has values of size [', num2str(g.dim), ...
	  '], the function returned is ', num2str(prod(g.dim)-1), ...
	  '-vector valued.'])
	 g.dim = prod(g.dim);
      end
      g.dim = g.dim-1;
      if g.dim<1
         error('SPLINES:FN2FM:ratneedsmoredim',...
               'A rational spline must have more than one component.')
      end
   end
case 'NP' % at present, this form only exists in the old version
   if isstruct(f)
      switch f.form(1:2)
      case 'NP', g = f;
      end
   else
      switch f(1)
      case 15, g = f;
      end
   end
case 'MA'     % convert univariate spline to old, nonstructure ppform.
   if isstruct(f)&&~iscell(f.order)
      switch f.form(1:2)
      case 'pp' 
         g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
      case {'B-','BB'}
         f = sp2pp(f);
         g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
      otherwise
         error('SPLINES:FN2FM:notintoMA',...
               'Cannot convert the given function into MATLAB''s ppform.')
      end
   else
      g = f;
   end
case 'ol'    % convert univariate structured form to corresponding 
             % formerly used array form.
   if ~isstruct(f)
      error('SPLINES:FN2FM:olneedsstruct',...
            'When FORM=''ol'', F must be a structure.')
   else
      if length(f.order)>1
         error('SPLINES:FN2FM:nonstructneedsuni',...
              ['Reversion to array is only possible for',...
                ' univariate functions.'])
      else
         switch f.form(1:2)
         case 'pp' 
            g = [10 f.dim f.pieces f.breaks(:).' f.order f.coefs(:).'];
         case {'B-','BB'}
            g = [11 f.dim f.number f.coefs(:).' f.order f.knots(:).'];
            if strcmp(f.form(1:2),'BB'), g(1) = 12; end
         otherwise
            error('SPLINES:FN2FM:unknownfn','Unknown function type encountered.')
         end
      end
   end
otherwise
   error('SPLINES:FN2FM:unknownfn',...
    ['The specified form, ''',form,''', is not (yet) recognized.']),
end

if ~exist('g','var')
   error('SPLINES:FN2FM:notintothatform',...
        ['Cannot convert the given function into ',form,'form.'])
end


%% Function 4
function pp = sp2pp(spline)
%SP2PP Convert from B-form to ppform.
%
%   SP2PP(SPLINE)  converts the B-form in SPLINE to the corresponding ppform
%   (on its basic interval).
%
%   For example,
%
%      p0 = ppmak([0 1],[3 0 0]); p1 = sp2pp(pp2sp(pprfn(p0,[.4 .6])));
%
%   gives p1 identical to p0 (up to round-off) since the spline has no
%   discontinuity in any derivative across the additional breaks introduced
%   by PPRFN, hence PP2SP ignores these additional breaks, and SP2PP does
%   not retain any knot multiplicities (like the knot multiplicities introduced
%   by PP2SP at the endpoints of the spline's basic interval).
%
%   See also PP2SP, SP2BB, FN2FM.

%   Copyright 1987-2008 The MathWorks, Inc.
%   $Revision: 1.17.4.3 $

if ~isstruct(spline), spline = fn2fm(spline); end

sizeval = fnbrk(spline,'dim');
if length(sizeval)>1, spline = fnchg(spline,'dz',prod(sizeval)); end

if iscell(spline.knots)   % we are dealing with a multivariate spline

   [t,a,n,k,d] = spbrk(spline);
   m = length(k);
   coefs = a; sizec = [prod(d),n]; % size(coefs);
   for i=m:-1:1
      ppi = sp2pp1(spmak(t{i},reshape(coefs,prod(sizec(1:m)),n(i))));
      breaks{i} = ppi.breaks;  sizec(m+1) = ppi.pieces*k(i);
      coefs = reshape(ppi.coefs,sizec);
      if m>1
         coefs = permute(coefs,[1,m+1,2:m]); sizec = sizec([1,m+1,2:m]);
      end
   end
   pp = ppmak(breaks,coefs,sizec);
      
else
   pp = sp2pp1(spline);
end

if length(sizeval)>1, pp = fnchg(pp,'dz',sizeval); end

function pp = sp2pp1(spline)
%  Take apart the  spline

[t,a,n,k,d] = spbrk(spline);

%  and augment the knot sequence so that first and last knot each have
%  multiplicity  k .

index = find(diff(t)>0); addl = k-index(1); addr = index(end)-n;
if (addl>0||addr>0)
   t = [repmat(t(1),1,addl) t(:).' repmat(t(n+k),1,addr)];
   a = [zeros(d,addl) a zeros(d,addr)];
end

%  From this, generate the pp description.

inter = find( diff(t)>0 ); l = length(inter);
if k>1
   temp = repmat(inter,d,1); dinter = temp(:);
   tx = repmat(2-k:k-1,d*l,1)+repmat(dinter,1,2*(k-1)); tx(:) = t(tx);
   tx = tx-repmat(t(dinter).',1,2*(k-1)); a = a(:);
   temp = repmat(d*inter,d,1)+repmat((1-d:0).',1,l); dinter(:) = temp(:);
   b = repmat(d*(1-k:0),d*l,1)+repmat(dinter,1,k); b(:) = a(b);
   c = sprpp(tx,b);
else temp = a(:,inter); c = temp(:);
end

%   put together the  pp

pp = ppmak([t(inter) t(inter(end)+1)],c,d);


%% Function 5
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
function [v,b] = sprpp(tx,a)
k = length(a(1,:)); km1 = k-1; b = a;
for r=1:km1
   for i=1:k-r
      b(:,i) =(tx(:,i+km1).*b(:,i)-tx(:,i+r-1).*b(:,i+1))./...
               (tx(:,i+km1)-tx(:,i+r-1));
   end
end

%  Use differentiation at  0  to generate the derivatives

v = b;
for r=2:k
   factor = (k-r+1)/(r-1);
   for i=k:-1:r
      v(:,i) = (v(:,i) - v(:,i-1))*factor./tx(:,i+k-r);
   end
end

v = v(:,k:-1:1);

%% Function 8
function pp = ppmak(breaks,coefs,d)

if nargin==0
    breaks=input('Give the (l+1)-vector of breaks  >');
    coefs=input('Give the (d by (k*l)) matrix of local pol. coefficients  >');
end

sizec = size(coefs);

if iscell(breaks)  % we are dealing with a tensor-product spline
    if nargin>2
        if prod(sizec)~=prod(d)
            error('SPLINES:PPMAK:coefsdontmatchsize', ...
                'The coefficient array is not of the explicitly specified size.')
        end
        sizec = d;
    end
    m = length(breaks);
    if length(sizec)<m
        error('SPLINES:PPMAK:coefsdontmatchbreaks', ...
            ['If BREAKS is a cell-array of length m, then COEFS must ',...
            'have at least m dimensions.'])
    end
    if length(sizec)==m,  % coefficients of a scalar-valued function
        sizec = [1 sizec];
    end
    sizeval = sizec(1:end-m);
    sizec = [prod(sizeval), sizec(end-m+(1:m))];
    coefs = reshape(coefs, sizec);
    
    for i=m:-1:1
        l(i) = length(breaks{i})-1;
        k(i) = fix(sizec(i+1)/l(i));
        if k(i)<=0||k(i)*l(i)~=sizec(i+1)
            error('SPLINES:PPMAK:piecesdontmatchcoefs', ...
                ['The specified number %g of polynomial pieces is', ...
                ' incompatible\nwith the total number %g of coefficients', ...
                ' supplied in variable %g.'], l(i),sizec(i+1),i)
        end
        breaks{i} = reshape(breaks{i},1,l(i)+1);
    end
else
    if nargin<3
        if isempty(coefs)
            error('SPLINES:PPMAK:emptycoefs','The coefficient sequence is empty.')
        end
        sizeval = sizec(1:end-1);
        d = prod(sizeval);
        kl = sizec(end);
        l=length(breaks)-1;
        k=fix(kl/l);
        if (k<=0)||(k*l~=kl)
            error('SPLINES:PPMAK:piecesdontmatchcoefs', ...
                ['The specified number %g of polynomial pieces is',...
                ' incompatible\nwith the total number %g of coefficients',...
                ' supplied.'],l,kl);
        elseif any(diff(breaks)<0)
            error('SPLINES:PPMAK:decreasingbreaks', ...
                'The break sequence should be nondecreasing.')
        elseif breaks(1)==breaks(l+1)
            error('SPLINES:PPMAK:extremebreakssame', ...
                'The extreme breaks should be different.')
        else
            % the ppformat expects coefs in array  (d*l) by k, while the standard
            % input supplies them in an array d by (k*l) . This requires the
            % following shuffling, from  D+d(-1+K + k(-1+L))=D-d +(K-k)d + dkL
            % to  D+d(-1+L + l(-1+K)=D-d +(L-l)d + dlK .
            % This used to be handled by the following:
            % c=coefs(:); temp = ([1-k:0].'*ones(1,l)+k*ones(k,1)*[1:l]).';
            % coefs=[1-d:0].'*ones(1,kl)+d*ones(d,1)*(temp(:).');
            % coefs(:)=c(coefs);
            % Thanks to multidimensional arrays, we can now simply say
            coefs = reshape(permute(reshape(coefs,[d,k,l]),[1,3,2]),d*l,k);
        end
    else % in the univariate case, a scalar D only specifies the dimension of
        % the target and COEFS must be a matrix (though that is not checked for);
        % but if D is a vector, then it is taken to be the intended size of
        % COEFS whatever the actual dimensions of COEFS might be.
        if length(d)==1
            k = sizec(end);
            l = prod(sizec(1:end-1))/d;
        else
            if prod(d)~=prod(sizec)
                error('SPLINES:PPMAK:coefsdontmatchsize', ...
                    ['The size of COEFS, [',num2str(sizec), ...
                    '], does not match the specified size, [',num2str(d),'].'])
            end
            k = d(end);
            l = d(end-1);
            d(end-1:end) = [];
            if isempty(d),
                d = 1;
            end
            % Interpret the coefficients, COEFS, to be of size SIZEC =: [d,l,k]
            coefs = reshape(coefs, prod(d)*l,k);
        end
        if l+1~=length(breaks)
            error('SPLINES:PPMAK:coefsdontmatchbreaks', ...
                'COEFS indicates %g piece(s) while BREAKS indicates %g.', ...
                l, length(breaks)-1)
        end
        sizeval = d;
    end
    breaks = reshape(breaks,1,l+1);
end
pp.form = 'pp';
pp.breaks = breaks;
pp.coefs = coefs;
pp.pieces = l;
pp.order = k;
pp.dim = sizeval;

%% Function 8
function varargout = ppbrk(pp,varargin)

if ~isstruct(pp)
   if pp(1)~=10
      error('SPLINES:PPBRK:unknownfn',...
      'The input array does not seem to describe a function in ppform.')
   else
      ppi = pp;
      di=ppi(2); li=ppi(3); ki=ppi(5+li);

      pp = struct('breaks',reshape(ppi(3+(1:li+1)),1,li+1), ...
                  'coefs',reshape(ppi(5+li+(1:di*li*ki)),di*li,ki), ...
                  'form','pp', 'dim',di, 'pieces',li, 'order',ki);
   end
end 

if ~strcmp(pp.form,'pp')
   error('SPLINES:PPBRK:notpp',...
   'The input does not seem to describe a function in ppform.')
end
if nargin>1 % we have to hand back one or more parts
   np = max(1,nargout);
   if np>length(varargin)
      error('SPLINES:PPBRK:moreoutthanin', ...
            'Too many output arguments for the given input.')
   end
   varargout = cell(1,np);
   for jp=1:np
      part = varargin{jp};

      if ischar(part)
         if isempty(part)
	    error('SPLINES:PPBRK:partemptystr',...
	    'Part specification should not be an empty string.')
	 end
         switch part(1)
         case 'f',       out1 = [pp.form,'form'];
         case 'd',       out1 = pp.dim;
         case {'l','p'}, out1 = pp.pieces;
         case 'b',       out1 = pp.breaks;
         case 'o',       out1 = pp.order;
         case 'c',       out1 = pp.coefs;
	 case 'v',       out1 = length(pp.order);
         case 'g',       % if the spline is univariate, scalar-valued,
                         % return the coefs in the form needed in the ppform
                         % used in PGS.
            if length(pp.dim)>1||pp.dim>1||iscell(pp.order)
               error('SPLINES:PPBRK:onlyuniscalar', ...
                     ['''%s'' is only available for scalar-valued',...
                              ' univariate pp functions.'],part)
            else
               k = pp.order;
               out1 = (pp.coefs(:,k:-1:1).').* ...
	                repmat(cumprod([1 1:k-1].'),1,pp.pieces);
            end
         case 'i'
            if iscell(pp.breaks)
               for i=length(pp.order):-1:1
                  out1{i} = pp.breaks{i}([1 end]); end
            else
               out1 = pp.breaks([1 end]);
            end
         otherwise
            error('SPLINES:PPBRK:unknownpart',...
	    '''%s'' is not part of a ppform.',part)
         end
      elseif isempty(part)
	 out1 = pp;
      else % we are to restrict PP to some interval or piece
	 sizeval = pp.dim; if length(sizeval)>1, pp.dim = prod(sizeval); end
         if iscell(part)  % we are dealing with a tensor-product spline
   
            [breaks,c,l,k,d] = ppbrk(pp); m = length(breaks);
            sizec = [d,l.*k]; %size(c);
            if length(sizec)~=m+1
	       error('SPLINES:PPBRK:inconsistentfn', ...
	       'Information in PP is inconsistent.'),
            end
            for i=m:-1:1
               dd = prod(sizec(1:m));
               ppi = ppbrk1(ppmak(breaks{i},reshape(c,dd*l(i),k(i)),dd),...
                           part{i}) ;
               breaks{i} = ppi.breaks; sizec(m+1) = ppi.pieces*k(i);
               c = reshape(ppi.coefs,sizec);
               if m>1
                  c = permute(c,[1,m+1,2:m]);
                  sizec(2:m+1) = sizec([m+1,2:m]);
               end
            end
            out1 = ppmak(breaks,c, sizec);
   
         else  % we are dealing with a univariate spline
   
            out1 = ppbrk1(pp,part);
         end
         if length(sizeval)>1, out1 = fnchg(out1,'dz',sizeval); end
      end
      varargout{jp} = out1;
   end
else
   if nargout==0
     if iscell(pp.breaks) % we have a multivariate spline and, at present,
                          % I can't think of anything clever to do; so...
       disp(pp)
     else
       disp('breaks(1:l+1)'),        disp(pp.breaks)
       disp('coefficients(d*l,k)'),  disp(pp.coefs)
       disp('pieces number l'),      disp(pp.pieces)
       disp('order k'),              disp(pp.order)
       disp('dimension d of target'),disp(pp.dim)
       % disp('dimension v of domain'),disp(length(pp.order))
     end
   else
      varargout = {pp.breaks, pp.coefs, pp.pieces, pp.order, pp.dim};
   end
end

function pppart = ppbrk1(pp,part)
%PPBRK1 restriction of pp to some piece or interval

if isempty(part)||ischar(part), pppart = pp; return, end

if size(part,2) > 1 , % extract the part relevant to the interval 
                      % specified by  part =: [a b]  
   pppart = ppcut(pp,part(1,1:2));
else                  % extract the part(1)-th polynomial piece of pp (if any)
   pppart = pppce(pp,part(1));
end

function ppcut = ppcut(pp,interv)
%PPCUT returns the part of pp  specified by the interval interv =: [a b]  

xl = interv(1); xr = interv(2); if xl>xr, xl = xr; xr = interv(1);  end
if xl==xr
   warning('SPLINES:PPBRK:PPCUT:trivialinterval', ...
           'No changes made since the given end points are equal.')
   ppcut = pp; return
end
 
%  the first pol. piece is  jl ,
% the one responsible for argument  xl
jl=pp.pieces; index=find(pp.breaks(2:jl)>xl); 
                                   % note that the resulting  index  ...
if (~isempty(index)), jl=index(1); % ... is shifted down by one  ...
end                                % ... because of  breaks(2: ...
%  if xl ~= breaks(jl), recenter the pol.coeffs.
x=xl-pp.breaks(jl);
di = pp.dim;
if x ~= 0
   a=pp.coefs(di*jl+(1-di:0),:);
   for ii=pp.order:-1:2
      for i=2:ii
         a(:,i)=x*a(:,i-1)+a(:,i);
      end
   end
   pp.coefs(di*jl+(1-di:0),:)=a;
end
 
%  the last pol. piece is  jr ,
% the one responsible for argument  xr .
jr=pp.pieces;index=find(pp.breaks(2:jr+1)>=xr); 
                                   % note that the resulting ...
if (~isempty(index)), jr=index(1); % index  is shifted down by
end                                % ... one because of  breaks(2: ...
 
%  put together the cut-down  pp
di = pp.dim;
ppcut = ppmak([xl pp.breaks(jl+1:jr) xr], ...
                        pp.coefs(di*(jl-1)+(1:di*(jr-jl+1)),:),di);

function pppce = pppce(pp,j)
%PPPCE returns the j-th polynomial piece of pp  (if any).

%  if  pp  has a  j-th  piece, ...
if (0<j)&&(j<=pp.pieces)  %             ...  extract it
   di = pp.dim;
   pppce = ppmak([pp.breaks(j) pp.breaks(j+1)], ...
              pp.coefs(di*j+(1-di:0),:),di);
else
   error('SPLINES:PPBRK:wrongpieceno', ...
   'The given pp function does not have %g pieces.',j);
end



%% Function 8
function v = ppual(pp,x,left)

if ~isstruct(pp), pp = fn2fm(pp); end

if iscell(pp.breaks)   % we are dealing with a multivariate spline

   [breaks,coefs,l,k,d] = ppbrk(pp); m = length(breaks);
   sizec = [d,l.*k]; % size(coefs)

   if nargin>2 % set up LEFT appropriately
      if ~iscell(left)
         temp = left; left = cell(1,m); [left{:}] = deal(temp);
      end
   else
      left = cell(1,m);
   end

   if iscell(x)  % evaluation on a mesh

      if length(x)~=m
         error('SPLINES:PPUAL:needgrid',...
	      ['X should specify a(n) ',num2str(m),'-dimensional grid.'])
      end

      v = coefs; sizev = sizec;
      nsizev = zeros(1,m);
      for i=m:-1:1
         nsizev(i) = length(x{i}(:)); dd = prod(sizev(1:m));
         v = reshape(ppual1(...
              ppmak(breaks{i},reshape(v,dd*l(i),k(i)),dd), x{i}, left{i} ),...
                [sizev(1:m),nsizev(i)]);
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

   else          % evaluation at scattered points

      % locate the scattered data in the break sequences:
      [mx,n] = size(x);
      if mx~=m, error('SPLINES:PPUAL:wrongx', ...
                     ['Each X(:,j) must be a ',num2str(m),'-vector.']), end

      ix = zeros(m,n);
      for i=1:m
         [ox,iindex] = sort(x(i,:));
         ix(i,iindex) = get_index(breaks{i}(2:end-1),ox,left{i});
      end

      % ... and now set up lockstep polynomial evaluation
      % %%First,  select the relevant portion of the coefficients array.
      % This has the additional pain that now there are k(i) coefficients
      % for the i-th univariate interval.
      % The coefficients sit in the (m+1)-dimensional array COEFS, with
      % the (i+1)st dimension containing the coefficients in the i-th
      % dimension, and organized to have first the highest coefficients
      % for each interval, then the next-highest, etc (i.e., as if coming
      % from an array of size [l(i),k(i)]).
      % ix(:,j) is the index vector for the lower corner of j-th point
      % The goal is to extract, for the j-th point, the requisite coefficients
      % from the equivalent one-dimensional array for COEFS, computing a
      % base index from ix(:,j), and adding to this the same set of offsets
      % computed from the l(i) and k(i).

      temp = l(1)*(0:k(1)-1)';
      for i=2:m
         lt = length(temp(:,1));
         temp = [repmat(temp,k(i),1), ...
                 reshape(repmat(l(i)*(0:k(i)-1),lt,1),k(i)*lt,1)];
      end
      % also take care of the possibility that the function in PP is
      % vector-valued:
      lt = length(temp(:,1));
      temp=[reshape(repmat((0:d-1).',1,lt),d*lt,1) temp(repmat(1:lt,d,1),:)];

      temp = num2cell(1+temp,1);
      offset = repmat(reshape(sub2ind(sizec,temp{:}),d*prod(k),1),1,n);

      temp = num2cell([ones(n,1) ix.'],1);
      base = repmat(sub2ind(sizec,temp{:}).',d*prod(k),1)-1;
      v = reshape(coefs(base+offset),[d,k,n]);

      % ... then do a version of local polynomial evaluation
      for i=m:-1:1
         s = reshape(x(i,:) - breaks{i}(ix(i,:)),[1,1,n]);
         otherk = d*prod(k(1:i-1));
         v = reshape(v,[otherk,k(i),n]);
         for j=2:k(i)
            v(:,1,:) = v(:,1,:).*repmat(s,[otherk,1,1])+v(:,j,:);
         end
         v(:,2:k(i),:) = [];
      end
      v = reshape(v,d,n);
   end
else
   if nargin<3, left = []; end
   v = ppual1(pp,x,left);
end

function v = ppual1(pp,x,left)
%PPUAL1 Evaluate univariate function in ppform.

[mx,nx] = size(x); lx = mx*nx; xs = reshape(x,1,lx);

%  take apart PP
[breaks,c,l,k,d] = ppbrk(pp);
%  if there are no points to evaluate at, return empty matrix of appropriate
%  size.
if lx==0, v = zeros(d,0); return, end

% for each data site, compute its break interval
[index,NaNx] = get_index(breaks(2:end-1),xs,left);
index(NaNx) = 1;

% now go to local coordinates ...
xs = xs-breaks(index);
if d>1 % ... replicate XS and INDEX in case PP is vector-valued ...
   xs = reshape(repmat(xs,d,1),1,d*lx);
   index = reshape(repmat(1+d*index,d,1)+repmat((-d:-1).',1,lx), d*lx, 1 );
end
% ... and apply nested multiplication:
v = c(index,1).';
for i=2:k
   v = xs.*v + c(index,i).';
end

if ~isempty(NaNx) && k==1 && l>1, v = reshape(v,d,lx); v(:,NaNx) = NaN; end
v = reshape(v,d*mx,nx);

function [index,NaNx] = get_index(mesh,sites,left)
%GET_INDEX appropriate mesh intervals for given ordered data sites

if isempty(left)||left(1)~='l'
   [ignored, index] = histc(sites,[-inf,mesh,inf]);
   NaNx = find(index==0); index = min(index,numel(mesh)+1);
else
   [ignored, index] = histc(-sites,[-inf,-fliplr(mesh),inf]);
   NaNx = find(index==0); index = max(numel(mesh)+2-index,1);
end
