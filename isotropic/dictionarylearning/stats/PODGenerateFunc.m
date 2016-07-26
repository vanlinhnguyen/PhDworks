%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% FUNCTION FOR PERFORMING POD ANALYSIS
%%%
%%% Function usage:
%%%    [N, PODzeta, ev1, PODsai] = PODGenerateFunc(Snapshots) ;
%%%    [N, PODzeta, ev1, PODsai] = PODGenerateFunc(Snapshots, M) ;
%%%
%%%  Input: 1. Snapshots -> The velocity snapshots. Has dim [P x Q] (Q: dim of input vector)
%%%         2. M (optional argument) -> If you want to extract only M                       %%%			eigenfunctions. 
%%%
%%%  Output: 1. N -> Number of POD modes extracted.
%%%          2. PODzeta -> The POD-coefficients. Dim [P x N]
%%%          3. ev1 -> Inverse of PODzeta. [N x P]
%%%          4. PODsai -> POD eigenfunctions. [N X Q]
%%%
%%%% Additional Notes:
%%%%
%%%%  In this function, the intgration is replace by Reimann sum. 
%%%%
%%%%  Here PODzeta*PODsai \approx Snapshots if N<P.
%%%%       PODzeta*PODsai = Snapshots if N==P.
%%%%
%%%%
%%%% The input Snapshots must have a specific form. Each snapshot must %%%% be laid out in the form of a 'row-vector'. If there are 'P'.
%%%% snapshots then there should be 'P' rows. The columns have a size %%%% 'Q' which represents the number of grid-points in the domain.
%%%%
%%%% If there are 2 components of the velocity field ,say [u,v], then %%%% each snapshot will be a row vector with column size '2Q'. This %%%% extends to 3D snapshots.
%%%%
%%%% Since the snapshots are originally "flattened" for POD analysis, %%%% the resulting eigenfunction will also have a flattened form. In %%%% order to get it back into the original form, you will have to use %%%% the reshape command and set it back into the original dimensions.
%%%%
%%%% For more information, contact me at mokhpar@iit.edu.
%%%%
%%%% Paritosh Mokhasi.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Nrsa1, PODzeta, e_vect1, PODsai, d] = PODGenerateFunc(Snapshots,Nrsa)

% function [Nrsa,PODzeta, PODsai] = PODGenerateFunc (Snapshots)

Nt = size(Snapshots,1);

%%%% ================================================
%%% ---- THIS PART IS THE POD ANALYSIS.

D1 = Snapshots*Snapshots'/Nt; 
D1 = triu(D1) + triu(D1,1)'; %%--makes it truely symmetric.
[v d] = eig(D1); 
[d,it] = sort(diag(d),'descend'); 
v = v(:,it); v = v';


if (nargin==1)
    Nrsa = d >0;
    Nrsa1 = sum(Nrsa);
else
    Nrsa = 1:Nrsa;
    Nrsa1 = max(Nrsa);
end

C1 = D1*Nt;
zp = C1*v';
pnorm = diag(sqrt(abs(v*zp))); 
e_vect1 = v./repmat(abs(pnorm),1,size(v,2)); 
 
e_vect1 = e_vect1(Nrsa,:);

PODsai = e_vect1*Snapshots;
PODzeta = C1*e_vect1';

clear v zp C1 D1
% %%% ==================================================
