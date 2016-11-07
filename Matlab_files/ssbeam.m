% example script file for BEAM
% A five element aluminum beam clamped at the left end
% with five elements.  The cross section is 5 cm square.
% No shear deformation is modeled.
% Extensional degrees of freedom are included (not zeroed).
% The force vector is zero since this is a dynamic analysis.
%node=[0 0;.1 0;.2 0;.3 0;.4 0;.5 0];
node=[0 0; .5 0; 1 0];
ncon=[1 2 1 1 1 0 1;
      2 3 1 1 1 0 1];
% ncon=[1 2 69e10 .0025 5.208e-7 0 2700;
%       2 3 69e10 .0025 5.208e-7 0 2700;
%       3 4 69e10 .0025 5.208e-7 0 2700;
%       4 5 69e10 .0025 5.208e-7 0 2700;
%       5 6 69e10 .0025 5.208e-7 0 2700];
zero=[1 1;%      1 2;
      2 1;
      3 1];
%      3 2];
conm=[];
force=[];
save ssbeam.con
