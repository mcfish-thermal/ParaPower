%Template: Describe what this test case consists of
%warning('This is a template only and will not result in a usable case')    

%Clear the main variables that are passed out from it.
clear Features ExternalConditions Params PottingMaterial Descr MatLib

MatLib=PPMatLib;
MatLib.AddMatl(PPMatPCM(  'name'  , 'Ga'  ...
                           ,'k_l'   , 24    ...
                           ,'rho_l' , 6093  ...
                           ,'cp_l'  , 397   ...
                           ,'lf'    , 80300 ...
                           ,'tmelt' , 29.8  ...
                           ,'cte'   , 0     ...
                           ,'E'     , 0     ...
                           ,'nu'    , 0     ...
                           ,'k'     , 33.7  ...
                           ,'rho'   , 5903  ...
                           ,'cp'     , 340   ...
                        )) ;
MatLib.AddMatl(PPMatSolid('name'  , 'cu'  ...
                           ,'cte'   , 2.4e-5...
                           ,'E'     , 1.1e11...
                           ,'nu'    , .37   ...
                           ,'k'     , 390   ...
                           ,'rho'   , 8900  ...
                           ,'cp'     , 390   ...
                        )) ;
                    
MatLib.AddMatl(PPMatSolid('name'  , 'SiC'  ...
                           ,'cte'   , 4e-6 ...
                           ,'E'     , 4.1e11...
                           ,'nu'    , .14   ...
                           ,'k'     , 120   ...
                           ,'rho'   , 3100  ...
                           ,'cp'     , 750   ...
                        )) ;
                    
MatLib.AddMatl(PPMatIBC('name'  , 'IBC_1'  ...
                           ,'h_ibc'   , 100 ...
                           ,'t_ibc'     , 0 ...
                        )) ;

                    MatLib.AddMatl(PPMatNull('name'  , 'No_Matl'  ...
                          ));
Features.x=[]; Features.y=[]; Features.z=[]; Features.Matl=[]; Features.Q=[]; Features.Matl=''; 
Features.dz=0; Features.dy=0; Features.dz=0;

Desc='Cyl_quarter_annulus_periodic';  %Description of the test case

ExternalConditions.h_Xminus=0;
ExternalConditions.h_Xplus =100;
ExternalConditions.h_Yminus=0;
ExternalConditions.h_Yplus =0;
ExternalConditions.h_Zminus=0;
ExternalConditions.h_Zplus =0;

ExternalConditions.Ta_Xminus=0;
ExternalConditions.Ta_Xplus =0;
ExternalConditions.Ta_Yminus=0;
ExternalConditions.Ta_Yplus =0;
ExternalConditions.Ta_Zminus=0;
ExternalConditions.Ta_Zplus =0;

ExternalConditions.Tproc=280; %Processing temperature, used for stress analysis

Params.Tinit     = 20;  %Initial temperature of all materials
Params.DeltaT    = 10; %Time step size, in seconds
Params.Tsteps    = 1000; %Number of time steps.

PottingMaterial  = 0;  %Material that surrounds features in each layer as defined by text strings in matlibfun. 
                       %If Material is 0, then the space is empty and not filled by any material.

%Each feature is defined separately.  There is no limit to the number of
%features that can be defined.  For each layer of the model, unless a
%feature exists, the material is defined as "potting material."  There is
%no checking to ensure that features do not overlap.  The behavior for
%overlapping features is not defined.

PottingMaterial  = 0;

Features(1).x    = [.005 .010];
Features(1).y    = [0 pi/2];
Features(1).z    = [0 .01];
Features(1).dx   = 10;
Features(1).dy   = 2;
Features(1).dz   = 2;
Features(1).Matl = 'Cu';
Features(1).Q    = 0;  %Total watts per features dissipated.

Features(2)      = Features(1);
Features(2).x    = [.005 .005];
Features(2).dx   = 1;
Features(2).Q    = 1.2;
Features(2).Matl = 'No_Matl';

Features(3)      = Features(1);
Features(3).x    = [.010 .020];
Features(3).dx   = 10;
Features(3).Matl = 'Ga';

Features(4)      = Features(1);
Features(4).x    = [.020 .025];
Features(4).dx   = 1;
Features(4).dy   = 6*4;
Features(4).Matl = 'IBC_1';

Features(5)      = Features(1);
Features(5).x    = [.020 .022];
Features(5).dx   = 2;
Features(5).y    = pi/2*[1/6 2/6];
Features(5).dy   = 2;
Features(5).z    = [.003 .006];
Features(5).dz   = 2;
Features(5).Matl = 'CU';
Features(5).Q    = 2;


%     Define Q as @(T,Q,Ti)interp1(T,Q,Ti,'spline');
%     if Q is cell then do interp, otherwise 
%         Q{ii,jj,kk}(Arg1 Arg2 Arg3)

TestCaseModel.Desc=Desc;
TestCaseModel.TCM=PPTCM(ExternalConditions,Features,Params,PottingMaterial,MatLib,'CS','cyl');
% TestCaseModel.TCM.Features=Features;
% TestCaseModel.TCM.Params=Params;
% TestCaseModel.TCM.PottingMaterial=PottingMaterial;
% TestCaseModel.TCM.ExternalConditions=ExternalConditions;
% TestCaseModel.TCM.MatLib=MatLib;
%TestCaseModel.TCM.VariableList=ParamList;

MFILE=mfilename('fullpath');

return

D=ParaPowerGUI_V2('GetResults');
temps=D.R.getState('Thermal');
temps_R=temps(:,1,1,2);

figure;
temps=NewResults.Tprnt(end-3,:,3,end);
thetas=[MI.OriginPoint(2) MI.OriginPoint(2)+cumsum(MI.Y)];
thetas=(thetas(1:end-1)+thetas(2:end))/2;
thetas=thetas*180/pi;
plot([thetas thetas(end)+MI.Y(1)*180/pi],[temps temps(1)])