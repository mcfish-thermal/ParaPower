classdef hookPPT < scPPT
    %Interface class between ParaPower Thermal solvers and Simulink System Object Block
    
    %To Do:  How to get I/O implementation to self organize according to obj.MI
    %MI.Matlib is not cultivated according to what is present in MI.Model,
    %just copied wholesale from the parent PPTCM.  This creates a lot of
    %error-prone handholding.
    
    
    properties (Access = protected)%(DiscreteState)
        bdry_watts  %output for PPMatIBC_SimInterface, one per material instance (row) per timestep (col)
        Q_cmd       %input command for heat rate to each Solid_SimInterface
        source_temps %output for PPMatSolid_SimInterface, one per material instance (row) per timestep (col) 
    end
    
    properties (Access = protected,Nontunable)
        Q_iface %cell array of masks of driven Solid_SimInterfaces of size(MI.Model(:))
        h_iface %-matnums of IBC_SimInterfaces for later comparison to fullheader
    end
    
    methods (Access = protected)
        function setupImpl(obj)
            setupImpl@scPPT(obj);
            obj.GlobalTime = -1;
            
           rollcall=unique(obj.MI.Model);
           rollcall=rollcall(rollcall>0); %cant index zero or negative mat numbers
           Types=obj.MI.MatLib.GetParam('Type'); 
           
           %need to map from list to vmatnum so we can pick the appropriate
           %columns from fullheader
           list_IBC_iface=rollcall(strcmp(Types(rollcall),'IBC_SimInterface'));
           vmatnum=rollcall(strcmp(Types(rollcall),'IBC')|strcmp(Types(rollcall),'IBC_SimInterface'));
           [~,h_iface]=ismember(list_IBC_iface,vmatnum);
           obj.h_iface=-h_iface;
           
           
           list_solid_iface=rollcall(strcmp(Types(rollcall),'Solid_SimInterface'));
           %load in complex VA matrix from formed model... will use to
           %create scaling masks for each solid interface material
           Q_iface=repmat({obj.MI.Vol_Area_Matrix},length(list_solid_iface),1);
           
           for si=1:length(list_solid_iface)
               if length(unique(obj.MI.FeatureMatrix(obj.MI.Model==list_solid_iface(si))))>1
                   warning('Solid_SimInterface material being used in multple features')
               end
               Q_iface{si}(obj.MI.Model~=list_solid_iface(si))=0; %zero out non iface mats
               sum_VA=sum(Q_iface{si},'all');
               if real(sum_VA)>0 %material presence has volume
                   Q_iface{si}=sparse(real(Q_iface{si}(:)))/real(sum_VA);  %scaling mask with entries v_el/v_tot
               elseif imag(sum_VA)>0 %material presence has area but no volume
                   Q_iface{si}=sparse(imag(Q_iface{si}(:)))/imag(sum_VA);  %scaling mask with entries a_el/a_tot
               else
                   error('Unexpected volume/area presence of Solid_SimInterface material')
               end
           end
           
           obj.Q_iface=Q_iface;
           
           
          
        end
        
        function [bdry_watts,source_temps,Tres,PHres] = stepImpl(obj,GlobalTime,htcs,Ta_vec,Q_cmd)        
            %Have to assume protocol that size htcs, Ta_vec == size h_iface
            if length(htcs)~=length(obj.h_iface) || length(Ta_vec)~=length(obj.h_iface)
                error('IBC interface variables htc or Ta_vec dimension does not match number of IBC_SimInterface materials present in MI.Model')
            end
            
            %fullheader=[header find(h)]; %fullheader is a rowvector of negative matnums and a subset of 1 thru 6
            %Ta_vec=Ta_void(-(header));
            %Ta_vec=[Ta_vec Ta(h~=0)];  %grab just those Ta corresponding to the active ext boundaries
            %from null_void_init
            
            %Filter by NaN, which indicate to use last known value
            %(possibly material default value)
            Ta_load=Ta_vec(isnan(Ta_vec));
            iface_load=obj.h_iface(isnan(Ta_vec));
            for hi=1:nnz(isnan(Ta_vec))
                obj.Ta_vec(obj.fullheader==iface_load(hi))=Ta_load(hi);
            end
            
            %obj.htcs=[hint(-fullheader(fullheader<0)) h(fullheader(fullheader>0))];  
            %built htcs from header and h lists from setup
            %from conduct_build
            htc_load=htc(isnan(htcs));
            iface_load=obj.h_iface(isnan(htcs));
            for hi=1:nnz(isnan(htcs))
                obj.htcs(obj.fullheader==iface_load(hi))=htc_load(hi);
            end

            
            [Tres, ~, PHres, ~]=stepImpl@scPPT(obj,GlobalTime);
            bdry_watts=obj.bdry_watts;
            source_temps=obj.source_temps;
        end
        
        function obj = pre_step_hook(obj)
            obj.bdry_watts = NaN(length(obj.Ta_vec),length(obj.GlobalTime)-1);
            
        end  
        
        function obj = post_step_hook(obj)
            %derive and overload me to insert postprocessing hook
            flux=obj.ExternalFlux(obj.T,obj.B,obj.Ta_vec);
            step=find(isnan(obj.bdry_watts(1,:)),1);
            obj.bdry_watts(:,step)=sum(flux,1)';
        end
        
        function resetImpl(obj)
            % Initialize / reset discrete-state properties
        end
    end
    
    methods (Static)
        % Overload Qinit to incorporate Solid_SimInterface materials
        function Qv=Qinit(Q,GlobalTime,Qmask)
    end
    
    
    %% Simulink I/O methods
    methods (Access = protected)
        function icon = getIconImpl(obj)
            % Define icon for System block
            icon = mfilename("sPPT"); % Use class name
            % icon = "My System"; % Example: text icon
            % icon = ["My","System"]; % Example: multi-line text icon
            % icon = matlab.system.display.Icon("myicon.jpg"); % Example: image file icon
        end     
        
        function [outsz_1,outsz_2,outsz_3] = getOutputSizeImpl(obj)
            numsteps = propagatedInputSize(obj,1);
            if isempty(obj.GlobalTime)
                outsz_1 = [10 1];  %Need to case-handle first iteration.
            else
                outsz_1 = [10 numsteps(2)];  %10 being the number of nodes in the input model
            end
            outsz_2=[3,10,4];
            outsz_3=outsz_2;
        end
        
        function [outtype_1,outtype_2,outtype_3] = isOutputFixedSizeImpl(obj)
            outtype_1 = propagatedInputFixedSize(obj,1);
            outtype_2 = true;
            outtype_3 = true;
        end
        
        function [type_1,type_2,type_3] = getOutputDataTypeImpl(obj)
            type_1 = 'double';
            type_2 = type_1;
            type_3 = type_1;
        end
        
        function [cx1,cx2,cx3] = isOutputComplexImpl(obj)
            cx1 = false;
            cx2 = false;
            cx3 = false;
        end
    end
    
    %% Simulink customization functions
    methods(Static, Access = protected)
        function header = getHeaderImpl
            % Define header panel for System block dialog
            header = matlab.system.display.Header(mfilename("class"));
        end

        function group = getPropertyGroupsImpl
            % Define property section(s) for System block dialog
            group = matlab.system.display.Section(mfilename("class"));
        end
    end
    end