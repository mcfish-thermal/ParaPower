classdef PPMatIBC_SimInterface < PPMatIBC
    %Derived IBC PPMat Object that is intended to interface with a simulink
    %model running a configured hookPPT System object block.  Has same
    %behavior as standard IBC under standalone use, but is flagged for 
    %override and I/O by the simulink model.
    
    %For reference, hookPPT writes h_IBC and T_ibc, and returns the aggregate 
    %heat rate rejected/accepted by the model at each timestep
    
    methods
        
        function [obj]=PPMatIBC_SimInterface(varargin)
            %Default Values
            Type='IBC_SimInterface';

            if nargin == 1 & ~iscell(varargin{1})
                Name=varargin{1};
                varargin={'name' Name};
            elseif nargin == 1 & iscell(varargin{1})
                varargin=varargin{1};
            end
            obj=obj@PPMatIBC([ 'type' Type varargin]);
            
            PropValPairs=obj.PropValPairs;
            obj.PropValPairs={};
            while ~isempty(PropValPairs) 
                [Prop, PropValPairs]=obj.Pop(PropValPairs);
                if ~ischar(Prop)
                    error('Property name must be a string.');
                end
                Pl=length(Prop);
                %disp(Prop)
                switch lower(Prop)
%                     case obj.strleft('name',Pl)
%                         [Value, PropValPairs]=obj.Pop(PropValPairs); 
%                         Name=Value;

                    otherwise
                        [Value, PropValPairs]=obj.Pop(PropValPairs); 
                        obj.PropValPairs=[Prop Value obj.PropValPairs ];
                end
            end
            obj.CheckProperties(mfilename('class'));  %Checks for left over properties.

%            obj.h_ibc=h_ibc;
%            obj.T_ibc=T_ibc;
        end
    end
end
