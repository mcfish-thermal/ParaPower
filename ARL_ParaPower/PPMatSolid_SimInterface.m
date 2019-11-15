classdef PPMatSolid_SimInterface < PPMatSolid
    %Derived Solid PPMat Object that is intended to interface with a simulink
    %model running a configured hookPPT System object block.  Has same
    %behavior as standard Solid under standalone use, but is flagged for 
    %override and I/O by the simulink model
    
    %For reference, hookPPT writes a Q(t) vector distributed across elements of this material
    %and returns the volume- (3D features) or area- (2D features) averaged
    %temperature.  Use of more than one feature per material instance not
    %supported.
    
    methods
        function [obj]=PPMatSolid_SimInterface(varargin)
            %Default Values
            Type='Solid_SimInterface';
            
            if nargin == 1 & ~iscell(varargin{1})
                Name=varargin{1};
                varargin={'name' Name};
            elseif nargin == 1 & iscell(varargin{1})
                varargin=varargin{1};
            end
            obj=obj@PPMatSolid([ 'type' Type varargin]);
            
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

        end
    end
end
