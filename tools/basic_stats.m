classdef basic_stats
    % basic_stats: calculates basic statistics for a matrix along a given
    % dimension ignoring nan's
    %   Calculates nanmean, nanmedian, nanstd
    
    properties (SetAccess = private, Hidden = false)
        mean
        median
        std
    end
    
    methods
        function obj = basic_stats(A,dim)
            % basic_stats: Construct an instance of this class
            %   calculate mean, median, std
            
            obj.mean    = nanmean(A,dim);
            obj.median  = nanmedian(A,dim);
            obj.std     = nanstd(A,[],dim);
        end
    end
end

