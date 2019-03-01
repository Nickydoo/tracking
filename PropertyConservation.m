function [varargout] = PropertyConservation(varargin)

% calcule la conservation d'une quantit� normalis�e
% Nicolas f�vrier 2019

    if nargout == nargin+1


        N = varargin{end};
        for idx = 1:nargin-1
            variable = varargin{idx};
            SUM = sum(sum([variable{:}]));
            for idxN = 1:N
                varargout{idx}(idxN) = mean(diff(variable{idxN})/SUM).^2;
            end
            varargout{nargout-1}(idx) = mean(varargout{idx});
            varargout{nargout}(idx) = sqrt(mean((varargout{idx}-varargout{nargout-1}(idx)).^2));
        end
        

    else
        error('pas le m�me nombre d''entr�e que de sortie')
    end
end


