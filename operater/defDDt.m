    function [D,Dt] = defDDt
        % defines finite difference operator D
        % and its transpose operator
        
        D = @(U) ForwardD(U);
        Dt = @(X,Y) Dive(X,Y);
        
        function [Dux,Duy] = ForwardD(U)
            % Forward finite difference operator
            Dux = [diff(U,1,2), U(:,1) - U(:,end)];
            Duy = [diff(U,1,1); U(1,:) - U(end,:)];
        end
        
        function DtXY = Dive(X,Y)
            % Transpose of the forward finite difference operator
            DtXY = [X(:,end) - X(:, 1), -diff(X,1,2)];
            DtXY = DtXY + [Y(end,:) - Y(1, :); -diff(Y,1,1)];
        end
    end