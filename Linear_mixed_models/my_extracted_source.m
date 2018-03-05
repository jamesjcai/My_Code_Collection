
        function fun = makeObjectiveFunctionForMinimization(slme)
%makeObjectiveFunctionForMinimization - Get objective function to minimize.
%   fun = makeObjectiveFunctionForMinimization(slme) takes an object slme
%   of type StandardLinearMixedModel and depending on slme.FitMethod makes
%   either the objective function for ML or REML.

            % (1) Make an objective function for minimization. 
            %     (a) If slme.FitMethod is 'ML' this is the negative 
            %     profiled log likelihood (beta and sigma profiled out).
            %
            %     (b) If slme.FitMethod is 'REML' this is the negative 
            %     profiled restricted log likelihood (sigma profiled out).
            switch lower(slme.FitMethod)
                case 'ml'
                    fun = makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta(slme);
                case 'reml'
                    fun = makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta(slme);
                otherwise 
                    % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                    error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
            end
            
        end % makeObjectiveFunctionForMinimization.
        

        function fun = makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta(slme)
%makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta - This returns the
%   negative profiled log likelihood with beta and sigma profiled out as a
%   function of theta.

            fun = @f0;            
            function y0 = f0(theta)                
                 L = BetaSigmaProfiledLogLikelihood(slme,theta);
                y0 = -1*L;
                % Don't let y0 be equal to -Inf.
                y0 = max(-realmax,y0);
            end
            
        end % end of makeNegBetaSigmaProfiledLogLikelihoodAsFunctionOfTheta.
       
        function fun = makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta(slme)
%makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta - This
%   returns the negative profiled restricted log likelihood with sigma
%   profiled out as a function of theta.

            fun = @f1;
            function y1 = f1(theta)
                 L = SigmaProfiledRestrictedLogLikelihood(slme,theta);
                y1 = -1*L;
                % Don't let y1 be equal to -Inf.
                y1 = max(-realmax,y1);
            end
            
        end % end of makeNegSigmaProfiledRestrictedLogLikelihoodAsFunctionOfTheta.
        
        function fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma(slme)
%makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma - This returns
%   the negative profiled log likelihood with beta profiled out as a
%   function of x = [theta;log(sigma)].
           
            fun = @f2;
            function y2 = f2(x)
                % (1) Extract theta and logsigma from x.
                   theta = x(1:end-1);
                logsigma = x(end);
                % (2) Evaluate BetaProfiledLogLikelihood at (theta,sigma)
                % and negate the result.
                   sigma = exp(logsigma);
                 L = BetaProfiledLogLikelihood(slme,theta,sigma);
                y2 = -1*L;
            end
            
        end % end of makeNegBetaProfiledLogLikelihoodAsFunctionOfThetaLogSigma.
        
        function fun = makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma(slme)
%makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma - This returns the
%   negative restricted log likelihood as a function of x =
%   [theta;log(sigma)].
            
            fun = @f3;
            function y3 = f3(x)
                % (1) Extract theta and logsigma from x.
                   theta = x(1:end-1);
                logsigma = x(end);
                % (2) Evaluate RestrictedLogLikelihood at (theta,sigma) and
                % negate the result.
                   sigma = exp(logsigma);
                   L = RestrictedLogLikelihood(slme,theta,sigma);
                  y3 = -1*L;                
            end
            
        end % end of makeNegRestrictedLogLikelihoodAsFunctionOfThetaLogSigma.
        
        function fun = makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma(slme)
%makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma - This returns
%   the negative profiled log likelihood with beta profiled out as a
%   function of x = [eta;log(sigma)] where eta is the Natural parameter
%   vector.

            fun = @f6;
            function y6 = f6(x)
                % (1) Extract eta, logsigma and sigma from x.
                     eta = x(1:end-1);
                logsigma = x(end);
                if slme.isSigmaFixed
                    sigma = slme.sigmaFixed;
                else
                    sigma = exp(logsigma);
                end
                % (2) Set sigma and eta in slme.Psi.
                Psi = slme.Psi;
                Psi = setSigma(Psi,sigma);
                Psi = setNaturalParameters(Psi,eta);
                % (3) Get theta from Psi.
                theta = getUnconstrainedParameters(Psi);
                % (4) Evaluate BetaProfiledLogLikelihood at (theta,sigma)
                % and negate the result.
                L = BetaProfiledLogLikelihood(slme,theta,sigma);
                y6 = -1*L;                
            end

        end % end of makeNegBetaProfiledLogLikelihoodAsFunctionOfEtaLogSigma.
        
        function fun = makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma(slme)
%makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma - This returns the
%   negative restricted log likelihood as a function of x =
%   [eta;log(sigma)] where eta is the Natural parameter vector.
            
            fun = @f7;
            function y7 = f7(x)
                % (1) Extract eta, logsigma and sigma from x.
                     eta = x(1:end-1);
                logsigma = x(end);
                if slme.isSigmaFixed
                    sigma = slme.sigmaFixed;
                else
                    sigma = exp(logsigma);
                end
                % (2) Set sigma and eta in slme.Psi.
                Psi = slme.Psi;
                Psi = setSigma(Psi,sigma);
                Psi = setNaturalParameters(Psi,eta);
                % (3) Get theta from Psi.
                theta = getUnconstrainedParameters(Psi);
                % (4) Evaluate RestrictedLogLikelihood at (theta,sigma) and
                % negate the result.
                L = RestrictedLogLikelihood(slme,theta,sigma);
                y7 = -1*L;
            end

            
        end % end of makeNegRestrictedLogLikelihoodAsFunctionOfEtaLogSigma.
        
        end
        
        
        function PRLogLik = SigmaProfiledRestrictedLogLikelihood(slme,theta)
%SigmaProfiledRestrictedLogLikelihood - Compute the profiled restricted log likelihood.
%   PRLogLik = SigmaProfiledRestrictedLogLikelihood(slme,theta) returns the
%   profiled restricted log likelihood of a StandardLinearMixedModel slme
%   at the specified value of theta with sigma profiled out.

            if slme.isSigmaFixed
                PRLogLik = RestrictedLogLikelihood(slme,theta,slme.sigmaFixed);
            else
                % (1) Using the given value of theta, solve mixed model
                % equations to get r2, R and R1.
                [~,~,~,r2,R,R1] = solveMixedModelEquations(slme,theta);

                % (2) Get log abs determinant of upper triangular R and lower
                % triangular R1.
                 logAbsDetR = slme.logAbsDetTriangular(R);
                logAbsDetR1 = slme.logAbsDetTriangular(R1);

                % (3) Compute profiled restricted log likelihood.
                N = slme.N;
                p = slme.p;
                PRLogLik = (-(N-p)/2) * ( 1 + log( 2*pi*r2/(N-p) ) ) ...
                    - logAbsDetR - logAbsDetR1;            
            end    
            


        function [betaHat,bHat,sigmaHat,r2,R,R1,Deltab] = solveMixedModelEquations(slme,theta)
%solveMixedModelEquations - Solve mixed model equations.
%   [betaHat,bHat,sigmaHat,r2,R,R1] = solveMixedModelEquations(slme,theta)
%   takes a StandardLinearMixedModel object slme and solves the mixed model
%   equations at the specified value of theta. First, we solve mixed model
%   equations using theta to get betaHat and bHat. Then we get sigmaHat
%   (for ML or REML), r2, R and R1. The meaning of r2, R and R1 is defined
%   in internal LME specs. Deltab is the estimated normalized posterior
%   mode of the random effects vector.
    
            % (1) Get X, y, Z, N, p and q from slme. Also get various
            % matrix and vector products to reuse precomputed values.
            X = slme.X;
            y = slme.y;
            Z = slme.Z;
            N = slme.N;
            p = slme.p;
            q = slme.q;            
            XtX = slme.XtX;
            Xty = slme.Xty;
            XtZ = slme.XtZ;
            Zty = slme.Zty;
            ZtZ = slme.ZtZ;
            
            % (2) Get slme.Psi, set current theta and set sigma = 1.
            Psi = slme.Psi;
            Psi = setUnconstrainedParameters(Psi,theta);
            Psi = setSigma(Psi,1);

            % (3) Get lower triangular Cholesky factor of D matrix in Psi.
            % Lambda has size q by q where q = size(Z,2).
            Lambda = getLowerTriangularCholeskyFactor(Psi);

            % (4) Compute the U matrix. U is N by q.
            %U = Z*Lambda;
            
            % TODO: Try to make steps 5 onward more efficient.
            % (5) Get R and S such that S*R'*R*S' = U'*U + eye(q) where R
            % and S are q by q. S is a permutation matrix: S*S' = eye(q).
            Iq = spdiags(ones(q,1),0,q,q);
            %[R,status,S] = chol(U'*U + Iq);
            [R,status,S] = chol(Lambda'*ZtZ*Lambda + Iq);
            
            % (6) Make sure status is zero.
            if (status ~= 0)
                % <entry key="ErrorSparseCholesky">An error occured during sparse Cholesky factorization.</entry>
                error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:ErrorSparseCholesky'));
            end
            
            % For steps (7) onwards, we may get a rank deficient X or a
            % singular factor R1. Hence we disable the nearly singular
            % matrix warning and then enable it at the end of this method
            % or on error.
            warnState = warning('query','all');
            warning('off','MATLAB:nearlySingularMatrix');
            warning('off','MATLAB:illConditionedMatrix');
            warning('off','MATLAB:singularMatrix');
            cleanupObj = onCleanup(@() warning(warnState));
            
            % (7) Compute P1 (q by q), Q1 (p by q) and R1 (p by p).
            %P1 = S*R';
            %Q1 = ((X'*U)*S) / R;
            Q1 = ((XtZ*Lambda)*S) / R;
            
            R1R1t = XtX - Q1*Q1';
            try
                %R1 = chol(X'*X - Q1*Q1','lower');
                %R1 = chol(XtX - Q1*Q1','lower');
                R1 = chol(R1R1t,'lower');
            catch ME %#ok<NASGU>
                % We know from theory that R1 must exist.
                %R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(X'*X - Q1*Q1');
                R1 = classreg.regr.lmeutils.covmats.CovarianceMatrix.singularLowerChol(R1R1t);
            end
                        
            % (8) Compute cDeltab (q by 1) and cbeta (p by 1).
            %cDeltab = P1 \ (U'*y);
            %cDeltab = R' \ (S'*(U'*y));
            cDeltab = R' \ (S'*((Lambda'*Zty)));
            %cbeta = R1 \ (X'*y - Q1*cDeltab);
            cbeta = R1 \ (Xty - Q1*cDeltab);
            
            % (9) Compute betaHat and Deltab.
            betaHat = R1' \ cbeta;
            %Deltab = P1' \ (cDeltab - Q1'*betaHat);
            Deltab = S*(R \ (cDeltab - Q1'*betaHat));

            % (10) Compute bHat.
            bHat = Lambda * Deltab;
            
            % (11) Compute r2.
            r2 = sum(Deltab.^2) + sum((y - X*betaHat - Z*bHat).^2);
            
            % (12) Compute sigmaHat.
            if slme.isSigmaFixed
                sigmaHat = slme.sigmaFixed;
            else
                switch lower(slme.FitMethod)
                    case 'ml'
                        sigma2 = r2/N;
                        sigmaHat = sqrt(sigma2);
                    case 'reml'
                        sigma2 = r2/(N - p);
                        sigmaHat = sqrt(sigma2);
                    otherwise
                        % <entry key="BadFitMethod">''FitMethod'' must be either ''ML'' or ''REML''.</entry>
                        error(message('stats:classreg:regr:lmeutils:StandardLinearMixedModel:BadFitMethod'));
                end
            end
            % NOTE: Restore the warning state to that saved in warnState.
            % This is automatically done by the cleanupObj going out of
            % scope. warning(warnState);
            
        end % solveMixedModelEquations.           
            