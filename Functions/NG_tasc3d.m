function icode = NG_tasc3d(model,c,nu,s, param, nbfam)

%%% Not working well with impricated model %%%

% Number of variable
nvar = size(model,1);

% Parameters tau that ensures a continuous transition of C(h,tau)
if nbfam ==1
    tau = 0;
else
    tau = 0 : 1/(nbfam-1) : 1;
end

% Compute nffam set of familly and applied transformation one step at a time.
icode = zeros(nbfam,1);
for nf = 1 : nbfam
    nf
    % According to param, we change the multivariace covariance model
    % Here, only range are program
    for i = 1 : nvar
        for j = i:nvar
            
            if isfield(param, 'rangeX')
                if ~isempty(param.rangeX{i,j})
                    if length(param.rangeX{i,j})==2
                        model{i,j}(:,2) = (param.rangeX{i,j}(2)-param.rangeX{i,j}(1))*tau(nf)+param.rangeX{i,j}(1);
                    elseif length(param.rangeX{i,j})==nbfam
                        model{i,j}(:,2) = param.rangeX{i,j}(nf);
                    end
                end
            end

            if isfield(param, 'rangeY')
                if length(param.rangeY{i,j})==2
                    model{i,j}(:,3) = (param.rangeY{i,j}(2)-param.rangeY{i,j}(1))*tau(nf)+param.rangeY{i,j}(1);
                elseif length(param.rangeY{i,j})==nbfam
                    model{i,j}(:,3) = param.rangeY{i,j}(nf);
                end
            end

            if isfield(param, 'rotX')
                if length(param.rangeY{i,j})==2
                    model{i,j}(:,5) = (param.rotX{i,j}(2)-param.rotX{i,j}(1))*tau(nf)+param.rotX{i,j}(1);
                elseif length(param.rotX{i,j})==nbfam
                    model{i,j}(:,5) = param.rotX{i,j}(nf);
                end
            end

            if isfield(param, 'nu')
                if length(param.nu{i,j})==2
                    nu{i,j} = (param.nu{i,j}(2)-param.nu{i,j}(1))*tau(nf)+param.nu{i,j}(1);
                elseif length(param.nu{i,j})==nbfam
                    nu{i,j} = param.nu{i,j}(nf);
                end
            end

            if isfield(param, 'c')
                if length(param.c{i,j})==2
                    c{i,j} = (param.c{i,j}(2)-param.c{i,j}(1))*tau(nf)+param.c{i,j}(1);
                elseif length(param.c{i,j})==nbfam
                    c{i,j} = param.c{i,j}(nf);
                end
            end

            if j ~= i
                model{j,i} = model{i,j}; nu{j,i} = nu{i,j}; c{j,i} = c{i,j};
            end
        end
    end
    [icode(nf),~]=tasc3d(model, c, nu, s);
end
