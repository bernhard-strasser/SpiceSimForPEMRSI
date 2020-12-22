function params = setDefaultField(params, fieldName, vDefault)
    if ~isfield(params,fieldName) || isempty(getfield(params,fieldName))
        params  = setfield(params,fieldName,vDefault);
    end
end