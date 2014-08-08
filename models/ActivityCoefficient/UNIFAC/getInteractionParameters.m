function a_mk   = getInteractionParameters()
%INTERACTION PARAMETERS
[~,interaction_parameters_chars]    = ...
    cargarCSV('interaction_parameters.csv');
a_mk                                = ...
    cellfun(@str2num,interaction_parameters_chars);
end