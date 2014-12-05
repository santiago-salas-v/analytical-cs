function a_mk   = getInteractionParameters()
%INTERACTION PARAMETERS returns interaction parameters between UNIFAC 
% subgroups
[~,interaction_parameters_chars]    = ...
    cargarCSV('interaction_parameters.csv');
a_mk                                = ...
    cellfun(@str2num,interaction_parameters_chars);
end