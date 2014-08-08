function a_mk   = getInteractionParameters()
%INTERACTION PARAMETERS
[~,a_mk_chars]  = cargarCSV('a_mk.csv');
a_mk            = cellfun(@str2num,a_mk_chars);
end