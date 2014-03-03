function [resultado,Datos]=cargarCSV(nombreCompleto_csv)
fid=fopen(nombreCompleto_csv);
if fid < 3
    ME=MException('FILE_READ:ERROR',['Cannot open file. ',...
        ' Is it open outside of Matlab?']);
    throw(ME);
end
i=1;
primeraLinea=deblank(fgetl(fid));
Datos=[{} {} {}];
% First line always has the separating character at the end
% If not by stating 'sep=|', the first line always has empty
% rows save for the first. Meaning a CSV of one column is not
% supported.
sepChar=primeraLinea(end);
if ~strncmp(primeraLinea,'sep=',4)    
    Datos(i,:)=regexp(primeraLinea,['\',sepChar],'split');
    i=i+1;
else    
    sepChar=primeraLinea(end);
end
while ~feof(fid)
    Datos(i,:)=regexp(fgetl(fid),['\',sepChar],'split');
    i=i+1;
end
fclose(fid);
exito=false;
resultado=exito;
end