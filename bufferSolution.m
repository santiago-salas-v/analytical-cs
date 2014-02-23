%% BUFFER SOLUTIONS
% Calcula disoluciones tamp�n / buffer, basado en el balance
% de materia y equilibrios �cido - base.
%% Reacciones
% $\begin{array}{ccccccc|cc} 
% HA & + & H_2O   & {\rightarrow \atop \leftarrow} &
% H_3O^{(+)} & + & A^{(-)} & 
% Ka = C_{H_3Oeq^{(+)}} \times C_{Aeq^{(-)}} / 
%      C_{HAeq}  & \\
%    & + & H_2O & {\rightarrow \atop \leftarrow} &
% H_3O^{(+)} & + &  OH^{(-)} & 
% Kw = C_{H_3Oeq^{(+)}} \times C_{OHeq^{(-)}}  = 
%      10^{-14} & = K_{ca} \times C_{H_2Oeq}^2 \\
% A^{(-)} & + & H_2O   & {\rightarrow \atop \leftarrow} &
% HA &  + & OH^{(-)}  & 
% Kb = C_{OHeq^{(-)}} \times C_{HAeq} / 
%      C_{Aeq^{(-)}}  & = K_{cb} \times C_{H_2Oeq} \\
% \end{array}$
%% Condiciones iniciales
% * Agregar buffer en forma de sal + �cido: 
% ${C_{HA_0}V_0 \over Vr} + {C_{NaA_0}V_0 \over Vr}$ .
% La sal se disocia
% completamente, el �cido de acuerdo con el equilibrio �cido-base.
% * Partir de pH 7. 
% * Par�metros:
% $Vr$, $C_{HA_inicial}=C_{HA_0}V_0/Vr$, 
% $C_{A^{(-)}_inicial}=C_{NaA_0}V_0/Vr$ (disociaci�n completa), 
% $C_{H_3O_0^{(+)}}=10^{-7}mol/L$, $C_{OH_0^{(-)}}=10^{-7}mol/L$.
%% Ecuaciones (por balance + equilibrio)
% Sistema algebr�ico no lineal de 4 ecuaciones *independientes*
% con 4 variables. 
% $\begin{array}{cccccccccc} 
% {C_{HA_0}V_0 \over Vr} & = & + & C_{HAeq} & + & 
% 0 & - & C_{H_3Oeq^{(+)}} & + & C_{OHeq^{(-)}}
% \\
% {C_{NaA_0}V_0 \over Vr} & = & + & 0 & + & 
% C_{Aeq^{(-)}} & + & C_{H_3Oeq^{(+)}} & - & C_{OHeq^{(-)}}
% \\
% 0 & = & + & Ka \times C_{Aeq^{(-)}} & - & 
% C_{Aeq^{(-)}} & \times & C_{H_3Oeq^{(+)}} & + & 0
% \\
% Kw & = & + & 0 & + & 0 & 
% + & C_{H_3Oeq^{(+)}} & \times & 
% C_{OHeq^{(-)}} \\
% \end{array}$
%% Soluci�n
% Reducci�n a 1 ecuaci�n de 3er orden.
% 
% $\begin{array}{cccccccccccc} 
% 0 & = & + & \left[ {C_{H_3Oeq^{(+)}} \over Ka} \right]^3 &
% + & \left( 1-{1 \over Ka} \times 
% {C_{NaA_0}V_0 \over Vr} \right) &   
% \left[ {C_{H_3Oeq^{(+)}} \over Ka} \right]^2 \\
% & & - & \left( {Kw \over Ka^2}-{1 \over Ka} \times 
% {C_{HA_0}V_0 \over Vr} \right) \times   
% \left[ {C_{H_3Oeq^{(+)}} \over Ka} \right] & 
% - & {Kw \over Ka^2 }
% \end{array}$
%
% * Fijar alimentaci�n.
% * Resolver ecuaci�n c�bica para 
% $C_{H_3Oeq^{(+)}}$ , calcular pH.
% * Calcular $C_{OHeq^{(-)}}$ por ec. 4.
% * Calcular $C_{Aeq^{(-)}}$ y $C_{HAeq}$ por ecs. 1,2.
function bufferSolution()
% Generar ventana
screenSize  = get(0,'MonitorPositions');
Ventana     = findobj('Tag','BufferSolutionGui');
pageicon    = imread(['utils',filesep,'help_ug.png'],...
    'BackgroundColor',[1,1,1]);
saveicon    = imread(['utils',filesep,'file_save.png'],...
    'BackgroundColor',[1,1,1]);
if isempty(Ventana) || ~ishandle(Ventana) || ~isscalar(Ventana)
    Ventana = figure('Menu','none','Toolbar','none',...
        'Color',[255,255,255]/255,'Name','Buffer Solution',...
        'WindowStyle','normal','NumberTitle','off','Resize','off',...
        'Units','pixels',...
        'Position',[screenSize(3:4)/2-[780,400]/2,780,400],...
        'Tag','BufferSolutionGui');
    Panel2          = uipanel('Parent',Ventana,...
        'Title','Buffer Solution',...
        'BackgroundColor',[255,255,255]/255);
    Tabla_Entrada   = uitable('Parent',Panel2,...
        'Units','normalized',...
        'Position',[0,1/4,30/100,3/4-1/8],...
        'RowName',{},...
        'ColumnName',{'VARIABLE','VALOR','UNIDADES'},...
        'ColumnWidth',{  80 'auto' 'auto' },...
        'ColumnEditable',[ false true false ]);
    Tabla_Salida   = uitable('Parent',Panel2,...
        'Units','normalized',...
        'Position',[1-30/100,1/4,30/100,3/4-1/8],...
        'RowName',{},...
        'ColumnName',{'VARIABLE','VALOR','UNIDADES'},...
        'ColumnWidth',{  80 'auto' 'auto' },...
        'ColumnEditable',[ false false false ]);
    texto1  = uicontrol('Parent',Panel2,...
        'Units','normalized',...
        'Style','text',...
        'String','PAR�METROS',...
        'BackgroundColor','white',...
        'Position',[0,1-1/8+1/4/6,1/6,1/4/6]);
    texto2  = uicontrol('Parent',Panel2,...
        'Units','normalized',...
        'Style','text',...
        'String','RESPUESTA',...
        'BackgroundColor','white',...
        'Position',[1-30/100+1/6/2,1-1/8+1/4/6,1/6,1/4/6]);
    Botones = uibuttongroup('Parent',Panel2,...
        'Units','normalized',...
        'Title','RAICES',...
        'BackgroundColor','white',...
        'Position',[1-30/100,1/4/10,30/100,9/4/10]);
    Boton_1 = uicontrol('Parent',Botones,...
        'Units','normalized',...
        'Style','radiobutton',...
        'String','pH=',...
        'BackgroundColor','white',...
        'Position',[1/22,1/11,1-1/22,3/11],...
        'Enable','off');
    Boton_2 = uicontrol('Parent',Botones,...
        'Units','normalized',...
        'Style','radiobutton',...
        'String','pH=',...
        'BackgroundColor','white',...
        'Position',[1/22,4/11,1-1/22,3/11],...
        'Enable','off');
    Boton_3 = uicontrol('Parent',Botones,...
        'Units','normalized',...
        'Style','radiobutton',...
        'String','pH=',...
        'BackgroundColor','white',...
        'Position',[1/22,7/11,1-1/22,3/11],...
        'Enable','off');
    DatosIniciales      =    ...
        {...
        'COND. IN./ALIM.' ''          '';...
        'C_HA0'         0.3         'gmol/L';...
        'C_A0'          0.2         'gmol/L';...
        'C_HAinicial'   0.3         'gmol/L';...
        'C_Ainicial'    0.2         'gmol/L';...
        'C_H3Oinicial'  10^-7       'gmol/L';...
        'C_OHinicial'   10^-7       'gmol/L';...
        'Ka'            5.70*10^-10 '';...
        'Kw'            10^-14      '';...
        'MEDIO REAC.'   ''          '';...
        'V0'            0.2         'L';...
        'Vr'            0.2         'L';...
        };
    DatosDeRespuesta      =    ...
        {...
        'pH'            []      '';...
        'pOH'           []      '';...
        'C_HA'          []      'gmol/L';...
        'C_A'           []      'gmol/L';...
        'C_H3O'         []      'gmol/L';...
        'C_OH'          []      'gmol/L';...
        'pKb'           []      '';...
        'pKa'           []      '';...
        };
    BotonCalcular = uicontrol('Parent',Panel2,...
        'Units','normalized',...
        'Style','pushbutton',...
        'String','CALCULAR',...
        'Position',[30/100/2-1/6/2,1/4/2-1/4/3/2,1/6,1/4/3]);
    Herramientas    = uitoolbar(Ventana);
    BotonDeReporte  = uipushtool(Herramientas,...
        'CData',            pageicon,...        
        'Tag',              'uipushtool1',...
        'TooltipString',    'Ver procedimiento de soluci�n');
    BotonDeGuardar  = uipushtool(Herramientas,...
        'CData',            uint8(saveicon),...        
        'Tag',              'uipushtool1',...
        'TooltipString',    'Ver procedimiento de soluci�n');
    handles = struct('Tabla_Entrada',Tabla_Entrada,...
        'Tabla_Salida',Tabla_Salida,'texto1',texto1,...
        'texto2',texto2,'Panel2',Panel2,...
        'Ventana',Ventana,'Botones',Botones,...
        'Boton_1',Boton_1,'Boton_2',Boton_2,...
        'Boton_3',Boton_3,'BotonDeReporte',BotonDeReporte,...
        'Herramientas',Herramientas,...
        'BotonDeGuardar',BotonDeGuardar);
    handles.BotonCalcular = BotonCalcular;
    set(Tabla_Entrada,'Data',DatosIniciales);
    set(Tabla_Salida,'Data',DatosDeRespuesta);
    set(BotonCalcular,...
        'Callback',{@funcionDeCalcular,handles});
    set(Tabla_Entrada,...
        'CellEditCallback',{@funcionDeCalcular,handles});
    set(Botones,...
        'SelectedObject',[],'SelectionChangeFcn',...
        {@CambioDeBoton,handles});
    set(BotonDeReporte,'ClickedCallback',...
        {@GenerarReporte,handles});
    funcionDeCalcular(Tabla_Entrada,{},handles);
else
    figure(Ventana);
end
end

function funcionDeCalcular(~,~,handles,varargin)
Datos_entrada   = get(handles.Tabla_Entrada,'Data');
Datos           = ...
    Datos_entrada(cellfun(@isvarname,Datos_entrada(:,1)),:);
Datos_s         = cell2struct(Datos(:,2),Datos(:,1),1);
Respuesta       = get(handles.Tabla_Salida,'Data');

C_HA0           = Datos_s.C_HA0;
C_A0            = Datos_s.C_A0;
Ka              = Datos_s.Ka;
Kw              = Datos_s.Kw;
V0              = Datos_s.V0;
Vr              = Datos_s.Vr;
C_HAinit        = C_HA0*V0/Vr;
C_Ainit         = C_A0*V0/Vr;
pKa             = -log10(Ka);
pKb             = -log10(Kw/Ka);

Datos_entrada{...
    strcmp(Datos_entrada(:,1),...
    'C_HAinicial'),2}...
    = C_HAinit;
Datos_entrada{...
    strcmp(Datos_entrada(:,1),...
    'C_Ainicial'),2}...
    = C_Ainit;

% Polinomio C(1)*X^N + ... + C(N)*X + C(N+1) = 0
% Se representa con el vector C = [C(1),C(2),...,C(N)];
Polinomio       = ...
    [...
    +1,...
    +(1-1/Ka*C_Ainit),...
    -(Kw/Ka^2-1/Ka*C_HAinit),...
    -Kw/Ka^2
    ];

raices          = roots(Polinomio);

C_H3O           = raices*Ka;
C_OH            = Kw./C_H3O;
C_HA            = C_HAinit  + C_H3O - C_OH;
C_A             = C_Ainit   - C_H3O + C_OH;
pH              = -log10(C_H3O);
pOH             = -log10(C_OH);

diferencias     = ....
    [(C_HA - C_HAinit)';(C_A - C_Ainit)']';

normas          = zeros(size(diferencias,1),1);

for i=1:size(normas,1)
    normas(i)   = norm(diferencias(i,:));
end

diferenciaMin   = find(normas==min(normas));

if      size(diferenciaMin,1) > 1
    diferenciaMin = diferenciaMin(1);
elseif  size(diferenciaMin,1) < 1
    diferenciaMin = 2;
end

if size(varargin,1) > 0 && ...
    strcmp(varargin{1},'indiceDepHImpuesto')
    diferenciaMin = varargin{2};
end

for i=1:size(Respuesta,1)
    if exist(Respuesta{i,1},'var')
        Respuesta{i,2}   = ...
            eval(Respuesta{i,1});
        if ~isscalar(Respuesta{i,2})
            Respuesta{i,2}   = ...
                eval([Respuesta{i,1},...
                '(',num2str(diferenciaMin),')']);
        elseif ~isscalar(Respuesta{i,2})
            Respuesta{i,2}   = ...
                num2str(eval(Respuesta{i,1}));
        end
    end
end

Botones_Matriz = get(handles.Botones,'Children');

for i=1:3
    set(Botones_Matriz(i),'String',...
        ['pH=',sprintf('%1.2e',pH(i)),...
        ' ||Cambio|| ',sprintf('%1.2e',normas(i)),...
        'gmol/L'] , 'UserData',i);
end

set(Botones_Matriz,'Enable','off');
set(handles.Botones,'SelectedObject',...
    Botones_Matriz(diferenciaMin));

set(Botones_Matriz,'Enable','on');
set(handles.Tabla_Salida,'Data',Respuesta);
set(handles.Tabla_Entrada,'Data',Datos_entrada);

end

function CambioDeBoton(hObject,eventData,handles,varargin)
funcionDeCalcular(hObject,eventData,handles,...
    'indiceDepHImpuesto',get(get(handles.Botones,...
    'SelectedObject'),'UserData'));
end

function GenerarReporte(~,~,~)
publish('bufferSolution.m','format','html','evalCode',true);
web('html/bufferSolution.html');
end

function GenerarCSV(hObject,eventData,handles)
[success,~] = mkdir('DATA');
if success 
    [FileName,PathName,~]=uigetfile('./DATA/*.mat;*.xlsx;*.xls;*.csv');    
else
    [FileName,PathName,~]=uigetfile('./*.mat;*.xlsx;*.xls;*.csv');
end
Datos={};
try
    if FileName~=0
        extension=regexp(FileName,'.mat$|.xls$|.xlsx$|.csv$','match');        
        if strcmp('.mat',extension)
            load([PathName filesep FileName],'Datos');
        elseif strcmp('.csv',extension)        
            [~,Datos]=cargarCSV([PathName FileName]);
        elseif strcmp(extension,'.xls') ||...
                strcmp(extension,'.xlsx')
            [~,~,Datos]=xlsread([PathName FileName]);
            Datos=quitarNaN(Datos);
        end
        %
        % Poner estos valores en la tabla uitable1 (a la izq.)}
        set(handles.Tabla_Entrada,'Data',Datos); %#ok<COLND>
        %
        % Correr el c�digo para actualizar ( o generar en dado caso) la gr�fica
        % solicitada.
        funcionDeCalcular(hObject,eventData,handles,varargin);
    end
catch exception
    msgbox([exception.identifier,'. ',...
        exception.message],'ERROR','error');
end
end