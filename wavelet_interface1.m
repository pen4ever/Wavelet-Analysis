

function varargout = wavelet_interface1(varargin)



set(0,'DefaultAxesFontSize', 8)
% WAVELET_INTERFACE1 MATLAB code for wavelet_interface1.fig
%      WAVELET_INTERFACE1, by itself, creates a new WAVELET_INTERFACE1 or raises the existing
%      singleton*.
%
%      H = WAVELET_INTERFACE1 returns the handle to a new WAVELET_INTERFACE1 or the handle to
%      the existing singleton*.
%
%      WAVELET_INTERFACE1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WAVELET_INTERFACE1.M with the given input arguments.
%
%      WAVELET_INTERFACE1('Property','Value',...) creates a new WAVELET_INTERFACE1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before wavelet_interface1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to wavelet_interface1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help wavelet_interface1

% Last Modified by GUIDE v2.5 05-Nov-2015 17:31:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @wavelet_interface1_OpeningFcn, ...
                   'gui_OutputFcn',  @wavelet_interface1_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% 
% function saveState(handles)
% state.loadButton = get(handles.pushbutton1, 'value');
% save('state.mat', 'state');

% --- Executes just before wavelet_interface1 is made visible.
function wavelet_interface1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to wavelet_interface1 (see VARARGIN)

% Choose default command line output for wavelet_interface1


handles.output1 = hObject;
xlabel(handles.Spectrum, 'Time','FontSize',14)
ylabel(handles.Spectrum, 'Frequency [Hz]','FontSize',14)
set(handles.Spectrum,'FontSize',14);
guidata(hObject, handles);

handles.output3 = hObject;
xlabel(handles.Original, 'Frame','FontSize',14)
ylabel(handles.Original, 'Data','FontSize',14)
set(handles.Original,'FontSize',14);
guidata(hObject, handles);





% Update handles structure


% UIWAIT makes wavelet_interface1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = wavelet_interface1_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output1;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% handles.file_1 is the name of the file 
[handles.file_1 handles.folderpath] = uigetfile(); %if using GUIDE will save it inside the handles structure.  
format = handles.file_1(length(handles.file_1)-3:length(handles.file_1));
clear handles.data
if  strcmp(format, '.txt') == 1
    data = load(handles.file_1);
    handles.data = data;
elseif strcmp(format, '.mat') == 1
    data = load(handles.file_1);
    handles.Fname = fieldnames(data)
    handles.data = data.(handles.Fname{1});
    size(handles.data)
end
guidata(hObject, handles)
plot(handles.Original,0:1:length(handles.data(:))-1, handles.data(:))
xlabel(handles.Original, 'Frame','FontSize',14)
ylabel(handles.Original, 'Data','FontSize',14)

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in WaveForm.
function WaveForm_Callback(hObject, eventdata, handles)
% hObject    handle to WaveForm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
contents = cellstr(get(hObject, 'String'));
wave = contents{get(hObject,'Value')};

if strcmpi(wave,'Morlet') == 1
    clearvars handles.wavename handles.powerr  handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'morl';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.8125./(handles.scales.*1/handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
   
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.mort;
  
    val_WAV = val_WAVs.valWAV.mort;
    xWAV = xval;
    xMaxWAV = 16;
    stepWAV = xWAV(2)-xWAV(1);

    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    stepSIG = 1;
    coefs = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
            f            = fliplr(val_WAV(j));
            coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
            ind          = ind+1;
    end
%     coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
    guidata(hObject, handles)
    
elseif strcmpi(wave,'Haar') == 1 
    clearvars handles.wavename handles.powerr  handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'haar';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.9961./(handles.scales.*1/handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.haar;
    val_WAV = val_WAVs.valWAV.haar;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
    guidata(hObject, handles)
      
  elseif strcmpi(wave,'Gaus1') == 1 
    clearvars handles.wavename handles.powerr  handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'gaus1';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.2./(handles.scales.*1/handles.Fs);
   % frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.gaus1;
    val_WAV = val_WAVs.valWAV.gaus1;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)    
      
   elseif strcmpi(wave,'Gaus2') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'gaus2';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.3./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
     xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.gaus2;
    val_WAV = val_WAVs.valWAV.gaus2;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)  
  
elseif strcmpi(wave,'Gaus3') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'gaus3';
    handles.wavename = wavename;
    Fs = 10000;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.4./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
     xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.gaus3;
    val_WAV = val_WAVs.valWAV.gaus3;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)  
      
 elseif strcmpi(wave,'Gaus4') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'gaus4';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.4./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.gaus4;
    val_WAV = val_WAVs.valWAV.gaus4;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
%    coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles) 
      
  elseif strcmpi(wave,'Gaus8') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'gaus8';
    handles.wavename = wavename;
    Fs = 10000;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.6./(handles.scales.*1/handles.Fs);
    
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.gaus8;
    val_WAV = val_WAVs.valWAV.gaus8;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end

    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles) 
      
      
      
  elseif strcmpi(wave,'Db1') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'db1';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.9961./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.db1;
    val_WAV = val_WAVs.valWAV.db1;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
    guidata(hObject, handles) 
      
      
  elseif strcmpi(wave,'Db2') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'db2';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.6667./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
     xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');
  
    xval = xvals.xval.db2;
    val_WAV = val_WAVs.valWAV.db2;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)  
      
      
   elseif strcmpi(wave,'Db3') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'db3';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.8./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    xval = xvals.xval.db3;
    val_WAV = val_WAVs.valWAV.db3;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end    
    
    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)     
      
      
 elseif strcmpi(wave,'Db4') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'db4';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.7143./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    xval = xvals.xval.db4;
    val_WAV = val_WAVs.valWAV.db4;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end     
    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles) 
      
      
  elseif strcmpi(wave,'Db8') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'db8';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.6667./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    xval = xvals.xval.db8;
    val_WAV = val_WAVs.valWAV.db8;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end       
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)       
      
      
      
    elseif strcmpi(wave,'Coif1') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'coif1';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.8./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    xval = xvals.xval.coif1;
    val_WAV = val_WAVs.valWAV.coif1;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end       

    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)
      
   elseif strcmpi(wave,'Coif2') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data;
    wavename = 'coif2';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.7273./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    xval = xvals.xval.coif2;
    val_WAV = val_WAVs.valWAV.coif2;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end       

    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)    
 
   elseif strcmpi(wave,'Coif3') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data; 
    wavename = 'coif3';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.7059./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');    
    
    xval = xvals.xval.coif3;
    val_WAV = val_WAVs.valWAV.coif3;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end       
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles) 
      
      
    elseif strcmpi(wave,'Coif4') == 1 
    clearvars handles.wavename handles.powerr handles.Spectrum
    cla(handles.Spectrum);
    cla(handles.Freq);
    data1 = handles.data; 
    wavename = 'coif4';
    handles.wavename = wavename;
    t1 = 0:1:length(handles.data)-1;
    t1 = t1.*(1/handles.Fs);
    frq = 0.6957./(handles.scales.*1/handles.Fs);
    %frq = scal2frq(handles.scales, handles.wavename, handles.Fs);
    handles.frq = frq;
    set(handles.MinFrq, 'String', num2str(min(handles.frq)));
    set(handles.MaxFrq, 'String', num2str(max(handles.frq)));
    plot(handles.Freq, handles.scales, handles.frq,'*');
    xlabel(handles.Freq, 'Scales','FontSize',14)
    ylabel(handles.Freq, 'Frequency [Hz]','FontSize',14)
    set(handles.Freq,'FontSize',14);
    xvals = load('xval.mat');
    val_WAVs = load('valWAV.mat');  
    xval = xvals.xval.coif4;
    val_WAV = val_WAVs.valWAV.coif4;
    xWAV = xval;
    stepWAV = xWAV(2)-xWAV(1);
    stepSIG = 1;
    xWAV = xWAV-xWAV(1);
    xMaxWAV = xWAV(end);
    ySIG    = handles.data;
    lenSIG = length(ySIG);
    nb_SCALES = length(handles.scales);
    coefs     = zeros(nb_SCALES,lenSIG);
    ind  = 1;
    for k = 1:nb_SCALES
        a = handles.scales(k);
        a_SIG = a/stepSIG;
        j = 1+floor((0:a_SIG*xMaxWAV)/(a_SIG*stepWAV));     
        if length(j)==1 , j = [1 1]; end
        f            = fliplr(val_WAV(j));
        coefs(ind,:) = -sqrt(a)*wkeep1(diff(wconv1(ySIG,f)),lenSIG);
        ind          = ind+1;
    end       
    
    %coefs = cwt(handles.data,handles.scales,handles.wavename);
    powerr = abs(coefs).^2;
    guidata(hObject, handles)
    [ff tt] = ndgrid(frq, t1);
     handles.powerr = powerr;
     handles.ff = ff;
     handles.tt = tt;
  
    axes(handles.Spectrum);
    ylim([min(frq), max(frq)]);
    xlim([min(t1), max(t1)]);
    hold on
    contourf(tt,ff,powerr)
    contour(tt,ff,powerr)
    hold off
      guidata(hObject, handles)     
end

% data1 = handles.data; 
% plot(handles.Fitted, data1(:,1), data1(:,2))


%coefs = cwt(data,scalesset,wave)
% Hints: contents = cellstr(get(hObject,'String')) returns WaveForm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from WaveForm


% --- Executes during object creation, after setting all properties.
function WaveForm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to WaveForm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Variable_Callback(hObject, eventdata, handles)
% hObject    handle to Variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Variable as text
%        str2double(get(hObject,'String')) returns contents of Variable as a double
handles.varname = get(hObject,'String');
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Variable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Variable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ScalesSet_Callback(hObject, eventdata, handles)
% hObject    handle to ScalesSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ScalesSet as text
%        str2double(get(hObject,'String')) returns contents of ScalesSet as a double
s= get(hObject, 'String');
[row, column] = size(s)
for i = 1:row                    
      eval(s(i,:));                  
end
handles.scales = scales;
guidata(hObject, handles)


% --- Executes during object creation, after setting all properties.
function ScalesSet_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScalesSet (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function SaveVarName_Callback(hObject, eventdata, handles)
% hObject    handle to SaveVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of SaveVarName as text
%        str2double(get(hObject,'String')) returns contents of SaveVarName as a double
% user_entry = str2double(get(hObject,'string'));
user_entry = get(hObject,'string');
handles.user_entry = user_entry;
guidata(hObject,handles)



% --- Executes during object creation, after setting all properties.
function SaveVarName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SaveVarName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in SaveVar.
function SaveVar_Callback(hObject, eventdata, handles)
% hObject    handle to SaveVar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

powerr = handles.powerr;
time = handles.tt;
frq = handles.ff;

filename = handles.user_entry;
uisave({'powerr', 'time', 'frq'},filename) 

% h = 365;
% g = 52;
% uisave({'h','g'},'var1');



function Fs_Callback(hObject, eventdata, handles)
% hObject    handle to Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Fs as text
%        str2double(get(hObject,'String')) returns contents of Fs as a double
user_entry_fs = get(hObject,'string');
handles.Fs = eval(user_entry_fs(:));
guidata(hObject,handles)

% --- Executes during object creation, after setting all properties.
function Fs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Fs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
