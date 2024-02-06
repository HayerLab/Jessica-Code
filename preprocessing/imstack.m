function hmf = imstack(IMGS, ud, MASK)
%im_stack_viewer Displays an image stack in a imtool-like window
% Allows for two-finger scroll to scrub through the images
% Right click to bring up image info and controls
% Requirements:
% ---- IMGS and MASK must be 3D or 4D matrices
% ---- Requires imopener, addRoiToolbar, orthoSlices, subtightplot

% USAGE

% imstack - brings up dialog for opening an image stack

% imstack(3D or 4D image_stack) - opens image stack. In a 4D image stack, the 3rd dimension
% holds channel intensity, and the 4th dimension is the section info

% imstack(image_stack, ud) - a ud structure must contain specific info
% about the image_stack (see init_ud)

% imstack(image_stack,ud,mask) - mask must match the same dimenions as
% image_stack

%INTERNAL INFO
% appdata:
% ---- DS - dataset created when using the count tool
% ---- IMGS - image stack
% ---- MASK - MASK stack, same size as IMGS
% ---- tool_handles - handle structure containing info for tools: count, threshold, group

% UserData:
% ---- UserData for the main figure is a structure that contains info for
% the image stack and mask.
% ---- Userdata for the children of the main figure may sometimes contain the
% handle to the main figure

% Script Folder:
% - Contains a thresh and count folder that contains scripts used  by the threshold and count tools
% - Scripts can be added to customize thresholding or pixel counting
% for a given project
% - Scripts added to the these folders will automatically be detected
% by the respective tool

% TODO: add 3 simultaneous orthagonal views for scrolling
% TODO: show isosurface
% TODO: UPDATE ADDROITOOLS (buggy)
% TODO: Update Bioformats Load - eliminate subfunctions?
% TODO: Change script when loading dataset
% TODO: Reimplement DICOM Reader
% TODO: Window Clean up on close


% Written by Ernesto Salcedo, PhD
% Senior Instructor, University of Colorado School of Medicine
% 2016

% Questions: Ernesto Salcedo @ ucdenver dot edu

%% -------- CREATE FIGURE / TOOL HANDLES -----------------------------------------------

hmf = figure(...
    'Tag','imstack_fig',...
    'Colormap', bone(255),...
    'NumberTitle', 'off',...
    'DeleteFcn',@delete_imstack,...
    'WindowScrollWheelFcn',@scroll_process,...
    'WindowKeyPressFcn', @key_press);

axes;
toolbar_cleanup(hmf);

% Add toolbar items
addRoiToolbar(hmf);

imstack_contents = what('imstack'); % find all MATLAB related files in imstack

% Check to make sure bioformats is up to date
% bfUpgradeCheck(true)

%% -------- CREATE MENUS -----------------------------------------------
%% FILE MENU
% Setup File menu: http://undocumentedmatlab.com/blog/customizing-standard-figure-toolbar-menubar/
hMenuFile = findall(gcf,'tag','figMenuFile');
delete(allchild(hMenuFile)); % clear all built-in menu items

% import_formats = {'bioformats', 'tiff stack', 'dicom stack', 'dicom virtual', 'image folder','workspace variable'};

import_formats = {'tiff stack', 'image folder','workspace variable', 'using custom script'};

if exist('bfmatlab','dir')
    import_formats = [{'bioformats'} import_formats];
else
    display ('bioformats toolbox not detected. Make certain bfmatlab is added to to the path')
end

for m = 1:numel(import_formats)
    uimenu('Label',sprintf('Open %s', import_formats{m}),...
        'Parent',hMenuFile,...
        'UserData',import_formats{m},...
        'Tag', 'ims_menu_open',...
        'Callback',@open_menu_callback);
end

uimenu(hMenuFile, 'Label','Revert Image Stack',...
    'Tag','revert_image',...
    'Separator','on',...
    'UserData',hmf,...
    'Accelerator','R',...
    'Callback', @menu_revert_STACK)

uimenu(hMenuFile, 'Label','Revert Mask Stack',...
    'Tag','revert_mask',...
    'Separator','off',...
    'UserData',hmf,...
    'Callback', @menu_revert_STACK)

uimenu(hMenuFile,'Label','Load dataset',...
    'Separator','on',...
    'UserData',hmf,...
    'Accelerator','l',...
    'Callback', @menu_load_save_dataset_callback );

uimenu(hMenuFile,'Label','Save Dataset',...
    'Accelerator','s',...
    'UserData',hmf,'Callback', @menu_load_save_dataset_callback )

uimenu(hMenuFile,'Label','Save Dataset as...',...
    'UserData',hmf,'Callback', @menu_load_save_dataset_callback )

uimenu(hMenuFile,'Label','Export Dataset as...',...
    'UserData',hmf,'Callback', @menu_load_save_dataset_callback )

uimenu(hMenuFile,'Label','ADD Dataset to WORKSPACE',...
    'UserData',hmf,'Callback', @menu_load_save_dataset_callback )

uimenu(hMenuFile,'Label', 'Add stack (and mask) to workspace', ...
    'Separator','on',...
    'Tag', 'fm_stack2workspace',...
    'Enable', 'off',...
    'Callback', @assign2workspace)
uimenu(hMenuFile,'Label', 'Add file paths to workspace', ...
    'Separator','off',...
    'Tag', 'fm_paths2workspace',...
    'Enable', 'off',...
    'Callback', @assign2workspace)
uimenu(hMenuFile,'Label', 'Add current image slice to workspace', ...
    'Tag','fm_slice2workspace',...
    'Callback', @assign2workspace);
uimenu(hMenuFile,'Label', 'Add imstack user data to workspace', ...
    'Tag','fm_ud2workspace',...
    'Callback', @assign2workspace);

uimenu(hMenuFile,'Label', 'Save MASK',...
    'Separator','off',...
    'Callback', @save_mask)

%% EDIT MENU
hMenuEdit = findall(gcf,'tag','figMenuEdit');
delete(allchild(hMenuEdit)); % clear all built-in menu items

uimenu(hMenuEdit,'Label', 'Crop Stack', 'Callback', @menu_crop_callback);
uimenu(hMenuEdit,'Label', 'Trim Stack', 'Callback', @menu_trim_stack_callback);
uimenu(hMenuEdit,'Label', 'Clear mask',...
    'Separator','on',...
    'Tag','hMenuEditMaskClear',...
    'Callback', @clear_MASK);

% uimenu(hMenuEdit,'Label', 'Color Map tool', ...
%         'Separator','on',...
%         'Callback', 'imcolormaptool');

%% REMOVE INSERT MENU
delete(findall(gcf,'tag','figMenuInsert'));

%% TOOL MENU

hMenuTools = findall(gcf,'tag','figMenuTools');
delete(allchild(hMenuTools)); % clear all built-in menu items

% initialize tool handles
ht = [];

% load function handles of scripts
% REQUIRES files in DEFAULT script folder to be name 'count_default', 'thresh_default', etc.
% cd(fullfile(imstack_contents.path,'scripts','DEFAULT'))

% cd(fullfile(imstack_contents.path,'scripts'))
load(fullfile(imstack_contents.path,'scripts','last_selected.mat'))
ht.tool_path = fullfile(imstack_contents.path,'scripts',script_folder);
ht.tool_names = {'count', 'thresh', 'segment','stats'};
ht = get_tool_scripts(ht);

setappdata(hmf,'tool_handles', ht);

% Find Folder Names in Scripts folder (ignore invisible folders)
scriptsf_content = dir(fullfile(imstack_contents.path,'scripts'));
scriptsf_content(~[scriptsf_content.isdir]) = [];
scriptsf_content(~cellfun(@isempty,regexp({scriptsf_content.name},'^\.')))=[];
folder_names = {scriptsf_content.name}';

uimenu(hMenuTools, 'Label','Threshold...', ...
    'Separator', 'off',...
    'Accelerator','T',...
    'UserData', hmf,...
    'Callback', @ims_threshold_tool);

uimenu(hMenuTools, 'Label','Segmentation Tool...', ...
    'Separator', 'off',...
    'Accelerator','G',...
    'UserData', hmf,...
    'Callback', @ims_segmentation_tool);

uimenu(hMenuTools, 'Label','Count Tool...', ...
    'Separator', 'off',...
    'Accelerator','Y',...
    'UserData', hmf,...
    'Callback', @ims_count_tool);

uimenu(hMenuTools, 'Label','Filter Image Stack...', ...
    'Tag', 'filter_image',...
    'Separator', 'on',...
    'Accelerator','F',...
    'UserData', hmf,...
    'Callback', @ims_filter_tool);

uimenu(hMenuTools, 'Label','Filter Mask Stack...', ...
    'Tag', 'filter_mask',...
    'Separator', 'off',...
    'UserData', hmf,...
    'Callback', @ims_filter_tool);

uimenu(hMenuTools, 'Label','Dataset Stats...', ...
    'Separator', 'on',...
    'UserData', hmf,...
    'Callback', @ims_stats_tool);

% Create  Set Scripts submenu
hMenuToolSet = uimenu(hMenuTools, 'Label','Set Tools Scripts', ...
    'Separator', 'on');
for fn = 1:numel(folder_names)
    hm = uimenu(hMenuToolSet, 'Label',folder_names{fn}, ...
        'Separator', 'off',...
        'Callback', @menu_tool_set_scripts);
    
    if strcmp(folder_names{fn}, script_folder)
        hm.Checked = 'on';
    end
end

uimenu(hMenuTools,'Label','cd to scripts folder',...
    'Separator','off',...
    'Callback', @(varargin) cd(fullfile(imstack_contents.path,'scripts')));


%% VIEW MENU
hViewMenu = findall(gcf,'tag','figMenuView');

uimenu(hViewMenu, 'Label', 'Show Montage',...
    'Separator', 'on',...
    'Tag', 'vm_montage',...
    'Callback', @call_montage);

uimenu(hViewMenu, 'Label', 'Show Max Projection',...
    'Tag', 'max_projection',...
    'Callback', @show_projection);
% uimenu(hViewMenu, 'Label', 'Projection - Split Channels',...
%     'Tag', 'vm_projection',...
%     'Callback', @show_projection);
uimenu(hViewMenu, 'Label', 'Isosurface - Mask',...
    'Separator', 'on',...
    'Tag', 'show_iso',...
    'Callback', @show_isosurface);

%% HELP MENU
hMenuHelp = findall(gcf,'tag','figMenuHelp');
url = 'https://dl.dropboxusercontent.com/u/10600990/site/MATLAB%20functions/html/imstack_doc.html';
uimenu(hMenuHelp, 'Label','IMSTACK User Guide', ...
    'Separator', 'on',...
    'Callback', @(varargin) web(url,'-browser'))

% -------- END MENUS -----------------------------------------------

%% PARSE INPUTS
% if no arguments passed in, go find an image stack and display
if nargin == 0
    
    [Selection,ok] = listdlg('ListString',import_formats,...
        'PromptString', 'SELECT SOURCE',...
        'SelectionMOde', 'single',...
        'Name','Open Image Stack');
    if ok
        ok = load_IMGS(hmf, import_formats{Selection});
    end
    
    if ~ok
        display('load canceled') 
        close(hmf)
        close(findall(0,'tag','bformat_load'))
        return
    end
else % if arguments passed in
    
    % If no structure passed in, get name of variable
    if nargin == 1 || isempty(ud)
        
        % handles 4D stacks that only have one channel
        sz = size(IMGS);
        if ndims(IMGS)>3 && sz(3) == 1
            IMGS = squeeze(IMGS);
        end
        
        % get name of input variable
        if ~isempty(inputname(1))
            ud.StudyDescription = inputname(1);
            ud.SeriesDescription = '';
        else
            ud.StudyDescription = 'IMGS';
            ud.SeriesDescription = '';
        end
    end
    
    setappdata(hmf,'IMGS', IMGS); %attach image stack to figure and display
    
    ud = init_ud(ud, hmf, size(IMGS),false); % also clears mask app data
    
    % if MASK passed in, attach to figure
    if nargin == 3
        setappdata(ud.hmf,'MASK', MASK)
    end
    % initialize settings
    init_im_display(ud); % menu items, context menu
end

% handles = guihandles(hmf);
% guidata(hmf,handles)

end
% -------- END imstack MAIN FUNCTION-------------------------------------------
%% Subfunctions------------------------------------------------------------
%% ------------------------------------------------------------------------

function open_menu_callback(obj,~)
% callback for menu items that start with 'Open' in the File menu

menu_items = findall(get(obj,'Parent'),'Tag','ims_menu_open');
set(menu_items,'Accelerator','')
set(obj,'Accelerator','o')

load_IMGS(imgcf, get(obj,'UserData'))

end

%% ---- LOAD STACK ----------------------------------------------------------
function  ok = load_IMGS(hmf, selection)
% find an image stack to load
% uses imopener - keep at this level
ok = true;
IMGS = [];

try
    switch selection
        case 'bioformats'
            bioformats_load(hmf);
            return
        case 'tiff stack'
            [IMGS, ud] = read_tiff_stack;
            
        case 'image folder'
            [IMGS, ud] = read_image_folder;
            
        case 'dicom stack'
            [ DCM_FULL_PATHS ] = dicom_get_paths; 
            [IMGS, ud] = dicom_read_engine(DCM_FULL_PATHS);
            
        case 'dicom virtual'
            [ DCM_FULL_PATHS ] = dicom_get_paths; 
            if isappdata(hmf, 'IMGS')
                rmappdata(hmf, 'IMGS');
            end
            [~,ud] = imopener(DCM_FULL_PATHS{1});
            setappdata(hmf,'PATHS', DCM_FULL_PATHS); % attach image stack to figure
            ud = init_ud(ud, hmf, numel(DCM_FULL_PATHS), true);
            
        case 'workspace variable'
            % load workspace variable - allow for structure load
            % warning allows for more than two variable selections
            ws_vars = evalin('base','whos'); % retrieve varibles from workspace
            names = {ws_vars.name};
            [Selection,ok] =  listdlg('ListString', names, ...
                'SelectionMode', 'multiple',...
                'Name', 'Workspace',...
                'PromptString', {'Select variable'});
            
            ws_vars = ws_vars(Selection); % reduces to only selected variables
            
            if numel(ws_vars) > 1 % more than one variable selected
                Selection = ~strcmp({ws_vars.class},'struct'); % change to logical array pointing to numeric array
                ud = evalin('base',ws_vars(~Selection).name); % not selection = structure
            else
                ud.StudyDescription = ws_vars(Selection).name;
                ud.SeriesDescription = '';
            end
            IMGS=evalin('base',ws_vars(Selection).name);
        case 'using custom script'
            [IMGS,ud] = read_script(hmf);
    end
    
    if ~isempty(IMGS) % successful open
        setappdata(hmf,'IMGS', IMGS); % attach image stack to figure
        ud = init_ud (ud, hmf, size(IMGS), false);
        %     ud = generate_adj_in(ud);
    end
    
    init_im_display(ud);
catch ME
    beep
    ok = false;
    fprintf('\nError Loading image stack: %s\n',ME.message)
    display(ME.stack)
    assignin('base','ME',ME)
    return
end
end

%% ------------------------------------------------------------------------
function call_montage(~,~)
% Needs to handle more than three channels
ud = get(imgcf,'UserData');
hfmontage = figure('Name',sprintf('MONTAGE: %s',ud.StudyDescription),'NumberTitle','off');
toolbar_cleanup(hfmontage);

IMGS = getappdata(ud.hmf,'IMGS');

if ud.ChannelCount < 3
    IMGS(:,:,3,:) = zeros(ud.StackSize(1), ud.StackSize(2),1,ud.StackSize(4));
end

montage(IMGS)
impixelinfo
end

%% ------------------------------------------------------------------------
function delete_imstack(~,~)
% figure deletefcn

ud = get(gcbo,'UserData');

if isempty(ud)
    return
end

ht = getappdata(gcbo,'tool_handles');

if isfield(ht,'tool_h')
    tags = fieldnames(ht.tool_h);
    % check for tool windows (thresh, count, etc) and delete
    for n=1:numel(tags)
        if ishandle(ht.tool_h.(tags{n}))
            delete(ht.tool_h.(tags{n}))
        end
    end
end

if isappdata(ud.hmf,'DS')
    save_data(ud.hmf,ht,'Save Dataset')
end
end

%% ------------------------------------------------------------------------
function menu_revert_STACK(obj,~)
% reverts stack to save original version (name preceded by a 'o')
% that was generated after using the filter tool

switch get(obj,'Tag')
    case 'revert_image'
        stack_name = 'IMGS';
    case 'revert_mask'
        stack_name = 'MASK';
end

hmf = get(obj,'UserData');
if isappdata(hmf,['o' stack_name])
    setappdata(hmf,stack_name, getappdata(hmf,['o' stack_name]));
else
    fprintf('original %s stack not found',stack_name);
    %     display('original stack not found')
end
end

%% ------------------------------------------------------------------------
function menu_load_save_dataset_callback(obj,~)
hmf = get(obj,'UserData');
ht = getappdata(hmf,'tool_handles');

label = get(obj,'Label');
switch label
    case 'Load dataset'
        % Load
        data_path = get_datapath(hmf,ht,label);
        
        if isempty(data_path)
            display ('load canceled')
            return
        end
        
        load(data_path); % should load a table variable named DS
        
        if istable(DS)
            setappdata(hmf,'DS', DS);
            if isfield(ht,'tool_h')
                structfun(@close,ht.tool_h); % close all tool windows
            end
            
            tool_path = ht.tool_path; % tool_path set at imstack start
            
            ht = DS.Properties.UserData; % replace with ht saved data
            
            ht.tool_path = tool_path;
            ht.data_path = data_path;
            ht.count.selected_row = size(DS,1); % set selection to end
            set_dataset_table(ht, DS);            
            group_table_settings(ht)
            
            ht = get_tool_scripts( ht);
            
            display('Dataset loaded');
        else
            display('wrong *.mat file. File must return a table variable')
        end
    case {'Save Dataset','Save Dataset as...'}
        ht = save_data(hmf,ht,label);
    case 'Export Dataset as...'
        DS = getappdata(hmf,'DS');
        [filename,pathname] = uiputfile(fullfile(pwd,'*.csv'));
        writetable(DS,fullfile(pathname,filename));
    case 'Add Dataset to workspace'
        DS = getappdata(hmf,'DS');
        assignin('base', 'DS', DS);
    otherwise
        display('Sorry, I got nothing')
        return
end
setappdata(hmf,'tool_handles',ht);
end

