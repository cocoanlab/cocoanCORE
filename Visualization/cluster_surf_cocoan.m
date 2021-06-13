function [p, colorbar, str] = cluster_surf_cocoan(r, varargin)

% Surface plot of clusters on a standard brain. 
%
% This was from canlab core toolbox, where this function is deprecated, but
% still very useful. Thus cocoanlab wanted to keep developing this in our
% github repository. The documentation below will be updated soon. We now
% made this substantially faster than before. 
%
% Some parts of this function requires the canlab core toolbox. 
%
% :Usage:
% ::
%
%    [suface_handle,colorchangestring] = cluster_surf_cocoan(r, P, depth, varargin)
%
% :Inputs: (to be updated)
%
%   **regions:**
%        region object
%
%   **COLORS:**
%       cell array of colors for each cluster: {[1 0 0] [0 1 0] [0 0 1]}
%       if number of colors specified is greater than number of clusters
%       structures entered, n+1 and n+2 colors are overlap of 2 and overlap
%       of all clusters, respectively.
%
%   **SURFACE MAT FILE:**
%       file name of mat file containing brain surface vertices and faces
%       as created with isosurface.
%
%   **SPECIAL SURFACE KEYWORDS:**
%        special string: 'bg' 'hipp' (hcmp,thal,amy)
%        number of mm to plot from surface (mmdeep)
%
%          - Special keywords for sets of surfaces are:
%            left, right, bg, limbic, cerebellum, brainstem
%
%          - Other keywords: 'left' 'right' 'amygdala' 'thalamus' 'hippocampus' '
%            'midbrain'  'caudate'   'globus pallidus'  'putamen'  'nucleus accumbens'
%             'hypothalamus' 'cerebellum'
%
%   **EXISTING SURFACE HANDLE(S):**
%        handles for surface patches, created,
%        e.g., with addbrain.m.  This lets you be very flexible in the
%        surfaces you image onto.
%
%   **'colorscale':**
%        This scales colors by Z-scores of voxels if used
%          - Uses input color, unless also used with 'heatmap'
%          - Z scores should be in ROW vector
%          - use with 'normalize' to scale Z-scores between -1 and 1
%          - will also create transparent effects, mixing
%            blob color with existing surface color in linear
%            proportion to Z-scores
%
%   **'heatmap':**
%        Map Z-scores to surface colors
%          - Used WITH or instead of 'colorscale'
%          - Blobs can have a range of colors
%          - Use with REFERENCE RANGE option below to control scale
%          - solid colors entered as input will be ignored
%          - use with 'colormaps' option below to be flexible
%            in which color maps you apply.
%          - if 'colorscale' is also used, will produce transparent blobs.
%
%   **REFERENCE RANGE:**
%        reference Z-scores range, [zmin_act zmax_act
%        zmax_negact zmin_negact], e.g., [0 5 -5 0], use only
%        with 'heatmap' option
%        to get refZ from clusters, try:
%        ::
%
%            clZ = cat(2,clusters.Z);
%            refZ = [min(clZ(clZ > 0)) max(clZ) min(clZ(clZ < 0)) min(clZ)];
%
%   **'colormaps':**
%        - followed by custom [colors x 3] matrices for positive colors
%          and negative colors.
%        - matlab can create some: e.g., colormap summer, jet, etc.
%          others can be created with colormap_tor.m
%
%        color [0 1 1] (cyan) is reserved for the overlap color btwn cluster sets.
%
% :Examples:
% ::
%
%    P = 'C:\tor_scripts\3DheadUtility\canonical_brains\surf_single_subj_T1_gray.mat';
%    cluster_surf(tcl,acl,P,10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
%
%    or P = h (surface handle) to use current surface in figure, and refZ
%    cluster_surf(tcl,acl,h,[3 5 -5 -3],10,{[0 1 0] [1 0 0]},'colorscale','heatmap')
%
% :More examples:
% ::
%
%    cluster_surf(cl,2,'heatmap');     % brain surface.  vertices colored @2 mm
%    cluster_surf(cl,2,'bg','heatmap');    % heatmap on basal ganglia
%    cluster_surf(cl,5,'left','heatmap');  % heatmap on left surface @5 mm
%    cluster_surf(cl,2,'right','heatmap');
%
%    % A multi-color, multi-threshold display on the cerebellum
%    colors = {[1 1 0] [1 .5 0] [1 .3 .3]};
%    tor_fig;
%    sh = cluster_surf(cl{3},colors(3),5,'cerebellum');
%    cluster_surf(cl{2},colors(2),5,sh);
%    cluster_surf(cl{1},colors(1),5,sh);
%
%    % Custom colormaps:
%    create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap','left');
%
%    poscm = colormap_tor([.2 .2 .4], [1 1 0], [.9 .6 .1]);  %slate to orange to yellow
%    negcm = colormap_tor([0 0 1], [0 .3 1]);  % light blue to dark blue
%    create_figure('Brain Surface'); cluster_surf(cl, 2, 'heatmap', 'colormaps', poscm, negcm, 'left');
%
%    % Single-color transparent map (green):
%    cluster_surf(cl, 2, {[0 1 0]}, 'colorscale', p3(2), 'normalize');
%
% :See Also: addbrain, img2surf.m, surface() methods for objects, cluster_cutaways
%
%    WARNING: SURFACE_CUTAWAY, TOR_3D, AND CLUSTER_SURF ARE DEPRECATED.
%    THEY USE AN OLDER STYLE OF SURFACE RENDERING THAT IS VERY SLOW. THEY
%    STILL WORK, BUT OBJECT-ORIENTED CANLAB TOOLS HAVE SHIFTED TO USING
%    ADDBRAIN.M COMBINED WITH RENDER_ON_SURFACE() METHOD, WHICH USES
%    MATLAB'S ISOCOLORS FOR DRAMATICALLY FASTER COLOR RENDERING.


% ..
%    Programmers' Notes
%    Created by Tor, a long time ago
%    updated Sept 2015 to keep up with matlab graphics and handle some weird
%    stuff with processing inputs.
%      - Figures are now scalars and we need to check for those first.
%      - Also changed default surface and colormap
%
%    Jan 2020: Updated to add 'noverbose' option, minor cosmetic code
%    cleanup
% ..

% -----------------------------------------------------------------------
%    set up input arguments and defaults
% -----------------------------------------------------------------------

mmdeep = 10;
cscale = 0;
heatm = false;
donormalize = false; % used with colorscale
actcolors = [];  % used with heatmap
adjust_var = [];
viewdeg = [135 30];
doverbose = true;
refZ = [];
colorbar = [];
prioritize_last = false;

if any(strcmp(varargin,'noverbose')), doverbose = false; end % do first

P = which('surf_spm2_brain_1mm.mat'); % default

% default color maps
% These match with brain_activations_display:
poscm = colormap_tor([0.96 0.41 0], [1 1 0]);  % warm
negcm = colormap_tor([.23 1 1], [0.11 0.46 1]);  % cools

% -----------------------------------------------------------------------
%    optional inputs
% -----------------------------------------------------------------------

for i = 1:length(varargin)
    if ischar(varargin{i})
        switch varargin{i}
            % functional commands
            case {'depth'}
                mmdeep = varargin{i+1};
            case {'underlay'}
                P = varargin{i+1};
            case {'noverbose'}
                doverbose = false;
            case {'colormaps'}
                if doverbose, disp('Using custom color maps.'); end
                poscm = varargin{i+1}; varargin{i+1} = [];
                negcm = varargin{i+2}; varargin{i+2} = [];
                heatm = 1;
            case {'colors', 'color'}
                mycolors = varargin{i+1};
            case {'colorscale'}
                cscale = 1;
            case {'normalize'}
                donormalize = 1;
            case {'heatmap'}
                heatm = 1;
            case {'refz', 'refZ'}
                refZ = varargin{i+1};
            case {'prioritize_last'}
                prioritize_last = varargin{i+1};
        end
    end
end

% read P
switch P
    case {'left'}
        P = which('surf_spm2_left.mat'); %which('surf_single_subj_grayL.mat');
        viewdeg = [90 0];
            
    case {'right'}
        P = which('surf_spm2_right.mat'); %which('surf_single_subj_grayR.mat');
        viewdeg = [270 0];
            
    case {'hires left'}
        P = which('surf_spm2_brain_left.mat'); %which('surf_single_subj_grayL.mat');
        viewdeg = [90 0];
            
    case {'hires right'}
        P = which('surf_spm2_brain_right.mat'); %which('surf_single_subj_grayR.mat');
        viewdeg = [270 0];
            
    case {'fsavg_right'} % uses freesurfer inflated brain
            % with Thomas Yeo group's RF_ANTs mapping
            % from MNI to Freesurfer. (https://doi.org/10.1002/hbm.24213)
        P = which('surf_freesurf_inflated_Right.mat');
        viewdeg = [270 0];
        adjust_var = 'fsavg_right'; % varargin for getVertexColors
            
    case {'fsavg_left'} % uses freesurfer inflated brain
        P = which('surf_freesurf_inflated_Left.mat');
        viewdeg = [90 0];
        adjust_var = 'fsavg_left'; % varargin for getVertexColors
        
    otherwise
        % do nothing

end  

% -------------------------------------------------------------------------
% * build xyz list
%
% also get cscale values for each coordinate, and alphascale values if both
% heatmap and colorscale options are entered
% -------------------------------------------------------------------------
    
xyz = cat(2,r(:).XYZmm)';
Z = cat(2,r(:).Z)';

% order voxels from lowest to highest, so that peak colors
% appear because they are plotted last
tmp_pos = [xyz(Z>=0,:) Z(Z>=0)];
tmp_neg = [xyz(Z<0,:) Z(Z<0)];

tmp_pos = sortrows(tmp_pos,4);
tmp_neg = sortrows(tmp_neg,4, 'descend');

xyz = [tmp_pos(:,1:3); tmp_neg(:,1:3)];
Z = [tmp_pos(:,4); tmp_neg(:,4)];
    
if cscale
    if donormalize
        Z = Z ./ max(Z);
    end
    
    if heatm
        % treat colorscale as alpha scaling to add transparent blobs
        % (preserve existing surface; good for isosurface objects)
        
        Za = abs(Z);
        a = prctile(Za, 85);  % midpoint of Z{i} defines transparency 0.5
        b = 1./mad(Za); % multiplier for Z for sigmoid
        
        sZ = 1 ./ (1 + exp(-b*(Za-a)));
        
        % fix, if all constant
        sZ(isnan(sZ)) = 1;
        
        alphascale = sZ;
        
        % this is further adusted based on radius
        alphascale = 5 * alphascale ./ mmdeep^3; % should be 3?
    end
    
else
    if heatm
        % if heat map only, set mycolor = [1 1 1]
        mycolors = [1 1 1];
    end
end

% ------------------------------------------------------------
% for heatmap option: get actcolors
% -------------------------------------------------------------
if heatm
    if doverbose, fprintf(' Getting heat-mapped colors\n'); end
    [actcolors, colorbar] = map_data_to_colormap_sub(Z, poscm, negcm, refZ);
end

% -------------------------------------------------------------------------
% * build function call
% -------------------------------------------------------------------------
if doverbose, fprintf(' Building color change function call\n'); end

str = ['[c,alld] = getVertexColors_cocoan(xyz, p, mycolors, [.5 .5 .5], ' num2str(mmdeep)];

if heatm
    str = [str ',''colorscale'',actcolors'];
    if cscale
        % treat colorscale as alpha scaling to add transparent blobs
        str = [str ',''alphascale'',alphascale'];
    end
elseif cscale
    % cscale alone - treat cscale as color-mapping index for single
    % colors in actcolors
    str = [str ',''colorscale'',Z'];
end

if ~isempty(adjust_var)
    str = [str ',''' adjust_var ''''];
end

if ~doverbose
    str = [str ', ''noverbose'''];
end

if prioritize_last
    str = [str ', ''prioritize_last'''];
end

str = [str ');'];

if exist('alphascale','var') % don't know what this is doing, but keep it from canlab version
    alphascale = alphascale * 12;
end

% -------------------------------------------------------------------------
% * run brain surface
% -------------------------------------------------------------------------
if ishandle(P)      % no input file, use existing handle
    if doverbose
        fprintf(' Using existing surface image\n');
        fprintf(' Running color change.\n');
    end
    for i = 1:length(P)
        p = P(i);
        if doverbose, disp([' eval: ' str]), end
        eval(str)
    end
else
    % we have either an input file or a special string ('bg')
    if doverbose, fprintf(' Loading surface image\n'); end
    [~, ~, etmp] = fileparts(P);
    
    if strcmp(etmp,'.mat')
        
        load(P);
        
        %%figure
        p = patch('Faces',faces,'Vertices',vertices,'FaceColor',[.5 .5 .5], ...
            'EdgeColor','none','SpecularStrength',.2,'FaceAlpha',1,'SpecularExponent',200);
        lighting gouraud;camlight right
        axis image;
        lightRestoreSingle(gca);
        %myLight = camlight(0,0);set(myLight,'Tag','myLight');
        %set(gcf, 'WindowButtonUpFcn', 'lightFollowView');lightFollowView
        
        view(viewdeg(1),viewdeg(2));
        drawnow
        
        
        % -------------------------------------------------------------------------
        % * run color change
        % -------------------------------------------------------------------------
        if doverbose
            fprintf(' Running color change.\n');
            disp([' eval: ' str])
        end
        eval(str);
        
        
        % this for subcortex stuff
%     elseif strcmp(P,'bg') % we can make some default structures here..
%                           % but comment out for now becuase they are
%                           % really old
%         p = [];
%         myp = addbrain('caudate');p = [p myp];
%         run_colorchange(myp,str,xyz,mycolors);
%         
%         myp = addbrain('globus pallidus');p = [p myp];
%         run_colorchange(myp,str,xyz,mycolors);
%         
%         myp = addbrain('putamen');p = [p myp];
%         run_colorchange(myp,str,xyz,mycolors);
%         
%         set(myp,'FaceAlpha',1);
%         
%         axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
        
%     elseif strcmp(P,'limbic') % we can make some default structures here..
%         p = [];
%         
%         myp = addbrain('amygdala');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         myp = addbrain('hypothalamus');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         myp = addbrain('hippocampus');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         myp = addbrain('thalamus');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         myp = addbrain('nucleus accumbens');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         
%         myp = addbrain('left');p = [p myp];
%         run_colorchange(myp,str,xyz, mycolors, actcolors);
%         set(myp,'FaceAlpha',1);
%         
%         axis image; axis vis3d; lighting gouraud; lightRestoreSingle(gca)
        
    else
        error('Must input mat surf file or img file to convert to surf')
    end
    
end % if ishandle

lighting gouraud
lightRestoreSingle(gca);
material dull
axis off
set(gcf,'Color','w')
%scn_export_papersetup(400);  % this will mess up movies!!

if doverbose
    disp('Finished!')
    disp('___________________________________________')
end

end % main function


function [actcolor, colorbar] = map_data_to_colormap_sub(datavalues, poscm, negcm, varargin)
% Usage
% ::
%
%    actcolor = map_data_to_colormap(datavaluesets, poscm, negcm, varargin)
%
% Given sets of data values (each cell is a row vector of data values,
% e.g., z-scores) and color maps for positive and negative values,
% returns mapped colors for each data value in order.
% These colors can be used for direct plotting.
%
% :Inputs:
%
%   input 1:
%        data values (e.g., z-scores). k data value sets, in cells.  Each cell contains row vector of data values
%
%   input 2/3:
%        color maps [n x 3] for positive and negative values
%
%   input 4:
%        optional: fixed range of data defining max and min colors
%
% :Examples:
% ::
%
%    poscm = colormap_tor([0 0 0], [1 1 0]);
%    negcm = colormap_tor([0 0 1], [0 0 0]);
%    Z = randn(40, 1)';
%    actcolors = map_data_to_colormap({Z}, poscm, negcm)
%    [Z' actcolors{1}]
%
% ..
% Tor Wager, Sept. 2007
% ..

    nposcolors = size(poscm, 1);
    nnegcolors = size(negcm, 1);

    % -------------------------------------------------------------
    % determine overall data range
    % -------------------------------------------------------------

    % exactly zero values will give wrong length, so make pos. small number
    datavalues(datavalues == 0) = 100 * eps;
    
    posvalues = datavalues(datavalues > 0);
    negvalues = datavalues(datavalues < 0);
   
    if ~isempty(posvalues)

        zrange = [min(posvalues) max(posvalues)];

        % input fixed data range for max and min colors
        if ~isempty(varargin) && ~isempty(varargin{1})
            zrange = varargin{1}(1:2);

        end

        zh =  linspace(zrange(1),zrange(2),nposcolors);

        if isempty(zh), zh = [1 1 0]; end   % only one element?
        
        colorbar.pos = [zh' poscm]; 
    end

    if ~isempty(negvalues)
        zrangec = [min(negvalues) max(negvalues)];

        % input fixed data range for max and min colors
        if ~isempty(varargin) && ~isempty(varargin{1})
            zrangec = varargin{1}(3:4);
        end

        zhc =  linspace(zrangec(1),zrangec(2),nnegcolors);

        if isempty(zhc), zhc = [0 0 1]; end   % only one element?
        
        colorbar.neg = [zhc' negcm]; 
    end
    

    % -------------------------------------------------------------
    % find color for each xyz
    % -------------------------------------------------------------
    
    if sum(isnan(datavalues))
        disp('Warning! NaNs in data values mapped to colors.  These will be mapped to black [0 0 0].');
    end
    
    actcolor = zeros(length(datavalues), 3);

    for i = 1:length(datavalues)

    	dv = datavalues(i);

        if dv < 0
            [mydistance, wh] = min(abs(zhc - dv), [], 2);
            actcolor(i,:) = negcm(wh, :);
            
        elseif dv >= 0
            [mydistance, wh] = min(abs(zh - dv), [], 2);
            actcolor(i,:) = poscm(wh, :);
            
        else
            % could be NaN; leave as 0
            
        end
    end

end




