%% Name Plate_local_coordination_date_angle.m
% This script measures the local coordination of neighboring migrating
% cells by calculating a velocity correlation function
% corr(r)=dot(u,v)/(norm(u)*norm(v)), which is the cos(theta), where theta
% is the angle between velocity vectors of neighboring cells. 
% 160715 Arnold Hayer

% Input: 
%     row_col_site_tracedata.mat - tracking data from Min's code
    
% Output: 
%     The input datastructure tracedata.mat, but extended by additional layers in the 3rd dimension for
%     motility parameters, for each cell, and each timepoint: 
%     tracedata(:,:,9)  - angle of movement
%     tracedata(:,:,10) - single cell velocity
%     tracedata(:,:,11) - coordination with neighboring cells in the front, located within a 60° sector ahead
%     tracedata(:,:,12) - coordination with lateral neighbors, located within lateral 120° sectors)
%     tracedata(:,:,13) - coordination with neighboring cells in the back, located within a 60° sector behind
%     tracedata(:,:,14) - coordination with all neighboring cells
%     tracedata(:,:,15) - persistence. col 1: net distance. Col2: total distance. Col3: net/total distance.
% 
% Uses non built-in functions: 
%     angdiff.m
% 

%% init
function tracedata = CoordinationSpeed(tracedata, rows,cols,SF,EF)

%% params
% coff for 4x no binning
% 25 pixels =40.5 µm radius
% 70 pixels = 113 µm
% 100 pixels = 162 µm
% 63 pixels=100 µm
% 32 pixels=50 µm
% for 10x : 100 µm = 154 pixels
coff=154; % distance within which cell coordination is measured. 
tstep=1; % time step
% SF=50; % startframe
% EF=70; % endframe, should be set to total frames -1
MinDist=1;
displ_coff=[];numcell=[];numcellcoor=[];

%% looping through different wells
vcell=[];
% loop through wells
% for rows=1:8
%     for cols=1:12
old_size = size(tracedata,3);
tracedata(:,:,old_size+1:old_size+6)=NaN(size(tracedata,1),size(tracedata,2),6); % Preallocate space for information added
% loop through frames
for frame=SF:EF
    disp(frame)
    tm1=[tracedata(:,frame,1) tracedata(:,frame,2)]; % x-y coordinates at time point 1
    tm4=[tracedata(:,frame+tstep,1) tracedata(:,frame+tstep,2)]; % x-y coordinates at time point 1+tstep
    len=length(tm1);

    % velocity and angle vectors (atot - vector angle, vtot - vector length);
    d=tm4-tm1; % Vectors connecting two successive points.
    [atot vtot]=cart2pol(d(:,1),d(:,2));
    tracedata(:,frame,9)=atot;
    tracedata(:,frame,10)=vtot;
    % pairwise distances
    % The distances are arranged in the order (2,1), (3,1), ...,
    % (n,1), (3,2), ..., (n,2), ..., (n,n?1))
    dtot=pdist(tm1);
    dtot=squareform(dtot);
    coordination=[];

    %looping through objects
    coor_obj=NaN(length(tracedata),1);
    neighbors=NaN(length(tracedata),1);
    disp(length(tracedata))
    for n=1:length(tracedata)
        
        %disp(n)
        %objects within distance (coff)
        if ~isempty(dtot)
            dff=dtot(n,:); % vector with distances between object n and all other objects
            dff=find(dff<coff & dff>0); % Selects cells within a ring with radius coff
        else
            dff=[];
        end

        % Calculate coordinaton for all neighbors
        if ~isempty(dff) %selects distances unequal zero
            %object vectors
            aval=atot(n); % angle of object n
            %Determine angle difference between object n and its
            angles=atot(dff,:); % angles of neighbors
            angles_diff=angdiff(angles,ones(size(dff')).*aval); % angular difference between query vector and its neighbors
            coor_all=cos(angles_diff); % equivalent to dot(u,v)/(norm(u)*norm(v))

            tracedata(n,frame,old_size+6)=mean(coor_all);
            % identify front, back, lateral neighbors
            % x-y coordinates of all neighbors
            LocNeighbors=[tracedata(dff,frame,1),tracedata(dff,frame,2)];
            % scatter(tracedata(dff,frame,1),tracedata(dff,frame,2)); hold on;
            % Generate vectors connecting cell being analyzed with its neighbors
            Connecting=LocNeighbors-repmat([tracedata(n,frame,1),tracedata(n,frame,2)],length(dff),1);
            %scatter(tracedata(n,frame,1),tracedata(n,frame,2),'r');
            %plot([tracedata(n,frame,1) tracedata(n,frame+1,1)],[tracedata(n,frame,2),tracedata(n,frame+1,2)]);
            [ConnectingAng,~]=cart2pol(Connecting(:,1),Connecting(:,2)); % angle of vectors connecting location of cell n and its neighbors
            % Absolute angular difference between movement of cell n and location of the neigbors in radians
            angles_diffLoc=abs(angdiff(ConnectingAng,ones(size(dff')).*aval));

            % Select cells satisfying "front" criterion (angdiff equal or smaller than 30°)
            frontcells=find(angles_diffLoc<=(pi/6));
            %scatter(tracedata(dff(frontcells),frame,1),tracedata(dff(frontcells),frame,2),'g');
            if ~isempty(frontcells) %
                tracedata(n,frame,old_size+3)=mean(coor_all(frontcells)); % equivalent to dot(u,v)/(norm(u)*norm(v))
            else
                tracedata(n,frame,old_size+3)=NaN;
            end
            % Select cells satisfying "lateral" criterion (angdiff greater than 30° and smaller than 150°)
            sidecells=find(angles_diffLoc>pi/6 & angles_diffLoc<(5*pi/6));
            %scatter(tracedata(dff(sidecells),frame,1),tracedata(dff(sidecells),frame,2),'y');
            if ~isempty(sidecells) %
                tracedata(n,frame,old_size+4)=mean(coor_all(sidecells)); % equivalent to dot(u,v)/(norm(u)*norm(v))
            else
                tracedata(n,frame,old_size+4)=NaN;
            end

            % Select cells satisfying "rear" criterion (angdiff equal or grater than 150° and equal or smaller than 180°)
            backcells=find(angles_diffLoc>=(5*pi/6));
            %scatter(tracedata(dff(backcells),frame,1),tracedata(dff(backcells),frame,2),'b');axis equal;
            if ~isempty(backcells) %
                tracedata(n,frame,old_size+5)=mean(coor_all(backcells)); % equivalent to dot(u,v)/(norm(u)*norm(v))
            else
                tracedata(n,frame,old_size+5)=NaN;
            end

        else
            tracedata(n,frame,old_size+3)=NaN;
            tracedata(n,frame,old_size+4)=NaN;
            tracedata(n,frame,old_size+5)=NaN;
        end
    end %objects
    disp([num2str(coff),'  ',num2str(rows),'  ', num2str(cols),'  ', num2str(frame)]);

end %frame
% compute persistence
% Net distance traveled between SF and EF
LocStart=[tracedata(:,SF,1) tracedata(:,SF,2)]; % x-y coordinates at time point SF
LocEnd=[tracedata(:,EF,1) tracedata(:,EF,2)]; % x-y coordinates at time point EF
StartEnd=LocEnd-LocStart;
[TotAng TotDist]=cart2pol(StartEnd(:,1),StartEnd(:,2));
tracedata(:,1,old_size+7)=TotDist;

% Sum of displacement per frame:
velocity=tracedata(:,:,10);
tracedata(:,2,old_size+7)=sum(velocity(:,SF:(EF-1)),2);

% Persistence = Net distance / total distance
tracedata(:,3,old_size+7)=tracedata(:,1,old_size+7)./tracedata(:,2,old_size+7);


disp('done!');
