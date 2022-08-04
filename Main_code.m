
clear 

% We need to provide the input files, which are 
% 1.A DEM.tiff file of the area (containg the grid data)
% 2.A KML file of conotur limits ( containing the boundary 

Dem_file='frasse_data.tif'  % put the name of Dem file in tiff formet 
KML_file='frasse_boundary.kml'  % put the name of kml file in kml formate 



% 3 approaches have been discussed here for getting the failure surface and
% volume. Based on the grid data and the kml file of counture limits of landslide
% Method 1-Constant average Dip angle for entire failure surface (assumed
%          that dip angle is constant for failure surface.
% Method 2-Varying Dip anle (weightaga is given based on-
%          distance of intersection points from elevation line ). Higher the
%          distance lower the weightage for Dip and strike angle variations
%          Elevation line is horizontal projection of steepest slope.
% Method 3-Varying Dip angle(The left and right profile has been
%          assumed as an arc of big Circle) . And then the Dip angle have
%          been linearly distributed to the intersection points
% Enter any of the 3 methods in the input option to proceed further
Option=input(['Option 1-Open a saved Project Enter-1 \n' ...
             'Option 2-Open a New project Enter-2\n' ...
             'Selct one of the above\n'])
if Option==1
    [file,path] = uigetfile('*.mat')
    pro=load(file);
    Method=pro.method
    if Method==1
        Right_Dip_angle =pro.r
        Left_Dip_angle  =pro.l
    else
        Rgt.top_dip= pro.r.top_dip
        Rgt.top_strik= pro.r.top_strik
        Rgt.botm_dip= pro.r.botm_dip
        Rgt.botm_strik= pro.r.botm_strik
        
        Lft.top_dip= pro.l.top_dip
        Lft.top_strik= pro.l.top_strik
        Lft.botm_dip= pro.l.botm_dip
        Lft.botm_strik= pro.l.botm_strik
        
    end
elseif Option==2
    
    Method = input (['Method 1-Constant average Dip angle  enter-1 \n' ...
                 'Method 2-Varying Dip and strik anle (waitege based on distance from elevation line) enter-2\n' ...
                 'Method 3-Varying Dip and strike angle(each profile is an arc of big Circle) enter-3\n' ...
                'Selct one of the above\n'])
        if (Method > 3 || Method <1)
           error('Wrong selection: Please select between 1 to 3')     
        end
    
    %For calculation the strike and dip angles are required.For right profile
    %of the failure surface as well as for the left profile. Provide the
    %suitable values for selected method
    
    if (Method==1)
       Constant_Dip = inputdlg({'Enter Right profile Dip angle','Enter Left profile Dip angle'},'Constant Dip', [1 50; 1 50]) 
       Right_Dip_angle = str2double(Constant_Dip{1})
       Left_Dip_angle = str2double(Constant_Dip{2})
       
       pro.r=Right_Dip_angle;
       pro.l=Left_Dip_angle;
       pro.method=Method
       ask=input(['Press-1 to Save data entered \n' ...
           'Or Press Enter to move farwar(without saving) \n' ...
           '(use the suitable formate to title to remember the name- like-dip_strik30_40_method1)\n'])
       if ask==1
           keep=uiputfile('*.mat');
           save(keep,'-struct','pro');
       end

             
    elseif(Method==2 || Method==3) 
        Right_profile = inputdlg({' Enter Dip angle(top side)','Enter Dip angle(bottom side)'}, ...
        'Right profile', [ 1 50; 1 50],{'50','20'})
        
        Rgt.top_dip= str2double(Right_profile{1});
        Rgt.top_strik= 120; %str2double(Right_profile{2});
        Rgt.botm_dip= str2double(Right_profile{2});
        Rgt.botm_strik=240;% str2double(Right_profile{4})
        
        if Rgt.top_strik > Rgt.botm_strik
            
            error("wrong input\n Strike at the top, should be less than at the bottom")           
        end
    
        Left_profile = inputdlg({' Enter Dip angle(top side)','Enter Dip angle(bottom side)'}, ...
        'Left Profile', [1 50; 1 50],{'50','20'})
         
        Lft.top_dip= str2double(Left_profile{1});
        Lft.top_strik= 240; %str2double(Left_profile{2});
        Lft.botm_dip= str2double(Left_profile{2});
        Lft.botm_strik=120; % str2double(Left_profile{4})
        
        if Lft.top_strik < Lft.botm_strik
            error("wrong input\n Strike ate the top should be greater than at the bottom")
        end
        pro.r=Rgt;
        pro.l=Lft;
        pro.method=Method

        ask=input(['Press-1 then Enter to Save data entered \n' ...
           'Or Press Enter to move farward (without saving)\n'])
       if ask==1
           keep=uiputfile('*.mat');
           save(keep,'-struct','pro');
       end
        
    end

end   
%####################################################################################
%Calculations start from here --- 
%Importing the Grid Data and Coordinate transformation 
%(Rotating the x-axis along cross-section lines and Y-axis along the elevation line  

[Z_mesh, Rf]= readgeoraster (Dem_file)

jh=isprop(Rf,'LatitudeLimits')
if jh==1
    compar=strcmp(Rf.CoordinateSystemType,'geographic')
    if compar==1
        [xc,yc,~] = deg2utm(Rf.LatitudeLimits(1), Rf.LongitudeLimits(1));
        [xc1,yc1,utmzone] = deg2utm(Rf.LatitudeLimits(2),Rf.LongitudeLimits(2));
        utmzone
    else
 
       xc = Rf.LongitudeLimits(1)
       yc = Rf.LatitudeLimits(1)
       xc1= Rf.LongitudeLimits(2)
       yc1= Rf.LatitudeLimits(2)
    end
else
    compar=strcmp(Rf.CoordinateSystemType,'geographic')
    if compar==1
       [xc,yc,~] =deg2utm(Rf.YWorldLimits(1), Rf.XWorldLimits(1));
       [xc1,yc1,utmzone] = deg2utm(Rf.YWorldLimits(2), Rf.XWorldLimits(2));
       utmzone
       %if error is like unrecognized method/property or field "...." for
       %class-____>> replace it with suitable name fo limits given in Rf.
       %structure  in the right side(as Rf.LatitudeLimits)
    else
        xc = Rf.XWorldLimits(1)
        yc = Rf.YWorldLimits(1)
        xc1= Rf.XWorldLimits(2)
        yc1= Rf.YWorldLimits(2)
    end
end
%xv = linspace(Rf.LongitudeLimits(1),Rf.LongitudeLimits(2),Rf.RasterSize(2))
%yv = linspace(Rf.LatitudeLimits(2), Rf.LatitudeLimits(1),Rf.RasterSize(1))
Grid.x = linspace(xc,xc1,Rf.RasterSize(2));
Grid.y = linspace(yc1,yc,Rf.RasterSize(1))

column=Rf.RasterSize(2)
  rows=Rf.RasterSize(1)
  Size_grid=size(Z_mesh)
[X_mesh,Y_mesh] = meshgrid(Grid.x, Grid.y)

Xcell_size=(max(Grid.x)-min(Grid.x))/(Rf.RasterSize(2)-1) % in meter
szzx=X_mesh(1,2)-X_mesh(1,1)
Ycell_size=(max(Grid.y)-min(Grid.y))/(Rf.RasterSize(1)-1) % in meter
szzy=Y_mesh(1,1)-Y_mesh(2,1)
surf(X_mesh, Y_mesh, Z_mesh);
title('Grid');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
if (rows*column > 3000 && rows*column <13000)
    l=0
    for i=1:2:column
        l=l+1;
        m=0;
        for j=1:2:rows
            m=m+1;
            xxx(m,l)=X_mesh(j,i);
            yyy(m,l)=Y_mesh(j,i);
            zzz(m,l)=Z_mesh(j,i);
            
        end
        
    end
    l
    m
    X_mesh=xxx
    Y_mesh=yyy
    Z_mesh=zzz
    rows=m
    column=l
    Xcell_size=(max(Grid.x)-min(Grid.x))/(column-1) % in meter
    Ycell_size=(max(Grid.y)-min(Grid.y))/(rows-1) % in meter
    szzx=2*szzx
    szzy=2*szzy
elseif (rows*column > 13000)
        l=0
    for i=1:3:column
        l=l+1;
        m=0;
        for j=1:3:rows
            m=m+1;
            xxx(m,l)=X_mesh(j,i);
            yyy(m,l)=Y_mesh(j,i);
            zzz(m,l)=Z_mesh(j,i);
            
        end
        
    end
    l
    m
    X_mesh=xxx
    Y_mesh=yyy
    Z_mesh=zzz
    rows=m
    column=l
    Xcell_size=(max(Grid.x)-min(Grid.x))/(column-1) % in meter
    Ycell_size=(max(Grid.y)-min(Grid.y))/(rows-1) % in meter
    szzx=3*szzx
    szzy=3*szzy
    
end
surf(X_mesh, Y_mesh, Z_mesh);

%Converting the grid into columns of x, y and z
x=zeros(1,column*rows);
y=zeros(1,column*rows);
z=zeros(1,column*rows);
l=1;
for i=1:column
    for j=1:rows
        x(l)=X_mesh(j,i);
        y(l)=Y_mesh(j,i);
        z(l)=Z_mesh(j,i);
        l=l+1;
   end
end
l
total.x=x;
total.y=y;
total.z=z

%Importing the contour curve(Limits of the failur boundaries) as a KML file

crv=kml2struct(KML_file) 
x0=(crv.Lon)  
y0=(crv.Lat)

%crv =load('border.dat')
%x0=crv(:,1)
%y0=crv(:,2)
%Conversion of Lat Lon into UTM projection(converting the
% Geographic coordinate system into WGS 1984 UTM projection) and plotting
if x0<90
    [x0,y0,utmzone] = deg2utm(y0,x0);
end
contor.x=x0
contor.y=y0

plot(total.x,total.y,'r+');
title('grid Points');
xlabel('Longitude');
ylabel("Latitude");
axis([(1-.00001)*min(total.x),(1+.00001)*max(total.x) (1-.00001)*min(total.y),(1+.00001)*max(total.y)])

rtm1=[1.001 0; 0 1.001];
m=[contor.x';contor.y']
%((max(x0)+min(x0))/2)*0.1
uj=rtm1*m;
quary.Xt=uj(1,:)-((max(contor.x)+min(contor.x))/2)*0.001;
quary.Yt=uj(2,:)-((max(contor.y)+min(contor.y))/2)*0.001;
plot(contor.x,contor.y);
title('Landslide Contour');
xlabel('Longitude');
ylabel("Latitude");
hold on
plot(quary.Xt,quary.Yt)
plot(total.x,total.y,'r+')

%Extracting the points lying in side or on the contour boundaries
[in,on] = inpolygon(total.x,total.y,quary.Xt,quary.Yt)
in_boundary=numel(total.x(in))
on_boundary=numel(total.x(on))
quary.x=[total.x(in), total.x(on)];
quary.y=[total.y(in), total.x(on)];
quary.z=[total.z(in), total.x(on)]

%Finding the line passing through heighest elevation point and lowest
%elevation point
[m,j]=max(quary.z)
[mi,k]=min(quary.z)
xma=quary.x(j)
yma=quary.y(j)
zma=quary.z(j)
xmi=quary.x(k)
ymi=quary.y(k)
zmi=quary.z(k)
if xma< xmi
   elevation_x=xma:xmi 
else
   elevation_x=xmi:xma 
end

elevation_y=(elevation_x-xmi)*(yma-ymi)/(xma-xmi)+ymi
plot(elevation_x,elevation_y)
plot(quary.x,quary.y,'b+') % points inside and on the contour
axis([(1-.00001)*min(total.x),(1+.00001)*max(total.x) (1-.00001)*min(total.y),(1+.00001)*max(total.y)])
hold off

%Angle of rotation(coordinate transformation angle)
mm=(yma-ymi)/(xma-xmi)
angle=atand(mm);
angle_from_yaxis=atand(1/mm)
rotation_angle=angle_from_yaxis  %rotation angle

%Avarage slope of elevation
delta_z=zma-zmi
distance_zmax_zmin=sqrt((xma-xmi)*(xma-xmi)+(yma-ymi)*(yma-ymi))
Avarage_slope=delta_z/distance_zmax_zmin
slop_Angle=atand(Avarage_slope)




%Coordinate transformation of quary data 
rtn_matrix=[cosd(rotation_angle) -sind(rotation_angle); sind(rotation_angle) cosd(rotation_angle)];
mq=[quary.x;quary.y];
uj=rtn_matrix*mq;
quary.Xt=uj(1,:)
quary.Yt=uj(2,:)

[quary.Yts, sortIndex] = sort(quary.Yt,"descend");
quary.Xts = quary.Xt(sortIndex);
quary.Zts = quary.z(sortIndex);


%Coordinate tronsformation(rotating the axis in direction of cross-section line and elevation line) of imported curve and plotting

rtn_matrix=[cosd(rotation_angle) -sind(rotation_angle); sind(rotation_angle) cosd(rotation_angle)]
mc=[contor.x';contor.y'];
Vfil=rtn_matrix*mc;
contor.xt=Vfil(1,:)
contor.yt=Vfil(2,:)

figure
plot(contor.xt,contor.yt);
hold on
title('Contour in Transformed coordinate');
xlabel('Longitude');
ylabel("Latitude");

%Elevation line rotation(Line joining the maximum and minimum elevation)
if xma< xmi
   elevation_x=xma:xmi 
else
   elevation_x=xmi:xma 
end
elevation_y=(elevation_x-xmi)*(yma-ymi)/(xma-xmi)+ymi

mcelv=[elevation_x;elevation_y];
Vfile=rtn_matrix*mcelv;
elevation_xt=Vfile(1,:)
elevation_yt=Vfile(2,:)
plot(elevation_xt,elevation_yt)
Xt_vertical=elevation_xt(1)
hold off

%Plotting of the crosssection lines
Yp=quary.Yts';
size(Yp);
cross_x=min(contor.xt)-2:1:max(contor.xt)+2;
cross_y=Yp+cross_x*0;
plot(cross_x,cross_y);
hold on
plot(quary.Xt,quary.Yt,'r+')
hold off
title('Cross Section lines in Transformed coordinate');
xlabel('Longitude');
ylabel("Latitude");



% Plotting the cross section and curve together in transformed coordinate
figure
plot(contor.xt,contor.yt);
hold on
plot(elevation_xt,elevation_yt)
plot(cross_x,cross_y);
title('Contour & line intersection');
xlabel('Longitude');
ylabel("Latitude");
hold off;


%Intersection of cross sections with curve and sorting(sorting the data
%based on y value for ease of calculation
Intersect_xy=zeros(length(quary.y),4)
for i=1:1:length(quary.y)
    Rt=cross_y(i,:);
    
   [xi,yi] = polyxpoly(contor.xt,contor.yt,cross_x,Rt);

    if length(xi)==2  
        Intersect_xy(i,1)=xi(1);
        Intersect_xy(i,3)=xi(2);
        Intersect_xy(i,2)=yi(1);
        Intersect_xy(i,4)=yi(2);
    elseif length(xi)==4
        if quary.Xts(i) >= xi(1) && quary.Xts(i) <= xi(2) 
           Intersect_xy(i,1)=xi(1);
           Intersect_xy(i,3)=xi(2);
           Intersect_xy(i,2)=yi(1);
           Intersect_xy(i,4)=yi(2);
        else
           Intersect_xy(i,1)=xi(3);
           Intersect_xy(i,3)=xi(4);
           Intersect_xy(i,2)=yi(3);
           Intersect_xy(i,4)=yi(4);
        end  
        
    elseif length(xi)==1
        Intersect_xy(i,1)=xi(1);
        Intersect_xy(i,3)=xi(1);
        Intersect_xy(i,2)=yi(1);
        Intersect_xy(i,4)=yi(1);
        
    elseif length(xi)==3
        if quary.Xts(i) >= xi(1) && quary.Xts(i) <=xi(2)
            Intersect_xy(i,1)=xi(1);
            Intersect_xy(i,3)=xi(2);
            Intersect_xy(i,2)=yi(1);
            Intersect_xy(i,4)=yi(2);
        else
            Intersect_xy(i,1)=xi(2);
            Intersect_xy(i,3)=xi(3);
            Intersect_xy(i,2)=yi(2);
            Intersect_xy(i,4)=yi(3);            
            
        end  
    elseif length(xi)==6
        if quary.Xts(i) >= xi(3) && quary.Xts(i) <=xi(4)
            Intersect_xy(i,1)=xi(3);
            Intersect_xy(i,3)=xi(4);
            Intersect_xy(i,2)=yi(3);
            Intersect_xy(i,4)=yi(4);            
        else
            Intersect_xy(i,1)=xi(5);
            Intersect_xy(i,3)=xi(6);
            Intersect_xy(i,2)=yi(5);
            Intersect_xy(i,4)=yi(6);           
        end
    else 
        Intersect_xy(i,1)=xi(1);
        Intersect_xy(i,3)=xi(2);
        Intersect_xy(i,2)=yi(1);
        Intersect_xy(i,4)=yi(2);       
    end
end

Intersect_xy
sect.Xlst=Intersect_xy(:,1);
sect.Ylst=Intersect_xy(:,2);
sect.Lt=[sect.Xlst,sect.Ylst]


sect.Xrst=Intersect_xy(:,3);
sect.Yrst=Intersect_xy(:,4);
sect.Rt=[sect.Xrst,sect.Yrst]

sect.LRt=[sect.Lt;sect.Rt];

quary.xs=quary.x(sortIndex);
quary.ys=quary.y(sortIndex);
quary.zs=quary.z(sortIndex);

%Transforming into original coordinate system(Retaing the original
%coordinates after calculations)
rtn_matrix=[cosd(-rotation_angle) -sind(-rotation_angle); sind(-rotation_angle) cosd(-rotation_angle)]
sect.L=rtn_matrix*sect.Lt';
sect.R=rtn_matrix*sect.Rt';
sect.LR=[sect.L sect.R]';
xy=[quary.xs' quary.ys']
xycomb=[xy;sect.LR]


% Calculating the Elevation(z values ) for Intersection Points. 

zi=[total.z';nan(length(sect.LR),1)]'
sect.LR(:,1)
xi=[total.x';sect.LR(1:length(sect.LR),1)]'
yi=[total.y';sect.LR(1:length(sect.LR),2)]'

[zgrid,xgrid,ygrid] = gridfit(xi,yi,zi,min(xi):10:max(xi),min(yi):10:max(yi));

figure
surf(xgrid,ygrid,zgrid);
hold on
plot3(xi,yi,zi,'ro');
zgrid
title('Grid for quary points');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
hold off

%Left Profile z value calculation for intersection points
znanL = isnan(zi(1,1:length(zi)-length(sect.L)))
ZfilL = interp2(xgrid,ygrid,zgrid,xi(znanL),yi(znanL))
voidL=isnan(ZfilL)

VfilL=find(voidL==1)

for a=1:1:length(VfilL)
    if VfilL(1,a)-3>0 && VfilL(1,a)+3<length(zi)
        if isnan(ZfilL(1,VfilL(1,a)-1)) && ~isnan(ZfilL(1,VfilL(1,a)+1))
            if isnan(ZfilL(1,VfilL(1,a)-2)) 
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-3)+ZfilL(1,VfilL(1,a)+1))/2;
            else, isnan(ZfilL(1,VfilL(1,a)-2)) 
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-2)+ZfilL(1,VfilL(1,a)+1))/2;
            end 
           
        elseif ~isnan(ZfilL(1,VfilL(1,a)-1)) &&  isnan(ZfilL(1,VfilL(1,a)+1))
            if isnan(ZfilL(1,VfilL(1,a)+2))
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)+3))/2;
            else
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)+2))/2; 
            end
        elseif isnan(ZfilL(1,VfilL(1,a)+1)) && isnan(ZfilL(1,VfilL(1,a)-1))
            if isnan(ZfilL(1,VfilL(1,a)+2)) && ~isnan(ZfilL(1,VfilL(1,a)-2))
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-2)+ZfilL(1,VfilL(1,a)+3))/2;
            elseif ~isnan(ZfilL(1,VfilL(1,a)+2)) && isnan(ZfilL(1,VfilL(1,a)-2))
                ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-3)+ZfilL(1,VfilL(1,a)+2))/2;
            end
        %elseif Vfil==length(Zfil)/2+1     
            %Zfil(1,Vfil(1,a))=(Zfil(1,Vfil(1,a)+1)+Zfil(1,Vfil(1,a)+2))/2;
       % elseif Vfil==length(Zfill)/2-1     
            %Zfil(1,Vfil(1,a))=(Zfil(1,Vfil(1,a)-11)+Zfil(1,Vfil(1,a)-2))/2;
        else
            
         ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)+1))/2;
         
        end 
    
    else
        if VfilL(1,a)-1==0
            ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)+1)+ZfilL(1,VfilL(1,a)+2))/2;
        elseif VfilL(1,a)-1>0 && VfilL(1,a)-3<=0
            ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)+1))/2; 
        elseif VfilL(1,a)+1==length(zi)
            ZfilL(1,VfilR(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)-2))/2; 
        elseif VfilL(1,a)+1<length(zi) && VfilL(1,a)+3>=length(zi)
            ZfilL(1,VfilL(1,a))=(ZfilL(1,VfilL(1,a)-1)+ZfilL(1,VfilL(1,a)+1))/2;        
           
        end
    end
end
ZfilL
check=isnan(ZfilL);
find(check==1)

%Right Profile z value calculation for intersection points
%zi=[total.z';nan(length(sect.R),1)]'
totl=[total.z';ZfilL';zi(1,length(zi)+1-length(sect.L):length(zi))']'
znanR = isnan(totl)
ZfilR = interp2(xgrid,ygrid,zgrid,xi(znanR),yi(znanR))
voidR=isnan(ZfilR)

VfilR=find(voidR==1)

for a=1:1:length(VfilR)
    if VfilR(1,a)-3>0 && VfilR(1,a)+3<length(zi) 
        if isnan(ZfilR(1,VfilR(1,a)-1)) && ~isnan(ZfilR(1,VfilR(1,a)+1))
            if isnan(ZfilR(1,VfilR(1,a)-2)) 
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-3)+ZfilR(1,VfilR(1,a)+1))/2;
            else, isnan(ZfilR(1,VfilR(1,a)-2)) 
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-2)+ZfilR(1,VfilR(1,a)+1))/2;
            end 
           
        elseif ~isnan(ZfilR(1,VfilR(1,a)-1)) &&  isnan(ZfilR(1,VfilR(1,a)+1))
            if isnan(ZfilR(1,VfilR(1,a)+2))
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)+3))/2;
            else
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)+2))/2; 
            end
        elseif isnan(ZfilR(1,VfilR(1,a)+1)) && isnan(ZfilR(1,VfilR(1,a)-1))
            if isnan(ZfilR(1,VfilR(1,a)+2)) && ~isnan(ZfilR(1,VfilR(1,a)-2))
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-2)+ZfilR(1,VfilR(1,a)+3))/2;
            elseif ~isnan(ZfilR(1,VfilR(1,a)+2)) && isnan(ZfilR(1,VfilR(1,a)-2))
                ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-3)+ZfilR(1,VfilR(1,a)+2))/2;
            end
        %elseif Vfil==length(Zfil)/2+1     
            %Zfil(1,Vfil(1,a))=(Zfil(1,Vfil(1,a)+1)+Zfil(1,Vfil(1,a)+2))/2;
       % elseif Vfil==length(Zfill)/2-1     
            %Zfil(1,Vfil(1,a))=(Zfil(1,Vfil(1,a)-11)+Zfil(1,Vfil(1,a)-2))/2;
        else
            
         ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)+1))/2;
        
        end   
    else
        if VfilR(1,a)-1==0
            ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)+1)+ZfilR(1,VfilR(1,a)+2))/2;
        elseif VfilR(1,a)-1>0 && VfilR(1,a)-3<=0
            ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)+1))/2; 
        elseif VfilR(1,a)+1==length(zi)
            ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)-2))/2; 
        elseif VfilR(1,a)+1<length(zi) && VfilR(1,a)+3>=length(zi)
            ZfilR(1,VfilR(1,a))=(ZfilR(1,VfilR(1,a)-1)+ZfilR(1,VfilR(1,a)+1))/2;        
           
        end
    end
end
ZfilR
check=isnan(ZfilR);
Ck=find(check==1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Combining the calculated elevation data ate intersected points into a
%single matrix for further proccessing

sect.Zl=ZfilL'
sect.Zr=ZfilR'
sect.Zlr=[sect.Zl;sect.Zr];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3 Different approches have been used to get the cross-section angles with
% failure surface for each cross-section lines, so that they can be used in
% calculation of elevation using spline method in the next step.

%Calculation of cross-section angles with failure surface for the right profile 

%Right profile-Right portion of the failure surface divided by the
% elevation line
% Elevation line-The line joining the max elevation to min within the contour boundaries


% For Right profile
Intersect_xy;
sect.Xrst=Intersect_xy(:,3);
sect.Yrst=Intersect_xy(:,4);
sect.Rt=[sect.Xrst,sect.Yrst]

sect.LRt=[sect.Lt;sect.Rt];

len=length(sect.Lt)
rrx=sect.Rt(:,1)';
rry=sect.Rt(:,2)';
plot(rrx,rry,'*');

daspect([1 1 1])
%axis square
hold on 
Yrn=rry';
size(Yrn);
xre=min(rrx)-2:1:max(rrx)+2;
cross_y=Yrn+xre*0;
plot(xre,cross_y);

xlabel('Longitude');
ylabel(" Latitude          ");
%Method-1 Right Profile
if Method==1
    % The true angle between cross section line and falure plane Method-1
    title('              Right Profile');
    hold off
    Csr=zeros(1,len)'; 
    Right_Dip_angle;
    for i=1:len
        Csr(i)=Right_Dip_angle;
    end
    Csr
%Method-2 Right Profile    
elseif(Method==2) 
    title('             Right Profile');
    hold off
   % Rgt.top_dip= str2double(Right_profile{1})
   % Rgt.top_strik= str2double(Right_profile{2})
   % Rgt.botm_dip= str2double(Right_profile{3})
   % Rgt.botm_strik= str2double(Right_profile{4})
    
    % Right Profile Dip angle Method-2
    sumr=0;
    for i=1:1:len
        sumr=1/abs(abs(sect.Xrst(i))-abs(Xt_vertical))+sumr;
    end
    Weightr=abs(sumr);
    kr=Rgt.top_dip;
    Dr=zeros(1,len)';
    Dr(1,1)=Rgt.top_dip;
    
    for i=1:1:len-1
        dip=kr-(abs(1/(abs(sect.Xrst(i))-abs(Xt_vertical)))/Weightr)*(Rgt.top_dip-Rgt.botm_dip);
        kr=dip;
        Dr(i+1,1)=dip;
    end
    Dr    
    
    %Right Profile strik angle Method-2
    Sr=zeros(1,len)';
    Sr(1,1)=Rgt.top_strik;
    sumr=0;
    for i=1:1:len
        sumr=1/abs(abs(sect.Xrst(i))-abs(Xt_vertical))+sumr;
    end
    lr=Rgt.top_strik;
    Weightr=abs(sumr);
    for i=1:1:len-1
        por=lr+(abs(1/(abs(sect.Xrst(i))-abs(Xt_vertical)))/Weightr)*(Rgt.botm_strik-Rgt.top_strik);
        lr=por;
        Sr(i+1,1)=por;
    end
    Sr
    
    % The true angle between cross section line and falure plane Method-2
    Csr=zeros(1,len)';    
    for a=1:1:len
       uu=tand(Dr(a,1))*cosd(180-Sr(a,1));
       fd=atand(uu);
       Csr(a,1)=fd;
    end 
    Csr
%Method-3 Right Profile    
elseif(Method==3)
    title('Intersection of normals from the approximated arc');

   % Rgt.top_dip= str2double(Right_profile{1})
   % Rgt.top_strik= str2double(Right_profile{2})
   % Rgt.botm_dip= str2double(Right_profile{3})
   % Rgt.botm_strik= str2double(Right_profile{4})
    
    normr1=Rgt.top_strik-90
    normr2=Rgt.botm_strik-90
    
    slopr1=tand(90-normr1);
    slopr2=tand(90-normr2);
    
     xrm=(sect.Yrst(1,1)+slopr2*sect.Xrst(length(sect.Xrst),1)- ...
        sect.Yrst(length(sect.Yrst),1)-slopr1*sect.Xrst(1,1))/(slopr2-slopr1)
    yrm=(slopr2*(sect.Yrst(1,1)+slopr1*sect.Xrst(length(sect.Xrst),1)) ...
        -slopr1*(sect.Yrst(length(sect.Yrst),1)+slopr2*sect.Xrst(1,1)))/(slopr2-slopr1)
    %[xrm,yrm] = polyxpoly(xr1,yr1,xr1,yr2)   
    
    if(sect.LRt(1,1)<0)
        eer1=xrm*1.00005;
        eer2=max(sect.LRt(:,1))*.99990;
    else
        eer1=xrm*.99995;
        eer2=max(sect.LRt(:,1))*1.00009;        
    end

    xr1=eer1:1:eer2;
    
    yr1=sect.Yrst(1,1)+(xr1-sect.Xrst(1,1))*slopr1;
    yr2=sect.Yrst(length(sect.Yrst),1)+(xr1-sect.Xrst(length(sect.Xrst),1))*slopr2;

    
    
    plot(xr1,yr1,xr1,yr2)
    axis([eer1 eer2 min(sect.LRt(:,2)) max(sect.LRt(:,2))])
    hold on
    plot(xrm, yrm,'*')
    hold off
    
    %Right Profile Dip angle Method-3
    mr=(Rgt.botm_dip-Rgt.top_dip)/(Rgt.botm_strik-Rgt.top_strik);
    kr=Rgt.top_dip;
    Dr=zeros(1,len)';
    Dr(1,1)=Rgt.top_dip;
    
    for a=1:1:len-1
        
        i=sect.Xrst(a)-sect.Xrst(a+1);
        j=sect.Yrst(a)-sect.Yrst(a+1);
        hh=sqrt(i*i+j*j);
        ii=sect.Xrst(a)-xrm;
        jj=sect.Yrst(a)-yrm;
        ln=sqrt(ii*ii+jj*jj);
        Sr=hh/ln;
        do=(180/pi)*Sr;
        gd=atand(hh/ln);
        dip=gd*mr+kr;
        kr=dip;
        Dr(a+1,1)=dip;
    end
    Dr

   %Right Profile strike angle Method-3
    Sr=zeros(1,len)';
    lr=Rgt.top_strik;
    Sr(1,1)=Rgt.top_strik;
    
    for a=1:1:len-1 
        i=sect.Xrst(a)-sect.Xrst(a+1);
        j=sect.Yrst(a)-sect.Yrst(a+1);
        hh=sqrt(i*i+j*j);
        ii=sect.Xrst(a)-xrm;
        jj=sect.Yrst(a)-yrm;
        ln=sqrt(ii*ii+jj*jj);
        u=hh/ln;
        gd=atand(u);
        strik=gd+lr;
        lr=strik;
        Sr(a+1,1)=strik;
    end
    Sr
    
    % The true angle between cross section line and falure plane Method-3
    Csr=zeros(1,len)';
    
    for a=1:1:len
       uu=tand(Dr(a,1))*cosd(180-Sr(a,1));
       fd=atand(uu);
       Csr(a,1)=fd;
    end 
    Csr

end    
    

%Calculation of cross-section angles with failure surface for the Left profile 
%Left profile-Right portion of the failure surface divided by the elevation line
% Elevation line- The line joining the max elevation to min within the contour boundaries


%Left profile

Intersect_xy;
sect.Xlst=Intersect_xy(:,1);
sect.Ylst=Intersect_xy(:,2);
sect.Lt=[sect.Xlst,sect.Ylst];

len=length(sect.Lt)
lld=sect.Lt(:,1)';
lle=sect.Lt(:,2)';
plot(lld,lle,'*')
hold on

daspect([1 1 1]);
%axis square

Yln=lle';
size(Yln);
xle=min(lld)-4:1:max(lld)+2;
Ylfn=Yln+xle*0;
plot(xle,Ylfn);

xlabel('Longitude               ');
ylabel("Latitude");

%Method-1 Left Profile
% The true angle between cross section line and plane Method-1

if (Method==1)
    hold off
    title('             Left Profile');
    Csl=zeros(1,len)';
    Left_Dip_angle;
    for i=1:len
        Csl(i)=Left_Dip_angle;
    end
    Csl
%Method-2 Left Profile 

elseif(Method==2)
    title('              Left Profile');
    hold off
   % Lft.top_dip= str2double(Left_profile{1})
   % Lft.top_strik= str2double(Left_profile{2})
   % Lft.botm_dip= str2double(Left_profile{3})
   % Lft.botm_strik= str2double(Left_profile{4})
    
    suml=0;
    for i=1:1:len
        suml=1/abs(abs(sect.Xlst(i))-abs(Xt_vertical))+suml;
    end
    %Left profile Strike angle Method-2
    Weightl=abs(suml);
    kl=Lft.top_dip;
    Dl=zeros(1,len)';
    Dl(1,1)=Lft.top_dip;
    
    for i=1:1:len-1
        dip=kl-(abs(1/(abs(sect.Xrst(i))-abs(Xt_vertical)))/Weightl)*(Lft.top_dip-Lft.botm_dip);
        kl=dip;
        Dl(i+1,1)=dip;
    end
    Dl   
    
    %Left profile Strike angle Method-2  
    Sl=zeros(1,len)'; 
    Sl(1,1)=Lft.top_strik;
    suml=0;
    for i=1:1:len
        suml=1/abs(abs(Xt_vertical)-abs(sect.Xlst(i)))+suml;
    end
    ll=Lft.top_strik;
    Weightl=abs(suml);
    for i=1:1:len
        po=ll-(abs(1/(abs(Xt_vertical)-abs(sect.Xlst(i))))/Weightl)*(Lft.top_strik-Lft.botm_strik);
        ll=po;
        Sl(i+1,1)=po;
    end
    Sl
    
    % The true angle between cross section line and plane Method-2
    Csl=zeros(1,len)';
    
    for a=1:1:len
        
       uu=tand(Dl(a,1))*cosd(180-Sl(a,1));
       fd=atand(uu);
       Csl(a,1)=fd;
    end 
    
    Csl
% Method-3 Left Profile  

elseif (Method==3)
    title('Intersection of normals from the approximated arc');
   % Lft.top_dip= str2double(Left_profile{1});
   % Lft.top_strik= str2double(Left_profile{2});
   % Lft.botm_dip= str2double(Left_profile{3});
   % Lft.botm_strik= str2double(Left_profile{4})
    
    norml1=Lft.top_strik-90
    norml2=Lft.botm_strik-90
    
    slopl1=tand(90-norml1);
    slopl2=tand(90-norml2);
    
    xlm=(sect.Ylst(1,1)+slopl2*sect.Xlst(length(sect.Xlst),1)- ...
        sect.Ylst(length(sect.Ylst),1)-slopl1*sect.Xlst(1,1))/(slopl2-slopl1)
    ylm=(slopl2*(sect.Ylst(1,1)+slopl1*sect.Xlst(length(sect.Xlst),1)) ...
        -slopl1*(sect.Ylst(length(sect.Ylst),1)+slopl2*sect.Xlst(1,1)))/(slopl2-slopl1)    
    
    if(sect.LRt(1,1)<0)
        eel1=min(sect.LRt(:,1))*1.00005;
        eel2=xlm*.99995;
    else
        eel1=min(sect.LRt(:,1))*.99995;
        eel2=xlm*1.00005;        
    end
    xl1=eel1:eel2;
    
    yl1=sect.Ylst(1,1)+(xl1-sect.Xlst(1,1))*slopl1;
    yl2=sect.Ylst(length(sect.Ylst),1)+(xl1-sect.Xlst(length(sect.Xlst),1))*slopl2;

    
    %[xlm,ylm] = polyxpoly(xl1,yl1,xl1,yl2)
    plot(xl1,yl1,xl1,yl2)
    hold on
    axis([eel1 eel2 min(sect.LRt(:,2)) max(sect.LRt(:,2))]);
    
    plot(xlm, ylm,'*');
    hold off  
    
    %Left profile Strike angle Method-3
    ml=(Lft.botm_dip-Lft.top_dip)/(Lft.botm_strik-Lft.top_strik);
    kl=Lft.top_dip;
    Dl=zeros(1,len)';
    Dl(1,1)=Lft.top_dip;
    
    for a=1:1:len-1     
        i=sect.Xlst(a)-sect.Xlst(a+1);
        j=sect.Ylst(a)-sect.Ylst(a+1);
        hh=sqrt(i*i+j*j);
        ii=sect.Xlst(a)-xlm;
        jj=sect.Ylst(a)-ylm;
        ln=sqrt(ii*ii+jj*jj);
        Sr=hh/ln;
        do=(180/pi)*Sr;
        gdl=atand(hh/ln);
        dipl=-gdl*ml+kl;
        kl=dipl;
        Dl(a+1,1)=dipl;
    end
    Dl

    %Left Profile strike angle Method-3
    Sl=zeros(1,len)';
    lr=Lft.top_strik;
    Sl(1,1)=Lft.top_strik;
    
    for a=1:1:len-1
        i=sect.Xlst(a)-sect.Xlst(a+1);
        j=sect.Ylst(a)-sect.Ylst(a+1);
        hh=sqrt(i*i+j*j);
        ii=sect.Xlst(a)-xlm;
        jj=sect.Ylst(a)-ylm;
        ln=sqrt(ii*ii+jj*jj);
        u=hh/ln;
        gdl=atand(u);
        strik=-gdl+lr;
        lr=strik;
        Sl(a+1,1)=strik;
    end
    Sl
    % The true angle between cross section line and plane Method-3
    Csl=zeros(1,len)';
    
    for a=1:1:len
        
       uu=tand(Dl(a,1))*cosd(180-Sl(a,1));
       fd=atand(uu);
       Csl(a,1)=fd;
    end 
    Csl

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sect.Zl=ZfilL'
sect.Zr=ZfilR'
sect.Zlr=[sect.Zl;sect.Zr];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%code for getting the Z values using spline and Poly3 functions.

zout=zeros(length(quary.y),1);
for a=1:1:length(quary.y)
    x=Intersect_xy(a,1):.01:Intersect_xy(a,3);
    Spm.Y_1	= sect.Zl(a,1);
    Spm.Y_2	=sect.Zr(a,1); 
    Spm.X_1	=sect.Xlst(a,1);
    Spm.X_2	=sect.Xrst(a);
    Spm.angle_left=-Csl(a,1)+180;
    Spm.angle_right=-Csr(a,1);
    Spm.slope_left=tand(Spm.angle_left)*(Spm.X_2-Spm.X_1);
    Spm.slope_right=tand(-Spm.angle_right)*(Spm.X_2-Spm.X_1);
    Spm.D_1	=Spm.slope_left;
    Spm.D_2	=Spm.slope_right; 
    Spm.shift=0;
    param=spline_param_new(Spm);
    b=poly3n(param,quary.Xts(a));
    zout(a,1)=b;
    spl=poly3n(param,x);
    plot(x,spl)  %Ploting the spline curve for each point
    hold on
    
end
hold off
title('Spline Curves for each cross-section');
xlabel('Longitude- transformed coordinate');
ylabel("Height");
       
%Elevation for each grid point- will be same for transformed
%original-coordinate system(Z value remains same -either in original
%coordinate system or in transformed coordinate system

zout;
bnew=zout'
zcomb=[zout;sect.Zlr]';
zcombd=[quary.zs';sect.Zlr]';


%Plotting the original surface height and failure surface in original coordinate system
figure
plot3((xycomb(:,1))',(xycomb(:,2))',zcomb,'*',"Color",'r')
title('Failure surface Z values for each point-spline all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on;

plot3((xycomb(:,1))',(xycomb(:,2))',zcombd,'o')
title('Origional Surface Z values for each point all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on;

figure
plot3((xycomb(:,1))',(xycomb(:,2))',zcombd,'o');
hold on 
plot3((xycomb(:,1))',(xycomb(:,2))',zcomb,'*',"Color",'r')
hold off
title('Failure Surface and Origional Surface together usin all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on;

% Plotting the height in Original coordinate system
figure
plot3(quary.xs,quary.ys,bnew,'*');
title('Failure Surface Height-using Spline');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on;

plot3(quary.xs,quary.ys,quary.zs,'*'); %surface height 
title(' original surface Height');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on;

figure
plot3(quary.xs,quary.ys,bnew,'*'); %failure surface height-using spline
title('Failure Surface& original surface together');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
hold on ;
grid on
plot3(quary.xs,quary.ys,quary.zs,'o');
hold off;




% Failure surface ploting with original grid and original Coordinate system

Grid.xv = linspace(min((xycomb(:,1))), max((xycomb(:,1))),column );
Grid.yv = linspace(min((xycomb(:,2))), max((xycomb(:,2))),rows);
[Grid.Xo1,Grid.Yo1] = meshgrid(Grid.xv, Grid.yv);
Grid.Zo1 = griddata((xycomb(:,1))',(xycomb(:,2))',zcombd,Grid.Xo1,Grid.Yo1)

surf(Grid.Xo1, Grid.Yo1, Grid.Zo1);
title('original surface with all data-Original grid');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])

%Failure surface with all data 

Grid.xv = linspace(min((xycomb(:,1))), max((xycomb(:,1))),column );
Grid.yv = linspace(min((xycomb(:,2))), max((xycomb(:,2))),rows);
[Grid.Xf1,Grid.Yf1] = meshgrid(Grid.xv, Grid.yv);
Grid.Zf1 = griddata((xycomb(:,1))',(xycomb(:,2))',zcomb,Grid.Xf1,Grid.Yf1)


surf(Grid.Xf1, Grid.Yf1, Grid.Zf1);
title('Failure Surface with all data-Original grid');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])
%Failure surface and origina surface with all data 
surf(Grid.Xo1, Grid.Yo1, Grid.Zo1);
title('Failure Surface & original surface together with all data-Original grid');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');

hold on 
surf(Grid.Xf1, Grid.Yf1, Grid.Zf1);
hold off  


%Ploting the whole grid 
%[X_mesh,Y_mesh] = meshgrid(Grid.x, Grid.y)
figure
grid on

surf(X_mesh, Y_mesh, Z_mesh);
title('Actual grid');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');

surf(Grid.Xf1, Grid.Yf1, Grid.Zf1);
title('Failure Surface and Whole actual gride together');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])
hold on 
surf(X_mesh, Y_mesh, Z_mesh);
hold off



% SurfacePlotting in Original coordinate system with data within contour
% lines

xv = linspace(min(quary.xs), max(quary.xs), 100);
yv = linspace(min(quary.ys), max(quary.ys), 100);
[Xo,Yo] = meshgrid(xv, yv);
Zo = griddata(quary.xs,quary.ys,quary.zs,Xo,Yo);


xv = linspace(min(quary.xs), max(quary.xs), 100);
yv = linspace(min(quary.ys), max(quary.ys), 100);
[Xf,Yf] = meshgrid(xv, yv);
Zf = griddata(quary.xs,quary.ys,bnew,Xf,Yf);

%Failure surface and orinal surface together 
surf(Xo, Yo, Zo);
title('Failure Surface & original surface together');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');

hold on
surf(Xf, Yf, Zf);
hold off
set(gca, 'ZLim',[min(zout) max(quary.Zts)]);

%Plotting after combining the whole data 
%original surface with all data 

figure 
stem3((xycomb(:,1))',(xycomb(:,2))' , zcombd);
title('Original Surface stem');
grid on

xv = linspace(min((xycomb(:,1))), max((xycomb(:,1))), 100);
yv = linspace(min((xycomb(:,2))), max((xycomb(:,2))), 100);
[Xo1,Yo1] = meshgrid(xv, yv);
Zo1 = griddata((xycomb(:,1))',(xycomb(:,2))',zcombd,Xo1,Yo1)

figure
surf(Xo1, Yo1, Zo1);
title('original surface with all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])

%Failure surface with all data 
figure
stem3((xycomb(:,1))',(xycomb(:,2))' , zcomb);
title('Failure Surface stem');
grid on

xv = linspace(min((xycomb(:,1))), max((xycomb(:,1))), 100)
yv = linspace(min((xycomb(:,2))), max((xycomb(:,2))), 100)
[Xf1,Yf1] = meshgrid(xv, yv);
Zf1 = griddata((xycomb(:,1))',(xycomb(:,2))',zcomb,Xf1,Yf1)


figure
surf(Xf1, Yf1, Zf1);
title('Failure Surface with all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])

%Failure surface and orinal surface together 
figure
surf(Xo1, Yo1, Zo1);
title('Failure Surface & original surface together with all data');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');

hold on
surf(Xf1, Yf1, Zf1);
hold off


%Ploting the whole grid 

figure
%stem3(xi1, yi1, zi);
grid on

xv = linspace(min(total.x), max(total.x), 100);
yv = linspace(min(total.y), max(total.y), 100);
[X1,Y1] = meshgrid(xv, yv);
Zo = griddata(total.x,total.y,total.z,X1,Y1)

%%Failure Surface and Whole gride together
surf(Xf1, Yf1, Zf1);
hold on
surf(X1, Y1, Zo);
title('Failure Surface and Whole gride together');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height');
grid on
set(gca, 'ZLim',[min(zout) max(quary.Zts)])

%%%%%Saving the required grid data in tiff formate
if ask==1
    keep=keep(1:end-4)
    %Save Original surface    
    orig =append(keep,'_Original_grid.tif')
    %save(orig,'Z_mesh','-mat')
    [Lat_min,Lon_min] = utm2deg(min(total.x),min(total.y),utmzone(1,:));
    [Lat_max,Lon_max] = utm2deg(max(total.x),max(total.y),utmzone(1,:));   
    R = georasterref('RasterSize',size(Z_mesh),'LatitudeLimits',[Lat_min,Lat_max], ...
        'LongitudeLimits',[Lon_min,Lon_max],'ColumnsStartFrom','north');
    geotiffwrite(orig,Z_mesh,R)
    
    
    %Save(filename,variables) save Original failure surface
    orig_fail =append(keep,'_Original_failure_grid.tif');
    gvb=Grid.Zf1;
    %save(orig_fail,'gvb', '-ASCII');
    [Lat_min,Lon_min] = utm2deg(min(Grid.Xf1,[],"all"),min(Grid.Yf1,[],"all"),utmzone(1,:));
    [Lat_max,Lon_max] = utm2deg(max(Grid.Xf1,[],"all"),max(Grid.Yf1,[],"all"),utmzone(1,:));   
    R = georasterref('RasterSize',size(gvb),'LatitudeLimits',[Lat_min,Lat_max], ...
        'LongitudeLimits',[Lon_min,Lon_max]);
    geotiffwrite(orig_fail,gvb,R)
    
    %Save refined whole gride
    refine_whole =append(keep,'_Refined_whole_grid.tif') ;
    %save(refine_whole,'Zo', '-ASCII');
    [Lat_min,Lon_min] = utm2deg(min(X1,[],"all"),min(Y1,[],"all"),utmzone(1,:));
    [Lat_max,Lon_max] = utm2deg(max(X1,[],"all"),max(Y1,[],"all"),utmzone(1,:));    
    R = georasterref('RasterSize',size(Zo),'LatitudeLimits',[Lat_min,Lat_max], ...
        'LongitudeLimits',[Lon_min,Lon_max]);    
    geotiffwrite(refine_whole,Zo,R)    
    
    
    %Save refined original surface
    refine_origin =append(keep,'_Refined_original_grid.tif') ;
    %save(refine_origin,'Zo1', '-ASCII');
    [Lat_min,Lon_min] = utm2deg(min(Xo1,[],"all"),min(Yo1,[],"all"),utmzone(1,:));
    [Lat_max,Lon_max] = utm2deg(max(Xo1,[],"all"),max(Yo1,[],"all"),utmzone(1,:));    
    R = georasterref('RasterSize',size(Zo1),'LatitudeLimits',[Lat_min,Lat_max], ...
        'LongitudeLimits',[Lon_min,Lon_max]);     
    geotiffwrite(refine_origin,Zo1,R)    
    
    %Save refined failure surface
    refine_failure =append(keep,'_Refined_failure_grid.tif') ;
    %save(refine_failure,'Zf1', '-ASCII');
    [Lat_min,Lon_min] = utm2deg(min(Xf1,[],"all"),min(Yf1,[],"all"),utmzone(1,:));
    [Lat_max,Lon_max] = utm2deg(max(Xf1,[],"all"),max(Yf1,[],"all"),utmzone(1,:));    
    R = georasterref('RasterSize',size(Zf1),'LatitudeLimits',[Lat_min,Lat_max], ...
        'LongitudeLimits',[Lon_min,Lon_max]);      
    geotiffwrite(refine_failure,Zf1,R)
    
end

isequal(quary.xs,quary.Xts)

dZ=zeros(1,length(quary.zs))
for k=1:length(quary.zs)
    dZ(k)=quary.zs(k)-zout(k);
end

dZ
figure
stem3(quary.xs', quary.ys' , dZ);
title('Height Difference stem');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height_difference b/w failure & original surf');
grid on


xv = linspace(min(quary.xs), max(quary.xs), 100);
yv = linspace(min(quary.ys), max(quary.ys), 100);
[Xd,Yd] = meshgrid(xv, yv);
Zd = griddata(quary.xs,quary.ys,dZ,Xd,Yd)

figure
surf(Xd, Yd, Zd);
title('Height difference b/w failure & original surf');
xlabel('Longitude');
ylabel("Latitude");
zlabel('Height ');
grid on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calcuation of volume
Xcell_size
Ycell_size
szzx;
szzy;
Volume=0;
for vi=1:length(zout)
    Volume=Volume+Xcell_size*Ycell_size*(quary.zs(vi)-zout(vi));
    Volume;
end
name1 = 'Volume'; 
fprintf('%s = %10.2f  Cubic Meter.\n',name1,Volume)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calculation of the Area of polygon
Area =polyarea(contor.x,contor.y);
name2= 'Area';
fprintf('%s = %10.2f  Square Meter.\n',name2,Area)










