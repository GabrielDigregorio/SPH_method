#GEOM1
%********************************************
% PLAYGROUND : Geometry example
%
% #GEOM
% #domain
%   lx, ly, lz ==> Lower coordinate of the domain
%   ux, uy, uz ==> Upper coordinate of the domain
%
% geometry : #brick, #cylin or #spher
% %param
%   c ==> 0=free, 1=moving, 2=fixed
%   s ==> spacing between particles
%   r ==> % randomness in particle's position
% %coord
%   x ==> x coordinate of the center of mass
%   y ==> y coordinate of the center of mass
%   z ==> z coordinate of the center of mass
% %dimen
%   L ==> Length of the geometry
%   W ==> Width  of the geometry
%   H ==> Height of the geometry
%   d ==> Small diameter of the geometry
%   D ==> Large diameter of the geometry
% %inclin
%   tetax ==> Angle of rotation around the x axis (in degre)
%   tetay ==> Angle of rotation around the y axis (in degre)
%   tetaz ==> Angle of rotation around the z axis (in degre)
% %movpara
%   movingDirection X, Y & Z =
%   speedLaw = dictate the behaviour of moving boundaries
%********************************************

#domsz
	ux = 0.16
	uy = 0.26
	uz = 0.31
	lx = -0.16
	ly = -0.31
	lz = -0.04
% Soup and container
#bathy
	batFile=Playgrounds/soup.txt
	file=1
	s=0.006
	r=0
	numberGroundParticles=1
	height0=0
	hFreeSurface=0.01

% Glass
% End wall
   #brick
       %param
           c=2
           s=0.006
           r=0
       %coord
           x=0
           y=-0.306
           z=0.125
       %dimen
           L=0.112
           W=0.012
           H=0.09
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=3
           angleLaw=2
           charactTime=3
           movingDirectionX=1
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=-0.10
           rotationCenterZ=0.08
           amplitude=-1.7
%Front wall
   #brick
       %param
           c=2
           s=0.006
           r=0
       %coord
           x=0
           y=-0.094
           z=0.125
       %dimen
           L=0.112
           W=0.012
           H=0.09
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=3
           angleLaw=2
           charactTime=3
           movingDirectionX=1
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=-0.10
           rotationCenterZ=0.08
           amplitude=-1.7
%Side wall
   #brick
       %param
           c=2
           s=0.006
           r=0
       %coord
           x=-0.056
           y=-0.2
           z=0.125
       %dimen
           L=0.012
           W=0.2
           H=0.09
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=3
           angleLaw=2
           charactTime=3
           movingDirectionX=1
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=-0.10
           rotationCenterZ=0.08
           amplitude=-1.7
%Side wall
   #brick
       %param
           c=2
           s=0.006
           r=0
       %coord
           x=0.056
           y=-0.2
           z=0.125
       %dimen
           L=0.012
           W=0.2
           H=0.09
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=3
           angleLaw=2
           charactTime=3
           movingDirectionX=1
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=-0.10
           rotationCenterZ=0.08
           amplitude=-1.7
%Bottom wall
   #brick
       %param
           c=2
           s=0.006
           r=0
       %coord
           x=0
           y=-0.2
           z=0.074
       %dimen
           L=0.112
           W=0.212
           H=0.012
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=3
           angleLaw=2
           charactTime=3
           movingDirectionX=1
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=-0.10
           rotationCenterZ=0.08
           amplitude=-1.7
%Fluid
   #brick
       %param
           c=0
           s=0.006
           r=0
       %coord
           x=0
           y=-0.2
           z=0.1206
       %dimen
           L=0.1
           W=0.2
           H=0.08
       %inclin
           tetax=0.0
           tetay=0.0
           tetaz=0.0
       %movpara
           speedLaw=0
           angleLaw=0
           charactTime=0
           movingDirectionX=0
           movingDirectionY=0
           movingDirectionZ=0
           rotationCenterX=0
           rotationCenterY=0
           rotationCenterZ=0
           amplitude=0
%Floor
% Soup and container
#bathy
	batFile=Playgrounds/soup2.txt
	file=1
	s=0.006
	r=0
	numberGroundParticles=1
	height0=-0.006
	hFreeSurface=-0.01
#END_G
