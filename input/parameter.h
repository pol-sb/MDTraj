
	!	Initialization parameters
	!---------------------------
	integer::natoms=32 !Number of particles
	!double precision::density
	double precision::L=10.0d0 !Lomgitude of the size box
	integer::dimension=3 !Dimension of the system.
	
	integer::structure=2 !For SC:1, for fcc:2, for diamond 3.  
	
	double precision::temp= 300.0d0  !Temperature (K)
	integer::vel_opt=1 ! if 1, implemented bimodal. Else, set up to 0
	
	!Simulation set up
	!-----------------
	integer::ntimes=100 !Steps of the simulation
	double precision:: dt= 0.005 !Time step (ps)
	
	
	!Force parameters
	!----------------
	double precision::rc= 5.0d0 !Cut-off
