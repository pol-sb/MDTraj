
	!	Initialization parameters
	!---------------------------
	integer::nc=4 !Number of particles per side
	double precision::density=0.4d0 !Lomgitude of the size box
	integer::dimension=3 !Dimension of the system.

	integer::structure=1 !For SC:1, for fcc:2, for diamond 3.

	double precision::temp= 300.0d0  !Temperature (K)
	integer::vel_opt=1 ! if 1, implemented bimodal. Else, set up to 0

	!Simulation set up
	!-----------------
	integer::ntimes=100 !Steps of the simulation
	integer::everyt=10
	double precision:: dt= 0.005 !Time step (ps)

	integer::thermo=1 !For no thermostat:0, for active thermostat:1

	!Force parameters
	!----------------
	double precision::rc= 0.05d0 !Cut-off

	!Lennard-Jones parameters
	!----------------
	double precision:: epsilon = 10.22, sigma = 2.556
