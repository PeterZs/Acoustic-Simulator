#Experiment_1.geo

algebraic3d

solid muffler=cylinder(0,0,0;0,6,0;1)
			  and plane(0,0,0;0,-1,0)
			  and plane(0,6,0;0,1,0);

solid cube=orthobrick(-3,3,-3;3,10,3);

solid shell=cylinder(0,0,0;0,6,0;1.05)
			  and plane(0,0,0;0,-1,0)
			  and plane(0,6,0;0,1,0);

solid room=cube and not shell;

solid model=room or muffler;

tlo model;
