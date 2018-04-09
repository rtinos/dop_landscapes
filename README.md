# dop_landscapes
This repository contains a benchmark dynamic optimization problem generator based on the analysis of 6 types of fitness landscape transformations. 

************** Class dop.h **************

The class dop.h is used to produce the following types of Dynamic Optimization Problems: 

0) Stationary problem (no change)
1) DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)	
2) DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)	
3) DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a set of templates)		
4) DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)		
5) DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates)							
6) DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	

Reference: Tinos, R. & Yang, S. (2014). "Analysis of fitness landscape modifications in evolutionary dynamic optimization", Information Sciences, 282: 214-236.

The class dop.h has 4 public methods:
1) dop::dop(int l_p): constructor. The parameter <l_p> is the dimension of the problem.
2) dop::~dop(void): destructor. 
3) double dop::transform( int *x ,  int *x_new): causes a transformation in binary solution vector <x>. It returns the transformed vector <x_new> and the fitness difference <df>.
4) void dop::change( int change_type_p, double rho, double f_range ): changes the enviroment. The parameters are: change type <change_type_p>; change severity <rho>; parameter <f_range> used to define the size of the fitness transformation in DOP Type 3.

Observation: the class dop.h uses some auxiliaty functions defined in util_functions.cpp

Example: Using the class dop.h.

	...
	#include "dop.h"
	// <tau>: change period
	// <l>: dimension of the problem
	...
	// create DOP
	dop *DOP = new dop(l);		
	...
	f_range=population.max_fitness-population.mean_fitness;
	if (f_range<0.1)
		f_range=population.max_fitness;	
	...	
	do {
		gen = gen + 1; 								
		if ( (gen-gen_init)>tau ){
			gen_init=gen;
			// Change landscape
		  	DOP->change(change_type_temp, rho, f_range);	
		}
		...
		//  Fitness computation for individual <x>
		delta_f=DOP->transform(x, x_transf);
		fitness = fitnessFunction( x_transf ) + delta_f;	
		...			...
	} while ( stop_criterion==0 );
	...

