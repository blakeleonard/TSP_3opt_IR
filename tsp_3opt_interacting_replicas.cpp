// An Interacting Replica approach incorporating 3-opt applied to the Traveling Salesman Problem.

// Blake Leonard 2011

// Physics Department
// Washington University in St. Louis


#include <iostream>
using namespace::std;

#include <fstream>
using namespace::std;

#include <cmath>
using namespace::std;

#include <ctime>
using namespace::std;

#include <algorithm>
using namespace::std;

#include <vector>
using namespace::std;

#include <functional>
using namespace::std;


int main()
{

	// Variable declaration and initial conditions	

	char *inname = "eil51.txt";         			// file of city positions

	const int N_max = 51;                			// number of cities
	const int L_optimum = 7542;

	const int R = 20;                    			// number of replicas
	const int m = 5000000;					// number of global runs

	const int optstep1 = 1000;              		// number of initial 3-opt steps
	const int optstep2 = 1000;           			// number of 3-opt steps after info step begins
	const int infostep = 1;              			// number of info steps
	const int rand_moves = 74; 
	const float allowable_L_change = 1.008;

	int L_best = 1000000000;
	int best_rep;
	int infoflag = 0;
	int cityposition [ 3*N_max ];
	int xplot [ N_max + 1 ];
	int info_moves;

	int yplot [ N_max + 1 ];
	float citydistance [ N_max + 1 ] [ N_max + 1 ];
	int mstep;
	int rand_city;
	int rand_city_vect [ N_max + 1 ];
	int file_counter = 0;
	float timer;
	float start_time;
	int order [ R + 1 ] [ N_max + 1 ];
	int L [ R + 1 ];
	int city;
	int nextcity;
	int shortestdistance;
	int flag;
	int choice;
	int temp_perm [ N_max + 1 ];

	int optstep = optstep1;
	int standard_city;
	int fstep;
	int sstep;
	int nstep;
	int city_location [ R + 1 ] [ N_max + 1 ];
	int old_order [ R + 1 ] [ N_max + 1 ];
	int average_position [ N_max + 1 ];
	int L_new [ R + 1 ];

	int randkvect [ 4 ];
	int temp_store;
	int selected_cities [ 7 ];
	int orig_edge1;
	int orig_edge2;
	int orig_edge3;
	int orig_total;
	int temp_edge1;
	int temp_edge2;
	int temp_edge3;

	int new_order [ N_max + 1 ];
	int counter1;
	int counter2;
	int counter3;
	int counter4;
	int temp_order1 [ N_max + 1 ];
	int temp_order2 [ N_max + 1 ];
	int temp_order3 [ N_max + 1 ];
	int temp_order4 [ N_max + 1 ];
	int lowest_total;
	int opt_flag;

	int istep;
	int jstep;
	int pstep;
	int kstep;
	int ostep;
	int count_max;
	int best_nn;
	int average_flag [ R + 1 ];
	int nn_counter [ N_max + 1];
	int average_counter;
	int rstep;


	// Define a template class vector of int

	typedef vector<int > IntVector ;


	//Define an iterator for template class vector of strings

	typedef IntVector::iterator IntVectorIt ;

	IntVector temp_total( 8 ) ;

	IntVectorIt start, end;

	temp_total [ 0 ] = 10000000;


	srand((unsigned)time(0));

    
	// Open file of city positions

	float i;

	ifstream infile(inname);

	if (!infile) 
	{
        
		cout << "There was a problem opening file "
             
			 << inname
             
			 << " for reading."
             
			 << endl;
        
	}
    
	while (infile >> i) 
	{
		file_counter = file_counter + 1;

		cityposition [ file_counter ] = i;
      
	}

	int N = file_counter / 3;

	for ( istep = 0; istep < N; istep ++ )
	{

		xplot [ istep + 1 ] = cityposition [ 3 * istep + 2 ];

		yplot [ istep + 1 ] = cityposition [ 3 * istep + 3 ]; 

	}


	// Create N * N matrix of city pairwise distances


	for ( istep = 1; istep <= N; istep ++ )
	{    

		for ( jstep = 1; jstep <= N; jstep ++ )
		{    
        
			citydistance [istep] [jstep] = sqrt ( ( xplot [istep] - xplot [jstep] ) * ( xplot [istep] - xplot [jstep] ) + ( yplot [istep] - yplot [jstep] ) * ( yplot [istep] - yplot [jstep] ) );
      
			citydistance [istep] [jstep] = floor ( 0.5 + citydistance [istep] [jstep] );	

		}
	}


	start_time = clock() * 0.001;


	for ( rstep = 1; rstep <= R; rstep ++ )          // Seed with random starting order
	{

		// Nearest Neighbor Algorithm to seed k-opt

		city = 1 + rand() % N;

		L [ rstep ] = 0;

		for ( jstep = 1; jstep <= N; jstep ++ )
		{

			shortestdistance = 10000000;

			for (istep = 1; istep <= N; istep ++ )
			{

				if ( ( citydistance[city][istep] < shortestdistance ) && ( city != istep ) ) 
				{  
 
					flag = 0;
            
					if (jstep != 1)
					{    
					
						for ( kstep = 1; kstep <= (jstep - 1); kstep ++ )
						{

							if (istep != order[rstep][kstep])
                    
								flag = flag + 1;
                    
						}
            
					}   
                
					if ( flag == (jstep - 1) || jstep == 1)
					{

						shortestdistance = citydistance[city][istep];
        
						nextcity = istep;

					}
        
				}
			
			}

  
			order[rstep][jstep] = city;

			city = nextcity;

		}
	

		// Determine Route Length
	
		for ( jstep = 1; jstep <= N; jstep ++ )
		{

			if ( jstep != N )

				L [ rstep ] = L [ rstep ] + citydistance [ order [ rstep ] [ jstep ] ] [ order [ rstep ] [ jstep + 1 ] ];

			else

				L [ rstep ] = L [ rstep ] + citydistance [ order [ rstep ] [ N ] ] [ order [ rstep ] [ 1 ] ];

		}
   
	}


	for ( kstep = 1; kstep <= m; kstep ++ )         // Global steps
	{   
    
		if ( infoflag == 1 )
        
			optstep = optstep2;


		// 3 - opt step
        
    
		for ( ostep = 1; ostep <= optstep; ostep ++ )         
		{       
    
			for ( rstep = 1; rstep <= R; rstep ++ )              // Cycle over replicas
			{

				// Pick random edges

				randkvect [ 1 ] = 1;

				randkvect [ 2 ] = 2;

				randkvect [ 3 ] = 3;


				while( ( randkvect [ 3 ] == randkvect [ 2 ] + 1 ) || ( randkvect [ 2 ] == randkvect [ 1 ] + 1 )  || ( randkvect [ 3 ] == N ) )       // Don't Allow adjacent edges, randkvect(3) = N arg more restrictive than neccessary, needs debugging
				{
                
					randkvect [ 1 ] = 1 + rand() % N;

					randkvect [ 2 ] = 1 + rand() % N;

					randkvect [ 3 ] = 1 + rand() % N;


					while ( randkvect [ 1 ] == randkvect [ 2 ] || randkvect [ 2 ] == randkvect [ 3 ] )

						randkvect [ 2 ] = 1 + rand() % N;                                  // Need to change this to run faster


					while ( randkvect [ 1 ] == randkvect [ 3 ] || randkvect [ 2 ] == randkvect [ 3 ] )

						randkvect [ 3 ] = 1 + rand() % N;


					if ( randkvect [ 3 ] < randkvect [ 2 ] )
					{

						temp_store = randkvect [ 2 ];

						randkvect [ 2 ] = randkvect [ 3 ];

						randkvect [ 3 ] = temp_store;

					}


					if ( randkvect [ 2 ] < randkvect [ 1 ] )
					{
                 
						temp_store = randkvect [ 1 ];

						randkvect [ 1 ] = randkvect [ 2 ];

						randkvect [ 2 ] = temp_store;

					}


					if ( randkvect [ 3 ] < randkvect [ 2 ] )
					{

						temp_store = randkvect [ 2 ];

						randkvect [ 2 ] = randkvect [ 3 ];

						randkvect [ 3 ] = temp_store;

					}

				}
				

				for ( istep = 1; istep <= 3; istep ++ )

					selected_cities [ istep ] = order [ rstep ] [ randkvect [ istep ] ];
  

				for ( istep = 1; istep <= 3; istep ++ )
				{

					if ( randkvect [ istep ] != N )

						selected_cities [ istep+3 ] = order [ rstep ] [ randkvect [ istep ] + 1 ];

					else

						selected_cities [ istep+3 ] = order [ rstep ] [ 1 ];
 
				}    


				// 8 possible combinations in 3 -opt

				orig_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 4 ] ];

				orig_edge2 = citydistance [ selected_cities [ 2 ] ] [ selected_cities [ 5 ] ];

				orig_edge3 = citydistance [ selected_cities [ 3 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 0 ] = orig_total = orig_edge1 + orig_edge2 + orig_edge3;


				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 4 ] ];

				temp_edge2 = citydistance [ selected_cities [ 2 ] ] [ selected_cities [ 3 ] ];

				temp_edge3 = citydistance [ selected_cities [ 5 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 1 ] = temp_edge1 + temp_edge2 + temp_edge3;


				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 2 ] ];

				temp_edge2 = citydistance [ selected_cities [ 4 ] ] [ selected_cities [ 5 ] ];

				temp_edge3 = citydistance [ selected_cities [ 3 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 2 ] = temp_edge1 + temp_edge2 + temp_edge3;

			
				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 2 ] ];

				temp_edge2 = citydistance [ selected_cities [ 4 ] ] [ selected_cities [ 3 ] ];

				temp_edge3 = citydistance [ selected_cities [ 5 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 3 ] = temp_edge1 + temp_edge2 + temp_edge3;


				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 5 ] ];

				temp_edge2 = citydistance [ selected_cities [ 3 ] ] [ selected_cities [ 4 ] ];

				temp_edge3 = citydistance [ selected_cities [ 2 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 4 ] = temp_edge1 + temp_edge2 + temp_edge3;


				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 5 ] ];

				temp_edge2 = citydistance [ selected_cities [ 3 ] ] [ selected_cities [ 2 ] ];

				temp_edge3 = citydistance [ selected_cities [ 4 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 5 ] = temp_edge1 + temp_edge2 + temp_edge3;

		
				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 3 ] ];

				temp_edge2 = citydistance [ selected_cities [ 5 ] ] [ selected_cities [ 4 ] ];

				temp_edge3 = citydistance [ selected_cities [ 2 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 6 ] = temp_edge1 + temp_edge2 + temp_edge3;


				temp_edge1 = citydistance [ selected_cities [ 1 ] ] [ selected_cities [ 3 ] ];

				temp_edge2 = citydistance [ selected_cities [ 5 ] ] [ selected_cities [ 2 ] ];

				temp_edge3 = citydistance [ selected_cities [ 4 ] ] [ selected_cities [ 6 ] ];

				temp_total [ 7 ] = temp_edge1 + temp_edge2 + temp_edge3;


				start = temp_total.begin() ;  

				end = temp_total.end() ;       

				lowest_total = *min_element( start, end );

				opt_flag = 0;


				if ( temp_total [ 1 ] == lowest_total && opt_flag == 0 )
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 1 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 4 ];


					counter1 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 4 & 2 in original order
					{
						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order 

						new_order[ istep ] = temp_order1 [ istep - 2 ];

				
					new_order [ counter1+3 ] = selected_cities [ 2 ];

					new_order [ counter1+4 ] = selected_cities [ 3 ];


					counter2 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{
						
						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )        // read into new order backwards

						new_order[ istep ] = temp_order2 [ counter1 + counter2 + 5 - istep ];

				
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 5 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{
						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];

				
					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];

                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{
                    
						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
             
					}


					opt_flag = 1;
			

				}


				if ( temp_total [ 2 ] == lowest_total && opt_flag == 0 )
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 2 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 2 ];


					counter1 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 4 & 2 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}


					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order1 [ counter1 - istep + 3 ];


					new_order [ counter1+3 ] = selected_cities [ 4 ];

					new_order [ counter1+4 ] = selected_cities [ 5 ];


					counter2 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )                               // read into new order

						new_order[ istep ] = temp_order2 [ istep - counter1 - 4 ];

			
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 3 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;

					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];


					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];

                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                				
				}

			
				if ( temp_total [ 3 ] == lowest_total  && opt_flag == 0 )
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 3 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 2 ];


					counter1 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 4 & 2 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order1 [ counter1 - istep + 3 ];

			
					new_order [ counter1+3 ] = selected_cities [ 4 ];

					new_order [ counter1+4 ] = selected_cities [ 3 ];


					counter2 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order2 [ counter1 + counter2 + 5 - istep ];

				
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 5 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];

	
					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];

                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                				
				}

			
				if ( temp_total [ 4 ] == lowest_total && opt_flag == 0)
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 4 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 5 ];


					counter1 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order 

						new_order[ istep ] = temp_order1 [ istep - 2 ];

				
					new_order [ counter1+3 ] = selected_cities [ 3 ];

					new_order [ counter1+4 ] = selected_cities [ 4 ];


					counter2 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 2 & 4 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )                               // read into new order

						new_order[ istep ] = temp_order2 [ istep - counter1 - 4 ];

			
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 2 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];

				
					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];

                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                			
				}

	
				if ( temp_total [ 5 ] == lowest_total && opt_flag == 0)
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 5 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 5 ];


					counter1 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order 

						new_order[ istep ] = temp_order1 [ istep - 2 ];

			
					new_order [ counter1+3 ] = selected_cities [ 3 ];

					new_order [ counter1+4 ] = selected_cities [ 2 ];


					counter2 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 2 & 4 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order2 [ counter1 + counter2 + 5 - istep  ];


					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 4 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];

				
					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];
          
                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                				
				}

	
				if ( temp_total [ 6 ] == lowest_total && opt_flag == 0)
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 6 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 3 ];


					counter1 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order1 [ counter1 + 3 - istep ];

		
					new_order [ counter1+3 ] = selected_cities [ 5 ];

					new_order [ counter1+4 ] = selected_cities [ 4 ];


					counter2 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 2 & 4 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )                               // read into new order

						new_order[ istep ] = temp_order2 [ istep - counter1 - 4 ];

			
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 2 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];
				 
				
					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];
             
                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                				
				}


				if ( temp_total [ 7 ] == lowest_total && opt_flag == 0)
				{    

					L [ rstep ] = L [ rstep ] - ( orig_total - temp_total [ 7 ] );
                
					new_order [ 1 ] = selected_cities [ 1 ];

					new_order [ 2 ] = selected_cities [ 3 ];


					counter1 = 0;

					for ( istep = (randkvect [ 2 ] + 2); istep <= ( randkvect [ 3 ] - 1); istep ++ )    // Capture order of cities between city 5 & 3 in original order
					{

						counter1 = counter1 + 1;

						temp_order1 [ counter1 ] = order [ rstep ] [ istep ];

					}

					for ( istep = 3; istep <= ( counter1 + 2 ); istep ++ )                               // read into new order backwards

						new_order[ istep ] = temp_order1 [ counter1 + 3 - istep ];

			
					new_order [ counter1+3 ] = selected_cities [ 5 ];

					new_order [ counter1+4 ] = selected_cities [ 2 ];


					counter2 = 0;

					for ( istep = (randkvect [ 1 ] + 2); istep <= ( randkvect [ 2 ] - 1); istep ++ )    // Capture order of cities between city 2 & 4 in original order
					{

						counter2 = counter2 + 1;

						temp_order2 [ counter2 ] = order [ rstep ] [ istep ];

					}

					for ( istep = (counter1+5); istep <= ( counter1 + counter2 + 4 ); istep ++ )            // read into new order backwards

						new_order[ istep ] = temp_order2 [ counter1 + counter2 + 5 - istep ];

		
					new_order [ counter1 + counter2 + 5 ] = selected_cities [ 4 ];

					new_order [ counter1 + counter2 + 6 ] = selected_cities [ 6 ];


					counter3 = 0;

					for ( istep = (randkvect [ 3 ] + 2); istep <= N; istep ++ )                     // Capture order of cities between city 6 and end in original order
					{

						counter3= counter3 + 1;

						temp_order3 [ counter3 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + 7); istep <= (counter1 + counter2 + counter3 + 6); istep ++ )      // read in

						new_order [ istep ] = temp_order3 [ istep - counter1 - counter2 - 6 ];


					counter4 = 0;


					for ( istep = 1; istep <= (randkvect [ 1 ] - 1); istep ++ )           // Capture order of cities between original 1st city and city 1
					{

						counter4 = counter4 + 1;

						temp_order4 [ counter4 ] = order [ rstep ] [ istep ];

					}


					for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )         // read in

						new_order [ istep ] = temp_order4 [ istep - counter1 - counter2 - counter3 - 6 ];


					for ( istep = 1; istep <= N; istep ++ )

						order [ rstep ] [ istep ] = new_order [ istep ];

                
					if ( L [ rstep ] < L_best )      // If current route is shorter than archived best, store
					{

						L_best = L [ rstep ];

						best_rep = rstep;

						timer = clock() * 0.001;

						if  ( L_best == L_optimum )

							cout << L_best << endl << timer - start_time << endl;
               
					}


					opt_flag = 1;

                				
				}


			}		


		}


		// Info step

		for ( istep = 1; istep <= infostep; istep ++ )
		{    
        
			standard_city = 1;
        
			while (standard_city == 1)  
			{ 
  
				standard_city = 1 + rand() % ( N - 1 );
             
			} 


			for ( nstep = 1; nstep <= N; nstep ++ )

				nn_counter [ nstep ] = 0;

		
			for ( rstep = 1; rstep <= R; rstep ++ )            
			{  
            
				// Record locations of cities in order in each replica
                     
				for ( fstep = 1; fstep <= N; fstep ++ )            
				{
   
					for ( sstep = 1; sstep <= N; sstep ++ )
					{   
  
						if ( order [ rstep ] [ fstep ] == sstep )
						{ 
    
							city_location [ rstep ] [ sstep ] = fstep;
                        
							break;
                        
						}                    
                  
					}

				}
             
            
				// Cycle order array so it begins with standard_city with 2nd
				// element the nearest neighbor
         
				for ( nstep = 1; nstep <= N; nstep ++ )        
                
					old_order [ rstep ] [ nstep ] = order [ rstep ] [ nstep ];
 
            
				if ( ( citydistance [ standard_city ] [ order [ rstep ] [ city_location [ rstep ] [ standard_city ] ] + 1 ] ) <= ( citydistance [ standard_city] [ order [ rstep ] [ city_location [ rstep ] [ standard_city ] ] - 1 ]  ) )
				{    

					for ( nstep = 1; nstep <= ( N + 1 - city_location [ rstep ] [ standard_city ] ); nstep ++ )

						order [ rstep ] [ nstep ] = old_order [ rstep ] [ city_location [ rstep ] [ standard_city ] - 1 + nstep ];

                
					for ( nstep = (N + 2 - city_location [ rstep ] [ standard_city ] ); nstep <= N; nstep ++ )

						order [ rstep ] [ nstep ] = old_order [ rstep ] [ city_location [ rstep ] [ standard_city ] - N - 1 + nstep ];

				}
				else
				{    
                
					for ( nstep = 1; nstep <= city_location [ rstep ] [ standard_city ]; nstep ++ )
                    
						order [ rstep ] [ nstep ] = old_order [ rstep ] [ city_location [ rstep ] [ standard_city ] - nstep + 1 ];
                    
                
					for ( nstep = ( city_location [ rstep ] [ standard_city] + 1 ); nstep <= N; nstep ++ )
                    
						order [ rstep ] [ nstep ] = old_order [ rstep ] [ city_location [ rstep ] [ standard_city ] - nstep + 1 + N ];
                    
				}            
            
			
				// Record locations of cities in order in each replica

				for ( fstep = 1; fstep <= N; fstep ++ )            
				{   

					for ( sstep = 1; sstep <= N; sstep ++ )
					{  
					
						if ( order [ rstep ] [ fstep ] == sstep )
						{    

							city_location [ rstep ] [ sstep ] = fstep;
                        
							break;
                        
						}                    
                  
					}
                
				}

				nn_counter [ order [ rstep ] [ 2 ] ] = nn_counter [ order [ rstep ] [ 2 ] ] + 1;           

			}


			// Determine city which is most common nearest neighbor to standard city

			count_max = 0;

			for ( nstep = 1; nstep <= N; nstep ++ )
			{

				if ( nn_counter [ nstep ] > count_max )
				{

					count_max = nn_counter [ nstep ];

					best_nn = nstep;

				}
			
			}


			// Flag Replicas which contain the most common nearest neighbor city in the correct position

			for ( rstep = 1; rstep <= R; rstep ++ )
			{

				if ( order [ rstep ] [ 2 ] == best_nn )

					average_flag [ rstep ] = 1;
			
				else

					average_flag [ rstep ] = 0;

			}


			// Compute average position over replicas for each city  
		
			info_moves = 0;
        
			for ( nstep = 1; nstep <= N; nstep ++ )
			{
			
				average_position [ nstep ] = 0;

				average_counter = 0;
            
				for ( rstep = 1; rstep <= R; rstep ++ )
				{
				
					if ( average_flag [ rstep ] == 1)
					{

						average_position [ nstep ] = average_position [ nstep ] + city_location [ rstep ] [ nstep ];

						average_counter = average_counter + 1;
                
					}

				}
            
				average_position [ nstep ] = average_position [ nstep ] / average_counter;
            
				average_position [ nstep ] = floor( average_position [ nstep ] + 0.5 );

			}
			
		
			// Make sequence of random moves placing cities in their average
			// positions
      
			// Make random permutation vector for cities

			for ( pstep = 1; pstep <= N; pstep ++ )

				temp_perm [ pstep ] = pstep;


			for ( pstep = 0; pstep < rand_moves; pstep ++ )
			{

				choice = 1 + rand() % ( N - pstep );

				rand_city_vect [ pstep + 1 ] = temp_perm [ choice ];

				for ( ostep = choice; ostep <= (N - 1 - pstep); ostep ++ )

					temp_perm [ ostep ] = temp_perm [ ostep + 1 ];

			}

        
			for ( mstep = 1; mstep <= rand_moves; mstep ++ )
			{ 
                   
				rand_city = rand_city_vect [ mstep ];
            
            
				// Make moves
            
				for ( rstep = 1; rstep <= R; rstep ++ )
				{   
                
					if ( average_flag [ rstep ] == 1 )
					{

						// Record locations of cities in order in each replica
                     
						for ( fstep = 1; fstep <= N; fstep ++ )            
						{

							for ( sstep = 1; sstep <= N; sstep ++ )
							{

								if ( order [ rstep ] [ fstep ] == sstep )
								{

									city_location [ rstep ] [ sstep ] = fstep;

									break;

								}                    

							}

						}

					 
						for ( nstep = 1; nstep <= N; nstep ++ ) 
                
							old_order [ rstep ] [ nstep ] = order [ rstep ] [ nstep ];

                
						if ( average_position [ rand_city ] < city_location [ rstep ] [ rand_city ] )
						{  
  
							order [ rstep ] [ average_position [ rand_city ] ] = rand_city;
                    
							for ( nstep = ( average_position [ rand_city ] + 1 ); nstep <= city_location [ rstep ] [ rand_city ]; nstep ++ )
                        
								order [ rstep ] [ nstep ] = old_order [ rstep ] [ nstep - 1 ];
                        
						}   
						else if ( average_position [ rand_city ] > city_location [ rstep ] [ rand_city ] )
						{   

							order [ rstep ] [ average_position [ rand_city ] ] = rand_city;
                    
							for ( nstep = city_location [ rstep ] [ rand_city ]; nstep <= ( average_position [ rand_city ] - 1 ); nstep ++ )
                        
								order [ rstep ] [ nstep ] = old_order [ rstep ] [ nstep + 1 ];
                
						}       
                   
                                        				
						//  Compute new path lengths
                
						L_new [ rstep ] = 0;
                
						for ( nstep = 1; nstep <= ( N - 1 ); nstep ++ )
                    
							L_new [ rstep ] = L_new [ rstep ] + citydistance [ order [ rstep ] [ nstep] ] [ order [ rstep ] [ nstep + 1 ] ];
                    
                
						L_new [ rstep ] = L_new [ rstep ] + citydistance [ order [ rstep ] [ N ] ] [ order [ rstep ] [ 1 ] ];
               
        
						// If new length is too long, switch back
        
						if ( L_new [ rstep ] > allowable_L_change * L [ rstep ] )
						{

							for ( nstep = 1; nstep <= N; nstep ++ )

								order [ rstep ] [ nstep ] = old_order [ rstep ] [ nstep];

						}
						else
						{

							L [ rstep ] = L_new [ rstep ];
						}
                
                
						if ( L [ rstep ] < L_best)      // If current route is shorter than archived best, store
						{

							L_best = L [ rstep ];

							best_rep = rstep;

							timer = clock() * 0.001;

							if  ( L_best == L_optimum )

								cout << L_best << " info step " << endl << timer - start_time << endl;

						}   
          
					}

				}

			}

		}

		infoflag = 1;

	}

	return 0;

}


