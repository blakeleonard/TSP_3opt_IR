// Standard 3-opt algorithm applied to the Traveling Salesman Problem (for comparison with Interacting Replica method).

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

	char *inname = "eil51.txt";        			// file of city positions

	const int N_max = 51;                    		// max number of cities
	const int L_optimum = 7542;

	const int R = 20;                    			// number of replicas
	const int optstep = 10000000;           		// number of 3-opt steps

	int L_best = 1000000000;
	int best_rep;
	int cityposition [ 3*N_max ];
	int xplot [ N_max + 1 ];

	int yplot [ N_max + 1 ];
	float citydistance [ N_max + 1 ] [ N_max + 1 ];
	int file_counter = 0;
	float timer;
	float start_time;
	int order [ R + 1 ] [ N_max + 1 ];
	int L [ R + 1 ];
	int city;
	int nextcity;
	int shortestdistance;
	int flag;

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
	int kstep;
	int ostep;
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
        
		cout << "There was a problem opening file " << inname << " for reading." << endl;
        
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


	// 3 - opt step
        
	for ( ostep = 1; ostep <= optstep; ostep ++ )         
	{       
    
		for ( rstep = 1; rstep <= R; rstep ++ )              // Cycle over replicas
		{

			// Pick random edges

			randkvect [ 1 ] = 1;

			randkvect [ 2 ] = 2;

			randkvect [ 3 ] = 3;


			while( ( randkvect [ 3 ] == randkvect [ 2 ] + 1 ) || ( randkvect [ 2 ] == randkvect [ 1 ] + 1 )  || ( randkvect [ 3 ] == N ) )		// Don't Allow adjacent edges, randkvect(3) = N arg more restrictive than neccessary, needs debugging
			{
                
				randkvect [ 1 ] = 1 + rand() % N;

				randkvect [ 2 ] = 1 + rand() % N;

				randkvect [ 3 ] = 1 + rand() % N;


				while ( randkvect [ 1 ] == randkvect [ 2 ] || randkvect [ 2 ] == randkvect [ 3 ] )

					randkvect [ 2 ] = 1 + rand() % N;                                  			// Need to change this to run faster


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


			// 8 possible combinations in 3 - opt

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

				for ( istep = (counter1 + counter2 + counter3 + 7); istep <= (counter1 + counter2 + counter3 + counter4 + 6); istep ++ )     // read in

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

					if  ( L_best == L_optimum )

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

	return 0;

}


