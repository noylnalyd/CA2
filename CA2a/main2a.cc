/*This is a skeleton code file for use with the Finite Element Method for Problems in Physics.
It uses the deal.II FEM library, dealii.org*/

//Include files
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "FEM2a.h"
#include "writeSolutions.h"

using namespace dealii;

//The main program, using the FEM class
int main (){
  try{
    deallog.depth_console (0);
		
    const int dimension = 2;
    
		//Define problem number: 1 or 2
		unsigned int problem = 2;
    FEM<dimension> problemObject(problem);
    
    //NOTE: This is where you define the number of elements in the mesh
    std::vector<unsigned int> num_of_elems(dimension);
    num_of_elems[0] = 5;
    num_of_elems[1] = 13; //For example, a 5 x 13 element mesh in 2D
    
    problemObject.generate_mesh(num_of_elems);
    problemObject.setup_system();
    problemObject.assemble_system();
    problemObject.solve();
    problemObject.output_results();
		if(problem == 2){
			std::cout << problemObject.l2norm_of_error() << std::endl;
		}

    //write solutions to h5 file
    char tag[21];
    sprintf(tag, "CA2a_Problem%d",problem);
		if(problem == 1){
	    writeSolutionsToFileCA2_1(problemObject.D, tag);
		}
		else if(problem == 2){
	    writeSolutionsToFileCA2_2(problemObject.D, problemObject.l2norm_of_error(), tag);
		}
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Exception on processing: " << std::endl
	      << exc.what() << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;

    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    std::cerr << "Unknown exception!" << std::endl
	      << "Aborting!" << std::endl
	      << "----------------------------------------------------"
	      << std::endl;
    return 1;
  }

  return 0;
}
