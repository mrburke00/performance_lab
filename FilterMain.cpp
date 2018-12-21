#include <stdio.h>
#include "cs1300bmp.h"
#include <iostream>
#include <fstream>
#include "Filter.h"
#include <omp.h>

using namespace std;

#include "rtdsc.h"

//
// Forward declare the functions
//
Filter * readFilter(string filename);
double applyFilter(Filter *filter, cs1300bmp *input, cs1300bmp *output);

int
main(int argc, char **argv)
{

  if ( argc < 2) {
    fprintf(stderr,"Usage: %s filter inputfile1 inputfile2 .... \n", argv[0]);
  }

  //
  // Convert to C++ strings to simplify manipulation
  //
  string filtername = argv[1];

  //
  // remove any ".filter" in the filtername
  //
  string filterOutputName = filtername;
  string::size_type loc = filterOutputName.find(".filter");
  if (loc != string::npos) {
    //
    // Remove the ".filter" name, which should occur on all the provided filters
    //
    filterOutputName = filtername.substr(0, loc);
  }

  Filter *filter = readFilter(filtername);

  double sum = 0.0;
  int samples = 0;

  for (int inNum = 2; inNum < argc; inNum++) {
    string inputFilename = argv[inNum];
    string outputFilename = "filtered-" + filterOutputName + "-" + inputFilename;
    struct cs1300bmp *input = new struct cs1300bmp;
    struct cs1300bmp *output = new struct cs1300bmp;
    int ok = cs1300bmp_readfile( (char *) inputFilename.c_str(), input);

    if ( ok ) {
      double sample = applyFilter(filter, input, output);
      sum += sample;
      samples++;
      cs1300bmp_writefile((char *) outputFilename.c_str(), output);
    }
    delete input;
    delete output;
  }
  fprintf(stdout, "Average cycles per sample is %f\n", sum / samples);

}

struct Filter *
readFilter(string filename)
{
  ifstream input(filename.c_str());

  if ( ! input.bad() ) {
    int size = 0;
    input >> size;
    Filter *filter = new Filter(size);
    int div;
    input >> div;
    filter -> setDivisor(div);
    for (int i=0; i < size; i++) {
      for (int j=0; j < size; j++) {
	int value;
	input >> value;
	filter -> set(i,j,value);
      }
    }
    return filter;
  }
}


double
applyFilter(struct Filter *filter, cs1300bmp *input, cs1300bmp *output){

  long long cycStart, cycStop;

  cycStart = rdtscll();

  short int hDim = input -> height - 1;
  short int wDim = input -> width - 1;
    
  output -> height = hDim + 1;
  output -> width = wDim+ 1;
    
  short int InputFilter[3][3];
  #pragma omp parallel for
  for(short int i = 0; i < 3; i++){
      InputFilter[i][0] = filter -> get(i, 0);
      InputFilter[i][1] = filter -> get(i, 1);
      InputFilter[i][2] = filter -> get(i, 2);
  }
  
  // Gauss
  if(InputFilter[0][1] == 4){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0; 
          
          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parTwo += input -> color[plane][addR][col] << 2; //0,1
          
          parOne += input -> color[plane][row][addC] << 2; //1,0
          parTwo += input -> color[plane][row][col] << 3;   //1,1
          parThree += input -> color[plane][row][addC] << 2; //1,2
            
          parTwo += input -> color[plane][addR][col] << 2; //2,1
            
          parOne += parTwo + parThree;
          parOne = parOne/24;
            
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
    
  //Edge
  else if(InputFilter[1][1] == -7){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0; 
          
          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += input -> color[plane][subR][subC]; //0,0
          parTwo += input -> color[plane][subR][col]; //0,1
          parThree += input -> color[plane][addR][col]; //0,2
          
          parOne += input -> color[plane][row][subC]; //1,0
          parTwo += -((input -> color[plane][row][col] << 2) + 3) ;   //1,1
          parThree += -(input -> color[plane][row][addC]); //1,2
           
          parOne += input -> color[plane][addR][subC]; //2,0
          parTwo += input -> color[plane][addR][col]; //2,1
          parThree += input -> color[plane][addR][addC]; //2,2
          
          parOne += parTwo + parThree;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
    
  //emboss
  else if(InputFilter[0][1] == 1 && InputFilter[0][2] == -1){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0;
          
          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += input -> color[plane][subR][subC]; //0,0
          parTwo += input -> color[plane][subR][col]; //0,1
          parThree += -(input -> color[plane][addR][col]); //0,2
          
          parOne += input -> color[plane][row][subC]; //1,0
          parTwo += input -> color[plane][row][col] ;   //1,1
          parThree += -(input -> color[plane][row][addC]); //1,2
           
          parOne += input -> color[plane][addR][subC]; //2,0
          parTwo += -(input -> color[plane][addR][col]); //2,1
          parThree += -(input -> color[plane][addR][addC]); //2,2
            
          parOne += parTwo + parThree;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
                    
  //Hline
  else if(InputFilter[0][1] == -2){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0;

          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += -(input -> color[plane][subR][subC]); //0,0
          parTwo += -(input -> color[plane][subR][col] << 1); //0,1
          parThree += -(input -> color[plane][addR][col]); //0,2
           
          parOne += input -> color[plane][addR][subC]; //2,0
          parTwo += input -> color[plane][addR][col] << 1; //2,1
          parThree += input -> color[plane][addR][addC]; //2,2
           
          parOne += parTwo + parThree;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
  
  //Vline
  else if(InputFilter[1][0] == -2){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0;

          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += -(input -> color[plane][subR][subC]); //0,0
          parThree += input -> color[plane][addR][col]; //0,2
          
          parOne += -(input -> color[plane][row][subC] << 1); //1,0
          parThree += -(input -> color[plane][row][addC] << 1); //1,2
           
          parOne += -(input -> color[plane][addR][subC]); //2,0
          parThree += input -> color[plane][addR][addC]; //2,2
            
          parOne += parOne + parThree;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
  
  //Sharpen
  else if(InputFilter[0][0] == 11){
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0;
          
          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += (input -> color[plane][subR][subC] << 3) + 3; //0,0
          parTwo += (input -> color[plane][subR][col] << 3) + 2; //0,1
          parThree += input -> color[plane][addR][col]; //0,2
          
          parOne += -(input -> color[plane][row][subC]); //1,0
          parTwo += -(input ->color[plane][row][col]);   //1,1
          parThree += -(input ->color[plane][row][addC]); //1,2
           
          parOne += -(input -> color[plane][addR][subC]); //2,0
          parTwo += -(input -> color[plane][addR][col]); //2,1
          parThree += -(input -> color[plane][addR][addC]); //2,2
          
          parOne += parTwo + parThree;      
          parOne = parOne / 20;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
        
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }
  
  //average
  else{
    #pragma omp parallel for
    for(short int plane = 0; plane < 3; plane++){
      for(short int row = 1; row < hDim; row++){
        const short int addR = row + 1;
        const short int subR = row - 1;
        for(short int col = 1; col < wDim; col++){
          short int parOne = 0;
          short int parTwo = 0;
          short int parThree = 0;
          
          const short int addC = col + 1;
          const short int subC = col - 1;
          
          parOne += input -> color[plane][subR][subC]; //0,0
          parTwo += input -> color[plane][subR][col]; //0,1
          parThree += input -> color[plane][addR][col]; //0,2
          
          parOne += input -> color[plane][row][subC]; //1,0
          parTwo += input ->color[plane][row][col];   //1,1
          parThree += input ->color[plane][row][addC]; //1,2
           
          parOne += input -> color[plane][addR][subC]; //2,0
          parTwo += input -> color[plane][addR][col]; //2,1
          parThree += input -> color[plane][addR][addC]; //2,2
            
          parOne += parTwo + parThree;
          parOne = parOne / 9;
          parOne = (parOne < 0) ? 0 : parOne;
          parOne = (parOne > 255) ? 255 : parOne;
          
          output -> color[plane][row][col] = parOne;
        }
      }      
    }
  }              

  cycStop = rdtscll();
  double diff = cycStop - cycStart;
  double diffPerPixel = diff / (output -> width * output -> height);
  fprintf(stderr, "Took %f cycles to process, or %f cycles per pixel\n",
	  diff, diff / (output -> width * output -> height));
  return diffPerPixel;
  }
