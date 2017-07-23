#include "Alignment.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

int main(int argc, char *argv[])
{
    double alpha = 0.4; 
    stringstream strm;

	int degree = 10;
    char* name1; char* name2;char *moduleFile;      
	
    try
	{
		
        if(argc < 2) {
            cout << "There should be 2 files as input!" <<endl;
            return -1;
        }
        else //input arguments
        {
            int i = 1; //counter for input parameters
            name1 = argv[ i++ ]; //first network
            name2 = argv[ i++ ]; //second network
            moduleFile = argv[ i++ ];
            
            while (i<argc) //check all the input parameters
            {
                if ( ( strlen(argv[i]) == 2 ) && ( argv[i][0]=='-' ) && ( argv[i][1]=='a') && ( i + 1 < argc) )
                {
                    i++;
                    if( argv[i-1][1]=='a' )
                    {
                        alpha = atof(argv[i]);
                        if( alpha <0 || alpha > 1)
                        {
                            cout << "Error : value of alpha must be between 0 and 1 :" << endl;
                            return -1;
                        }
                    }
                    i++;
                }
                else
                {
                    cout << "Error in argument : " << argv[i] << endl;
                    return -1;
                }
            } //end while
           
            cout << "  _   _           _           _           ___   _" <<endl;
            cout << "/  \\/  \\   __   _| |  _   _  | |    ___  / _ \\ | |    _   __   _  ___" <<endl;
            cout << "| /\\/\\ | / _ \\ / _ | | |_| | | |/| / -_\\ |'-'| | |/| | | / _ \\ |\\/_  |" << endl;
            cout << "|/    \\| \\__ / \\___| |_____| |__ / \\_\\_  |/ \\| |__ / |_| \\__ / |_/ |_|" << endl;
            cout << "                                                          _//"<<endl;
            cout << "                                                         \\_/"<<endl;
            cout << endl <<"is running ... " << endl << endl;

            Network network1(name1);
            Network network2(name2);
            
            network1.makeSkeleton(degree);
            network2.makeSkeleton(degree);

            //Initializes the alignment class
            Alignment alignment( network1, network2 );
            string khar = moduleFile;
            alignment.readAlignedCluster(khar);
            
            strm.str("");strm.clear();
            strm  << name1 << "-" << name2 <<".blast";
            alignment.setSimilarities(strm.str());

            //making the name for output file
            strm.str("");strm.clear();
            strm  << name1 << "-" << name2 <<"-a" << alpha;

            alignment.Dynamic_hungarian(strm.str(), alpha);
        }
    }
    catch(exception &e)
    {
        cout << "Error in arguments or input files!" << endl;
        e.what();
        return -1;
    }
        
    return 0;
}



