#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include "math.h"
#include <time.h>
#include <map>
#include <list>

using namespace std;

class Network
{
    typedef map<string, int , less<string> > MapString2Int;
	typedef map<int, list<int>, less<int> > MapInt2List;
	typedef map<int, string, less<int> > MapInt2Name;

	MapInt2Name mNames;
    MapInt2Name cNames;
    
public:
    MapString2Int mapName;
    MapString2Int clusterName;
    
	int **neighbor; //neighbor of each node of the network
	int N;
	int **cluster; //clusters of each node of the network
    int **clusterMember;
    int size; //number of nodes
	int clusterSize;
    int maxDeg;//maximum degree of the network
    int *deg; // degree of each node
    int *cs; // number of clusters a nod is in.
    int *ms; //number of memebers of a acluster.
    int *remNode;  //removed nodes
    int *newDeg;    //new degree of each node after each removing 
    float *nodeWeight; //weight of each node 
	float **edgeWeight; //weight of each edge
    int **remEdge;
    int numOfEdge; //number of edges 
	char* name; //name of the network
    
	Network(char *nam); 
	
    Network(void);
    
	//destructor
    ~Network(void); 
    
    void skeletonInitialValue();
    void makeSkeleton(int t);
    void removeDegOne();
    void removeDeg(int degree);
    string getName(int id)
	{
		return mNames[ id ];
	};
  
    string getClusterName(int id)
	{
		return cNames[ id ];
	};
    
	int getID(string name)
	{
		return mapName[ name ];
	};
    int getCluster(string name)
	{
		return clusterName[ name ];
	};

private:
};

