#include "Network.h"

using namespace std;
class Alignment
{
public:
	float **similarity; // The similarity values of nodes
	float **similarity_hung;float **similarity_first;  
    int **M;
    int *best;
    int *xy;
    int *update;
    int *yx;
    int N, max_match;
    float *lx, *ly;
    float *slack;
    int *slackx;
    bool *S; bool *T;
    int *prev;
    
    float **clusterSim;
    float **clusterNum;	
   
	int *alignment; 
	int maxDeg;
	int *alignedCluster;
    bool reverse; 
	Network network1;
	Network network2;
    
	Alignment( Network net1, Network net2 );
    Alignment(void);
    ~Alignment(void);
    void align(double a);
    void computeClusterSimilarity(void);
    void readAlignedCluster(string name);
	void setSimilarities(string name);
	float hungarian(double aa,double kk);
	void Dynamic_hungarian(string name,double kk);

private:
    void init_labels(double aa,double kk);
    void update_labels();
    void add_to_tree(int x, int prevx);
    void augment();
    int findBestAligned();
};
