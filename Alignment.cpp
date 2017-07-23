#include "Alignment.h"
#include <fstream>

using namespace std;

//constructor
Alignment::Alignment( Network net1, Network net2)
{
	
    //compare networks to find the biggest one
    if( net1.size > net2.size )
	{
		N = net2.size;
        reverse = true;
		network1 = net2;
		network2 = net1;
	}
	else
	{
		N = net1.size;
        reverse = false;
		network1 = net1;
		network2 = net2;
	}
    
	//maximum degree of the network
    if(network1.maxDeg > network2.maxDeg)
		maxDeg = network1.maxDeg;
	else
		maxDeg = network2.maxDeg;

	similarity = new float*[network1.size];
	for(int i=0; i<network1.size; i++)
	{
		similarity[i] = new float[network2.size];
	}

	M = new int*[network1.size];    
	for(int i=0; i<network1.size; i++)
	{
		M[i] = new int[network2.size];
	}

    xy = new int[N];    
    yx = new int[N];
    update = new int[N];
    slack = new float[N];
    slackx = new int[N];
    lx = new float[N];
    ly = new float[N];
    S = new bool[N];
    T = new bool[N];
    prev = new int[N];
    
    clusterSim = new float*[network1.size];
 	for(int i=0; i<network1.size; i++)
	{
		clusterSim[i] = new float[network2.size];
	}
    
    for(int i=0; i<network1.size; i++)
		for(int j=0; j<network2.size; j++)
        {
			clusterSim[i][j] = 0;
		}
    
    clusterNum = new float*[network1.size];
 	for(int i=0; i<network1.size; i++)
		clusterNum[i] = new float[network2.size];
	
    
    for(int i=0; i<network1.size; i++)
		for(int j=0; j<network2.size; j++)
			clusterNum[i][j] = 0;
}

//read the file of aligned clusters
void Alignment::readAlignedCluster(string name)
{
    string line;
    ifstream input(name.c_str());
    string term1,term2,token;
    int id1,id2;
    float score;

    if(!input) //there is not such file 
    {
        cout << "Error: Module similarity file doesn't exist"<<endl;
        exit (1);
    }
 
    while (getline(input, line))  //reading network file
    {
        istringstream tokenizer(line);
        getline(tokenizer, term1, '\t');
        getline(tokenizer, term2, '\t');
        
        if (reverse) {
            id2 = network2.getCluster( term1 );
            id1 = network1.getCluster( term2 );            
        }
        else {
            id1 = network1.getCluster( term1 );
            id2 = network2.getCluster( term2 );
        }
        
        getline(tokenizer, token, '\t');
        
        score = atof(token.c_str());	
        
        for (int i=0; i<network1.ms[id1]; i++) {            
            for (int j=0; j<network2.ms[id2]; j++) {              
                clusterSim[network1.clusterMember[id1][i]][network2.clusterMember[id2][j]]+=score;
                clusterNum[network1.clusterMember[id1][i]][network2.clusterMember[id2][j]]++;
            }
        }
    }
    
    for (int i=0; i<network1.size; i++) 
        for (int j=0; j<network2.size; j++)                
            if(clusterNum[i][j]!=0)
                clusterSim[i][j]=clusterSim[i][j]/clusterNum[i][j];
    
    
    float min=0,max=0; int counter=0;
    for(int i=0; i<network1.size; i++)
		for(int j=0; j<network2.size; j++)
            if(max < clusterSim[i][j]) 
                max = clusterSim[i][j];
            
		
    
    //==== Normalize =====
    for(int i=0; i<network1.size; i++)
        for(int j=0; j<network2.size; j++){
            clusterSim[i][j] = clusterSim[i][j]/max;
        }
}

void Alignment::setSimilarities(string name)
{
    float bb=0.01,counter=0;
    float max=0;
    string line;
    string term1,term2,token;
    int id1, id2;
    istringstream Stream;
    
    for(int i=0; i<network1.size; i++)
		for(int j=0; j<network2.size; j++)
		{
			similarity[i][j] = 0;
		}
	
    if( bb != 1 )
	{
        bb=1;
        ifstream input(name.c_str());
        
        if(!input) //there is not such file 
        {
            cout << "Error: Blast file doesn't exist"<<endl;
            exit (1);
        }

		while (getline(input, line))
		{
            Stream.clear();
            Stream.str(line);
			getline(Stream, term1, '\t');
			getline(Stream, term2, '\t');
			getline(Stream, token, '\t');

           
			if(term1.at(term1.length()-1)=='\n')
				term1 = term1.substr(0,term1.length()-1);
            
			if(term2.at(term2.length()-1)=='\n')
				term2 = term2.substr(0,term2.length()-1);
            
			if(reverse) {
                id2 = network2.getID( term1 );
                id1 =  network1.getID( term2 );
            }
            else {
                id1 = network1.getID( term1 );
                id2 =  network2.getID( term2 );
            }
            similarity[ id1 ][ id2 ] =  atof(token.c_str()) ;
            if(max < similarity[ id1 ][ id2 ])
                max=similarity[ id1 ][ id2 ];
		}
	    
        input.close();
        
        for(int i=0; i<network1.size; i++)
            for(int j=0; j<network2.size; j++)
                similarity[i][j] = similarity[i][j]/max; 
        
        for(int i=0; i<network1.size; i++)
            for(int j=0; j<network2.size; j++)
                similarity[i][j] = (float)( ( 1 - bb) *  similarity[i][j] + bb * clusterSim[i][j] );

    }
   
}


void Alignment::init_labels(double aa,double kk)
{
        for(int i=0; i<N; i++)
        for(int j=0; j<N; j++)
            if(i < network1.size && j<network2.size) {
                similarity_hung[i][j]=  (1-aa) * similarity_hung[i][j] +  aa * similarity[i][j];
                similarity_first[i][j]= (1-kk) * similarity_hung[i][j] + kk * similarity[i][j];
            }
        else
            similarity_hung[i][j]=-10000000;
    
        for(int i=0;i<N;i++) {lx[i]=0;ly[i]=0;}
        for (int x = 0; x < N; x++) 
            for (int y = 0; y < N; y++)
                if (similarity_hung[x][y]>lx[x]) 
                    lx[x] = similarity_hung[x][y];							
    }
    void Alignment::update_labels()
    {
        int x, y;float  delta = 1000000;
        for (y = 0; y < N; y++) //calculate delta using slack
            if (!T[y] && slack[y] < delta)
                delta = slack[y];
        for (x = 0; x < N; x++) //update X labels
            if (S[x]) lx[x] -= delta;
        for (y = 0; y < N; y++) //update Y labels
            if (T[y]) ly[y] += delta;
        for (y = 0; y < N; y++) //update slack array
            if (!T[y])
                slack[y] -= delta;

    }
    void Alignment::add_to_tree(int x, int prevx) 
    //x - current vertex,prevx - vertex from X before x in the alternating path,
    //so we add edges (prevx, xy[x]), (xy[x], x)
    {
        S[x] = true; //add x to S
        prev[x] = prevx; //we need this when augmenting
        for (int y = 0; y < N; y++) //update slacks, because we add new vertex to S
            if (lx[x] + ly[y] - similarity_hung[x][y] < slack[y])
            {
                slack[y] = lx[x] + ly[y] - similarity_hung[x][y];
                slackx[y] = x;
            }
    }
    void Alignment::augment() //main function of the akgorithm
    {
        if (max_match == N) return;
        int x,y,root;
        int  wr = 0, rd = 0; //q - queue for bfs, wr,rd - write and read pos in queue
        int *q;
        q = new int[N];
        for(int i=0;i<N;i++) {S[i]=false;T[i]=false;prev[i]=-1;}
        for (x = 0; x < N; x++) //finding root of the tree
            if (xy[x] == -1)
            {
                q[wr++] = root = x;
                prev[x] = -2;
                S[x] = true;
                break;
            }
        for (y = 0; y < N; y++) //initializing slack array
        {
            slack[y] = lx[root] + ly[y] - similarity_hung[root][y];
            slackx[y] = root;
        }
        while (true) { //main cycle
            while (rd < wr) { //building tree with bfs cycle
                x = q[rd++]; //current vertex from X part
                for (y = 0; y < N; y++) //iterate through all edges in equality graph 
                    if (similarity_hung[x][y] == lx[x] + ly[y] && !T[y]) {
                        if (yx[y] == -1) break; //an exposed vertex in Y found, so
                        //augmenting path exist!
                        T[y] = true; //else just add y to T
                        q[wr++] = yx[y]; //add vertex yx[y], which is matched
                        //with y, to the queue
                        add_to_tree(yx[y], x); //add edges (x,y) and (y,yx[y]) to the tree
                    }
                if (y < N) break; //augmenting path found
            }
            if (y < N) break; //augmenting path found
            update_labels(); //augmenting path not found, so improve labeling
            wr = rd = 0;
            for (y = 0; y < N; y++) 
            //in this cycle we add edges that were added to the equality graph as a
            //result of improving the labeling, we add edge (slackx[y], y) to the tree if
            //and only if !T[y] && slack[y] == 0, also with this edge we add another one
            //(y, yx[y]) or augment the matching, if y was exposed 
                if (!T[y] && slack[y] == 0) {
                    if (yx[y] == -1) {
                        x = slackx[y];
                        break;
                    }
                    else {
                        T[y] = true; //else just add y to T
                        if (!S[yx[y]]) {
                            q[wr++] = yx[y]; //add vertex yx[y], which is matched with
                            //y, to the queue
                            add_to_tree(yx[y], slackx[y]); ////and add edges (x,y) and (y,
                            //yx[y]) to the tree
                        }
                    }        
                }
            if (y < N) break; //augmenting path found!            
        }
        delete [] q;
        if (y < N) {
            max_match++; //increment matching
            //in this cycle we inverse edges along augmenting path
            for (int cx = x, cy = y, ty; cx != -2; cx = prev[cx], cy = ty)
            {
                ty = xy[cx];
                yx[cy] = cx;
                xy[cx] = cy;
            }
            augment(); //recall function, go to step 1 of the algorithm
        }
    }

    float Alignment::hungarian(double aa,double kk)
    {
        float ret = 0; //weight of the optimal matching
        max_match = 0; //number of vertices in current matching
        for(int i=0;i < N;i++) {xy[i]=-1;yx[i]=-1;update[i]=0;};
        init_labels(aa,kk); //step 0
        augment(); //steps 1-3
        for (int x = 0; x < network1.size; x++) //forming answer there
            ret += similarity_hung[x][xy[x]];
        return ret;
    }

    int Alignment::findBestAligned()
    {
        float temp = -100000000;
        int ba;

        for(int i=0; i<network1.size; i++)
            for(int j=0; j<network2.size; j++)
                if(xy[i] == j && !update[i]) 
                    if(similarity_first[i][j]>temp) {
                        temp = similarity_first[i][j];
                        ba = i;		
                    }
                
        update[ba] = 1;
        return ba;
    }
    
    void Alignment::Dynamic_hungarian(string name, double kk)
    {
        float lambda = 0.2;
        double aa;
        if(kk > 0.2) aa = kk - 0.2; else aa = kk; 
        int totalScore;
        int coeff;
        if(network2.numOfEdge>network1.numOfEdge) 
            coeff = network2.numOfEdge/network1.numOfEdge;
        else 
            coeff = network1.numOfEdge/network2.numOfEdge;

        best = new int[network1.size];
        similarity_hung = new float*[N];
        similarity_first = new float*[N];
        
        for(int i=0; i<N; i++){
            similarity_hung[i] = new float[network2.size];
            similarity_first[i] = new float[network2.size];
        }
        for(int i=0; i<N; i++)
            for(int j=0; j<network2.size; j++){
                similarity_hung[i][j] = 0;
                similarity_first[i][j] = 0;
            }
        
        // ------------------------------
        //hubalign score
        
        float *nodeScore1 = new float[network1.size]; //scores of nodes of smaller network
        float *nodeScore2 = new float[network2.size]; //scores of nodes of bigger network

        for(int c1=0; c1< network1.size; c1++)
            nodeScore1[c1]=(1-lambda)*network1.nodeWeight[c1];
        for(int c1=0; c1< network2.size; c1++)
            nodeScore2[c1]=(1-lambda)*network2.nodeWeight[c1];
        
        //find max score 
        //finding the nodescore
        for (int c1=0; c1<network1.size; c1++){
            for (int c2=0; c2<network1.size; c2++) 
                nodeScore1[c1]+= lambda*network1.edgeWeight[c1][c2];
        }
        for (int c1=0; c1<network2.size; c1++){
            for (int c2=0; c2<network2.size; c2++)
                nodeScore2[c1] += lambda*network2.edgeWeight[c1][c2];
        } 
        
        //======first network
        float max = -10000;
        for (int c1=0; c1<network1.size; c1++) {
            if (max < nodeScore1[c1]) {
                max = nodeScore1[c1];
            }
        }
        
        //====== second network
        for (int c1=0; c1<network2.size; c1++) 
            if (max < nodeScore2[c1]) 
                max = nodeScore2[c1];
        for (int c1=0; c1<network1.size; c1++) 
            nodeScore1[c1] = nodeScore1[c1]/max;
            
		for (int c1=0; c1<network2.size; c1++) 
            nodeScore2[c1] = nodeScore2[c1]/max;
            
        float min;
        for(int i=0; i<network1.size; i++) 
            for(int j=0; j<network2.size; j++) { 
                if(nodeScore1[i]>nodeScore2[j])
                    min = nodeScore2[j];
                else
                    min = nodeScore1[i];
                similarity_hung[i][j] = min;//just the topology
            }
        
        hungarian(aa,kk);

        int bestAligned,tempNode,candid;
        float ret = 0;
        
        for (int k = 0; k < network1.size; k++) {
            
            bestAligned=findBestAligned();
            tempNode = xy[bestAligned];
            best[k] = bestAligned; //to save then for next step
            similarity_hung[bestAligned][tempNode]+=1000; 
            for (int i=0; i<network1.deg[bestAligned]; i++) {
                candid = network1.neighbor[bestAligned][i];
                for (int j=0; j<network2.deg[tempNode]; j++) {
                    similarity_hung[candid][network2.neighbor[tempNode][j]]=similarity_hung[candid][network2.neighbor[tempNode][j]]+(float)(coeff/max);
                }
                yx[xy[candid]] = -1;
                xy[candid] = -1;
                lx[candid] = 1000000;
                for (int y = 0; y < N; y++) //update Y labels
                    if(lx[candid] > ly[y] - similarity_hung[candid][y])
                        lx[candid] = ly[y] - similarity_hung[candid][y];
                
                max_match--;
                augment();                
            }
        
        }

        string alignFile = name;
        //insert a definite suffix for the alignment file
        alignFile.append(".alignment");
        ofstream alignmentFile( alignFile.c_str());        
/*        if(reverse)
            for (int x = 0; x < network1.size; x++)
                alignmentFile <<network2.getName( xy[x] )<< "\t" <<network1.getName( x ) << endl;            
        else
*/            for (int x = 0; x < network1.size; x++)
                alignmentFile <<network1.getName( x )  << "\t" <<network2.getName( xy[x] ) << endl;
        
        alignmentFile.close();							
	cout <<  "Finished !!" << endl << endl;            
    }

//instructor
Alignment::Alignment(void)
{
}
//destructor
Alignment::~Alignment(void)
{
}
