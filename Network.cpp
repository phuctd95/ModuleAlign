#include "Network.h"

Network::Network(char *nam)
{
    stringstream strm1, strm2;
    ifstream input;
    try
	{
		name = nam;
        string line;
		string token;
        int id1, id2;

		MapInt2List mapNeighbor;
        MapInt2List mapCluster;
        MapInt2List mapMember;

        strm1 << name << ".net";
        strm2 << name << ".cluster"; 
        
		input.open((strm1.str()).c_str());
        if(!input) //there is not such file 
        {
            cout << "Error: Network file doesn't exist"<<endl;
            exit (1);
        }
        
        numOfEdge = 0; //number of edges of the network
		maxDeg = 0;    //max degree of the network
        
		while (getline(input, line))  //reading network file
		{
			istringstream tokenizer(line);
			string token;
			getline(tokenizer, token, '\t');
            if(token.length()==0) //the input node is incorrect
             {
                 cout << "Error: No node in first column" <<endl;
                 exit (1);
             }

            //assign first node of the neywork an ID
            mapName.insert(MapString2Int::value_type(token,(int)mapName.size()));            
            
            id1 = mapName[ token ];
			mNames.insert(MapInt2Name::value_type(id1, token));
			list<int> neighbs1;
			mapNeighbor.insert(MapInt2List::value_type( id1, neighbs1));

            //read the second node in a row
			getline(tokenizer, token, '\t');
            if(token.length()==0)//the input node is incorrect
            {
                cout << "Error: No node in second column" <<endl;
                exit (1);
            }
            if(token.at(token.length()-1)==13)
			{
				token = token.substr(0,token.length()-1);
			}
            
            //assign second node of the neywork an ID			
            mapName.insert(MapString2Int::value_type(token,(int)mapName.size()));
            id2 = mapName[ token ];
            mNames.insert(MapInt2Name::value_type(id2, token));
			list<int> neighbs2;
			mapNeighbor.insert(MapInt2List::value_type( id2, neighbs2));
            
            //insert first node in neigbor list of the second nod and vise versa.
			mapNeighbor[ id1 ].push_front( id2 );
			mapNeighbor[ id2 ].push_front( id1 );

		}
        input.close();
        
        //=========== reading the cluster file =========
        
        input.open((strm2.str()).c_str());
        if(!input) //there is not such file 
        {
            cout << "Error: Cluster file doesn't exist"<<endl;
            exit (1);
        }

        while (getline(input, line))  //reading network file
		{
			istringstream tokenizer(line);
            
            //read the cluster 
			getline(tokenizer, token, '\t');
            
            if(token.length()==0) //the input cluster is incorrect
            {
                cout << "Error: No cluster in first column" <<endl;
                exit (1);
            }
            
            //assign the cluster an ID
            clusterName.insert(MapString2Int::value_type(token,(int)clusterName.size()));            
            
            id1 = clusterName[ token ];
            cNames.insert(MapInt2Name::value_type(id1, token));
                        
            //read the second node in a row
            getline(tokenizer, token, '\t');
            
            //find the nodes in the cluster
            istringstream temptokenizer(token);
            while (getline(temptokenizer, token, ',')) {
            
                //assign the node an ID			
                if(mapName.find(token)!=mapName.end()){
                    id2 = mapName[ token ];
            
                	//creat the corresponding list for its clusters
                	list<int> clusters;
                	mapCluster.insert(MapInt2List::value_type( id1, clusters));

                	//creat the corresponding list for the cluster members
                	list<int> members;
                	mapMember.insert(MapInt2List::value_type( id1, members));

                	//insert cluster in cluster list of the node and vise versa.
                	mapCluster[ id2 ].push_front( id1 );
                	//insert node in member list of a cluster.
                	mapMember[ id1 ].push_front( id2 );
								}
            }

		}
        input.close();
        //======================
        
        size = mapName.size(); //number of nodes of the network
        clusterSize = clusterName.size();
        //finding the neighbors and edges of the network with respect to mapNeigbor
        neighbor = new int*[size];
        cluster = new int*[size];
	    clusterMember = new int*[clusterSize];
        deg = new int[size];
        cs = new int [size]; //cs : cluster size: number of clusters that a node is in.    
        ms = new int [clusterSize];
        
		list<int> tempList;
        list<int> tempClusterList;
		list<int>::iterator it;

		//makeing the members of each cluster
        for(int i=0; i<clusterSize; i++)
		{
            tempList = mapMember[i];
            tempList.sort();
            tempList.unique();
            
            ms[i] = tempList.size();
            clusterMember[i] = new int[ms[i]];
            int j;
			for( j=0,it=tempList.begin(); it!= tempList.end() ; it++, j++)
				clusterMember[i][j] = *it;

        }        
        
		for(int i=0; i<size; i++)
		{

            tempList = mapNeighbor[ i ];
            tempClusterList = mapCluster[ i ];
			tempList.sort();
			tempList.unique();
            tempClusterList.sort();
            tempClusterList.unique();


			deg[i] = tempList.size();
            cs[i] = tempClusterList.size();
			neighbor[i] = new int[deg[i]];
            cluster[i] = new int[cs[i]];
            
			numOfEdge += deg[i];

			if(deg[i] > maxDeg)
				maxDeg = deg[i];
            int j;
			for(j=0,it=tempList.begin(); it!= tempList.end() ; it++, j++)
				neighbor[i][j] = *it;
            for(j=0,it=tempClusterList.begin(); it!= tempClusterList.end() ; it++, j++)
				cluster[i][j] = *it;

		}
        
		numOfEdge = numOfEdge / 2;
        //initialize some variables after findin the size of the network
        nodeWeight = new float[size]; //weight of each node that changes during making the skeletone
        newDeg = new int[size];  //new deg of each node during making the skeletone
        edgeWeight = new float*[size]; //weight of each edge that changes during making the skeletone
        remEdge = new int*[size];
        remNode = new int[size];  //removed nodes
        for(int c1=0; c1<size; c1++)
            edgeWeight[c1]=new float[size];
        for(int c1=0; c1<size; c1++)
            remEdge[c1]=new int[size];

        
        //initial value for weights of nodes and edges
        skeletonInitialValue();
	}
	catch (exception &e)
	{
        cerr << "Error in input file" << endl;
        exit(1);
    }
}

void Network::skeletonInitialValue()
{
    for (int c1=0; c1<size; c1++)
        nodeWeight[c1]=0;
    for (int c1=0; c1<size; c1++)
        for (int c2=0; c2<size; c2++)
            edgeWeight[c1][c2]=0;
    
    for (int c1=0; c1<size; c1++)
        for (int c2=0; c2<deg[c1]; c2++)
            edgeWeight[c1][neighbor[c1][c2]]=1;
    for (int c1=0; c1<size; c1++)
        remNode[c1]=0;
    for (int c1=0; c1<size; c1++)
        newDeg[c1]=deg[c1];
}

void Network::makeSkeleton(int t)
{

    removeDegOne();

    for(int c1=2; c1<t; c1++)
        removeDeg(c1);

    int max =0 ;
    for(int c=0; c<size; c++)
        if(remNode[c]==0)
            if(max < newDeg[c])
                max=newDeg[c];    
}

void Network::removeDegOne()
{
    bool flag = true;
    while(flag){ //while there is a node with degree one
        flag=false;
        //for each node of the network
        for (int c1=0; c1<size; c1++) 
        {
            if(newDeg[c1]==1) {//remove
                for(int c2=0; c2<deg[c1]; c2++)
                    if (remNode[neighbor[c1][c2]]==0){ //loop
                        nodeWeight[neighbor[c1][c2]]+=nodeWeight[c1]+edgeWeight[c1][neighbor[c1][c2]]; //neighbor of deleted node
                        newDeg[neighbor[c1][c2]]--;
                        newDeg[c1]--;
                    }
                remNode[c1]=1; //node is removed
            }
            
        }
        //check wether or not any node with degree one has produced in teh network after removing
        for (int c1=0; c1<size; c1++)
            if (newDeg[c1]==1)
                flag = true;
    }
    //removes nodes with no edges after removing the nodes with degree one
    for(int c1=0; c1<size; c1++)
        if(newDeg[c1]<1 && remNode[c1]==0){
            remNode[c1]=1;
        }
}

void Network::removeDeg(int degree)
{
    int *ngh; //neighbors of vertex with degree d
    int nghSize; //number of neigbors
    bool flag=true;
    int d = degree;
    float temp;
    
    while(flag){ //while there is such node
        flag=false;
        ngh = new int[d];
        for (int c1=0; c1<size; c1++) 
            if(newDeg[c1]==d) {//remove
                temp = nodeWeight[c1];
                nghSize = 0;
                for(int c2=0; c2<deg[c1]; c2++) //find two unremoved neighbors
                    if (remNode[neighbor[c1][c2]]==0){
                        ngh[nghSize]= neighbor[c1][c2];
                        nghSize++;
                        temp += edgeWeight[c1][neighbor[c1][c2]];
                    }
                
                newDeg[c1]=0;
                for(int c2=0; c2<nghSize; c2++){
                    newDeg[ngh[c2]]--;
                    for(int c3=c2+1; c3<nghSize; c3++){
                        edgeWeight[ngh[c2]][ngh[c3]]+= temp/(d*(d-1)/2);
                        edgeWeight[ngh[c3]][ngh[c2]]= edgeWeight[ngh[c2]][ngh[c3]]; 
                    }  
                }    
                remNode[c1]=1;
            }
        //check if any node with degree equal or less than the input degree has produced after removing
        for(int c2=degree; c2>1; c2--)
            if(!flag){
                for (int c1=0; c1<size; c1++)
                    if (newDeg[c1]==c2){
                        flag = true;
                        d=c2;
                        break;
                    }
            }
            else
                break;
        if(!flag)                
            for (int c1=0; c1<size; c1++)
                if (newDeg[c1]==1){
                    flag = true;
                    d=degree;
                    removeDegOne();
                    break;
                }        
    }
    
    for(int c1=0; c1<size; c1++)
        if(newDeg[c1]<1 && remNode[c1]==0){
            remNode[c1]=1;
        }
    delete [] ngh; 
}


//constructor
Network::Network(void)
{
}
//destructor
Network::~Network(void)
{
}

