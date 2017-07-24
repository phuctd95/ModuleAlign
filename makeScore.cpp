
#include <iostream>
#include <map>
#include <sstream>
#include <fstream>
#include <string>
#include <cstdlib>

using namespace std;

void  get_gene_root_number(string tree, int& geneNum, int& rootNum) {

    string line;
    int term1;
    ifstream input;
    istringstream Stream;

    // ------------------------------
    input.open(tree.c_str());
    //find geneNum
    getline(input, line);getline(input, line);
    Stream.clear();
    Stream.str(line);
    Stream >> term1;
    geneNum = term1-1;
    //find rootNum
    input.seekg(0,ios_base::end);
    char ch = ' ';                        //Init ch not equal to '\n'
    while(ch != '\n'){
    input.seekg(-2,std::ios_base::cur); //Two steps back, this means we
                                              //will NOT check the last character
        if((int)input.tellg() <= 0){        //If passed the start of the file,
            input.seekg(0);                 //this is the start of the line
            break;
        }
        input.get(ch);                      //Check the next character
    }

    getline(input, line);
    Stream.clear();
    Stream.str(line);
    Stream >> term1;
    rootNum = term1;
		input.close();
}

void  make_cluster_members(string tree, string gene_mapping, int geneNum, int rootNum, int* cluster_size, string* cluster_children, string* cluster_member) //cluster_size shows the members of each clusters
{
    int* lnodes;int* rnodes;
		string line,term2,term3,temp;
    int term1;
    ifstream input;
		istringstream Stream;
    map<int, string> geneMap;
    map<int, string>::iterator git;

		// ------------------------------
    git=geneMap.begin();
    input.open(gene_mapping.c_str());
		while (!input.eof()) {
        getline(input, line);
				Stream.clear();        
				Stream.str(line);
        Stream >> term1;
        Stream >> term2;
        if(geneMap.find(term1)==geneMap.end()){
            geneMap.insert(git, pair<int, string>(term1,term2));
            git++;
        }
    }
    input.close();

		// ------------------------------
    input.open(tree.c_str());
		
    for (int i=1; i<=geneNum; i++) { //geneNum may be ?????????????
        stringstream ss;
        ss << i;
        cluster_children[i]=ss.str();
				cluster_member[i]=geneMap.find(i)->second;
    }

    for (int i=1; i<=geneNum; i++) 
        cluster_size[i]=1;

    lnodes=new int[rootNum+1];
    rnodes=new int[rootNum+1];

		// ------------------------------

    input.seekg(0,ios_base::beg);
		getline(input, line);		
    while (!input.eof()) {
        getline(input, line);
		    Stream.clear();
  		  Stream.str(line);
        Stream >> term1;
        Stream >> term2;
        Stream >> term3;	
        lnodes[term1]=atoi(term2.c_str());
        rnodes[term1]=atoi(term3.c_str());
				temp = term2 + "," + term3;
				cluster_children[term1]=temp;
        cluster_member[term1]=cluster_member[lnodes[term1]]+","+cluster_member[rnodes[term1]];
        cluster_size[term1]=cluster_size[lnodes[term1]]+cluster_size[rnodes[term1]];
		}
    input.close();
}     


//////////////////////////////////////////////////////////////////////

int main(int argc, const char * argv[])
{      
    string line,outstring;
    int term1,row,col, geneNum1, geneNum2, rootNum1, rootNum2;
    string term2,term4,term5;
    double term3, sum;
    ifstream input;
    ofstream output;
		istringstream Stream1, Stream2;
		stringstream strm1,strm2;    

    string name1 = argv[1];
    string name2 = argv[2];
		string similarity = argv[3];
    string outFile = argv[4];
		
		// -----------------------------
		// get number of genes and root
		strm1  << "Tree/" << name1 << "-ML.tree";
		get_gene_root_number(strm1.str(), geneNum1, rootNum1);
		strm2  << "Tree/" << name2 << "-ML.tree";
		get_gene_root_number(strm2.str(), geneNum2, rootNum2);
				
		// -----------------------------
		//initialization
    string* C1_children = new string[rootNum1+1];
		string* C2_children = new string[rootNum2+1];
    string* C1_member = new string[rootNum1+1];
    string* C2_member = new string[rootNum2+1];    
		int* C1_size = new int[rootNum1+1];
		int* C2_size = new int[rootNum2+1];

		// ------------------------------
		// compute members and size of the clusters 
		strm2.str("");strm2.clear();
		strm2  << "Node/" << name1 << ".node";
		make_cluster_members(strm1.str(),strm2.str(),geneNum1,rootNum1,C1_size,C1_children,C1_member);
		strm1.str("");strm1.clear();
    strm2.str("");strm2.clear();
		strm1 << "Tree/"<< name2 << "-ML.tree";
    strm2  << "Node/" <<name2 << ".node";
		make_cluster_members(strm1.str(),strm2.str(),geneNum2,rootNum2,C2_size,C2_children,C2_member);		
    
		// ------------------------------ 
		// make the blast matrix
    
    double** blast = new double*[geneNum1+1];
    
    for (int c=1; c<=geneNum1; c++) {
        blast[c]=new double[geneNum2+1];
    }
    
    for (int c1=1; c1<=geneNum1; c1++) {
        for (int c2=1; c2<=geneNum2; c2++) {
            blast[c1][c2]=0;
        }
    }  

    input.open(similarity.c_str());
    getline(input, line);
    while (!input.eof()) {
				Stream1.clear();
        Stream1.str(line);
        Stream1 >> row;
        Stream1 >> col;
        Stream1 >> term3;
        blast[row][col] = term3;
        getline(input, line); 
    }
    input.close();

    // ------------------------------ 
    double** score = new double*[rootNum1+1];
    int** visit = new int*[rootNum1+1]; 
    double** sumArray = new double*[rootNum1+1];
		for (int c=1; c<=rootNum1; c++) {
        score[c]=new double[rootNum2+1];
    }
    
    for (int c1=1; c1<=rootNum1; c1++) {
        for (int c2=1; c2<=rootNum2; c2++) {
            score[c1][c2]=0;
        }
    }
    for (int c=1; c<=rootNum1; c++) {
        visit[c]=new int[rootNum2+1];
    }

    for (int c1=1; c1<=rootNum1; c1++) {
        for (int c2=1; c2<=rootNum2; c2++) {
            visit[c1][c2]=0;
        }
    }
    for (int c=1; c<=rootNum1; c++) {
        sumArray[c]=new double[rootNum2+1];
    }

    for (int c1=1; c1<=rootNum1; c1++) {
        for (int c2=1; c2<=rootNum2; c2++) {
            sumArray[c1][c2]=0;
        }
    }

    for (int c1=1; c1<=geneNum1; c1++) {
        for (int c2=1; c2<=geneNum2; c2++) {
            score[c1][c2]=blast[c1][c2];
						sumArray[c1][c2]=blast[c1][c2];
						visit[c1][c2]=1;
				}
    }

    // ------------------------------  
    for (int i=1; i<rootNum1; i++) {
        for (int j=1; j<rootNum2; j++) {
            if(visit[i][j]==0) {
                sum=0;
								Stream1.clear();
                Stream1.str(C1_children[i]);
                while(getline(Stream1,term2, ',')) {
										Stream2.clear();
										Stream2.str(C2_children[j]);
                    while (getline(Stream2,term4, ','))
                        sum = sum + sumArray[atoi(term2.c_str())][atoi(term4.c_str())];
                }
								score[i][j]=sum/(C1_size[i]*C2_size[j]);
								sumArray[i][j]=sum;
            		visit[i][j]=1;
						}
        }
    }
    // ------------------------------  
		output.open(outFile.c_str());
    for (int i=1; i<rootNum1; i++) {
        for (int j=1; j<rootNum2; j++) {
            if(score[i][j]!=0 && (i>geneNum1 || j>geneNum2))
							 output << i << "\t" << j  << "\t" <<score[i][j]<<"\n";									 
        }       
    }
		output.close();

		strm1.str("");strm1.clear();
		strm1 << name1 << ".cluster";
    output.open((strm1.str()).c_str());
    for (int i=1; i<rootNum1; i++) 
    		output << i << "\t" << C1_member[i]<<"\n";
    output.close();

    strm1.str("");strm1.clear();
    strm1 << name2 << ".cluster";
    output.open((strm1.str()).c_str());
    for (int i=1; i<rootNum2; i++) 
    		output << i << "\t" << C2_member[i]<<"\n";
    output.close();
    return 0;
}        
