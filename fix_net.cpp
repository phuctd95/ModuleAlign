#include <bits/stdc++.h>

using namespace std;

string in_file;
vector <int> n1;
vector <int> n2;
vector <bool> check;
int main()
{
	cin >> in_file;
    freopen((in_file+".net").c_str(),"r",stdin);
    int u,v;
    while (cin >> u >> v)
    {
        if (u == v) continue;
        n1.push_back(u);
        n2.push_back(v);
        while (check.size() <= max(u,v)) check.push_back(true);
        check[u] = false;
        check[v] = false;
    }
    cerr << check.size()<< endl;
    freopen((in_file+".net").c_str(),"w",stdout);
    for (int i = 0; i < n1.size(); ++i)
        cout << n1[i] << "\t" << n2[i] << endl;
    for (int i = 0; i < check.size(); ++i)
        if (check[i])
        {
            cout << i << "\t" << i-1 << endl;
            cerr << i << endl;
        }
    return 0;
}
