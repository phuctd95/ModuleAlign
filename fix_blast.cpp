#include <bits/stdc++.h>

using namespace std;

string in_file;
vector <int> n1;
vector <int> n2;
vector <double> n3;
int main()
{
    cin >> in_file;
    in_file += ".blast";
    freopen(in_file.c_str(),"r",stdin);
    int u,v;
    double s;
    while (cin >> u >> v >> s)
    {
        n1.push_back(u);
        n2.push_back(v);
        n3.push_back(s);
    }
    freopen(in_file.c_str(),"w",stdout);
    for (int i = 0; i < n1.size(); ++i)
        cout << n1[i] << "\t" << n2[i] << "\t" <<n3[i] << endl;
    return 0;
}
