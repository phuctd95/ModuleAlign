#include <bits/stdc++.h>

using namespace std;

string in_file;
vector <int> n1;
vector <int> n2;
vector <bool> check;
int main()
{
    freopen("scere.net","r",stdin);
    int u,v;
    while (cin >> u >> v)
    {
        n1.push_back(u);
        n2.push_back(v);
        while (check.size() <= max(u,v)) check.push_back(true);
        check[u] = false;
        check[v] = false;
    }
    cerr << check.size()<< endl;
    freopen("scere.net","w",stdout);
    for (int i = 0; i < n1.size(); ++i)
        cout << n1[i] << " " << n2[i] << endl;
    for (int i = 0; i < check.size(); ++i)
        if (check[i])
        {
            cout << i << " " << i << endl;
            cerr << i << endl;
        }
    return 0;
}
