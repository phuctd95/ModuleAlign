#include <bits/stdc++.h>

using namespace std;

string s[4] = {"celeg","scere","dmela","hsapi"};

int main()
{
    freopen("run.sh","w",stdout);

        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j)
                for (double a = 0.3; a <= 0.7; a+= 0.1)
                cout << "./script.sh "<< s[i] << " "<< s[j] << " "<<a << endl;
    return 0;
}
