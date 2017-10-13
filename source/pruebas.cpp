#include<iostream>
#include<cstdlib>
#include<vector>

#include "SparseMatrix.cpp"

int main(void)
{
    int i;

    SparseMatrix<bool> s = SparseMatrix<bool>(4);
    SparseMatrix<bool> t = SparseMatrix<bool>(4);

    s.push_back(data<bool>(0,1,true));
    s.push_back(data<bool>(0,2,true));
    s.push_back(data<bool>(0,3,true));
    s.push_back(data<bool>(1,0,true));
    s.push_back(data<bool>(1,3,true));
    s.push_back(data<bool>(2,2,true));
    s.push_back(data<bool>(3,1,true));

    t.push_back(data<bool>(0,0,true));
    t.push_back(data<bool>(1,1,true));
    t.push_back(data<bool>(2,2,true));
    t.push_back(data<bool>(3,3,true));


    /*vector<double> a = vector<double>(4,1);
    vector<double> b = s*a;

    for (i=0; i < b.size(); i++) cout << b[i] << " ";*/

    //SparseMatrix<double> u = s.pow(3);

    //cout << u.m.size() << endl;
    //for (i=0; i < u.m.size(); i++) cout << u.m[i].x << " " << u.m[i].y << " " << u.m[i].value << endl;

    vector<double> eigens = s.dom_eigen(0.001);

    cout << eigens[eigens.size() - 1] << endl;

    return 0;
}
