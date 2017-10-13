#include<vector>
#include<iostream>
#include<algorithm>
#include<vector>
#include<cmath>

using namespace std;

template<class T = bool>
struct data
{
    unsigned int x, y;
    T value;
    data(unsigned int _x, unsigned int _y, T val) {x=_x; y=_y; value=val;};
    data() {x=0; y=0; value=NULL;};
};


template<class T = bool>
class SparseMatrix
{
public:

    SparseMatrix(int msize);
    SparseMatrix(int type, int msize);

    void push_back(data<T> d); //Adds an element

    vector<double> dom_eigen(double epsilon);


    //Multiplication by vectors or matrices
    vector<double> operator *(vector<double> &v);
    SparseMatrix<double> operator *(SparseMatrix<T> &s);
    //Matrix exponentiation
    SparseMatrix<double> pow(int n);

    //Modify private variable msize
    void set_size(int msize);
    int get_size();

    static const int SM_DIAGONAL = 0;

    int EIGEN_MAX_ITS = 20;

    vector<data<T>> m;

private:

    unsigned int m_dim;

    SparseMatrix<double> convert_double();

    void normalize_vector(vector<double> &v);
    double scalar_product(vector<double> &u, vector<double> &v);

};



template<class T>
SparseMatrix<T>::SparseMatrix(int msize)
{
    m = vector<data<T>>();
    m_dim = msize;
}

template<class T>
SparseMatrix<T>::SparseMatrix(int type, int msize)
{
    m = vector<data<T>>();
    m_dim = msize;

    if (type == SM_DIAGONAL)
    {
        int i;
        if (typeid(T) == typeid(bool))  for (i=0; i < m_dim; i++)  m.push_back(data<bool>(i,i,true));
        else if (typeid(T) == typeid(double)) for (i=0; i < m_dim; i++)  m.push_back(data<double>(i,i,1.0));
        else cout << "ERROR_TYPE: please use bool or double as SparseMatrix type" << endl;
    }
}


template<class T>
void SparseMatrix<T>::push_back(data<T> d)
{
    m.push_back(d);
}

template<class T>
void SparseMatrix<T>::set_size(int msize)
{
    m_dim = msize;
}
template<class T>
int SparseMatrix<T>::get_size()
{
    return m_dim;
}

template<class T>
SparseMatrix<double> SparseMatrix<T>::convert_double()
{
    int i;
    SparseMatrix<double> s(m_dim);

    for (i=0; i < m.size(); i++)
    {
        s.push_back(data<double>(m[i].x, m[i].y, (double)m[i].value)); //Same info, but in double
    }

    return s;
}


template<class T>
vector<double> SparseMatrix<T>::operator *(vector<double> &v)
{
    int i;
    vector<double> u = vector<double>(v.size(), 0.0);
    data<T> aux = data<T>();

    //Make the product
    for (i=0; i < m.size(); i++)
    {
        aux = m[i];
        u[aux.x] += aux.value * v[aux.y];
    }

    return u;
}



template<class T>
SparseMatrix<double> SparseMatrix<T>::operator *(SparseMatrix<T> &s)
{
    int i,j;
    int msize = m.size();

    vector<double> proxy = vector<double>(msize*msize, 0.0); //To sum the different columns
    vector<int> cx = vector<int>(msize, 0); //To store where the thing goes in the new matrix
    vector<int> cy = vector<int>(msize, 0);

    SparseMatrix<double> u = SparseMatrix<double>(m_dim); //New matrix

    data<T> aux = data<T>();

    //Multiplication needs only msize elements!
    for (i=0; i < msize; i++)
    {
        for (j=0; j < msize; j++)
        {
            //Only these elements survive
            if (m[i].y == s.m[j].x)
            {
                cx[i] = m[i].x;
                cy[j] = s.m[j].y;
                proxy[j + i*msize] += m[i].value * s.m[j].value; //Add therm
            }
        }
    }

    //Once we have computed positions, we add them to the matrix
    for (i=0; i < msize; i++)
    {
        for (j=0; j < msize; j++)
        {
            if ( proxy[j + i*msize] != 0.0)
            {
                u.push_back( data<double>(cx[i],cy[j],proxy[j + i*msize]) );
            }
        }
    }
    return u;
}

template<class T>
SparseMatrix<double> SparseMatrix<T>::pow(int n)
{
    if (n == 2)
    {
        return *(this) * *(this);
    }
    else
    {
        int i=0;
        SparseMatrix<double> this_double = this->convert_double();
        SparseMatrix<double> s = this_double * this_double;

        while (i < n-2) //So if n=3 we get this*this*this
        {
            s = this_double * s;
            i++;
        }
        return s;
    }

}

template<class T>
vector<double> SparseMatrix<T>::dom_eigen(double epsilon)
{
    int i;

    //Init random initializers to get random vector
    random_device rnd_device;
    mt19937 gen(rnd_device());
    uniform_real_distribution<double> ran_u(0.0, 1.0);

    vector<double> b = vector<double>(m_dim);
    ector<double> bm = vector<double>(m_dim);
    double eigen, old_eigen;

    //Make the vector completely random
    i = 0;
    while (i < m_dim)
    {
        b[i] = ran_u(gen);
        i++;
    }

    normalize_vector(b);

    //Very different values to avoid not-entering in loop
    eigen = 0.0;
    old_eigen = 1000.0;
    while (abs(eigen-old_eigen) > epsilon)
    {
        bm = *(this) * b;
        normalize_vector(bm);
        old_eigen = eigen;
        eigen = scalar_product(b, bm);
    }

    bm.push_back(eigen);

    return bm;
}

template<class T>
void SparseMatrix::normalize_vector(vector<double> &v)
{
    int i;
    double sq = 0.0;

    i = 0;
    while(i < v.size())
    {
        sq += v[i]*v[i];
        i++;
    }
    sq = 1.0/sqrt(sq);

    i = 0;
    while(i < v.size())
    {
        v[i] *= sq;
        i++;
    }
    return;
}


double scalar_product(vector<double> &u, vector<double> &v)
{
    int i;
    double p = 0.0;

    i = 0;
    while(i < v.size())
    {
        p += v[i] * u[i];
        i++;
    }
    return p;
}

