#include<vector>
#include<iostream>
#include<algorithm>
#include<random>
#include<vector>
#include<cmath>

using namespace std;

template<typename T = bool>
struct data
{
    unsigned int x, y;
    T value;
    data(unsigned int _x, unsigned int _y, T val) {x=_x; y=_y; value=val;};
    data() {x=0; y=0; value=NULL;};
};


template<typename T = bool>
class SparseMatrix
{
public:

    SparseMatrix();
    SparseMatrix(int msize, bool symmetric);
    SparseMatrix(int type, int msize);

    void push_back(data<T> d); //Adds an element
    void erase(int index);

    double trace();

    //Multiplication by vectors or matrices
    vector<double> operator *(vector<double> &v);
    template<typename R>
    SparseMatrix<double> operator *(SparseMatrix<R> &s);
    data<T> operator [](int index);

    //Matrix exponentiation
    SparseMatrix<double> pow(int n);
    //Compute eigenvector
    vector<double> dom_eigen(double epsilon, int max_it);


    //Modify private variable msize
    void set_max_size(int msize);
    int get_max_size();
    int get_number_elements();

    static const int SM_DIAGONAL = 0;

    bool is_symmetric;


    vector<data<T>> m;

private:

    unsigned int m_dim;


    SparseMatrix<double> convert_double();

    void normalize_vector(vector<double> &v);
    double scalar_product(vector<double> &u, vector<double> &v);

};


template<typename T>
SparseMatrix<T>::SparseMatrix()
{
    ;
}



template<typename T>
SparseMatrix<T>::SparseMatrix(int msize, bool symmetric)
{
    m = vector<data<T>>();
    m_dim = msize;

    is_symmetric = symmetric;
}

template<typename T>
SparseMatrix<T>::SparseMatrix(int type, int msize)
{
    m = vector<data<T>>();
    m_dim = msize;

    if (type == SM_DIAGONAL)
    {
        int i;
        is_symmetric = true;
        if (typeid(T) == typeid(bool))  for (i=0; i < m_dim; i++)  m.push_back(data<bool>(i,i,true));
        else if (typeid(T) == typeid(double)) for (i=0; i < m_dim; i++)  m.push_back(data<double>(i,i,1.0));
        else cout << "ERROR [SparseMatrix]: please use bool or double as SparseMatrix type" << endl;
    }
}


template<typename T>
void SparseMatrix<T>::push_back(data<T> d)
{
    m.push_back(d);
}

template<typename T>
void SparseMatrix<T>::erase(int index)
{
    m.erase(m.begin() + index);
}

template<typename T>
void SparseMatrix<T>::set_max_size(int msize)
{
    m_dim = msize;
}
template<typename T>
int SparseMatrix<T>::get_max_size()
{
    return m_dim;
}
template<typename T>
int SparseMatrix<T>::get_number_elements()
{
    return m.size();
}

template<typename T>
SparseMatrix<double> SparseMatrix<T>::convert_double()
{
    int i;
    SparseMatrix<double> s(m_dim, is_symmetric);

    for (i=0; i < m.size(); i++)
    {
        s.push_back(data<double>(m[i].x, m[i].y, (double)m[i].value)); //Same info, but in double
    }

    return s;
}


template<typename T>
vector<double> SparseMatrix<T>::operator *(vector<double> &v)
{
    int i;
    vector<double> u = vector<double>(v.size(), 0.0);
    data<T> aux = data<T>();

    //Make the product
    if (is_symmetric)
    {
        for (i=0; i < m.size(); i++)
        {
            aux = m[i];
            u[aux.x] += aux.value * v[aux.y];
            //If it is symmetric, then we have to invert also this
            if (aux.x != aux.y)
            {
                u[aux.y] += aux.value * v[aux.x];
            }
        }
    }
    else
    {
        for (i=0; i < m.size(); i++)
        {
            aux = m[i];
            u[aux.x] += aux.value * v[aux.y];
        }
    }


    return u;
}



template<typename T> template <typename R>
SparseMatrix<double> SparseMatrix<T>::operator *(SparseMatrix<R> &s)
{
    int i,j;
    int m_index;
    int msize = m.size();
    int smsize = s.m.size();

    int aux_dim = max(msize, smsize); //The new matrix is not going to have more than this number of elements

    vector< vector<double> > proxy = vector< vector<double> >(aux_dim, vector<double>(aux_dim, 0.0));


    bool mdiag, smdiag;

    //New matrix. A*B is symmetric if [A,B] = 0. If we are computing A^n, this holds, so check it
    SparseMatrix<double> u = SparseMatrix<double>(m_dim, false);

    data<T> aux = data<T>();

    //Multiplication needs only msize elements!
    if (is_symmetric && s.is_symmetric)
    {
        for (i=0; i < msize; i++)
        {
            mdiag = m[i].x == m[i].y;

            for (j=0; j < smsize; j++)
            {
                smdiag = s.m[j].x == s.m[j].y;
                //Only these elements survive...
                if (mdiag)
                {
                    if (smdiag)
                    {
                        //If both are in the diagonal, then we only have to compare any two...
                        if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                    }
                    else
                    {
                        //If this is a diagonal element of my matrix, I have to compare one element in m with the two of s
                        if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                        if (m[i].y == s.m[j].y) proxy[m[i].x][s.m[j].x] += m[i].value * s.m[j].value;
                    }
                }
                else if (smdiag)
                {
                    //In this case we always have mdiag = false, so only s matrix is symmetric.
                    //Then compare one element in s with the two in m
                    if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                    if (m[i].x == s.m[j].x) proxy[m[i].y][s.m[j].y] += m[i].value * s.m[j].value;
                }
                else
                {
                    //Then we have to evaluate all the possible permutations
                    if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                    //If s matrix is symmetric, then swapping  x and y in m also works.
                    if (m[i].y == s.m[j].y) proxy[m[i].x][s.m[j].x] += m[i].value * s.m[j].value;
                    //If my matrix is symmetric, then swapping  x and y in m also works
                    if (m[i].x == s.m[j].x) proxy[m[i].y][s.m[j].y] += m[i].value * s.m[j].value;
                    //Swap on both!
                    if (m[i].x == s.m[j].y) proxy[m[i].y][s.m[j].x] += m[i].value * s.m[j].value;
                }

            }
        }
    }
    else if (is_symmetric && !s.is_symmetric)
    {
        for (i=0; i < msize; i++)
        {
            for (j=0; j < smsize; j++)
            {
                //Only these elements survive
                if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                //If my matrix is symmetric, then swapping  x and y in m also works
                if (m[i].x == s.m[j].x) proxy[m[i].y][s.m[j].y] += m[i].value * s.m[j].value;
            }
        }
    }
    else if (!is_symmetric && s.is_symmetric)
    {
        for (i=0; i < msize; i++)
        {
            for (j=0; j < smsize; j++)
            {
                //Only these elements survive
                if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
                //If s matrix is symmetric, then swapping  x and y in m also works
                else if (m[i].y == s.m[j].y)  proxy[m[i].x][s.m[j].x] += m[i].value * s.m[j].value;
            }
        }
    }
    else
    {
        for (i=0; i < msize; i++)
        {
            for (j=0; j < smsize; j++)
            {
                if (m[i].y == s.m[j].x) proxy[m[i].x][s.m[j].y] += m[i].value * s.m[j].value;
            }
        }
    }

    //Once we have computed positions, we add them to the matrix
    for (i=0; i < aux_dim; i++)
    {
        for (j=0; j < aux_dim; j++)
        {
            if (proxy[i][j] != 0.0)  u.push_back( data<double>(i, j, proxy[i][j]) );
        }
    }

    return u;
}



template<typename T>
double SparseMatrix<T>::trace()
{

    int i;
    double sum = 0.0;

    i = m.size() - 1;

    if (typeid(T) == typeid(bool))
    {
        //In this case, we know that the value is one, so we can sum directly over the condition
        while (i > 0)
        {
            sum += m[i].x == m[i].y;
            i--;
        }
    }
    else
    {
        //In other case, first check condition and then sum value...
        while (i > 0)
        {
            sum += m[i].x == m[i].y ? m[i].value : 0.0;
            i--;
        }
    }

    return sum;
}

template<typename T>
data<T> SparseMatrix<T>::operator [](int index)
{
    return m[index];
}

template<typename T>
SparseMatrix<double> SparseMatrix<T>::pow(int n)
{

    int i=0;
    SparseMatrix<double> s = *(this) * *(this);
    while (i < n-2) //So if n=3 we get this*this*this
    {
        s = *(this) * s;
        i++;
    }


    return s;

}

//Note: it return a vector double where the last element is the eigenvalue
template<typename T>
vector<double> SparseMatrix<T>::dom_eigen(double epsilon, int max_it)
{
    int i,j;

    double scp1, scp2; //Auxiliary variables to do scalar products fast


    //Init random initializers to get random vector
    random_device rnd_device;
    mt19937 gen(rnd_device());
    uniform_real_distribution<double> ran_u(0.0, 1.0);

    //Declare random vector of size m_dim
    vector<double> b = vector<double>(m_dim);
    vector<double> bm = vector<double>(m_dim);
    double eigen, old_eigen;
    double calc_error;

    //Make the vector completely random to avoid being orthogonal to eigenvector
    i = 0;
    while (i < m_dim)
    {
        b[i] = ran_u(gen);
        i++;
    }
    normalize_vector(b); //Normalize the stuff

    //Very different values to avoid not-entering in loop
    eigen = 0.0;
    old_eigen = 1000.0;
    calc_error = 1.0;
    i = 0;
    while (calc_error > epsilon && i < max_it)
    {
        bm = *(this) * b;

        old_eigen = eigen; //Update value
        j = 0;
        //Compute eigenvalue using Rayleigh's quotient.
        //Both scalar products are evaluated at the same time in order to be fast
        scp1 = scp2 = 0.0;
        while (j < b.size())
        {
            scp1 += b[j] * bm[j];
            scp2 += b[j] * b[j];
            j++;
        }
        eigen = scp1 / scp2; //Get new eigenvalue


        normalize_vector(bm); //Normalize this stuff

        b = bm; //Update

        if (old_eigen != 0) calc_error = abs((eigen - old_eigen) / old_eigen);

        i++;
    }

    if (i >= max_it) cout << "WARNING [SparseMatrix]: eigenvalue max number of iterations reached. Computation probably did NOT converge." << endl;

    bm.push_back(eigen); //Add the eigenvalue as the last element of the vector

    return bm;
}

template<typename T>
void SparseMatrix<T>::normalize_vector(vector<double> &v)
{
    int i;
    double sq = 0.0;

    //Compute sum of squares
    i = 0;
    while(i < v.size())
    {
        sq += v[i]*v[i];
        i++;
    }
    sq = 1.0/sqrt(sq); //Inverse square root directly for speed

    //Then apply it
    i = 0;
    while(i < v.size())
    {
        v[i] *= sq;
        i++;
    }
    return;
}

//Computes scalar product between two vectors
template<typename T>
double SparseMatrix<T>::scalar_product(vector<double> &u, vector<double> &v)
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

