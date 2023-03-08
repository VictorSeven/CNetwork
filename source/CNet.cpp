#include "CNet_Dir.cpp"

/* =======================================================================================================

Victor Buendía's CNetwork Class

This software is Open Source and redistributed with a MIT LICENSE (Read LICENSE file for more details).
Please follow license terms when using this code.

========================================================================================================== */

// ========================================================================================================
// ========================================================================================================
// ========================================================================================================


/** \brief CNetwork base class
*
*   CNetwork is the core class for a weighted, undirected network.
*   Needs two template arguments: class associated to nodes and links
*
*/
template <class T = bool, class B = bool>
class CNetwork: public DirectedCNetwork<T,B>
{
    public:


        void add_nodes(int n);
        bool remove_node(int index);
        void add_link(int from, int to);
        void add_link(int from, int to, B w);
        bool remove_link(int from, int to);
        void remove_link(int index_link);


        double mean_degree() const;


        double clustering_coef(int node_index) const;
        double mean_clustering_coef() const;


        void degree_distribution(vector<int> &distribution, bool normalized = false) const;
        void degree_correlation(vector<int> &distribution, vector<double> &correlation, bool normalized = false) const;


        int in_degree(const int node_index) const;
        int out_degree(const int node_index) const;
        int degree(const int node_index) const;


        vector<unsigned int> get_neighs_out(const int node_index) const;
        vector<unsigned int> get_neighs_in(const int node_index) const;
        int get_out(const int node_index, const int k) const;
        int get_in(const int node_index, const int k) const;

        int get_link_index(const int from, const int to) const;


        vector<unsigned int> get_neighs(const int node_index) const;
        int get_neigh_at(const int node_index, const int k) const;

        void create_2d_lattice(const int L, bool eight=false, bool periodic=true);
        void create_erdos_renyi(int nodes, double mean_k, unsigned int random_seed=123456789, unsigned int n0=0, unsigned int nf=0);
        void create_albert_barabasi(int n, int m0, int m, unsigned int random_seed = 123456789);

        bool read_mtx(string filename);

        CNetwork();
        CNetwork(int max_size);
        CNetwork(vector<CNetwork<T,B> > &submodules, const bool freemem=true);

        void clear_network();


};

using CNb = CNetwork<bool, bool>;
using CNi = CNetwork<int, bool>;
using CNd = CNetwork<double, bool>;

using WCN = CNetwork<void, double>;
using WCNb = CNetwork<bool, double>;
using WCNi = CNetwork<int, double>;
using WCNd = CNetwork<double, double>;


// ========================================================================================================
// ========================================================================================================
// ========================================================================================================


/** \brief CNetwork void constructor
*
* Creates a new CNetwork with no size.
*/
template <class T, typename B>
CNetwork<T,B>::CNetwork() : DirectedCNetwork<T,B>()
{
    this->directed = false;
    return;
}


/** \brief CNetwork standard constructor
*  \param max_size: maximum size of the network
*
* Creates a new CNetwork with a limit to the number of nodes. The maximum is fixed and
* the memory for the nodes is not allocated.
*/
template <class T, typename B>
CNetwork<T,B>::CNetwork(int max_size) : DirectedCNetwork<T,B>(max_size)
{
    this->directed = false;
    return;
}


/** \brief DirectedCNetwork constructor via submodules
*  \param submodules: a list of networks to build this network from
*
* Creates a new DirectedCNetwork which contains all the previous submodules, generating
* a graph with disconnected sub-graphs.
*/
template <class T, typename B>
CNetwork<T,B>::CNetwork(vector<CNetwork<T,B> > &submodules, const bool freemem) 
{
    //Side node: implementation cannot be done just calling the parent method first
    //and then setting directed false, since we have to use the add_link function, which needs
    //to have directed false from the beginning.

    int i,j,k;
    int size_added;
    this->directed = false; 

    //Get size and clear network
    this->max_net_size = 0;
    for (i=0; i < submodules.size(); i++) this->max_net_size += submodules[i].get_node_count();
    this->clear_network();

    //Add all the nodes we will need
    this->add_nodes(this->max_net_size);

    size_added = 0;
    for (i=0; i < submodules.size(); i++) 
    {
        //Generate network links using the information from the neighbours
        for (j=0; j < submodules[i].get_node_count(); j++)
        {
            for (k=0; k < submodules[i].out_degree(j); k++)
            {
                this->add_link(j + size_added, submodules[i].get_out(j,k) + size_added);
            }
        }
        size_added += submodules[i].get_node_count();


        //Free submodule memory to avoid high RAM allocations
        if (freemem) submodules[i].clear_network();
    }
    return;
}



/** \brief Delete all the data stored by the network
*  \param max_size: maximum size of the network
*
* Delete everything stored by the network.
*/
template <class T, typename B>
void CNetwork<T,B>::clear_network()
{
    this->current_size = 0; //Init size to 0
    this->link_count = 0; //Init link count

    //Create vectors in order to add things
    this->adjm = SparseMatrix<B>(this->max_net_size, true);

    this->neighs = vector< vector<unsigned int> >(0, vector<unsigned int>(0));


    this->prop_d = map<string, vector<double> >();
    this->prop_i = map<string, vector<int> >();
    this->prop_b = map<string, vector<bool> >();
    this->prop_s = map<string, vector<string> >();

    this->value = vector<T>();


    return;
}


// ========================================================================================================
// ========================================================================================================
// ========================================================================================================


/** \brief Add new nodes to the network
*  \param n: number of nodes to add
*
* Add new nodes to the network
*/
template <class T, typename B>
void CNetwork<T,B>::add_nodes(int n)
{
    int old_size = this->current_size; //Store old size
    int next_size = old_size + n; //New size of the network

    //Set the size to the maximum if it exceeds it
    this->current_size = next_size > this->max_net_size ? this->max_net_size : next_size;

    this->value.resize(this->current_size); //Increase the number of value.size without adding any element to it

    for (int i = old_size; i < this->current_size; i++)
    {
        this->neighs.push_back(vector<unsigned int>()); //Add a new container for neighbours
    }
    return;
}




/** \brief Remove a node from the network
*  \param index: index of the node to remove
*
* Remove the selected node from the network. All links related to this node will be erased.
*/
template <class T, typename B>
bool CNetwork<T,B>::remove_node(int index)
{
    int i,j,k;

    if (index >= 0 and index < this->current_size)
    {
        this->neighs.erase(this->neighs.begin() + index); //Delete its entry in the neighbours
        this->value.erase(this->value.begin() + index); //Delete its value

        //Re-scale all the neighbours index to fit the net, and also erase
        //the entry of the removed node as neighbour of others
        for (i=0; i < this->neighs.size(); i++)
        {
            for (j=0; j < degree(i); j++)
            {
                k = get_neigh_at(i,j);
                if (k == index) //If it is the one, I want to erase, remove
                {
                    this->neighs[i].erase(this->neighs[i].begin() + j);
                }
                else if (k > index) //If it is higher, move 1 to left
                {
                    this->neighs[i][j] -= 1;
                }
            }
        }

        vector<int> who_to_erase = vector<int>(); //Track who you want to delete

        for (i=0; i < this->link_count; i++)
        {
            //Separe cases from-to in order to correctly reduce degree
            if (this->adjm[i].x == index || this->adjm[i].y == index)
            {
                //Our node is detected, erase the entry to that node
                who_to_erase.push_back(i);
            }

            //Also reduce index of nodes higher than our
            if (this->adjm[i].x > index)
            {
                this->adjm[i].x -= 1;
            }
            if (this->adjm[i].y > index)
            {
                this->adjm[i].y -= 1;
            }
        }

        //Perform the erase process
        for (i=0; i < who_to_erase.size(); i++)
        {
            this->adjm.erase(who_to_erase[i]); //Delete all links
            //Index was recorded before, we do 2*who in order to get as even-odd again
        }

        this->link_count -= who_to_erase.size() / 2; //Reduce the number of links
        this->current_size -= 1; //Reduce the current size of the network in one unit

        return true;

    }
    else return false;
}


/** \brief Add a link to the network
*  \param from: index of the origin node
*  \param to: index of the target node
*
* Adds a link from the nodes from and to
*/
template <class T, typename B>
void CNetwork<T,B>::add_link(int from, int to)
{
    this->adjm.push_back(data<bool>(from, to, true)); //Assume this method is for bools

    this->neighs[from].push_back(to); //Add the node to the neighbours
    this->neighs[to].push_back(from); //And do it in the other sense also

    this->link_count += 1; //Create one link more
    return;
}


/** \brief Add a link to the network
*  \param from: index of the origin node
*  \param to: index of the target node
*  \param w: weight of the link.
*
* Add a weighted link between nodes from and to
*/
template <class T, typename B>
void CNetwork<T,B>::add_link(int from, int to, B w)
{
    this->adjm.push_back(data<B>(from, to, w)); //Assume this method is for weighted things

    this->neighs[from].push_back(to); //Add the node to the neighbours
    this->neighs[to].push_back(from); //And do it in the other sense also

    this->link_count += 1; //Create one link more

    return;
}


/** \brief Remove a link from the network
*  \param from: Index of the origin node
*  \param to: Index of the target node
*  \return false if there is no link between from and to.
*
* Remove the a link between nodes from and to, if it exist.
*/
template <class T, typename B>
bool CNetwork<T,B>::remove_link(int from, int to)
{
    //Reduce the degree of the nodes
    auto index_it = find(this->neighs[from].begin(), this->neighs[from].end(), to); //Relative index of TO in terms of FROM
    int index_neigh = distance(this->neighs[from].begin(), index_it); //Get the relative index as an int


    if (index_neigh >= 0 and index_neigh < this->neighs[from].size())
    {
        //int index_neigh = distance(neighs[from].begin(), index_it); //Get the relative index as an int
        int index_link = get_link_index(from, to);

        this->adjm.erase(index_link); //Delete the link

        this->link_count -= 1;


        this->neighs[from].erase(this->neighs[from].begin()+index_neigh); //Erase the node it pointed in the neighbours list

        //Do the same process, but now in the other node, the to one.
        //Since this node was in the neigh list of FROM, we know that FROM has to be in the list of TO
        //That's why we don't check again if index_it < neighs end
        index_it = find(this->neighs[to].begin(), this->neighs[to].end(), from);
        index_neigh = distance(this->neighs[to].begin(), index_it);

        this->neighs[to].erase(this->neighs[to].begin()+index_neigh);

        return true;
    }
    else return false;

}


/** \brief Remove a link from the network
*  \param index: index of the link to erase
*  \return false if there is no link between from and to.
*
* Remove the a link between nodes from and to, if it exist.
*/
template <class T, typename B>
void CNetwork<T,B>::remove_link(int index_link)
{
    unsigned int from, to;
    vector<int> aux;

    //Get who are from and to nodes
    aux = this->get_link(index_link);
    from = aux[0];
    to = aux[1];

    this->adjm.erase(index_link); //Delete the link
    this->link_count -= 1;


    auto index_it = find(this->neighs[from].begin(), this->neighs[from].end(), to); //Relative index of TO in terms of FROM
    int index_neigh = distance(this->neighs[from].begin(), index_it); //Get the relative index as an int

    this->neighs[from].erase(this->neighs[from].begin()+index_neigh); //Erase the node it pointed in the neighbours list

    //Do the same process, but now in the other node, the to one.
    //Since this node was in the neigh list of FROM, we know that FROM has to be in the pointing_in list of TO
    //That's why we don't check again if index_it < neighs end
    index_it = find(this->neighs[to].begin(), this->neighs[to].end(), from); //find(neighs[to].begin(), neighs[to].end(), from);
    index_neigh = distance(this->neighs[to].begin(), index_it);

    this->neighs[to].erase(this->neighs[to].begin()+index_neigh);

    return;

}


// ========================================================================================================
// ========================================================================================================
// ========================================================================================================

/** \brief Compute mean degree of the network
*  \return Mean degree
*
* Computes the mean degree of the network
*/
template <class T, typename B>
double CNetwork<T,B>::mean_degree() const
{
    int i;
    double sum = 0.0;
    for (i=0; i < this->current_size; i++) sum += degree(i);

    return sum / (1.0 * this->current_size);
}


/** \brief Computes the clustering coefficient of a node
*  \param index: target node
*  \return clustering coefficient of target node
*
* Compute a clustering coefficient of a target node.
*/
template <class T, typename B>
double CNetwork<T,B>::clustering_coef(int node_index) const
{
    if (degree(node_index) > 1) //If we have more than one neighbour...
    {
        int i,j;
        int counter; //Count of pairs


        vector<unsigned int> nodes_neigh = get_neighs(node_index);
        vector<unsigned int> nodes_check;

        counter = 0;
        for (i=0; i < nodes_neigh.size(); i++) //Get neighbours of our node
        {
            nodes_check = this->get_neighs(nodes_neigh[i]); //Get neighbours of node i

            //For the next nodes, (start in j=i+1 to avoid double-count a pair)
            for (j=i+1; j < nodes_neigh.size(); j++)
            {
                //Use the  find function to see if this node is connected to any other neighbours of our node.
                //In that case, increase the counter
                if (find(nodes_check.begin(), nodes_check.end(), nodes_neigh[j]) != nodes_check.end()) counter += 1;
            }
        }

        return 2.0 * counter / (degree(node_index) * (degree(node_index) - 1)); //Finish computation and return clustering coefficient
    }
    else //... in other case, we cannot have common neighbours
    {
        return 0.0;
    }

}

/** \brief Computes verage clustering coefficient
*  \return clustering coefficient of the network
*
* Computes the clustering coefficient of each element, and takes the average.
* It is not the same as the one computed counting triangles.
*/
template <class T, typename B>
double CNetwork<T,B>::mean_clustering_coef() const
{
    int i;
    double sum = 0.0; //Get the sum,

    //Sum over the network
    for (i=0; i < this->current_size; i++)
    {
        sum += clustering_coef(i);
    }
    //Divide by current size
    return sum / (this->current_size * 1.0);
}

/** \brief Computes the degree distribution of the network
*  \param[out] distribution: index j contains number of nodes with degree j. It is the degree distribution
*  \param normalized: optional. If set to true, returns a normalized distribution. False by default.
*
* Compute the degree distribution of the network. If you also need the correlations, please use instead degree_correlation.
*/
template <class T, typename B>
void CNetwork<T,B>::degree_distribution(vector<int> &distribution, bool normalized) const
{
    int i;
    distribution = vector<int>(this->current_size, 0);

    //Select in, out, or full degree distribution and compute it:
    for (i=0; i < this->current_size; i++) distribution[degree(i)] += 1;

    //Erase the 0s at the end of the array.
    i = this->current_size - 1; //Start counter
    while (distribution[i] == 0)
    {
        //distribution.erase(distribution.begin() + i);
        distribution.pop_back();
        i -= 1;
    }

    //Normalize the distribution if it has been indicated
    if (normalized)
    {
        for (i=0; i < distribution.size(); i++)
        {
            distribution[i] /= 1.0*this->link_count;
        }
    }

    return;
}


/** \brief Computes the degree distribution and correlations
*  \param[out] distribution: index j contains number of nodes with degree j. It is the degree distribution
*  \param[out] correlation: index j contains average degree of neighbours of a node with degree j
*  \param normalized: optional. If set to true, returns a normalized distribution. False by default.
*
* Computes the average number of neighbours that a node with degree j has. It also computes and stores the degree distribution,
* since both quantities are usually needed.
*/
template <class T, typename B>
void CNetwork<T,B>::degree_correlation(vector<int> &distribution, vector<double> &correlation, bool normalized) const
{
    DirectedCNetwork<T,B>::degree_correlation(distribution, this->TOTAL_DEGREE, normalized);
    return;
}



// ========================================================================================================
// ========================================================================================================
// ========================================================================================================


/** \brief Creates a 2d lattice.
*
* This method creates a regular 2d lattice with 4 neighbours. It is useful to
* compare networks to lattices.
*/
template<class T, typename B>
void CNetwork<T,B>::create_2d_lattice(const int L, bool eight, bool periodic)
{
    int i,j,x,y, xm, ym, xp, yp;
    int u,d,l,r;

    add_nodes(L*L); //Create the nodes

    //To set coordinates in the desired order
    const int UP = 2;
    const int DOWN = 3;
    const int LEFT = 1;
    const int RIGHT = 0;


    if (!eight)
    {
        vector<int> coor(4);

        if (periodic)
        {
            //Store links for each node
            for (i=0; i < this->current_size; i++)
            {
                //Get x,y coords of node in this lattice
                x = i % L;
                y = i / L;

                //Coordinates of right, left, up and down, using that
                //index = x + L * y
                coor[RIGHT] = (x+1)%L + L*y;
                coor[LEFT] = x > 0 ? x-1 + L*y : L-1 + L*y;
                coor[UP] = x + L * ((y+1) % L);
                coor[DOWN] = y > 0 ? x + L*(y-1) : x + L*(L-1);


                //Create the links from this node. This differs from add_link in the fact we
                //do NOT create explicitly the link "to->from". It is created later when the "to" node is visited.
                for (j=0; j < 4; j++)
                {
                    this->adjm.push_back(data<bool>(i, coor[j], true));
                    this->neighs[i].push_back(coor[j]);
                }
            }

            this->link_count = 2*this->current_size; //Set manually the number of links
        }
        else
        {
            for (x = 0; x < L-1; x++)
            {
                for (y=0; y < L-1; y++)
                {
                    this->adjm.push_back(data<bool>(x+y*L, x+1+y*L, true));
                    this->neighs[x+y*L].push_back( x+1+y*L);
                    this->adjm.push_back(data<bool>(x+1+y*L, x+y*L, true));
                    this->neighs[x+1+y*L].push_back( x+y*L);

                    this->adjm.push_back(data<bool>(x+y*L, x+(y+1)*L, true));
                    this->neighs[x+y*L].push_back(x+(y+1)*L);
                    this->adjm.push_back(data<bool>(x+(y+1)*L, x+y*L, true));
                    this->neighs[x+(y+1)*L].push_back(x+y*L);
                }
            }

            for (i=0; i < L-1; i++)
            {
                this->adjm.push_back(data<bool>(L-1+i*L, L-1+(i+1)*L, true));
                this->neighs[L-1+i*L].push_back(L-1+(i+1)*L);
                this->adjm.push_back(data<bool>(L-1+(i+1)*L, L-1+i*L, true));
                this->neighs[L-1+(i+1)*L].push_back(L-1+i*L);

                this->adjm.push_back(data<bool>(i+(L-1)*L, i+1+(L-1)*L, true));
                this->neighs[i+(L-1)*L].push_back(i+1+(L-1)*L);
                this->adjm.push_back(data<bool>(i+1+(L-1)*L, i+(L-1)*L, true));
                this->neighs[i+1+(L-1)*L].push_back(i+(L-1)*L);

            }
            this->link_count = 2*L*(L-1); //Set manually the number of links
        }

    }
    else //8-neighbour lattice
    {
        vector<int> coor(8);

        const int UPRI = 4;
        const int UPLE = 5;
        const int DORI = 6;
        const int DOLE = 7;


        //Store links for each node
        for (i=0; i < this->current_size; i++)
        {
            //Get x,y coords of node in this lattice
            x = i % L;
            y = i / L;

            //Coordinates of right, left, up and down, using that
            //index = x + L * y
            xm = x > 0 ? x-1 : L-1;
            ym = y > 0 ? y-1 : L-1;
            xp = (x+1) % L;
            yp = (y+1) % L;

            coor[RIGHT] = xp + L*y;
            coor[LEFT] = xm + L*y;
            coor[UP] = x + L*yp;
            coor[DOWN] = x + L*ym;

            coor[UPRI] = xp+ L*yp;
            coor[UPLE] =  xm + L*yp;
            coor[DORI] = xp + L*ym;
            coor[DOLE] = xm + L*ym;


            //Create the links from this node. This differs from add_link in the fact we
            //do NOT create explicitly the link "to->from". It is created later when the "to" node is visited.
            for (j=0; j < 8; j++)
            {
                this->adjm.push_back(data<bool>(i, coor[j], true));
                this->neighs[i].push_back(coor[j]);
            }
        }

        this->link_count = 4*this->current_size; //Set manually the number of links
    }



    return;
}


/** \brief Generates a random Erdos-Renyi network
*  \param nodes: nodes of the network
*  \param mean_k: average degree
*  \param random_seed: optional, default 123456789. Same seed gives the same network.
*
* Generates an Erdos-Renyi network. The random seed should be specified for obtaining different networks each iteration.
*/
template <class T, typename B>
void CNetwork<T,B>::create_erdos_renyi(int n, double mean_k, unsigned int random_seed, unsigned int n0, unsigned int nf)
{
    //TODO: make faster. MAKE SAFER

    int i,j,k;
    double p = mean_k / (n - 1.0);


    //A negative value for nf (default) gives the entire network
    //Any nf > 0 will be used to crerate a network for the subset of nodes [n0, nf]
    nf = nf == 0 ? this->max_net_size : nf;

    mt19937 gen(random_seed);; //Create the generator
    uniform_real_distribution<double> ran_u(0.0,1.0); //Uniform number distribution
    uniform_int_distribution<int>  index(n0,nf-1); //-1 because closed interval for ints

    add_nodes(nf-n0); //Create the nodes

    //For the n (n-1) / 2 pairs, link them with probability p
    for (i=n0; i < nf; i++)
    {
        for (j=i+1; j < nf; j++)
        {
            //With probability p, add link. 
            if (ran_u(gen) <= p) this->add_link(i, j);
        }
    }

    return;

}

/** \brief Generates an Albert-Barabasi network
*  \param n: nodes of the network
*  \param m0: initial number of fully-connected nodes
*  \param m: new links added for each node
*  \param random_seed: optional, default 123456789. Same seed gives the same network.
*
* Generates an Albert-Barabasi network based in the algorithm given by Newman. The random seed should be
* specified for obtaining different networks each iteration.
*/
template <class T, typename B>
void CNetwork<T,B>::create_albert_barabasi(int n, int m0, int m, unsigned int random_seed)
{
    int i,j,k,l;
    double r;

    mt19937 gen(random_seed);; //Create the generator
    uniform_real_distribution<double> random(0.0,1.0); //Uniform number distribution
    uniform_int_distribution<int>  index(0,1); //Useful to select indices

    int index_add;
    vector<int> yet_linked(m-1, -1); //Stores with which nodes I have visited. I only need to remember m-1 -I don't have to store the last
    bool found;

    //Create fully connected network with m0 nodes
    this->add_nodes(n);
    for (i=0; i < m0; i++)
    {
        for (j=i+1; j < m0; j++)
        {
            this->add_link(i,j);
        }
    }

    int nodes_in_net = m0;

    //Add then more nodes...
    for (i = m0; i < n; i++)
    {
        nodes_in_net++;

        yet_linked = vector<int>(m, -1);
        k = 0;
        //For every link we want to do,
        for (j=0; j < m; j++)
        {
            //Generate a random number
            if (random(gen) <= 0.5)
            {
                //With half probability, add it to highly connected node, selecting randomly an edge.
                index = uniform_int_distribution<int>(0, this->link_count-1);//-1 because interval is closed
                index_add = this->adjm[index(gen)].y;
            }
            else
            {
                //If not, connect to a random node
                index = uniform_int_distribution<int>(0, nodes_in_net-2); //We don't want current_size-1 which is the actual node, so up to -2
                index_add = index(gen);
            }

            //Check there that we don't have still a link between the two...
            found = false;
            l = 0;
            while (l < m && !found && yet_linked[l] != -1) 
            {
                found = index_add == yet_linked[l];
                l++;
            }
            //If there is no previous link between both, then add it.
            if (!found)
            {
                this->add_link(nodes_in_net-1, index_add); //Add the link
                yet_linked[k] = index_add; //Say we explored it
                k++; //To fill next vector element
            }

        }
    }

    return;
}

/** \brief Read network data from plain MTX format
*  \param filename: name of the file (including extension)
*
* The MTX format is defined by a simple plaintext representation of the adjancency matrix.
* This function reads any MTX-like format
*/
template <class T, typename B>
bool CNetwork<T,B>::read_mtx(string filename)
{
    //Destroy this object and create new network
    clear_network();

    //To see if file was correctly opened
    bool file_opened_ok;

    bool read_header = false; //To see if we have read the dim1xdim2 links line
    string line; //Store a line
    ifstream input; //File
    int from, to, w; //Auxiliary

    //Open the file and checks avaiable
    input.open(filename + ".mtx");

    file_opened_ok = input.is_open();

    if (file_opened_ok)
    {
        while(getline(input, line)) //File we have not finished,
        {
            if (line[0] != '%') //If first char is not a comment
            {
                istringstream iss(line); //Transform into stream
                if (read_header) //Have we read the header?
                {
                    //Then check if network is weighted
                    if (typeid(B) == typeid(bool))
                    {
                        iss >> from >> to; //Fill vars with data from stream
                        add_link(from, to);
                    }
                    else
                    {
                        iss >> from >> to >> w; //Fill vars with data from stream
                        add_link(from, to, w);
                    }

                }
                else
                {
                    read_header = true; //Then mark header as read
                    iss >> from >> from >> to; //Matrix NxN so first two numbers are the same. "From" is N and "To" is the number of links
                    add_nodes(from); //Add N nodes to the network
                }

            }

        }
    }
    else cout << "[CNetwork Error]: The requested MTX file " << filename << ", could not be opened." << endl;
    input.close();

    return file_opened_ok;
}


// ========================================================================================================
// ========================================================================================================
// ========================================================================================================


/** \brief Gets the degree of the target node
*  \param node_index: target node
*  \return degree of target node
*
* Returns the degree of the target node
*/
template <class T, typename B>
int CNetwork<T,B>::degree(const int node_index) const
{
    return this->neighs[node_index].size();
}

/** \brief Gets the in-degree of the target node
*  \param node_index: target node
*  \return in-degree of target node
*
* Returns the in-degree of the target node
*/
template <class T, typename B>
int CNetwork<T,B>::in_degree(const int node_index) const
{
    return this->neighs[node_index].size();
}

/** \brief Gets the out-degree of the target node
*  \param node_index: target node
*  \return out-degree of target node
*
* Returns the out-degree of the target node
*/
template <class T, typename B>
int CNetwork<T,B>::out_degree(const int node_index) const
{
    return this->neighs[node_index].size();
}


/** \brief Get neighbours of given node
*  \param node_index: target node
*  \return vector containing the neighbours of node_index
*
* Returns the a vector with the indices of the neighbours of node_index
*/
template <class T, typename B>
vector<unsigned int> CNetwork<T,B>::get_neighs_out(const int node_index) const
{
    return this->neighs[node_index];
}


/** \brief Get neighbours of given node
*  \param node_index: target node
*  \return vector containing the neighbours of node_index
*
* Returns the a vector with the indices of the neighbours of node_index
*/
template <class T, typename B>
vector<unsigned int> CNetwork<T,B>::get_neighs_in(const int node_index) const
{
    return this->neighs[node_index];
}


/** \brief Selects the k-th neighbour of node_index
*  \param node_index: target node
*  \param k: node to be selected
*  \return index of the k-th neighbour of the list
*
* Returns the index of the k-th neighbour of the given node.
*/
template <class T, typename B>
int CNetwork<T,B>::get_out(const int node_index, const int k) const
{
    return this->neighs[node_index][k];
}


/** \brief Selects the k-th neighbour of node_index
*  \param node_index: target node
*  \param k: node to be selected
*  \return index of the k-th neighbour of the list
*
* Returns the index of the k-th neighbour of the given node.
*/
template <class T, typename B>
int CNetwork<T,B>::get_in(const int node_index, const int k) const
{
    return this->neighs[node_index][k];
}



/** \brief Gets a link between two specified nodes
*  \param from: origin node
*  \param to: destination node
*  \return index of the link between from and destination.
*
* Returns the index of the link that connects nodes from and to. If there is no link, then
* it returns -1.
*/
template <class T, typename B>
int CNetwork<T,B>::get_link_index(const int from, const int to) const
{
    int i,even,odd;
    bool found = false;
    i = 0;

    while (i < this->link_count and not found)
    {
        found = (this->adjm[i].x == from and this->adjm[i].y == to) or (this->adjm[i].x == to and this->adjm[i].y == from);
        i += 1;
    }

    return found ? i-1 : -1; //Remember we have just summed i
}


/** \brief Get the neighbours of a given node
*  \param node_index: target node
*  \return vector containing the indices of the neighbour nodes of the target node
*
* Returns the a vector with the indices of the neighbours of the specified node.
*/
template <class T, typename B>
vector<unsigned int> CNetwork<T,B>::get_neighs(const int node_index) const
{
    return this->neighs[node_index];
}


/** \brief Selects a neighbour of a given node
*  \param node_index: target node
*  \param k: neighbour to be selected
*  \return index of the k-th neighbour of the target node
*
* Returns the index of the k-th neighbour of the target node. Neighbours are unsorted
*/
template <class T, typename B>
int CNetwork<T,B>::get_neigh_at(const int node_index, const int k) const
{
    return this->neighs[node_index][k];
}

