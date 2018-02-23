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

        bool remove_node(int index);
        void add_link(int from, int to);
        void add_link(int from, int to, B w);
        bool remove_link(int from, int to);


        double mean_degree() const;


        double clustering_coef(int node_index) const;
        double mean_clustering_coef() const;


        void degree_distribution(vector<int> &distribution, bool normalized = false) const;
        void degree_correlation(vector<int> &distribution, vector<double> &correlation, bool normalized = false) const;


        int in_degree(int node_index) const;
        int out_degree(int node_index) const;
        int degree(int node_index) const;


        vector<unsigned int> get_neighs_out(int node_index) const;
        vector<unsigned int> get_neighs_in(int node_index) const;
        int get_out(int node_index, int k) const;
        int get_in(int node_index, int k) const;

        int get_link_index(int from, int to) const;


        vector<unsigned int> get_neighs(int node_index) const;
        int get_neigh_at(int node_index, int k) const;

        CNetwork(int max_size);

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
    return DirectedCNetwork<T,B>::mean_degree(this->TOTAL_DEGREE);
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
    DirectedCNetwork<T,B>::degree_distribution(distribution, this->TOTAL_DEGREE, normalized);
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
int CNetwork<T,B>::degree(int node_index) const
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
int CNetwork<T,B>::in_degree(int node_index) const
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
int CNetwork<T,B>::out_degree(int node_index) const
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
vector<unsigned int> CNetwork<T,B>::get_neighs_out(int node_index) const
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
vector<unsigned int> CNetwork<T,B>::get_neighs_in(int node_index) const
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
int CNetwork<T,B>::get_out(int node_index, int k) const
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
int CNetwork<T,B>::get_in(int node_index, int k) const
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
int CNetwork<T,B>::get_link_index(int from, int to) const
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
vector<unsigned int> CNetwork<T,B>::get_neighs(int node_index) const
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
int CNetwork<T,B>::get_neigh_at(int node_index, int k) const
{
    return this->neighs[node_index][k];
}

