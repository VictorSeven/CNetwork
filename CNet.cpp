#include<iostream>
#include<cstdlib>
#include<vector>
#include<random>
#include<algorithm>
#include<fstream>
#include<sstream>

#define MAX_SIZE 1000000
#define IN 0
#define OUT 1

using namespace std;


/// ==================================== Class definition ==================================
class CNetwork
{
    public:

        void add_nodes(int n);
        bool remove_node(int index);
        void add_link(int from, int to);
        void add_link(int from, int to, double w);
        bool remove_link(int from, int to);
        void create_adjacency_matrix();

        double mean_degree();
        double clustering_coef(int node_index);
        double mean_clustering_coef();
        void bread_first_search(int node, vector<int> &node_indices, vector<int> &dist);
        int component_size(vector<int> &node_in_this_component, vector<int> &size_of_components);
        vector<int> degree_distribution();
        vector<double> degree_correlation(vector<int> &distribution);
        double average_pathlenght();

        void create_albert_barabasi(int m0, int m, unsigned int random_seed);
        void create_configurational(int nodes, int kmin, double gamma, unsigned int random_seed);
        void create_wats_strogatz(int nodes, int regular_connections, double p, unsigned int random_seed);
        void create_erdos_renyi(int nodes, double mean_k, unsigned int random_seed);


        int degree(int node_index);
        vector<int> get_link(int link_index);
        double get_weight(int link_index);
        int get_link_index(int from, int to);
        int get_node_count();
        int get_link_count();
        vector<unsigned int> get_neighs(int node_index);
        int get_num_neighs(int node_index);
        int get_neigh_at(int node_index, int k);
        int get_a(int i, int j);


        void write_graphml(string filename, vector<string> labels = vector<string>(),  vector<string> info = vector<string>());
        void write_mtx(string filename);
        void read_mtx(string filename);

        CNetwork(int max_size, bool weight);

    private:
        void clear_network(int max_size, bool weight);

        bool weighted_net;

        int max_net_size;
        int current_size;
        int link_count;
        vector<unsigned int> links;
        vector<double> weight;
        vector<string> labels;
        vector< vector<unsigned int> > neighs;
        vector< vector<bool> > a;
        vector< vector<double> > a_w;
        vector<int> s;

        void matrixDotVector(vector< vector<double> >a, vector<double> v, vector<double>& r,  int n);
        double contract_vector(vector< vector<double> > a, vector<double> v);
        double calculateLambda(vector< vector<double> > a, vector<double> v, int n);
        double vectorNorm(vector<double>& v, int n);
        double largest_eigenvalue(vector< vector<double> > matrix, vector<double>& v, int n, double approx_error, int max_it);



};
/// ========================================================================================

/// ==================================== Constructor =======================================

///Create a newtork with a maximum size max_size
CNetwork::CNetwork(int max_size, bool weight)
{
    clear_network(max_size, weight);
    return;
}


void CNetwork::clear_network(int max_size, bool wght)
{
    weighted_net = wght; //Assign if we want a weighted network or not

    current_size = 0; //Init size to 0
    max_net_size = max_size; //Set the max size of the network
    link_count = 0; //Init link count
    //Create vectors in order to add things
    links = vector<unsigned int>();
    weight = vector<double>();
    neighs = vector< vector<unsigned int> >(0, vector<unsigned int>(0));

    a = vector< vector<bool> >();
    a_w = vector< vector<double> >();
    s = vector<int>();
    labels = vector<string>();

    return;
}



/// ========================================================================================



/// ================================= Add info functions ===================================

///Add n nodes to the network
void CNetwork::add_nodes(int n)
{
    int next_size = current_size + n; //New size of the network

    //Set the size to the maximum if it exceeds it
    current_size = next_size > max_net_size ? max_net_size : next_size;
    for (int i = 0; i < n; i++)
    {
        neighs.push_back(vector<unsigned int>()); //Add a new container for neighbours
    }
    return;
}

bool CNetwork::remove_node(int index)
{
    int i,j,k;

    if (index >= 0 and index < current_size)
    {
        neighs.erase(neighs.begin() + index); //Delete its entry in the neighbours

        //Re-scale all the neighbours index to fit the net, and also erase
        //the entry of the removed node as neighbour of others
        for (i=0; i < neighs.size(); i++)
        {
            for (j=0; j < get_num_neighs(i); j++)
            {
                k = get_neigh_at(i,j);
                if (k == index) //If it is the one, I want to erase, remove
                {
                    neighs[i].erase(neighs[i].begin() + j);
                }
                else if (k > index) //If it is higher, move 1 to left
                {
                    neighs[i][j] -= 1;
                }
            }
        }

        vector<int> who_to_erase = vector<int>(); //Track who you want to delete

        for (i=0; i < link_count; i++)
        {
            //Separe cases 2*i and 2*i+1 in order to correctly reduce degree
            if (links[2*i] == index)
            {
                //Our node is detected, erase the entry to that node
                who_to_erase.push_back(i);
            }
            if (links[2*i+1] == index)
            {
                who_to_erase.push_back(i);
            }

            //Also reduce index of nodes higher than our
            if (links[2*i] > index)
            {
                links[2*i] -= 1;
            }
            if (links[2*i+1] > index)
            {
                links[2*i+1] -= 1;
            }
        }

        //Perform the erase process
        for (i=0; i < who_to_erase.size(); i++)
        {
            //Index was recorded before, we do 2*who in order to get as even-odd again
            links.erase(links.begin() + 2*who_to_erase[i] - 2*i); //-2*i to take in account how list moves
            links.erase(links.begin() + 2*who_to_erase[i] - 2*i);  //This will eliminate the 2*i+1 since list moved to left
            //This uses directly the index
            weight.erase(weight.begin() + who_to_erase[i] - i);
        }

        link_count -= who_to_erase.size() / 2; //Reduce the number of links
        current_size -= 1; //Reduce the current size of the network in one unit

        return true;

    }
    else return false;
}

void CNetwork::add_link(int from, int to)
{
    links.push_back(from); //Even to origin,
    links.push_back(to); //Odd to destiny
    weight.push_back(1.0); //Unit weight
    neighs[from].push_back(to); //Add the node to the neighbours
    neighs[to].push_back(from); //And do it in the other sense also
    link_count += 1; //Create one link more


    return;
}

///Create a link between nodes from and to and weight w
void CNetwork::add_link(int from, int to, double w)
{

    int i = 0;
    bool existe = false;
    int rep_index;

    while (i < link_count and !existe)
    {
        if ((links[2*i] == from and links[2*i+1] == to) or (links[2*i] == to and links[2*i+1] == from))
        {
            rep_index = i;

            existe = true;
        }
        i+=1;
    }

    if (!existe)
    {
        links.push_back(from); //Even to origin,
        links.push_back(to); //Odd to destiny

        ///TODO maybe change
        weight.push_back(w); //Unit weight

        neighs[from].push_back(to); //Add the node to the neighbours
        neighs[to].push_back(from); //And do it in the other sense also

        link_count += 1; //Create one link more
    }
    else
    {
        weight[rep_index] += w;
    }
    return;
}

///Remove a link between from and to
bool CNetwork::remove_link(int from, int to)
{
    //Reduce the degree of the nodes
    auto index_it = find(neighs[from].begin(), neighs[from].end(), to); //Relative index of TO in terms of FROM
    int index_neigh = distance(neighs[from].begin(), index_it); //Get the relative index as an int

    if (index_neigh >= 0 and index_neigh < neighs[from].size())
    {
        int index_link = get_link_index(from, to);

        weight.erase(weight.begin() + index_link); //Erase this link weight
        link_count -= 1;

        links.erase(links.begin()+2*index_link);
        links.erase(links.begin()+2*index_link);//Use this last index obtained to erase from link array

        neighs[from].erase(neighs[from].begin()+index_neigh); //Erase the node it pointed in the neighbours list

        //Do the same process, but now in the other node, the to one.
        //Since this node was in the neigh list of FROM, we know that FROM has to be in the list of TO
        //That's why we don't check again if index_it < neighs end
        index_it = find(neighs[to].begin(), neighs[to].end(), from);
        index_neigh = distance(neighs[to].begin(), index_it);

        neighs[to].erase(neighs[to].begin()+index_neigh);

        return true;
    }
    else return false;

}


///Writes the adjacency matrix
void CNetwork::create_adjacency_matrix()
{
    int i,j;
    int aux;

    if (!weighted_net)
    {
        a = vector< vector<bool> >(current_size, vector<bool>(current_size, false) );
        //cout << "Adj. Matrix Byte Size = " << sizeof(a) + sizeof(vector<bool>) * current_size * current_size;

        for (i=0; i < current_size; i++)
        {
            for (j=0; j < get_num_neighs(i); j++)
            {
                a[i][get_neigh_at(i,j)] = true;
            }
        }
    }
    else
    {
        a_w = vector< vector<double> >(current_size, vector<double>(current_size, 0.0) );
        //cout << "Adj. Matrix Byte Size = " << sizeof(a) + sizeof(vector<bool>) * current_size * current_size;

        for (i=0; i < current_size; i++)
        {
            for (j=0; j < get_num_neighs(i); j++)
            {
                a_w[i][get_neigh_at(i,j)] = weight[get_link_index(i, get_neigh_at(i,j))];
            }
        }
    }

    return;
}

/// ========================================================================================


/// ================================= Topology functions ===================================

///Compute the mean degree of the network and returns it
double CNetwork::mean_degree()
{
    int i;
    double sum = 0.0; //Get the sum,

    //Sum over the network
    for (i=0; i < current_size; i++)
    {
        sum += degree(i);
    }
    //Divide by current size
    return sum / (current_size * 1.0);
}

///Computes the clustering coefficient of a particular node
double CNetwork::clustering_coef(int node_index)
{
    int i,j;
    int counter; //Count of pairs


    vector<unsigned int> nodes_neigh = get_neighs(node_index);
    vector<unsigned int> nodes_check;

    counter = 0;
    for (i=0; i < nodes_neigh.size(); i++) //Get neighbours of our node
    {
        nodes_check = get_neighs(nodes_neigh[i]); //Get neighbours of node i

        //For the next nodes, (start in j=i+1 to avoid double-count a pair)
        for (j=i+1; j < nodes_neigh.size(); j++)
        {
            //Use the fast find function from algorithm lib to see if this node is connected to any other neighbour of our node.
            //In that case, increase the counter
            if (find(nodes_check.begin(), nodes_check.end(), nodes_neigh[j]) != nodes_check.end()) counter += 1;
        }
    }

    if (degree(node_index) > 1)
    {
        return 2.0 * counter / (degree(node_index) * (degree(node_index) - 1)); //Finish computation and return clustering coefficient
    }
    else
    {
        return 0.0;
    }

}

///Computes the average clustering coefficient of the network
double CNetwork::mean_clustering_coef()
{
    int i;
    double sum = 0.0; //Get the sum,

    //Sum over the network
    for (i=0; i < current_size; i++)
    {
        sum += clustering_coef(i);
    }
    //Divide by current size
    return sum / (current_size * 1.0);
}


///Compute all the pathlenghts from node using the optimized version from Newman's book
void CNetwork::bread_first_search(int node, vector<int> &node_indices, vector<int> &dist)
{
    int w, r; //Write and read pointers;
    int d; //Current distance
    int i,j; //Counter
    int cur_node; //Current node we are evaluating
    int neigh; //Current neighbour we are seeing

    //Init everything
    for (i=0; i < current_size; i++)
    {
        node_indices[i] = -1;
        dist[i] = -1;
    }

    //Get first node
    node_indices[0] = node;
    dist[node] = 0;
    //Quantities,
    r = 0; //Reading first element,
    w = 1; //Writing second

    while(w != r)
    {
        cur_node = node_indices[r]; //The node we want to evaluate now
        r += 1; //We will want to take the next
        d = dist[cur_node]; //Get the distance this node is at,
        for (j=0; j < get_num_neighs(cur_node); j++)
        {
            neigh = get_neigh_at(cur_node, j); //Get the neighbour
            if (dist[neigh] == -1) //If distance is unknown,
            {
                dist[neigh] = d+1; //Put its distance,
                node_indices[w] = neigh; //Add it into the list of known nodes
                w += 1; //and increase write pointer
            }
        }
    }

    return;
}

///Use the breadth first search to get size of bigger component of the network
int CNetwork::component_size(vector<int> &node_in_this_component, vector<int> &size_of_components)
{
    int i,j,k;

    int remaining = current_size;

    node_in_this_component = vector<int>();
    size_of_components = vector<int>();

    vector<bool> visited(current_size, false); //List of visited nodes
    vector<int> node_list(current_size); //List of nodes
    vector<int> dist(current_size); //Distance to nodes
    int result;

    i = 0; //Init counters
    k = 0;
    while (i < current_size and remaining > 0) //While we have work to do,
    {
        if (!visited[i]) //If we have not visited this node before,
        {
            //cout << i << endl;
            remaining -= 1; //We have one less to see
            visited[i] = true; //Mark it
            bread_first_search(i, node_list, dist); //Compute paths from i to all other nodes and store them

            size_of_components.push_back(0); //Size is equal to zero when we start
            node_in_this_component.push_back(i); //The node we are now visiting is in this component - store it

            for (j=0; j < current_size; j++)
            {
                if (dist[j] >= 0) //If distance is positive, this is inside my cluster
                {
                    size_of_components[k] += 1;
                    visited[j] = true;
                    remaining -= 1;
                }
            }

            k += 1; //Increase the counter for size of components

        }

        i += 1;
    }

    result = 0;
    for (i=0; i < size_of_components.size(); i++)
    {
        if (size_of_components[i] > result)
        {
            result = size_of_components[i];
        }
    }


    return result;
}

///Computes the average path lenght of the network
double CNetwork::average_pathlenght()
{
    int i,j,k; //Counters
    //double sum;
    //int connected;
    int remaining = current_size; //How many nodes we have yet to evaluate
    int node; //auxiliar variable

    vector<bool> visited(current_size, false); //List of visited nodes
    vector<int> cluster_index = vector<int>(); //What nodes are reachable from me
    double pathlenght, maxpathlenght; //To account for pathlenght
    vector<int> node_list(current_size); //List of nodes
    vector<int> dist(current_size); //Distance to nodes

    maxpathlenght = 0.9; //Start with value > 0 so not every value will override this
    i = 0; //Init counters
    while (i < current_size and remaining > 0) //While we have work to do,
    {
        pathlenght = 0.0; //Start sum for average
        if (!visited[i]) //If we have not visited this node before,
        {
            //cout << i << endl;
            remaining -= 1; //We have one less to see
            visited[i] = true; //Mark it
            bread_first_search(i, node_list, dist); //Compute paths from i to all other nodes and store them

            //See which nodes are reachable from i
            for (j=0; j < current_size; j++)
            {
                if (dist[j] >= 0) //The condition is that distance is positive
                {
                    //cout << "-- " << j << endl;
                    cluster_index.push_back(j); //Store reachable nodes
                    pathlenght += dist[j]; //Add the average distance to them
                }
            }

            //Add the contribution of this marked nodes to the average pathlenght of this cluster
            for (j=0; j < cluster_index.size(); j++)
            {
                node = cluster_index[j]; //Get the selected node
                //If we have not visited it yet (to avoid counting first node found twice)
                if (!visited[node])
                {
                    bread_first_search(node, node_list, dist); //Compute paths from i to all other nodes and store them
                    visited[node] = true; //Mark as visited now
                    remaining -= 1; //Eliminate from remaining list
                    //We know beforehand which nodes we have to sum, because this are the reachable ones.
                    for (k=0; k < cluster_index.size(); k++)
                    {
                        pathlenght += dist[cluster_index[k]];
                    }
                }
            }

            //Finish pathlenght. Avoid single node divergence
            if (cluster_index.size() > 1)
            {
                pathlenght /= cluster_index.size() * (cluster_index.size() - 1);
            }
            else
            {
                pathlenght = 0.0;
            }

            //Check if this is max
            if (pathlenght > maxpathlenght)
            {
                maxpathlenght = pathlenght;
            }

            cluster_index.clear(); //Clear this var for another run
        }

        i += 1;
    }
    //cout << "end" << endl;

    //int node_list[current_size], dist[current_size]; //List of nodes and distances to them

    /*double sum = 0.0;
    for (i=0; i < current_size; i++)
    {
        bread_first_search(i, node_list, dist); //Compute paths from i to all other nodes and store them
        for (j=i+1; j < current_size; j++)
        {
            if (dist[j] > 0) //In other case, there is no path from i to j
            {
                sum += dist[j];
            }
        }
    }

    return 2.0*sum/((current_size-1)*current_size); //Finsh the mean */

    return maxpathlenght; //Return path lenght of largest component.
}

vector<int> CNetwork::degree_distribution()
{
    int i;
    vector<int> result(current_size, 0);
    //Compute the degree distribution
    for (i=0; i < current_size; i++)
    {
        result[degree(i)] += 1;
    }
    //Erase the 0 at the end of the array
    i = current_size - 1; //Start counter
    while (result[i] == 0)
    {
        result.erase(result.begin() + i);
        i -= 1;
    }
    return result;
}

///Compute the degree correlation (also return degree distribution)
vector<double> CNetwork::degree_correlation(vector<int> &distribution)
{
    int i,j,k;
    int index, numneighs;
    double mean_neigh_degree;
    vector<double> result;

    int maxdegree = 0;
    for (i=0; i < current_size; i++)
    {
        if (degree(i) > maxdegree)
        {
            maxdegree = degree(i);
        }
    }

    result = vector<double>(maxdegree, 0.0);
    distribution = vector<int>(maxdegree, 0);

    for (i=0; i < current_size; i++)
    {
        index = degree(i);
        distribution[index] += 1;
        numneighs = get_num_neighs(i);
        mean_neigh_degree = 0.0;
        for (j=0; j < numneighs; j++)
        {
            k = get_neigh_at(i,j);
            mean_neigh_degree += degree(k);
        }
        if (numneighs != 0) mean_neigh_degree /= 1.0 * numneighs;
        result[index] += mean_neigh_degree;
    }

    //To finish average over nodes, divide by the number of nodes with degree k
    //This are in distribution. Check that we are not dividing by 0
    for (i=0; i < result.size(); i++)
    {
        if (distribution[i] != 0) result[i] /= 1.0*distribution[i];
    }


    return result;
}

/// ========================================================================================


/// =================================== Network creation ===================================


void CNetwork::create_erdos_renyi(int nodes, double mean_k, unsigned int random_seed)
{
    int i,j,k;
    double p = mean_k / (nodes - 1.0);

    mt19937 gen(random_seed);; //Create the generator
    uniform_real_distribution<double> ran_u(0.0,1.0); //Uniform number distribution
    uniform_int_distribution<int>  index(0,nodes-1); //-1 because closed interval for ints

    add_nodes(nodes); //Create the nodes

    //For the n (n-1) / 2 pairs, link them with probability p
    for (i=0; i < current_size; i++)
    {
        for (j=i+1; j < current_size; j++)
        {
            if (ran_u(gen) <= p)
            {
                add_link(i, j);
            }
        }
    }

    return;

}


void CNetwork::create_configurational(int nodes, int mink, double gamma, unsigned int random_seed)
{
    int i,j;
    int n_links;
    int max_size; //Maximum size if we want an uncorrelated network

    vector<int> node_degree;
    vector<int> link_vector;

    mt19937 gen(random_seed); //Create the generator
    uniform_real_distribution<double> ran_u(0.0,1.0); //Uniform number distribution

    add_nodes(nodes);
    node_degree = vector<int>(current_size);
    max_size = sqrt(current_size);

    n_links = 0;

    double kmax = pow(max_size, 1.0-gamma);
    double kmin = pow(mink, 1.0-gamma);
    double invgamma = 1.0 / (1.0 - gamma);


    for (i=0; i < current_size; i++)
    {
        node_degree[i] = floor( pow( ran_u(gen)*(kmax - kmin) + kmin, invgamma ) );
        n_links += node_degree[i];
    }

    if (n_links % 2 == 1)
    {
        node_degree[0] += 1;
        n_links += 1;
    }

    link_vector = vector<int>(n_links);

    int k=0;
    for (i=0; i < current_size; i++)
    {
        for (j=0; j < node_degree[i]; j++)
        {
            link_vector[k] = i;
            k += 1;
        }
    }



    random_shuffle(link_vector.begin(), link_vector.end());


    for (i=0; i < link_vector.size()/2; i++)
    {
        if (link_vector[2*i] != link_vector[2*i+1])
        {
            add_link(link_vector[2*i], link_vector[2*i+1]);
        }
    }

    return;

}

void CNetwork::create_wats_strogatz(int nodes, int num_forward_edges, double p, unsigned int random_seed)
{
    int i,j;
    int to;
    mt19937 gen(random_seed);; //Create the generator
    uniform_real_distribution<double> ran_u(0.0,1.0); //Uniform number distribution
    uniform_int_distribution<int>  index(0,nodes-1); //-1 because closed interval for ints
    vector<unsigned int> aux;
    bool eliminated;

    //Add the nodes
    add_nodes(nodes);

    for (i=0; i < current_size; i++)
    {
        for (j=1; j <= num_forward_edges; j++) //At least one node.
        {
            add_link(i, (i+j)%current_size); //Add the link
        }
    }


    vector<int> who_to_erase = vector<int>(0);
    vector<int> who_to_erase_index = vector<int>(0);
    for (i=0; i < current_size; i++) //To -1 because the last is already connected to the first
    {
        //for (j=0; j < get_neighs(i).size(); j++)
        for (j = num_forward_edges; j < 2 * num_forward_edges; j++)
        {
            if (ran_u(gen) <= p)
            {
                who_to_erase.push_back(i);
                who_to_erase_index.push_back(get_neigh_at(i,j));
            }
        }
    }


    for (i=0; i < who_to_erase.size(); i++)
    {
        eliminated = remove_link(who_to_erase[i], who_to_erase_index[i]);
        do
        {
            to = index(gen);
            aux = get_neighs(to);

        }
        while (to == who_to_erase[i] or to == who_to_erase_index[i] or (find(aux.begin(), aux.end(), who_to_erase[i]) != aux.end()) );
        if (eliminated) add_link(who_to_erase[i], to);
    }

    return;
}


///Creates an Albert-Barabasi free scale network
void CNetwork::create_albert_barabasi(int m0, int m, unsigned int random_seed)
{
    int i,j;
    double r;

    mt19937 gen(random_seed);; //Create the generator
    uniform_real_distribution<double> random(0.0,1.0); //Uniform number distribution
    uniform_int_distribution<int>  index(0,1); //Useful to select indices

    //Create fully connected network with m0 nodes
    add_nodes(m0);
    for (i=0; i < m0; i++)
    {
        for (j=i+1; j < m0; j++)
        {
            add_link(i,j);
        }
    }

    //Add then more nodes...
    for (i = m0; i < max_net_size; i++)
    {
        add_nodes(1);
        //For every link we want to do,
        for (j=0; j < m; j++)
        {
            //Generate a random number
            if (random(gen) <= 0.5)
            {
                //With half probability, add it to highly connected node,
                //selecting randomly an edge.
                index = uniform_int_distribution<int>(0, link_count-1);//-1 because interval is closed
                add_link(current_size-1, links[2*index(gen)+1]); //Current_size-1 because we've added a node we don't have to account for; BEFORE WAS [index][OUT]
            }
            else
            {
                //If not, connect to a random node
                index = uniform_int_distribution<int>(0, current_size-2); //We don't want current_size-1 which is the actual node, so up to -2
                add_link(current_size-1, index(gen)); //Add the link
            }
        }
    }

    return;
}

/// ========================================================================================



/// ===================================== Getting info =====================================


///Get the degree of provided index
int CNetwork::degree(int node_index)
{
    //return degree[node_index];
    return neighs[node_index].size();
}


///Get the degree of provided index
int CNetwork::get_link_index(int from, int to)
{
    int i,even,odd;
    bool found = false;
    i = 0;

    while (i < link_count and not found)
    {
        found = (links[2*i] == from and links[2*i+1] == to) or (links[2*i] == to and links[2*i+1] == from);
        i += 1;
    }

    return found ? i-1 : -1; //Remember we have just summed i
}

///Returns the link in format (from, to) as a two-component vector
vector<int> CNetwork::get_link(int link_index)
{
    return {links[2*link_index], links[2*link_index+1]};
}

///Get the weight
double CNetwork::get_weight(int link_index)
{
    return weight[link_index];
}

///Get how many nodes we have
int CNetwork::get_node_count()
{
    return current_size;
}

///Get how many links we have
int CNetwork::get_link_count()
{
    return link_count;
}

///Returns a vector of neighbors of node_index
vector<unsigned int> CNetwork::get_neighs(int node_index)
{
    return neighs[node_index];
}

///Returns how meany neighbours a node has
int CNetwork::get_num_neighs(int node_index)
{
    return neighs[node_index].size();
}

///Get the k neighbour of the node node_index
int CNetwork::get_neigh_at(int node_index, int k)
{
    return neighs[node_index][k];
}

///Returns element of the adjacency matrix, depending if
///we have selected or not a weighted network
int CNetwork::get_a(int i, int j)
{
    return weighted_net ? a_w[i][j] : int(a[i][j]);
}


///This function creates an standard graphml file format, which is able to store all the data of CNetwork:
///nodes, node labels, links and weights.
void CNetwork::write_graphml(string filename, vector<string> labels, vector<string> info)
{
    int i;
    ofstream output;

    output.open(filename + ".graphml"); //Open the file


    //Write the XML header
    output << "<?xml version='1.0' encoding='utf-8'?>" << endl; //Encoding
    //All the URL info for graphml format
    output << "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" ";
    output << "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns ";
    output << "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\"> " << endl;
    //Define the weight as an attribute
    output << "<key attr.name=\"weight\" attr.type=\"double\" for=\"edge\" id=\"w\" />" << endl;
    output << "<key attr.name=\"label\" attr.type=\"string\" for=\"node\" id=\"s\" />" << endl;
    //Create graph
    output << "<graph edgedefault=\"undirected\">" << endl;

    //If we have labels
    if (labels.size() > 0)
    {
        if (info.size() > 0)
        {
            for (i=0; i < current_size; i++)
            {
                output << "<node id=\"" << labels[i] << "\">" << endl; //Then define nodes as labels
                output << "<data key=\"s\">" << info[i] << "</data>" << endl;
                output << "</node>" << endl;
            }
        }
        else
        {
            for (i=0; i < current_size; i++)
            {
                output << "<node id=\"" << labels[i] << "\"/>" << endl; //Then define nodes as labels
            }
        }

        for (i=0; i < link_count; i++)
        {
            //And link using the labels
            output << "<edge source=\"" << labels[links[2*i]] << "\" target=\"" << labels[links[2*i+1]] << "\">" << endl;
            output << "<data key=\"w\">" << weight[i] << "</data>" << endl;
            output << "</edge>" << endl;
        }
    }
    else
    {
        //Use node indices as labels
        if (info.size() > 0)
        {
            for (i=0; i < current_size; i++)
            {
                output << "<node id=\"" << i << "\">" << endl; //Then define nodes as labels
                output << "<data key=\"s\">" << info[i] << "</data>" << endl;
                output << "</node>" << endl;
            }
        }
        else
        {
            for (i=0; i < current_size; i++)
            {
                output << "<node id=\"" << i << "\"/>" << endl; //Then define nodes as labels
            }
        }
        for (i=0; i < link_count; i++)
        {
            output << "<edge source=\"" << links[2*i] << "\" target=\"" << links[2*i+1] << "\">" << endl;
            output << "<data key=\"w\">" << weight[i] << "</data>" << endl;
            output << "</edge>" << endl;
        }
    }


    //Close the open tags and the file
    output << "</graph>" << endl << "</graphml>" << endl;
    output.close();
}

///Write the mtx format
void CNetwork::write_mtx(string filename)
{
    ofstream output;
    int i;

    output.open(filename + ".mtx");
    //Write that this is a NxN matrix with link_count links
    output << current_size << " " << current_size << " " << link_count << endl;
    //Check if the used decided to use a weighted net
    if (!weighted_net)
    {
        for (i=0; i < link_count; i++)
        {
            output << links[2*i] << " " << links[2*i+1] << endl; //Using a weight 1
        }
    }
    else
    {
        for (i=0; i < link_count; i++)
        {
            output << links[2*i] << " " << links[2*i+1] << " " << weight[i] << endl; //Arbitrary weight
        }
    }
    output.close();
}


void CNetwork::read_mtx(string filename)
{
    //Destroy this object and create new network
    clear_network(max_net_size, weighted_net);

    bool read_header = false; //To see if we have read the dim1xdim2 links line
    string line; //Store a line
    ifstream input; //File
    int from, to, w; //Auxiliary

    //Open the file and checks avaiable
    input.open(filename);
    if (input.is_open())
    {
        while(getline(input, line)) //File we have not finished,
        {
            if (line[0] != '%') //If first char is not a comment
            {
                istringstream iss(line); //Transform into stream
                if (read_header) //Have we read the header?
                {
                    //Then check if network is weighted
                    if (!weighted_net)
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
            else
            {
                 cout << line[0] << endl; //Ignore line
            }

        }
    }
    input.close();
}

/// ========================================================================================


/// ============================= AUXILIARY FUNCTIONS ======================================


//Producto de una matriz por un vector. Almacena el resultado de A*v en el vector r
void CNetwork::matrixDotVector(vector< vector<double> > a, vector<double> v, vector<double>& r,  int n)
{
	int i,j;
	double sum;

	for (i=0; i < n; i++)
	{
		r[i] = v[i];
	}

	sum = 0.0;
	for (i=0; i < n; i++)
	{
		sum = 0.0;
		for (j=0; j < n; j++)
		{
			sum += a[i][j]*v[j];
		}
		r[i] = sum;
	}


	return;
}

//Calcula el autovalor de un vector según la fórmula de Rayleigh y la devuelve
double CNetwork::calculateLambda(vector< vector<double> > a, vector<double> v, int n)
{
    int i;
    double lambda;
    double sum1,sum2;
    vector<double> aux(v.size());

    matrixDotVector(a,v,aux,n);

    sum1=sum2=0.0;
    for (i=0; i < n;i++)
    {
        sum1 += aux[i]*v[i];
        sum2 += v[i]*v[i];
    }

    lambda = sum1 / sum2;

    return lambda;
}

double CNetwork::vectorNorm(vector<double>& v, int n)
{
    int i;

    double sum =0.0;
    for (i=0; i < n; i++)
    {
        sum += v[i]*v[i];
    }

    return sqrt(sum);
}

///Get the largest eigen_value of a matrix
double CNetwork::largest_eigenvalue(vector< vector<double> > matrix, vector<double>& v, int n, double approx_error, int max_it)
{
    double calc_error = 1e10; //Start with high value so loop starts
    double lamb1, lamb2;
    double norma;
    vector<double> aux(n);
    int i,j;

    lamb1 = calculateLambda(matrix,v,n);
    while(abs(calc_error) > approx_error && i < max_it)
    {
        matrixDotVector(matrix,v,aux,n); //Iteration A^k * v_k

        //Compute lambda_k y lambda_(k+1)
        //lamb1 = calculateLambda(matrix,v,n);
        lamb2 = calculateLambda(matrix,aux,n);

        //Calcular el error
        if (lamb2 != 0.0)
        {
           calc_error = abs(lamb1-lamb2)/lamb2 * 100.0; //Compute the error
        }

        //Normalize the vector
        norma = vectorNorm(aux,n);

        for (j=0; j < n; j++)
        {
            v[j] = aux[j]/norma; //Set to v_k the v_k+1 vector
        }
        i++;

        lamb1 = lamb2; //Store old value
    }

    return lamb2;
}

///Contracts vector with a matrix
double contract_vector(vector< vector<double> > a, vector<double> v)
{
    int i, j;
    int n = v.size();
    double sum = 0.0;

    for (i=0; i < n; i++)
    {
        sum += a[i][i] * v[i] * v[i];
        for (j=i+1; j < n; j++)
        {
            sum += 2.0 * a[i][j] * v[i] * v[j];
        }
    }

    return sum;
}

