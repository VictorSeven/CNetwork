# CNetwork


**CNetwork is a simple C++ class used to do network analysis.** The aim of this class is to have a single .cpp file which can be imported in any project, **without having to compile any source code.** The code runs moderately fast and allows the user to do analysis with networks of decent sizes (~10<sup>6</sup>). 

Although is not as complete as other libraries, such as NetworkX or igraph, it is faster than NetworkX (see benchmark below) and it is easier to code than igraph for C++. 
This makes CNetwork a good choice for making simple, fast code for network analysis.

Main features:

- Functions to create well-known network models, such as Erdos-Renyi, Watts-Strogatz, Barabási-Albert or any scale-free distribution using the configuration model.
- Functions to manually add or remove nodes and links, as well as accessing them. 
- Get node properties, such as degree, clustering, or distance between nodes.
-  Import mtx format.
- Allows to have custom node properties to run dynamics over networks. 
- Export your network to Graphml / Gephi compatible format to visualize it. It is possible to define custom tags or properties such as color or position.

## An example of implementation

Here you have a very simple example of code that creates an Albert-Barabási network and computes the average degree. Simply download the source from the repositories, and do something like:

    #include<iostream>
    #include<cstdlib>
    
    #include"CNet.cpp"
    
    using namespace std;
    
    int main(void)
    {
        CNetwork<> net(20000); //Create a network of max size 20000 nodes
        net.create_albert_barabasi(2, 2, 21354647); //Fill this size with AB model
        cout << net.mean_degree() << endl; //Compute mean degree
    }


See the wiki for more details on the use of CNetwork.

## Simple benchmark

To measure the speed of CNetwork, I used the following example: create an Albert-Barabási network of one million nodes and get the largest eigenvalue.

    #include<iostream>
    #include<cstdlib>
    
    #include"CNet.cpp"
    
    using namespace std;
    
    int main(void)
    {
        CNetwork<> net(20000); //Create a network of max size 20000 nodes
        net.create_albert_barabasi(2, 2, 5464531); //Fill this size with AB model
        //compute_eigenv returns a (N+1) vector where the largest eigenvalue is the last element
        vector<double> eigenv = net.compute_eigenv(0.01);
        cout << eigenv[eigenv.size() -1] << endl; //Get the eigenvalue
    }


The code above only needs **1.953** seconds (less than 2s!) to run on a moderately old computer (i3-370M processor, 2.4 GHz, from 2010). For comparison, I executed a similar code on Python's NetworkX. Only the creation of the network via

    %timeit nx.barabasi_albert_graph(1000000, 2, 5464531) 

gives 16.1 seconds. Saving the network into a variable a computing its Laplacian, 

    %timeit nx.normalized_laplacian_matrix(barabasi)

gives 15.4 seconds, making a total of ~30 seconds.  

Probably other libraries, as igraph, are able to run faster than CNetwork. However, **the real advantage of CNetwork is being as simple to code as NetworkX but several times faster.** 
