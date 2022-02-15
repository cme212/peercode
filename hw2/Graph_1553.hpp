#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <sstream>
#include <string>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V=int,typename E=int>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
 
 /* Forward declaration of a node linked list class (used for incident iterators)*/
  struct node_list;
   
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  
 
  /** Predeclaration of Node type. */
  using node_value_type= V ;
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  /** Predeclaration of Edge type. */
  using edge_value_type=E;
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    size_=0;
    num_removed_node_=0;
    num_edge_=0;
    num_removed_edge_=0;
    max_size_=600;
    nodes_=std::vector<Node> ();
    map_small_big_nodes_=std::vector<size_type> ();
    map_big_small_nodes_=std::vector<size_type> ();
    valid_nodes_=std::vector<bool> ();
    points_=std::vector<Point> ();
    edges_=std::vector<Edge> ();
    valid_edges_=std::vector<bool> ();
    incidents_=std::vector<node_list> ();
    values_=std::vector <node_value_type>();
    crosses_=init_crosses(max_size_);
    edge_values_=std::vector<edge_value_type> ();
    degrees_=std::vector<size_type>();
    null_point=Point(NAN);
    //empty_list=node_list();
    // HW0: YOUR CODE HERE
  }

/* Above please note that incidents is a vector of node_lists storing for each node its incident list
 Crosses is an adjacent matrix where crosses[i,j]=True if and only if edge(i,j) in graph, it is initiated
 to be of size 600 and is updated if the size of the graph becomes too big*/

  /** Default destructor */
  ~Graph() = default;



  //
  // NODES
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node: private totally_ordered<Node> {
   public:
    /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::node_type x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     */

    /* A node has two attributes: a pointer to its graph and its index in the graph*/
    Node() :graph_(nullptr),index_(size_type(-1)) {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      /* The loop belows check if the graph wasn't cleared since the creation of the node, 
       in that case the node is invalid*/
      if ((*graph_).size_<=index_){
       return (*graph_).null_point;}
      return (*graph_).points_[index_];
    }


   Point& position() {
      if ((*graph_).size_<=index_){
       return (*graph_).null_point;}
      return (*graph_).points_[index_];
    }



    /** Return this node's index, a number in the range [0, graph_size). */
    /* The index_ attribute stores the index of the node at its creation (ie: it is 
       not modified if other nodes are deleted). If we want its current index we need to apply
       the mapping map_big_small_nodes*/
    size_type index() const {
      // HW0: YOUR CODE HERE
      if (((*graph_).size_<=index_)){return size_type(-1);}
      return (*graph_).map_big_small_nodes_[index_];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      bool test1=(index()==n.index());
      bool test2=(graph_==n.graph_);
      return test1&test2;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& n) const {
      /* If the 2 nodes are on the same graph we compare the index, if not
       we compare the memory adresses of their graph as strings*/
      bool test1=(graph_==n.graph_);
      if (test1) {return index()<n.index();}
      else { 
      std::stringstream ss1;
      std::stringstream ss2;
      ss1<<graph_;
      ss2<<n.graph_;
      return ss1.str()<ss2.str();}    
    }
   

   size_type degree () const {
   if (index_>(*graph_).size()-1){return 0;}
   return (*graph_).degrees_[index_];
    }


  incident_iterator edge_begin () const {
  return IncidentIterator(this,&(*graph_).incidents_[index_]);
}

  incident_iterator edge_end () const {
  return IncidentIterator(this,&(*graph_).empty_list);

}
  
const node_value_type& value() const {return (*graph_).values_[index_];}

node_value_type& value() {return (*graph_).values_[index_];}
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph * graph_;
    size_type index_;
    Node (Graph * graph, size_type index)
    :graph_(graph),index_(index){}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_-num_removed_node_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */


/* Please note that all the operations below are in O(1) except the reshaping of the crosses
 matrix if the graph becomes too big, since in that case the max_size is multiplied by 4 the 
 amortized complexity is in O(1)*/
  Node add_node(const Point& position, const node_value_type & val=node_value_type ()) {
    size_type index=size_;
    Graph *graph=this;
    Node n=Node(graph,index);
    nodes_.push_back(n);
    values_.push_back(val);
    points_.push_back(position);
    node_list incident_n=empty_list;
    incidents_.push_back(incident_n);
    degrees_.push_back(0);
    if (size_>=1){
    //We also need to update the mappings, og_index-> current_index
    // and current_index -> og_index
    map_small_big_nodes_.push_back(map_small_big_nodes_[size_-1]+1);
    size_type last=map_big_small_nodes_[size_-1];
    map_big_small_nodes_.push_back(last+1);}
    else {
      map_small_big_nodes_.push_back(0);
      map_big_small_nodes_.push_back(0);}
    valid_nodes_.push_back(true);
    size_++;
    if (size_>max_size_){
       crosses_=update_crosses();}
    return n;        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
 /* We just check the value of the pointer graph_ and we check that the graph wasn't cleared in between, ie:
   that the index is <graph.size*/
  bool has_node(const Node& n) const {
    Graph * graph=const_cast<Graph *>(this);
    return (graph==n.graph_)&(n.index_<size_)&(valid_nodes_[n.index_]);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return nodes_[map_small_big_nodes_[i]];        // Invalid node
  }

  /** Removes the node  from the graph
   * @return false e if the node e was already removed true otherwise
   * @post new size()= old size()-1 if the node was removed
   * @post new num_edge()= old num_edge()-@a n.degree()
   * @post new node(@a n.index())=old node(@a n.index()+1)
   * Complexity: O(n_nodes)
   */

size_type remove_node (const Node & n){
   if (valid_nodes_[n.index_]==false) {return 0;}
   num_removed_node_++;
   valid_nodes_[n.index_]=false;
   //We remove the incident edges
   for (IncidentIterator ei=n.edge_begin();ei!=n.edge_end();++ei){remove_edge(*ei);} 
   //We need to update the mapping
  for (size_type i=0;i<size_-1;++i){
        if (map_small_big_nodes_[i]>=n.index_){map_small_big_nodes_[i]=map_small_big_nodes_[i+1];}
        if (i>=n.index_){--map_big_small_nodes_[i];}}
   map_big_small_nodes_[size_-1]--;
   if (size_==1){ map_small_big_nodes_[size_-1]=0;}
   else {map_small_big_nodes_[size_-1]=map_small_big_nodes_[size_-2]+1;}
   return 1;}

  /** Removes the node  * @a n_it from the graph
   * @return @a n_it
   * @post new size()= old size()-1 if the node was removed
   * @post new num_edge()= old num_edge()-*@a n_it.degree() if the node was removed
   * @post new node(*@a n_it.index())=old node(*@a n_it.index()+1)
   * Complexity: O(n_nodes)
   */


node_iterator remove_node(node_iterator n_it){
       remove_node(*n_it);
       return n_it;}


  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */

 /* An Edge has 2 attributes */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() :index_1_(size_type(-1)),index_2_(size_type(-1)),index_(size_type(-1)),graph_(nullptr){}

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return (*graph_).nodes_[index_1_];      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return (*graph_).nodes_[index_2_];      // Invalid Node
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
                 // Quiet compiler warning
      //HW0: YOUR CODE HERE
      bool test1=(node1()==e.node1())&(node2()==e.node2());
      bool test2=(node1()==e.node2())&(node2()==e.node1());
      return test1|test2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

/* To ensure tricodomoy we first check that e1!=e2, then if e1.n1==e2.n1 we compare the second nodes
  else we compare the first nodes*/
    bool operator<(const Edge& e) const {
            // Quiet compiler warning
      //HW0: YOUR CODE HERE    
      bool test1=(node1()==e.node1())&(node2()==e.node2());
      bool test2=(node2()==e.node1())&(node1()==e.node2());
      if (test1|test2){return false;}
      if (node1()==e.node1()){return node2()<e.node2();}
      return node1()<e.node1();
    }


double length() const {
  Point p1=node1().position();
  Point p2=node2().position();
  return norm(p1-p2);}

edge_value_type& value() {
  return (*graph_).edge_values_[index_];}

const edge_value_type& value() const {
  return (*graph_).edge_values_[index_];}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type index_1_;
    size_type index_2_;
    size_type index_;
    Graph * graph_;
    Edge(size_type index_1, size_type index_2,size_type index,Graph * graph):index_1_(index_1), index_2_(index_2), index_(index), graph_(graph){}      
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HER
    return num_edge_-num_removed_edge_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  /* We just return the edge at the correct index*/
  Edge edge(size_type i) const {
    
    if (num_removed_edge_==0) {return edges_[i];} 
    size_type index=0;
    for (size_type k=0;k<num_edge_;++k){
         index+=valid_edges_[k];
         if (index==i+1){return edges_[k];}}
    return Edge();
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */


/* In that implementation we just need to check the crosses matrix*/

bool has_edge (const Node& a, const Node& b) const {
   return crosses_[a.index_][b.index_]; 
}

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

/* We just need to update the different data structures*/

Edge add_edge(const Node& a, const Node&b, const edge_value_type & val=edge_value_type ()){
    size_type index_a=a.index_;
    size_type index_b=b.index_;     
    const Edge test_e=Edge(index_a,index_b,num_edge_,this);
    if (crosses_[a.index()][b.index()]==true){return test_e;}
    edges_.push_back(test_e);
    incidents_[index_a].push_back(b,num_edge_);
    incidents_[index_b].push_back(a,num_edge_);
    crosses_[index_a][index_b]=true;
    crosses_[index_b][index_a]=true;
    degrees_[index_a]++;
    degrees_[index_b]++;
    edge_values_.push_back(val);
    valid_edges_.push_back(true);
    num_edge_++;
    return test_e;}


  /** Removes the edge to the graph
   * @return false if the edge was already removed true otherwise
   * @post has_edge(@a e.node1(), @a e.node2())==false
   * @post new num_edge()= old num_edge()-1
   * @post num @a e.node1().degree() = old @a e.node1().degree()-1 if edge was present
   * @post num @a e.node2().degree() = old @a e.node2().degree()-1 if edge was present
   * Complexity: O(1), constant operations are done
   */

size_type remove_edge(const Edge& e){
    size_type index_1=e.index_1_;
    size_type index_2=e.index_2_;
    //Checking if the edge was already removed
    if (crosses_[index_1][index_2]==false){return 0;}
    //Updating the different variables
    size_type index=e.index_;
    crosses_[index_1][index_2]=false;
    crosses_[index_2][index_1]=false;
    valid_edges_[index]=false;
    num_removed_edge_++;
    degrees_[index_1]--;
    degrees_[index_2]--;
    return 1;}

  /** Removes edge(@a, @b ) to the graph
   * @return True if the edge was already removed false otherwise
   * @post has_edge(@ a, @ b)==false
   * @post new num_edge()= old num_edge()-1 if edge was present
   * @post num @& e.node1().degree() = old @a e.node1().degree()-1 if edge was present
   * @post num @a e.node2().degree() = old @a e.node2().degree()-1 if edge was present
   * Complexity: O(num_edge), we need to retrieve the index
   */


size_type remove_edge(const Node & a, const Node & b){
  Edge e=Edge(a.index_,b.index_,0,this);
  //We need to retrieve the index by iterating through the edges
  //This is why this version is in O(n_edges)
  for (EdgeIterator ei=edge_begin();ei!=edge_end();++ei){
         if (e==(*ei)){e=(*ei);}}
  const Edge ce=e;
  return remove_edge(ce);}

  /** Removes the edge *@a e_it to the graph
   * @return @a e_it
   * @post has_edge( *@a e_it.node1(), *@a e_it.node2())==false
   * @post new num_edge()= old num_edge()-1
   * @post num *@a e_it.node1().degree() = old *@a e_it.node1().degree()-1 if edge was present
   * @post num *@a e_it.node2().degree() = old *@a e_it.node2().degree()-1 if edge was present
   * Complexity: O(1), constant operations are done
   */



edge_iterator remove_edge(edge_iterator e_it){
    const Edge e=(*e_it);
    size_type i=remove_edge(e);
    return e_it;}

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    /*for (size_type i=0;i<size_;i++){
    Node * n=new &nodes_[i];
    delete &nodes_[i]; delete &points_[i];}
    for (size_type i=0;i<num_edge_;i++){delete &edges_[i];}*/
    nodes_.clear();
    map_small_big_nodes_.clear();
    map_big_small_nodes_.clear();
    valid_nodes_.clear();
    points_.clear();
    edges_.clear();
    incidents_.clear();
    values_.clear();
    size_=0;
    num_removed_node_=0;
    max_size_=600;
    crosses_=init_crosses(max_size_);
    edge_values_.clear();
    valid_edges_.clear();
    num_edge_=0;
    num_removed_edge_=0;
    }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */



/* Our node iterator is simply conssited by a pointer to the graph and the index of the node
 in its nodes_ vector*/
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() :graph_(nullptr),index_(size_type(-1)){}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
    //We separate the case where the index is bigger then the size of the graph
    if (index_>=(*graph_).size()){return Node();}
    return (*graph_).node(index_);}
    NodeIterator& operator++(){
       size_type size=(*graph_).size();
       if (index_>=size){
       index_=size_type(-1);
       }
       else {
       index_=index_+1;}
       return *this;

     }
    bool operator==(const NodeIterator& ni) const{
         return (graph_==ni.graph_)&(index_==ni.index_);}

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph * graph_;
    size_type index_;
    NodeIterator (const Graph * graph,size_type index) : graph_(graph),index_(index){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
  return NodeIterator(this,0);}
  node_iterator node_end() const { return NodeIterator(this,size());}

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */


/* Our incident iterator contains  a pointer to the genrator node and a pointer to a node_list
 whose head is the neighbor currently visited */
  class IncidentIterator: private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() :node_(nullptr),nl_(nullptr){
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Edge operator*() const { 
     //We separate the case where the list is empty
     if ((*nl_).size()==0){return Edge();}
     Node n=*node_;
     node_type n2=(*nl_).front();
     size_type edge_index=(*nl_).index_front();
     return Edge(n.index_,n2.index_,edge_index,n.graph_);}
    IncidentIterator& operator++(){
        Node n = *node_;
        if ((*nl_).size()==0){return *this;}
        // We simply replace the list by its tail
        nl_=(*nl_).next();
        if ((*nl_).size()==0){return *this;}
        if ((*(n.graph_)).valid_edges_[(*nl_).index_front()]==false){return ++(*this);}
        return *this;}
    bool operator==(const IncidentIterator& ii) const {
    return ((*nl_)==(*ii.nl_))&(*node_==*(ii.node_));}

   private:
    friend class Graph;
       // HW1 #3: YOUR CODE HERE
    const Node* node_;
    node_list * nl_;
    IncidentIterator(const Node * node,node_list * nl) 
    :node_(node),nl_(nl) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */


/* Here, to iterate over edges we simpy need to store a pointer to the graph and the index in to the edges vector*/
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), index_(size_type(-1)){
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
     Edge operator*() const { 
     /* We separate the case where we have passed the number of edges*/
     if (index_>=(*graph_).num_edge_){return Edge();}
     return (*graph_).edges_[index_];}
     EdgeIterator& operator++(){
        size_type num_edge=(*graph_).num_edge_;
        if (index_>=num_edge){ index_=size_type(-1); return *this;}
        else {index_=index_+1;
             if (index_>=num_edge){return *this;}
             //We need to skip invalid edges
             if ((*graph_).valid_edges_[index_]==false){return ++(*this);}
             else {return *this;}}
         return *this;}
     bool operator==(const EdgeIterator& ei) const{ return (graph_==ei.graph_)&(index_==ei.index_);}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph * graph_;
    size_type index_; 
    EdgeIterator (const Graph * graph, size_type index) :graph_(graph), index_(index) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const { 
  for (size_type i=0; i<num_edge_;++i){if (valid_edges_[i]==true){return EdgeIterator(this,i);}}
return EdgeIterator(this,0);}
  edge_iterator edge_end() const { return EdgeIterator(this,num_edge_);}

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  /* We store the points and the nodes in two different vectors*/
  /* We also store the values in a vector*/
  std::vector <Node> nodes_;
  std::vector <size_type> map_small_big_nodes_;
  std::vector <size_type> map_big_small_nodes_;
  std::vector <size_type> degrees_;
  std::vector <bool> valid_nodes_;
  std::vector <Point> points_;
  std::vector <Edge> edges_;
  std::vector <bool> valid_edges_;
  std::vector <node_list> incidents_;
  std::vector <std::vector <bool>> crosses_;
  std::vector <node_value_type> values_;
  std::vector <edge_value_type> edge_values_;
  size_type size_;
  size_type num_removed_node_;
  size_type max_size_;
  size_type num_edge_;
  size_type num_removed_edge_;
  Point null_point;
  node_value_type empty_value;
  edge_value_type empty_edge_value;
  //node_list empty_list=node_list();


/* The function below is meant to initiate the crosses matrix as false only */
  std::vector<std::vector <bool>> init_crosses(size_type n){
     std::vector<std::vector <bool>> output=std::vector<std::vector<bool>> ();
     for (size_type i=0;i<n;i++){
             std::vector<bool> local_output=std::vector<bool>();
     for (size_type j=0; j<n;j++){local_output.push_back(false);}
             output.push_back(local_output);}
     return output;}

/* The function below copies the crosses vector into a bigger one, and update the max_size, note that we
 make it voluntary very big to not have to run this update enough, this ensure a O(1) amortized complexity for 
 add_node*/
   std::vector<std::vector <bool>> update_crosses(){
   std::vector<std::vector <bool>> new_crosses=init_crosses(4*max_size_);
   for (size_type i=0; i< max_size_;i++){
   for (size_type j=0;j<max_size_;j++){new_crosses[i][j]=crosses_[i][j];}}
   max_size_*=4;
   return new_crosses;
   }


/* Implementation of the node list: the list contains its size, its head node, the associated edge index  and a pointer to the next list*/
 struct node_list {
   node_list() :n_(Node()), edge_index_(size_type(-1)),next_(nullptr),size_(0) {}
   Node n_;
   size_type edge_index_;
   node_list * next_;
   size_type size_;
   
   Node front(){return n_;}
   size_type index_front() {return edge_index_;}
   node_list  * next() {return next_;}
   size_type size() {return size_;}

   void push_back(Node  new_n,size_type new_edge_index){
          /* We create a new pointer*/
          node_list * new_node_list = (node_list * ) malloc(sizeof(struct node_list));
          (*new_node_list).set_n_(n_);
          (*new_node_list).set_edge_index_(edge_index_);
          (*new_node_list).set_next_(next_);
          (*new_node_list).set_size(size_);
          next_=new_node_list;
          n_=new_n;
          edge_index_=new_edge_index;
          size_++;}
                   
   bool operator == (const node_list nl){ 
        if (nl.size_==0){return 0==size_;}
        return (n_==nl.n_)&(next_==nl.next_)&(edge_index_==nl.edge_index_);}
   node_list (const Node n, const size_type edge_index, node_list * next) :n_(n), edge_index_(edge_index),next_(next), size_(1+(*next).size_){}
   
  private:
   void set_n_(Node n) {n_=n;}
   void set_edge_index_(size_type edge_index) {edge_index_=edge_index;}
   void set_next_(node_list * next) {next_=next;}
   void set_size(size_type size) {size_=size;}

};

node_list empty_list = node_list();

};

#endif // CME212_GRAPH_HPP
