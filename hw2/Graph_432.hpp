#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

// inspired mainly from peercode/hw0/Graph_276.hpp but also picked up a couple of minute things from other peercode files

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  /** Predeclartion of internal_node type */
  struct internal_node; 
  /** Predeclaration of internal_edge type */
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //


  /** Type of value to be stored in node */
  using node_value_type = V;

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
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
    // HW0: YOUR CODE HERE
  }

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
    Node(): graph_(nullptr), uid_(0) {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch_node().position_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return uid_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return Attribute stored in node:modifiable reference*/
    node_value_type& value(){
      return fetch_node().attribute_;
    }

    /** Return Attribute stored in node:constant reference */
    const node_value_type& value() const{
      return fetch_node().attribute_;
    }

    /** Return number of incident edges */
    size_type degree() const{
      return fetch_node().incident_edges_.size();
    }

    /** Starting pointer of incident iterator */ 
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uid_, fetch_node().incident_edges_.begin());
    }
    /** Ending pointer of incident iterator  */ 
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uid_, fetch_node().incident_edges_.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((this->graph_ == n.graph_) && (this->index() ==  n.index()));
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
      // HW0: YOUR CODE HERE
       if(this->graph_ != n.graph_)
         {
             return(this->graph_ < n.graph_);
         }
         else
         {
             return(this->uid_ < n.uid_);
         }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* const graph_; //graph_ is a pointer to GRAPH
    size_type uid_; //Index of current node

    // Private Constructor
    Node(const Graph* graph, size_type uid): graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

    // Use for fetching node stored in graph
    internal_node& fetch_node() const
    {
      return graph_->nodes_[uid_];
    }
 
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();
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
  Node add_node(const Point& position, const node_value_type& attribute = node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type nnodes = num_nodes();
    nodes_.push_back(internal_node(position, attribute));
    return node(nnodes);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (this == n.graph_ && n.index() < num_nodes());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this,i);
  }

  //
  // EDGES
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge(): graph_(nullptr), uid_(0) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if(reverse_flag_)
      {
        return (graph_->node(fetch_edge().uid_node2_));
      }
      else
      {
        return (graph_->node(fetch_edge().uid_node1_));
      }

    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if(reverse_flag_)
      {
        return (graph_->node(fetch_edge().uid_node1_));
      }
      else
      {
        return (graph_->node(fetch_edge().uid_node2_));
      }

    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (fetch_edge() == e.fetch_edge()); //Operator Defined in Struct internal_edge
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (fetch_edge() < e.fetch_edge()); //Operator Defined in Struct internal_edge
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

   Graph* const graph_; //pointer to GRAPH
   size_type uid_; //Index of Current Edge
   bool reverse_flag_; //If true then order of nodes is reversed (Required for incident edge iterator
   
   //Depreceated Code
   //size_type uidNode1_; 
   //size_type uidNode2_;
   
   //Private Constructor
   Edge(const Graph* graph, size_type uid, bool reverse_flag = false): graph_(const_cast<Graph*>(graph)), uid_(uid), reverse_flag_(reverse_flag){
   }

   // Use for fetching edge stored in graph
    internal_edge& fetch_edge() const
    {
      return graph_->edges_[uid_];
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const{
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, false);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if(a.graph_ == b.graph_ && this == a.graph_)
    {
        size_type node1_i = std::min(a.index(),b.index());
        size_type node2_i = std::max(a.index(),b.index());
        for(typename std::vector<internal_edge>::const_iterator i = edges_.begin(); i!=edges_.end(); ++i)
        {
            if(i->uid_node1_ == node1_i && i->uid_node2_ == node2_i)
            {
                return true;
            }
        }

    } 
    return false;
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == b.graph_ && this == a.graph_);
    
    size_type node1_i = std::min(a.index(),b.index());
    size_type node2_i = std::max(a.index(),b.index());
    size_type nedges = num_edges();
    
    bool reverse_flag = a.index() > b.index();
    //Check if edge exists already
    for(std::vector<size_type>::const_iterator i = a.fetch_node().incident_edges_.begin(); i!=a.fetch_node().incident_edges_.end(); ++i)
    {
        if(edges_[*i].uid_node1_ == node1_i && edges_[*i].uid_node2_ == node2_i)
        {
            return Edge(this, *i, reverse_flag);
        }
    }

   //Add Edge
   edges_.push_back(internal_edge(a.index(), b.index()));
   a.fetch_node().incident_edges_.push_back(nedges);
   b.fetch_node().incident_edges_.push_back(nedges); 
   return Edge(this, nedges, reverse_flag);

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator(): graph_(nullptr), uid_(0) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Dereference pointer to get node */
    Node operator*() const{
      return graph_->node(uid_);
    }

    /** Increment pointer to next stored node */
    NodeIterator& operator++() {
      uid_++;
      return *this;
    }

    /** Are two iterators equal? */ 
    bool operator==(const NodeIterator& n) const{
      return ((this->graph_ == n.graph_) && (this->uid_ ==  n.uid_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_; //Pointer to parent graph
    size_type uid_; //id of node

    //Private Constructor
    NodeIterator(const Graph* graph, size_type uid): graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Pointer to first node*/ 
  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }
 
  /** Pointer to last node*/ 
  node_iterator node_end() const{
    return NodeIterator(this,this->nodes_.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(): graph_(nullptr), uid_node_(0), edge_itr_() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference iterator to get edge*/ 
    Edge operator*() const{
      Edge temp = graph_->edge(*edge_itr_);
      temp.reverse_flag_ = graph_->node(uid_node_) == temp.node2(); //Reverse if node spawning is stored in node2
      return temp;
    }

    /** Increment iterator to next stored edge*/
    IncidentIterator& operator++(){
      edge_itr_++;
      return *this;
    }

    /** Are two iterators equal? */ 
    bool operator==(const IncidentIterator& n) const{
      return ((this->graph_ == n.graph_) && (this->uid_node_ == n.uid_node_) && (this->edge_itr_ == n.edge_itr_));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    Graph* const graph_; //Pointer to parent graph
    size_type uid_node_; //id of node
    typename std::vector<size_type>::const_iterator edge_itr_; //incident edge iterator used in has_edge()

    //Private Constructor
    IncidentIterator(const Graph* graph, size_type uid_node, typename std::vector<size_type>::const_iterator edge_itr): graph_(const_cast<Graph*>(graph)), uid_node_(uid_node), edge_itr_(edge_itr){
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(): graph_(nullptr), uid_(0) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Dereference to get edge*/ 
    Edge operator*() const{
      return graph_->edge(uid_);
    }
    
    /** Increment to next stored edge  */ 
    EdgeIterator& operator++(){
      uid_++;
      return *this;
    }
    
    /** Are two iterators equal*/ 
    bool operator==(const EdgeIterator& n) const{
      return ((this->graph_ == n.graph_) && (this->uid_ ==  n.uid_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* const graph_; //Pointer to parent graph
    size_type uid_; //id of edge

    //Private Constructor
    EdgeIterator(const Graph* graph, size_type uid): graph_(const_cast<Graph*>(graph)), uid_(uid){
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Pointer to first edge*/ 
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0);
  }

  /** */ edge_iterator edge_end() const{
    return EdgeIterator(this, this->edges_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  //internal type for node elements
  struct internal_node{
    //const size_type uid_node_; //Node id //Depreceated Code
    const Point position_; //Position of Node
    node_value_type attribute_; //Attribute stored within node
    std::vector<size_type> incident_edges_; //ids of incident edges

    //constructor
    internal_node(const Point& position, node_value_type attribute): position_(position), attribute_(attribute), incident_edges_(){
    }
      
  };

 //internal type for edge elements
 struct internal_edge : private totally_ordered<internal_edge> {
     //const size_type uid_edge_; //Edge id //Depreceated Code
     const size_type uid_node1_; //Corresponding Node 1 id
     const size_type uid_node2_; //Corresponding Node 2 id
     
     //constructor
     internal_edge(const size_type uid_node1, const size_type uid_node2): uid_node1_(std::min(uid_node1,uid_node2)), uid_node2_(std::max(uid_node1,uid_node2)){
     } //We want node1_ to have smaller index and node2_ to have larger index as this is easy to code and leads to an unique representation of an edge
     
     //Overloading Operators for Iterator Purposes
     bool operator==(const internal_edge& e) const
     {
         return (this->uid_node1_ == e.uid_node1_ && this->uid_node2_ == e.uid_node2_);
     }

     bool operator<(const internal_edge &e) const
     {
         if(this->uid_node1_ != e.uid_node1_)
         {
             return(this->uid_node1_ < e.uid_node1_);
         }
         else
         {
             return(this->uid_node2_ < e.uid_node2_);
         }
     }
 };

std::vector<internal_node> nodes_; //Stores Nodes
std::vector<internal_edge> edges_; // Stores Edges
};

#endif // CME212_GRAPH_HPP
