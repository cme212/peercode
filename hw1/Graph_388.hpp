#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <unordered_map>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  using node_value_type = V;

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
  Graph(): node_vec(), edge_vec() {
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
  class Node : private totally_ordered<Node> {
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
    Node() {
      // HW0: YOUR CODE HERE
      graph_ptr = nullptr;
      idx_ = 0; 
    }


    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // dereference graph pointer,
      // find internal node with the node's idx
      // and return the Point object associated with it
      return graph_ptr->node_vec[idx_].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      internal_node& n = (*this->graph_ptr).node_vec[this->idx_]; 
      return n.node_value;
    }

    const node_value_type& value() const {
      const internal_node& n = (*this->graph_ptr).node_vec[this->idx_]; 
      return n.node_value;
    }
    
    size_type degree() const {
      return (*this->graph_ptr).node_vec[idx_].adj_nodes.size(); 
    }
    
    incident_iterator edge_begin() const {
      return incident_iterator(graph_ptr, idx_, 0);
    }

    incident_iterator edge_end() const{
      return incident_iterator(graph_ptr, idx_,degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //(void) n;          // Quiet compiler warning
      return (graph_ptr == 
        n.graph_ptr) && (idx_ == n.idx_);
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
      //(void) n;           // Quiet compiler warning
      // only compare nodes in the same graph
      return (n.graph_ptr == graph_ptr) && (idx_ < n.idx_) ;  
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_ptr; // pointer to Graph class that contains Nodes
    size_type idx_; // unique identification/index for a node

    /*Private Constructor*/
    Node(const Graph* g, size_type idx)
        : graph_ptr(const_cast<Graph*>(g)), idx_(idx) {
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_vec.size();
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
  Node add_node(const Point& position, const node_value_type& v = node_value_type() ) {
    // HW0: YOUR CODE HERE
    //(void) position;      // Quiet compiler warning
    
    // Insert node in node vec
    internal_node new_node;
    new_node.position = position;  
    std::vector<adjacent_nodes> a_nodes; // Initialize vector of adjacent nodes
    new_node.adj_nodes = a_nodes; 
    new_node.node_value = v;
    node_vec.push_back(new_node);
    
    //(void) position;      // Quiet compiler warning
    // Return Node that points to new node
    return Node(this, num_nodes()-1);      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    //(void) n;            // Quiet compiler warning
    // check  if Node's idx is smaller than the size
    return (n.graph_ptr == this) && (n.idx_ < size());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    // create new Node obj with idx i;
    return Node(this, i);
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ptr = nullptr;
      idx_ = 0;
      d_idx_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_ptr, d_idx_);       
    }

    /** Return the other node of this Edge */
    Node node2() const {   
      // HW0: YOUR CODE HERE
      if(graph_ptr->edge_vec[idx_].first == d_idx_){
        return Node(graph_ptr, graph_ptr->edge_vec[idx_].second);
      }else{
        return Node(graph_ptr, graph_ptr->edge_vec[idx_].first);
      }    
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return (graph_ptr == e.graph_ptr) && 
        (idx_ == e.idx_); 
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return (graph_ptr == e.graph_ptr) && (idx_ < e.idx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graph_ptr; 
    size_type idx_;
    size_type d_idx_; // Index for dominant node

    /* Private Constructor */
    Edge(const Graph* g, size_type idx, size_type d_n)
        : graph_ptr(const_cast<Graph*>(g)), idx_(idx), d_idx_(d_n) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //(void) i;             // Quiet compiler warning
    return Edge(this, i, edge_vec[i].first);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //(void) a; (void) b;   // Quiet compiler warning
     // HW0: YOUR CODE HERE
    std::vector<adjacent_nodes> a_nodes = a.graph_ptr->node_vec[a.idx_].adj_nodes; 

    for (size_type i = 0; i < a_nodes.size(); ++i){
      // Check if b is an adjacent node to a 
        if(a_nodes[i].connected_node == b.idx_){
          return true;
        }
    }

    //(void) a; (void) b;   // Quiet compiler warning
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
    //(void) a, (void) b;   // Quiet compiler warning
    
   // If edge exists, then find and return it 
    assert(a.graph_ptr == b.graph_ptr);

    std::vector<adjacent_nodes> a_nodes = a.graph_ptr->node_vec[a.idx_].adj_nodes; 

    for (size_type i = 0; i < a_nodes.size(); ++i){
      // Check if b is an adjacent node to a 
        if(a_nodes[i].connected_node == b.idx_){
          return Edge(this, a_nodes[i].edge_id, a.idx_);
        }
    }

    // else add new edge 
    std::pair<size_type, size_type> new_edge (a.idx_, b.idx_);

    edge_vec.push_back(new_edge);
    
    // create adjacent node objects for a and b
    adjacent_nodes a1, a2;

    a1.edge_id = num_edges() - 1;
    a1.connected_node = b.idx_;
    a.graph_ptr->node_vec[a.idx_].adj_nodes.push_back(a1);

    a2.edge_id = num_edges() - 1;
    a2.connected_node = a.idx_;
    b.graph_ptr->node_vec[b.idx_].adj_nodes.push_back(a2);

    //(void) a, (void) b;   // Quiet compiler warning
    return Edge(this, num_edges()-1, a.idx_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_vec.clear();
    edge_vec.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      graph_ptr = nullptr;
      nidx_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const{
      return Node(graph_ptr, nidx_);
    }

    NodeIterator& operator++(){
      ++nidx_;
      return *this;
    }

    bool operator==(const NodeIterator& it) const {
      return graph_ptr == it.graph_ptr and 
        nidx_ == it.nidx_;
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_ptr; 
    size_type nidx_;

    /*Private Constructor*/
    NodeIterator(const Graph* g, size_type nidx)
        : graph_ptr(const_cast<Graph*>(g)), nidx_(nidx) {
    }
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    unsigned int begin_i = 0;
    return NodeIterator(this, begin_i);
  }
  
  node_iterator node_end() const {
    return NodeIterator(this,size()); 
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graph_ptr = nullptr;
      n_id_ = 0;
      a_id_ = 0; 
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      std::vector<adjacent_nodes> a_nodes = graph_ptr->node_vec[n_id_].adj_nodes; 
      size_type e_id = a_nodes[a_id_].edge_id; 
      return(Edge(graph_ptr,e_id, n_id_));
    }

    IncidentIterator& operator++(){
      a_id_++; 
      return(*this);
    }

    bool operator==(const IncidentIterator& a) const{
       // Edges must belong to the same graph
      return (a.graph_ptr == graph_ptr) && 
      (a.n_id_ == n_id_) && (a.a_id_ == a_id_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_ptr;
    size_type n_id_;
    size_type a_id_; 
    
    IncidentIterator(const Graph* g_ptr, size_type n_id, size_type a_id)
        : graph_ptr(const_cast<Graph*>(g_ptr)), n_id_(n_id), a_id_(a_id) {
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graph_ptr = nullptr;
      eidx_ = 0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(graph_ptr, eidx_, graph_ptr->edge_vec[eidx_].first);
    }

    EdgeIterator& operator++(){
      ++eidx_;
      return *this;
    }

    bool operator==(const EdgeIterator& it) const {
      return graph_ptr == it.graph_ptr && 
        eidx_ == it.eidx_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_ptr;
    size_type eidx_;

    /*Private Constructor*/
    EdgeIterator(const Graph* g, size_type eidx)
        : graph_ptr(const_cast<Graph*>(g)), eidx_(eidx) {
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
   edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
   }
  
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }


 private:
 // struct to hold adjacent nodes to a particular node
  struct adjacent_nodes{
     size_type edge_id;//edge connecting the node
     size_type connected_node;//adjacent node
  };

  // struct to hold internal node 
  struct internal_node {
    Point position; 
    node_value_type node_value;
    std::vector<adjacent_nodes> adj_nodes; //Store adjacent nodes 
  };

  // vector containing node indices as pairs
  std::vector<internal_node> node_vec;

  // Vector with node index pairs
  std::vector<std::pair<size_type,size_type>> edge_vec; 

  // HW0 YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
