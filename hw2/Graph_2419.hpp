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

template <typename V, typename E> 

class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  typedef V node_value_type;
  typedef E edge_value_type;
  //
  // PUBLIC TYPE DEFINITIONS
  //

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
    Node(): graph_ptr(nullptr), idx_ (0){
    }

    Point& position(){
      // dereference graph pointer,
      // find internal node with the node's idx
      // and return the Point object associated with it
      return graph_ptr->node_vec[graph_ptr->valid_nodes[idx_]].position;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // dereference graph pointer,
      // find internal node with the node's idx
      // and return the Point object associated with it
      return graph_ptr->node_vec[graph_ptr->valid_nodes[idx_]].position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){

      internal_node& n = (*this->graph_ptr).node_vec[graph_ptr->valid_nodes[idx_]]; 
      return n.node_value;
    }

    const node_value_type& value() const {
      const internal_node& n = (*this->graph_ptr).node_vec[graph_ptr->valid_nodes[idx_]]; 
      return n.node_value;
    }
    
    size_type degree() const {
      return graph_ptr->node_vec[graph_ptr->valid_nodes[idx_]].adj_nodes.size(); 
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
      // only compare nodes in the same graph
      return ((n.graph_ptr < graph_ptr) || (idx_ < n.idx_));  
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_ptr; // pointer to Graph class that contains Nodes
    size_type idx_; // user facing index for a node

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
    return valid_nodes.size();
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
    
    // Insert node in node vec
    internal_node new_node;
    new_node.position = position;  
    std::vector<adjacent_nodes> a_nodes; // Initialize vector of adjacent nodes
    new_node.adj_nodes = a_nodes; 
    new_node.node_value = v;
    new_node.n_idx_ = valid_nodes.size(); // Push user facing idx
    node_vec.push_back(new_node);
    
    // Add uid of new node to valid nodes
    valid_nodes.push_back(node_vec.size() - 1);

    // Return Node that points to new node
    return Node(this, valid_nodes.size()-1);      
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // check  if Node's idx is smaller than the size
    return ((n.graph_ptr == this) && (n.idx_ < size()));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
    Edge(): graph_ptr(nullptr), idx_(0), n1_idx_(0), n2_idx_(0) {
    }

    size_type get_edge_idx() const{
      return idx_; // Return user-facing idx of this edge
    }
    
    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_ptr, n1_idx_);       
    }

    /** Return the other node of this Edge */
    Node node2() const {   
      // HW0: YOUR CODE HERE
      return Node(graph_ptr, n2_idx_);   
    }

    edge_value_type& value(){
      internal_edge& e = (*this->graph_ptr).edge_vec[((*this->graph_ptr).valid_edges)[idx_]]; 
      return e.edge_value;
    }

    const edge_value_type& value() const {
      const internal_edge& e = (*this->graph_ptr).edge_vec[((*this->graph_ptr).valid_edges)[idx_]]; 
      return e.edge_value;
    }

    double length() const{
      // Retrieve nodes associated with edge 
      auto n1 = node1();
      auto n2 = node2();
      // Return L2 norm of edge
      return norm(n1.position() - n2.position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
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
      //HW0: YOUR CODE HERE
      return ((graph_ptr < e.graph_ptr) || (idx_ < e.idx_));
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
    size_type n1_idx_; // Index for node1
    size_type n2_idx_; // Index for node2

    /* Private Constructor */
    Edge(const Graph* g, size_type idx, size_type d_n, size_type s_n)
        : graph_ptr(const_cast<Graph*>(g)), idx_(idx), n1_idx_(d_n), n2_idx_(s_n){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return valid_edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, i, 
      edge_vec[valid_edges[i]].node_pair.first, 
      edge_vec[valid_edges[i]].node_pair.second);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
     // HW0: YOUR CODE HERE
    std::vector<adjacent_nodes> a_nodes = node_vec[valid_nodes[a.idx_]].adj_nodes; 

    for(unsigned int i =0; i < a.degree(); ++i){
      if(valid_nodes[b.idx_] == a_nodes[i].connected_node_uid){
        return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& e = edge_value_type()) {
    // HW0: YOUR CODE HERE
    
   // If edge exists, then find and return it 
    assert(a.graph_ptr == b.graph_ptr);

    std::vector<adjacent_nodes> a_nodes = node_vec[valid_nodes[a.idx_]].adj_nodes; 
    assert(a.degree() == a_nodes.size());

    for (size_type i = 0; i < a.degree(); ++i){
      // Check if b is an adjacent node to a 
        if(a_nodes[i].connected_node_uid == valid_nodes[b.idx_]){
          return Edge(this, edge_vec[a_nodes[i].edge_uid].e_idx_, a.idx_, b.idx_);
        }
    }

    // add edge uid to valid edges
    valid_edges.push_back(edge_vec.size());

    // else add new edge 
    internal_edge new_edge; 
    
    new_edge.node_pair.first = a.idx_;
    new_edge.node_pair.second = b.idx_;
    new_edge.edge_value = e;
    // Add user facing index
    new_edge.e_idx_ = valid_edges.size()-1;

    edge_vec.push_back(new_edge);
    
    // create adjacent node objects for a and b
    adjacent_nodes a1, a2;

    // Push UID of edge 
    a1.edge_uid = edge_vec.size()-1;
    a1.connected_node_uid = (b.graph_ptr->valid_nodes)[b.idx_];
    a.graph_ptr->node_vec[a.graph_ptr->valid_nodes[a.idx_]].adj_nodes.push_back(a1);

    a2.edge_uid = edge_vec.size()-1;
    a2.connected_node_uid = (a.graph_ptr->valid_nodes)[a.idx_];
    b.graph_ptr->node_vec[b.graph_ptr->valid_nodes[b.idx_]].adj_nodes.push_back(a2);

    return Edge(this, valid_edges.size() - 1, a.idx_, b.idx_);
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
    valid_edges.clear();
    valid_nodes.clear();
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
    Node operator*() const{
      return Node(graph_ptr, nidx_);
    }

    NodeIterator& operator++(){
      ++nidx_;
      return *this;
    }

    bool operator==(const NodeIterator& it) const {
      return (graph_ptr == it.graph_ptr) && 
        (nidx_ == it.nidx_);
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
    return NodeIterator(this, 0);
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
      n_id_ = 0; // user facing node id
      a_id_ = 0; // Idx in adj nodes vector
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const{
      std::vector<adjacent_nodes> a_nodes = graph_ptr->node_vec[graph_ptr->valid_nodes[n_id_]].adj_nodes; 
      size_type e_uid = a_nodes[a_id_].edge_uid; 
      size_type n2_id = a_nodes[a_id_].connected_node_uid; 
      return(Edge(graph_ptr,(graph_ptr->edge_vec[e_uid]).e_idx_, n_id_, (graph_ptr->node_vec)[n2_id].n_idx_));
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
      eidx_ = 0; // User facing ID
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(graph_ptr, 
                  eidx_, 
                  graph_ptr->edge_vec[(graph_ptr->valid_edges)[eidx_]].node_pair.first,
                  graph_ptr->edge_vec[(graph_ptr->valid_edges)[eidx_]].node_pair.second);
    }

    EdgeIterator& operator++(){
      ++ eidx_;
      return *this;
    }

    bool operator==(const EdgeIterator& it) const {
      return (graph_ptr == it.graph_ptr) && 
        (eidx_ == it.eidx_);
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

 public:

  /** Remove an edge of the graph given an Edge object.
   *  returning 1 if successful, 0 if not.
   * @pre Edge @e e is a valid edge of this graph
   * @return size_type of 1 if Edge @e was 
   * succesfully removed, 0 if Edge @e did not exist. 
   * @post has_edge( @e e ) == false
   * @post If old has_edge( @e e ),
   *           new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   * 
   * Invalidates Edge @e, pops element from valid_edges[e.get_edge-idx()] and 
   * removes corresponding adjacent node objects from e.node1() and e.node2().
   * 
   * Complexity: No more than O(num_edges()). More like O(e.node1().degree() + e.node2().degree())
   */
  size_type remove_edge(const Edge& e){
    // Check if the graph has edge
    if(!has_edge(e.node1(), e.node2())){
      return 0;
    }

    size_type idx = e.get_edge_idx(); // Get user facing idx
    size_type e_uid = valid_edges[idx]; 

    // Remove edge from adjacent nodes data structures
  
    size_type n1_idx = e.node1().index();
    
    size_type n2_idx = e.node2().index();

    for(size_type i = 0; i < e.node1().degree(); ++i){
      // Check if this is the edge we're trying to remove

      // Check if UIDs match
      if(node_vec[valid_nodes[n1_idx]].adj_nodes[i].edge_uid == e_uid){
        std::iter_swap(node_vec[valid_nodes[n1_idx]].adj_nodes.begin() + i, 
                        node_vec[valid_nodes[n1_idx]].adj_nodes.end() - 1);   
        node_vec[valid_nodes[n1_idx]].adj_nodes.pop_back(); // Remove last element
        break;
      }
    } 

    for(size_type i = 0; i < e.node2().degree(); ++i){
      // Check if this is the edge we're trying to remove
      if(node_vec[valid_nodes[n2_idx]].adj_nodes[i].edge_uid == e_uid){
        // Swap adj node to be removed with last element
        std::iter_swap(node_vec[valid_nodes[n2_idx]].adj_nodes.begin() + i, 
                        node_vec[valid_nodes[n2_idx]].adj_nodes.end() - 1);   
        
        node_vec[valid_nodes[n2_idx]].adj_nodes.pop_back(); // Remove last element
        break;
      } 
    }

    // Swap element to be deleted with last element
    std::iter_swap(valid_edges.begin() + idx, valid_edges.end() - 1); 
    valid_edges.pop_back(); // Remove last element
    // Change idx of corresponding internal edge
    edge_vec[valid_edges[idx]].e_idx_ = idx;

    return 1; 
  }

  /** Remove an edge of the graph given two adjacent nodes @n n1 and @n n2
   *  returning 1 if successful, 0 if not. Calls remove_edge(@e e) where
   * @e e is the incident edge connecting the 2 nodes.
   * @pre Nodes @n n1 and @n2 are valid nodes of this graph
   * @return size_type of 1 if edge 
   * succesfully removed, 0 if nodes are not adjacent.
   * @post has_edge( @e e ) == false
   * @post If old has_edge( @e e ),
   *           new num_edges() == old num_edges() - 1.
   *       Else, new num_edges() == old num_edges().
   * 
   * Invalidates Edge @e, pops element from valid_edges[e.get_edge-idx()] and 
   * removes corresponding adjacent node objects from @n n1 and @n n2.
   * 
   * Complexity: No more than O(num_edges()). More like O(n1.degree() + n2.degree())
   */

  size_type remove_edge(const Node& n1, const Node& n2){
      for(auto i = n1.edge_begin(); i != n1.edge_end(); ++i){
        // Find edge that is incident with n2
        if((*i).node2() == n2){
          size_type result = remove_edge(*i); // Remove incident edge
          return result;
        }
      }
      // else return 0
      return 0;
  }
  
   /** Remove an edge from the graph using an edge iterator e_it returning 1 if successful, 
   * 0 if not. Calls remove_edge(@e e) where @e e = (*e_it.)
   * @return edge_iterator pointing to first valid edge in graph
   * @post new num_edges() == old num_edges() - 1.
   * @post reduces degree of nodes incident to @e by 1. 
   * Invalidates Edge @e, pops element from valid_edges[@e.get_edge-idx()] and 
   * removes corresponding adjacent node objects from e.node1() and e.node2().
   *
   * Complexity: No more than O(num_edges()). More like O(e.node1().degree() + e.node2().degree())
   */

  edge_iterator remove_edge(edge_iterator e_it){
    // Dereference edge iterator to return edge and remove it
    remove_edge(*e_it);
    return this->edge_begin();
  }

  /** Remove a node of the graph given a Node object.
   *  returning 1 if successful, 0 if not.
   * @pre Node @n n is a valid node of this graph
   * @return size_type of 1 if Node @n was 
   * succesfully removed, 0 if Node @n did not exist. 
   * @post has_node( @n n ) == false
   * @post If old has_node( @n n ),
   *           new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   * 
   * Invalidates Node @n n, pops element from valid_nodes[n.index()]
   * and removes all incident edges to @n n by calling remove_edge().
   * 
   * Complexity: No more than O(n.degree() ^ 2), which by sparsity is less than
   * O(num_nodes()).
   */

  size_type remove_node(const Node& n){
    // Check if node exists
    if(!has_node(n)){
      return 0;
    }

    size_type idx = n.index();

    // Remove all incident edges
    auto iit = n.edge_begin();
    while(n.degree() > 0){
      remove_edge(*iit);
    }

    // Swap element to be deleted with last element
    std::iter_swap(valid_nodes.begin() + n.index(), valid_nodes.end() - 1); 
    valid_nodes.pop_back(); // Remove last element
    // Change idx of corresponding internal node
    node_vec[valid_nodes[idx]].n_idx_ = idx;

    return 1;

  }

 /** Remove a node of the graph given a Node iterator object.
   *  returning 1 if successful, 0 if not.
   * @pre Node iterator n_it is a valid iterator to a node
   *  of this graph
   * @return node_iterator pointing to first node in graph
   * succesfully removed, 0 if not.
   * @post has_node( @n (*n_it) ) == false
   * @post If old has_node( @n (*n_it) ),
   *           new num_nodes() == old num_nodes() - 1.
   *       Else, new num_nodes() == old num_nodes().
   * 
   * Invalidates Node @n (*n_it), pops element from valid_nodes[(*n_it).index()]
   * and removes all incident edges to @n (*n_it) by calling remove_edge().
   * 
   * Complexity: No more than O(n.degree() ^ 2), which by sparsity is less than
   * O(num_nodes()).
   */
  node_iterator remove_node(node_iterator n_it){
    remove_node(*n_it);
    return this->node_begin();
  }

 private:
 // struct to hold adjacent nodes to a particular node
  struct adjacent_nodes{
     size_type edge_uid;//uid of edge connecting the node
     size_type connected_node_uid;//uid of adjacent node
  };

  // struct to hold internal node 
  struct internal_node {
    Point position; 
    node_value_type node_value;
    std::vector<adjacent_nodes> adj_nodes; //Store adjacent nodes 
    size_type n_idx_; // User facing idx
  };

  // Struct to hold internal edge
  struct internal_edge{
    std::pair<size_type,size_type> node_pair;
    edge_value_type edge_value;
    size_type e_idx_; // User facing idx
  };

  // vector containing node indices as pairs
  std::vector<internal_node> node_vec;

  // Vector with valid nodes
  std::vector<size_type> valid_nodes;

  // Vector with node index pairs
  std::vector<internal_edge> edge_vec; 

  // Vector with valid edges
  std::vector<size_type> valid_edges;

};

#endif // CME212_GRAPH_HPP
