#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <iostream>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph {
 private:
  
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct internal_node;
  struct adjacent_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Allow the Nodes to support a user-specified value. */
  using node_value_type = V;
  using edge_value_type = E;

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

  /** Construct an empty graph */
  Graph() :
    nodes_(), edges_(), edge_values_() {
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
      graph_ = nullptr;
      id_ = 0; 
    }

    /** Return this node's position */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_-> nodes_[id_].position_;
    }
    
    /** Return a reference to this node's position */
    Point& position() {
      return graph_->nodes_[id_].position_;
    }

    /** Return this node's index, a number in the range [0, graph_size) */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's user-specified value */
    node_value_type& value(){
      return graph_ -> nodes_[id_].value_;
    }

    /** Return a constant reference to this node's user-specified value */
    const node_value_type& value() const{
      return graph_ -> nodes_[id_].value_;
    }

    /** Return the number of incident edges to this node */
    size_type degree() const{
      return graph_ -> nodes_[id_].adj_nodes_.size();
    }

    /** Return begin incident iterator for this node */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, id_, 0);
    }

    /** Return end incident iterator for this node */
    incident_iterator edge_end () const{
      return IncidentIterator(graph_, id_, degree()); 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& node) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == node.graph_) && (id_ == node.id_));
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node& node) const {
      // HW0: YOUR CODE HERE
      assert(graph_ == node.graph_);
      return id_ < node.id_;
    }
  

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type id_;

    Node(const graph_type* graph, size_type id)
    : graph_(const_cast<graph_type*>(graph)), id_(id) {    
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    //Create New node, add its attributes and push to vector
    internal_node new_node;
    new_node.position_ = position; 
    new_node.value_ = value;
    nodes_.push_back(new_node);

    return Node(this, size() - 1);  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return (n.id_ < size() && n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = nullptr; 
      index_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_);    
    }

    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position()-node2().position());
    }

    /** Return current edge's value */
    edge_value_type& value(){
      return graph_ -> edge_values_[index_];
    }

    /** Return constant reference to current node's value */
    const edge_value_type& value() const{
      return graph_ -> edge_values_[index_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& edge) const {
      //HW0: YOUR CODE HERE
      return ((index_ == edge.index_) && (graph_ == edge.graph_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& edge) const {
      //HW0: YOUR CODE HERE
      assert(graph_ == edge.graph_);
      return index_ < edge.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_; 
    size_type index_;
    size_type node1_;
    size_type node2_; 

    Edge(const graph_type* graph, size_type index, size_type node1, size_type node2)
    : graph_(const_cast<graph_type*>(graph)), index_(index), node1_(node1), node2_(node2){
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
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
    size_type node1_ = edges_[i].first;
    size_type node2_ = edges_[i].second;
    return Edge(this, i, node1_, node2_);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //Pull the adjacent node list for node a and check if node b is in the list
    std::vector<adjacent_node> adj_nodes_ = nodes_[a.id_].adj_nodes_;
    for (unsigned i = 0; i < adj_nodes_.size(); ++i){
      if (b.id_ == adj_nodes_[i].adj_node_){
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    // If edge exists, find edge id and return edge
    std::vector<adjacent_node> adj_nodes_ = nodes_[a.id_].adj_nodes_;
    for (unsigned i = 0; i < adj_nodes_.size(); ++i){
      if (b.id_ == adj_nodes_[i].adj_node_){
        size_type index_ = adj_nodes_[i].edge_index_;
        return edge(index_);
      }
    }
    // Add new edge to edge list and create initialize adjacent nodes
    edges_.push_back(std::make_pair(a.id_, b.id_));
    adjacent_node adj1, adj2; 

    // Assign adj node id and edge id and add to vector of adj nodes
    adj1.adj_node_ = b.id_;
    adj1.edge_index_ = num_edges() - 1;
    adj2.adj_node_ = a.id_;
    adj2.edge_index_ = num_edges() - 1;
    nodes_[a.id_].adj_nodes_.push_back(adj1);
    nodes_[b.id_].adj_nodes_.push_back(adj2);

    return edge(num_edges() - 1);

    }


  /** Remove node from this graph.
   * @post num_nodes() == previous num_nodes() - 1
   * @post g.node(i).index() == i for all i with 0 <= i < g.num_nodes()
   * @post  g.node(n.index()) == n
   * 
   * Invalidates an outstanding Node object and all incident Edges.
   */
  
  size_type remove_node(const Node& node) {
    (void) node;
    return 0;
  }

  node_iterator remove_node(node_iterator n_it) {
    (void) n_it;
    return NodeIterator();
  }

  /** Remove edge from this graph.
   * @post num_edges() == previous num_edges() - 1
   * @post g.edge(i).index < num_edges()
   * @post return number of unique edges
   *
   * Invalidates an outstanding Edge object.
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (!has_edge(a, b)) {
      return 0;
    }
    // Find edge index
    size_type index_ = 0;
    std::vector<adjacent_node> adj_nodes_ = nodes_[a.id_].adj_nodes_;
    for (unsigned i = 0; i < adj_nodes_.size(); ++i){
      if (b.id_ == adj_nodes_[i].adj_node_){
        index_ = adj_nodes_[i].edge_index_;
        // Remove adjacent node indicating presence of Edge
        adj_nodes_.erase(adj_nodes_.begin()+i);
      }
    }
    // Also go to adjacent node and remove the node indicating presence of Edge
    std::vector<adjacent_node> adj_adj_nodes_ = nodes_[b.id_].adj_nodes_;
    for (unsigned i = 0; i < adj_adj_nodes_.size(); ++i){
      if (a.id_ == adj_adj_nodes_[i].adj_node_){
        // Remove adjacent node indicating presence of Edge
        adj_adj_nodes_.erase(adj_adj_nodes_.begin()+i);
      }
    }
      
    
    // Swap this edge with the edge at the end of the vector of edges
    std::swap(edges_[index_], edges_[num_edges()-1]);
    // Reassign edge index
    //edges_[index_].index() = index_;
    // Pop the the last edge from the vector of edges
    edges_.pop_back();

    return 1;
  }

  /** Two other functions to remove edge, each calling the preceding
    * remove_edge function which takes in two nodes and arguments.
    */  
  size_type remove_edge(const Edge& edge) {
    node_type a = edge.node1();
    node_type b = edge.node2();
    return remove_edge(a,b);
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    node_type a = (*e_it).node1();
    node_type b = (*e_it).node2();
    return remove_edge(a,b);
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator */
    NodeIterator() {
      graph_ = nullptr;
      node_id_ = 0;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    Node operator*() const {
      return Node(graph_, node_id_);
    }

    NodeIterator& operator++(){
      node_id_++; 
      return *this;
    }

    bool operator==(const NodeIterator& iter) const{
      return ((graph_ == iter.graph_) && (node_id_ == iter.node_id_));
    } 

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type node_id_;

    NodeIterator(const graph_type* graph, size_type node_id)
    : graph_(const_cast<graph_type*>(graph)), node_id_(node_id) {    
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }

  node_iterator node_end() const{
    return NodeIterator(this, num_nodes());
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

    /** Construct an invalid IncidentIterator */
    IncidentIterator() {
      graph_ = nullptr;
      node_id_ = 0;
      adj_id_ = 0;

    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const{
      //Retrieve edge id and adj node id
      size_type edge_id_ = graph_ -> nodes_[node_id_].adj_nodes_[adj_id_].edge_index_;
      size_type adj_node_ = graph_ -> nodes_[node_id_].adj_nodes_[adj_id_].adj_node_;
      return Edge(graph_, edge_id_, node_id_, adj_node_);
    }

    IncidentIterator& operator++(){
      adj_id_++; 
      return *this;
    }

    bool operator==(const IncidentIterator& iter) const{
      return ((graph_ == iter.graph_) &&
              (node_id_ == iter.node_id_) &&
              (adj_id_ == iter.adj_id_));
    } 

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node_id_;
    size_type adj_id_;

    IncidentIterator(const graph_type* graph, size_type node_id, size_type adj_id)
    : graph_(const_cast<graph_type*>(graph)), node_id_(node_id), adj_id_(adj_id) {    
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

    /** Construct an invalid EdgeIterator */
    EdgeIterator() {
      graph_ = nullptr;
      index_ = 0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      //Grab node1 and node2 id
      size_type node1_id_ = graph_ -> edges_[index_].first;
      size_type node2_id_ = graph_ -> edges_[index_].second;
      return Edge(graph_, index_, node1_id_, node2_id_);
    }

    EdgeIterator& operator++(){
      index_++; 
      return *this;
    }

    bool operator==(const EdgeIterator& iter) const{
      return ((graph_ == iter.graph_) && (index_ == iter.index_));
    } 


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    EdgeIterator(const graph_type* graph, size_type index)
    : graph_(const_cast<graph_type*>(graph)), index_(index) {    
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }

  edge_iterator edge_end() const{
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  struct internal_node {
    Point position_;
    node_value_type value_;
    std::vector<adjacent_node> adj_nodes_; 
  };

  struct adjacent_node{
    size_type adj_node_;
    size_type edge_index_;
  };

  // Vector of internal nodes
  std::vector<internal_node> nodes_;
  // Vector of node pairs, each pair corresponding to an edge
  std::vector<std::pair<size_type, size_type>> edges_;
  // Map of edge id's to edge values
  std::map<size_type, edge_value_type> edge_values_;
};

#endif // CME212_GRAPH_HPP
