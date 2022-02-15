#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <map>
#include <tuple>

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

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_value_type = E;

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
  using edge_map_type = std::map<std::tuple<size_type, size_type>, edge_type>;
  using edge_value_map_type = std::map<std::tuple<size_type, size_type>, edge_value_type>;
  // typedef typename edge_map_type::iterator edge_map_iterator;
  // typedef EdgeMap::iterator EdgeMapIterator; 

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
   : nodes_(){ //HW0: YOUR CODE HERE
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
  class Node : private totally_ordered<Node>{
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
    Node()
        : graph_(nullptr), index_(0) { 
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->positions_[index_];
    }

    /** Return this node's position reference. */
    Point& position() {
      return graph_->positions_[index_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // return size_type(-1);
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      return graph_->values_[index_];
    }
    const node_value_type& value() const {
      return graph_->values_.at[index_];
    }

    size_type degree() const {
      return graph_->adjacency_list_[index_].size();
    }
    
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, index_, 0);
    }
    
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, index_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (graph_ == n.graph_) && (index_ == n.index_);
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
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Point* position_;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type index_;
    Node(const Graph* graph, size_type index)
        : graph_(const_cast<Graph*>(graph)), index_(index) {
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    node_type n = Node(this, nodes_.size());
    positions_.push_back(position);
    values_.push_back(val);
    nodes_.push_back(n);
    return n;   
  }

  //remove node from graph
  size_type remove_node(const Node& n){
    if (has_node(n)) {
      // first remove all edges to node
      for (auto it = n.edge_begin(); it != n.edge_end(); it = n.edge_begin()){
        remove_edge(*it);
      }
      if (n.index() == nodes_.size() - 1) {
        nodes_.pop_back();
        positions_.pop_back();
        values_.pop_back();
      } else {
        // update position to be last position
        positions_[n.index()] = positions_.back();
        // update value to be last value
        values_[n.index()] = values_.back();

        //update adjacency list so final node is replaced with this node
        for (auto it = nodes_.back().edge_begin(); it != nodes_.back().edge_end(); it = nodes_.back().edge_begin()){
          add_edge(n, (*it).node2());
          remove_edge(*it);
        }
        //remove last node
        nodes_.pop_back();
        //remove last position and value
        positions_.pop_back();
        values_.pop_back();
      }
      return 1;
    }
    return 0;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (size() > n.index()) && (n.graph_ == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (size() > i) {
      return nodes_.at(i);
    } else {
      std::cout << "Node index out of range" << std::endl;
      return Node();
    }
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
    Edge()
    : graph_(nullptr), node1_index_(0), node2_index_(0) {
    }

    /* Return edge length */
    double length() const {
      return norm_2(graph_->positions_[node1_index_] - graph_->positions_[node2_index_]);
    }

    edge_value_type& value() {
      if (node1_index_ > node2_index_) {
        return graph_->edge_values_[{node1_index_, node2_index_}];
      } else {
        return graph_->edge_values_[{node2_index_, node1_index_}];
      }
    }
    
    const edge_value_type& value() const {
      if (node1_index_ > node2_index_) {
        return graph_->edge_values_.at({node1_index_, node2_index_});
      } else {
        return graph_->edge_values_.at({node2_index_, node1_index_});
      }
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->nodes_.at(node1_index_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->nodes_.at(node2_index_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1() && node2() == e.node2()) ||
             (node1() == e.node2() && node2() == e.node1());
      // (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      // return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //this is the same logic that order tuples
      if (node1() < e.node1()) {
        return true;
      } else if (node1() == e.node1()) {
        return node2() < e.node2();
      } else {
        return false;
      }
      // (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      // return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_index_;
    size_type node2_index_;
    Edge(const Graph* graph, size_type node1_index, size_type node2_index)
        : graph_(const_cast<Graph*>(graph)), node1_index_(node1_index),
          node2_index_(node2_index) {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
    // HW0: YOUR CODE HERE
    // return 0;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i < edges_.size()) {
      auto it = edges_.begin();
      std::advance(it, i);
      return it->second;
    } else {
      std::cout << "Edge index out of range" << std::endl;
      return Edge();
    }
  }

  //remove edge from graph, needs to update adjacency_list_, edges_ and edge_values_
  size_type remove_edge(const Node& n1, const Node& n2){
    if (has_edge(n1, n2)){
      size_type index1 = n1.index();
      size_type index2 = n2.index();
      if (index1 < index2) {
        std::swap(index1, index2);
      }
      edges_.erase({index1, index2});
      edge_values_.erase({index1, index2});
      for (auto it = adjacency_list_[index1].begin(); it != adjacency_list_[index1].end(); ++it) {
        if ((*it) == index2) {
          adjacency_list_[index1].erase(it);
          break;
        }
      }
      for (auto it = adjacency_list_[index2].begin(); it != adjacency_list_[index2].end(); ++it) {
        if ((*it) == index1) {
          adjacency_list_[index2].erase(it);
          break;
        }
      }
      return 1;
    } else {
      return 0;
    }
  }
  
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    if ((a.graph_ != this) || (b.graph_ != this)) {
      std::cout << "Nodes are not in this graph" << std::endl;
      return false;
    }
    //Check if edge in edges_ map
    if (edges_.count({a.index(), b.index()}) || edges_.count({b.index(), a.index()})) {
      return true;
    }
    return false;
  }

    bool has_edge_directed(const Node& a, const Node& b) const {
    if ((a.graph_ != this) || (b.graph_ != this)) {
      std::cout << "Nodes are not in this graph" << std::endl;
      return false;
    }
    //Check if edge in edges_ map
    if (edges_.count({a.index(), b.index()})){
      return true;
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    if ((a.graph_ != this) || (b.graph_ != this)) {
      std::cout << "Node a or b is not in this graph" << std::endl;
      return Edge();
    } else if (a == b){
      std::cout << "Node a and b are the same" << std::endl;
      return Edge();
    }
    
    edge_type e; 
    if (a.index() > b.index()) {
      if (has_edge_directed(a, b)) {
        return edges_.at({a.index(),b.index()});
      }
      e = Edge(this, a.index(), b.index());
      edges_.insert({{a.index(),b.index()}, e});
      edge_values_.insert({{a.index(),b.index()}, val});
    } else {
      if (has_edge_directed(b, a)) {
        return edges_.at({b.index(),a.index()});
      }
      e = Edge(this, b.index(), a.index());
      edges_.insert({{b.index(),a.index()}, e});
      edge_values_.insert({{b.index(),a.index()}, val});
    }
    adjacency_list_[a.index()].push_back(b.index());
    adjacency_list_[b.index()].push_back(a.index());
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_.clear();
    edges_.clear();
    positions_.clear();
    // HW0: YOUR CODE HERE
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    Node operator*() const {
      return *it_;
    }

    NodeIterator& operator++() {
      ++it_;
      return *this;
    }

    bool operator==(const NodeIterator& n) const {
      return it_ == n.it_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    typename std::vector<node_type>::iterator it_; //iterator for nodes_
    NodeIterator(const Graph* graph, size_type begin)
        : graph_(const_cast<Graph*>(graph)) { if (begin){it_ = graph_->nodes_.begin();}
        else {it_ = graph_->nodes_.end();}
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 1);
  }
  node_iterator node_end() const {
    return NodeIterator(this, 0);
  }

  node_iterator remove_node(node_iterator n_it){
    if (n_it == node_end()) {
      return n_it;
    }
    //Remove node
    remove_node(*n_it);

    return node_begin();
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : graph_(nullptr), node_index_(0), edge_index_(0) { 
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return Edge(graph_, node_index_, graph_->adjacency_list_[node_index_][edge_index_]);
    }

    IncidentIterator& operator++() {
      if (edge_index_ < graph_->adjacency_list_[node_index_].size()) {
        ++edge_index_;
      }
      return *this;
    }

    bool operator==(const IncidentIterator& it) const {
      return (node_index_ == it.node_index_) && (edge_index_ == it.edge_index_);
    }

   private:
    friend class Graph;
    Graph* graph_;
    size_type node_index_;
    size_type edge_index_; //index of current edge
    IncidentIterator(const Graph* graph, size_type node, size_type index)
        : graph_(const_cast<Graph*>(graph)), node_index_(node), 
        edge_index_(index) {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return (*it_).second;
    }
    EdgeIterator& operator++() {
      ++it_;
      return *this;
    }
    bool operator==(const EdgeIterator& e) const {
      return it_ == e.it_;
    }

   private:
    friend class Graph;
    Graph* graph_;
    typename edge_map_type::iterator it_; //iterator for edges_
    EdgeIterator(const Graph* graph, size_type begin)
        : graph_(const_cast<Graph*>(graph)) { if (begin){it_ = graph_->edges_.begin();}
        else {it_ = graph_->edges_.end();}
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 1);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, 0);
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it == edge_end()) {
      return e_it;
    }
    remove_edge(*e_it);
    return edge_begin();
  }

 private:
  std::vector<node_type> nodes_;
  std::vector<node_value_type> values_;
  std::vector<Point> positions_;
  edge_map_type edges_;
  edge_value_map_type edge_values_;
  std::unordered_map<size_type, std::vector<size_type> > adjacency_list_;
};

#endif // CME212_GRAPH_HPP
