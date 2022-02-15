#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

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
template <typename V = int, typename E = int>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  // struct internal_element;
  
  struct node_struct;
  struct edge_struct;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Templated type */
  using node_value_type = V;
  using edge_value_type = E;

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

    //nodes - default no nodes
    //edges - default no edges
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
      // uid_ = -1; // should be out-of-bounds 
    }

    /** Return this node's position for editing */
    Point& position () { return fetch().point; }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return fetch().idx;
    }

    /** Return the value of a node in the graph.
     *
     * Complexity: O(1).
     */
    node_value_type& value() { return fetch().value; }
    
    /** Return the value of a node in the graph as const.
     *
     * Complexity: O(1).
     */
    const node_value_type& value() const { return fetch().value; }
    
    /** Return degree of node */
    size_type degree() const {
      return graph_->neighbor_map_[uid_].size();
    }

    /** Initialize an incident iterator at first neighboring edge */
    incident_iterator edge_begin() const { return IncidentIterator(graph_, uid_, 0); }
    
    /** Return an incident iterator at one past final neighboring edge */
    incident_iterator edge_end() const { return IncidentIterator(graph_, uid_, degree()); }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ && uid_ == n.uid_) {
        return true;
      }
      return false;
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
      // Compare pointers to node parent graphs
      if (graph_ < n.graph_) { 
        return true;
      }
      // If referencing the same graph, compare uids
      if (graph_ == n.graph_ && uid_ < n.uid_) {
        return true;
      }
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    graph_type* graph_;
    // This node's unique identification number
    size_type uid_;


    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

    /** Helper methods to return the appropriate node in the graph.
     * This method accesses the element from the appropriate vector
     * when the node is or is not const
     */
    node_struct& fetch() const {
      assert(uid_ >= 0 && uid_ < graph_->nodes_.size()+1); // Allow to access one past
      return graph_->nodes_[uid_];
    }
    node_struct& fetch() {
      assert(uid_ >= 0 && uid_ < graph_->nodes_.size()+1); // Allow to access one past
      return graph_->nodes_[uid_];
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // Return the node vector length
    return node_ids_.size();
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
    size_type uid = nodes_.size();
    size_type idx = node_ids_.size();
    nodes_.emplace_back(position, idx, value);
    node_ids_.push_back(uid);

    // initialize empty neighbor set
    neighbor_map_[uid] = std::vector<size_type> {};

    return Node(this, uid); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1). 
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (n.graph_ == this){
      if (n.index() >= 0 && n.index() < node_ids_.size())
        return true;

    }
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, node_ids_[i]);
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

    /** Return the value of an edge in the graph.
     *
     * Complexity: O(1).
     */
    edge_value_type& value() { 
      size_type lower = std::min(node1_, node2_);
      size_type upper = std::max(node1_, node2_);
      return graph_->edges_[std::make_tuple(lower,upper)]; 
    }
    
    /** Return the value of an edge in the graph as const.
     *
     * Complexity: O(1).
     */
    const edge_value_type& value() const { 
      size_type lower = std::min(node1_, node2_);
      size_type upper = std::max(node1_, node2_);
      return graph_->edges_[std::make_tuple(lower,upper)]; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ == e.graph_ && 
        ((node1_ == e.node1_ && node2_ == e.node2_) ||
        (node1_ == e.node2_ && node2_ == e.node1_))) {
        return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ < e.graph_) { 
        return true;
      }
      // If referencing the same graph, compare node numberings
      if (graph_ == e.graph_) {
        if (std::min(node1_, node2_) < std::min(e.node1_, e.node2_)) {
          return true;
        } 
        if (std::min(node1_, node2_) == std::min(e.node1_, e.node2_) &&
          std::max(node1_, node2_) < std::max(e.node1_, e.node2_)) {
          return true;
        } 
      }
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    graph_type* graph_;
    // The nodes in this edge
    size_type node1_;
    size_type node2_;

    /** Private Constructor */
    Edge(const Graph* graph, size_type node1, size_type node2)
        : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_(node2) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_ids_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    std::tuple<size_type, size_type> nodes = edge_ids_[i];
    return Edge(this, std::get<0>(nodes), std::get<1>(nodes));   
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    if (edges_.count(std::make_tuple(a.uid_, b.uid_)) > 0 || 
      edges_.count(std::make_tuple(b.uid_, a.uid_)) > 0) {
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
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    size_type lower = std::min(a.uid_, b.uid_);
    size_type upper = std::max(a.uid_, b.uid_);

    std::tuple<size_type, size_type> check_edge = std::make_tuple(lower,upper);
    if (edges_.count(check_edge) < 1) {
      
      // add edge tuple to list
      edge_ids_.push_back(check_edge);

      // make edge tuple to map to default edge_value_type true
      edges_[check_edge] = edge_value_type {};

      // update neighbor list
      add_neighbors(lower, upper);
    }
    return Edge(this, a.uid_, b.uid_); 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_ids_.clear();
    nodes_.clear();
    edges_.clear();
    edge_ids_.clear();
    neighbor_map_.clear();
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
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    /** Return node object corresponding to node iterator*/
    Node operator*() const {
      return Node(graph_, graph_->node_ids_[node_idx_]);
    }
    
    /** Increment node iterator */
    NodeIterator& operator++() {
      node_idx_++;
      return *this;
    }
    
    /** Check equality between node iterators */
    bool operator==(const NodeIterator& it) const {
      return (graph_ == it.graph_ && node_idx_ == it.node_idx_);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // store parent graph and node index
    graph_type* graph_;
    size_type node_idx_;

    // Hidden constructor
    NodeIterator(const graph_type* graph, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), node_idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Return iterator pointing to first node in node list */
  node_iterator node_begin() const { return NodeIterator(this, 0); }
  
  /** Return iterator pointing to one past final node in node list */
  node_iterator node_end() const { return NodeIterator(this, node_ids_.size()); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return edge iterator points to */
    Edge operator*() const { 
      return Edge(graph_, node_id_, 
        graph_->neighbor_map_[node_id_][neighbor_idx_]); 
    }
    
    /** Increment incident iterator */
    IncidentIterator& operator++() {
      neighbor_idx_++;
      return *this;
    }
    
    /** Check equality between incident iterators. 
     * Makes use of equality between edges. */
    bool operator==(const IncidentIterator& it) const {
      return (graph_== it.graph_ && 
        node_id_ == it.node_id_ &&
        neighbor_idx_ == it.neighbor_idx_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    
    // Store parent graph, node id, and index in neighbor list
    graph_type* graph_;
    size_type node_id_;
    size_type neighbor_idx_;

    // Hidden constructor
    IncidentIterator(const graph_type* graph, size_type id, size_type neighbor)
      : graph_(const_cast<graph_type*>(graph)), node_id_(id), neighbor_idx_(neighbor) {}
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    
    /** Return edge to which edge iterator points */
    Edge operator*() const {
      return Edge(graph_, std::get<0>(graph_->edge_ids_[edge_idx_]), 
        std::get<1>(graph_->edge_ids_[edge_idx_]));
    }
    
    /** Increment edge iterator */
    EdgeIterator& operator++() {
      edge_idx_++;
      return *this;
    }
    
    /** Check equality between two edge iterators*/
    bool operator==(const EdgeIterator& it) const {
      return (graph_ == it.graph_ && edge_idx_==it.edge_idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    
    // Store parent graph and edge index
    graph_type* graph_;
    size_type edge_idx_;

    // Hidden initializer
    EdgeIterator(const graph_type* graph, size_type idx)
      : graph_(const_cast<graph_type*>(graph)), edge_idx_(idx) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  /** Iterator pointing to beginning of graph edges */
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  /** Iterator pointing to one past graph edges */
  edge_iterator edge_end() const { return EdgeIterator(this, edge_ids_.size()); }

  // Node and Edge removal functions


  /** Remove an edge from the graph and return whether successful
   * @param[in] a node1 of the edge
   * @param[in] b node2 of the edge
   * @return a size_type which is 1 if the edge removal 
   *         was successful, else 0
   * @post new num_edges() == old num_edges() - 1 if successful
   * @post new edges_ has removed the edge as a key if successful
   * @post new neighbor_map_ has removed the edge from neighbor sets of 
   *       both lower and upper node indices if successful
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(num_edges) amortized operations (in call to remove_edge).
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type lower = std::min(a.index(), b.index());
    size_type upper = std::max(a.index(), b.index());
    std::tuple<size_type, size_type> idx = std::make_tuple(lower,upper);
    if (edges_.count(idx)) {
      remove_edge(lower, upper);
      return 1;
    } else{
      return 0;
    }
  }

  /** Remove an edge from the graph and return whether successful
   * @param[in] e reference to the Edge object
   * @return a size_type which is 1 if the edge removal 
   *         was successful, else 0
   * @post new num_edges() == old num_edges() - 1 if successful
   * @post new edges_ has removed the edge as a key if successful
   * @post new neighbor_map_ has removed the edge from neighbor sets of 
   *       both lower and upper node indices if successful
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(num_edges) amortized operations (in call to remove_edge).
   */
  size_type remove_edge(const Edge& e) {
    size_type lower = std::min(e.node1_, e.node2_);
    size_type upper = std::max(e.node1_, e.node2_);
    std::tuple<size_type, size_type> idx = std::make_tuple(lower,upper);
    if (edges_.count(idx)) {
      remove_edge(lower, upper);
      return 1;
    } else{
      return 0;
    }
  }

  /** Remove an edge from the graph and return an edge iterator
   * @param[in] e_it iterator pointing to edge
   * @return an iterator pointing to the start of the edge list if
   *         successful, else ++e_it
   * @post new num_edges() == old num_edges() - 1 if successful
   * @post new edges_ has removed the edge as a key if successful
   * @post new neighbor_map_ has removed the edge from neighbor sets of 
   *       both lower and upper node indices if successful
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(num_edges) amortized operations (in call to remove_edge).
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge e = *e_it;
    size_type lower = std::min(e.node1_, e.node2_);
    size_type upper = std::max(e.node1_, e.node2_);
    std::tuple<size_type, size_type> idx = std::make_tuple(lower,upper);
    if (edges_.count(idx)) {
      remove_edge(lower, upper);
      return this->edge_begin();
    } else {
      return ++e_it;
    }
  }

  /** Remove a node from the graph and return whether successful
   * @param[in] a node to be removed
   * @return a size_type which is 1 if the node removal 
   *         was successful (if the node exists), else 0
   * @post new num_nodes() == old num_nodes() - 1 if successful
   * @post all edges incident to @a a have been removed if successful
   *
   * Complexity: O(1 + degree * num_edges) amortized operations 
   *             (note we make degree calls to remove_edge).
   */
  size_type remove_node(const Node& a) {
    // if not a valid node, return 0
    if (!has_node(a))
      return 0;
    // remove all incident edges
    for (incident_iterator it = a.edge_begin(); it!= a.edge_end();) {
      Edge e = *it;
      size_type success = remove_edge(e);
      if (success) {
        it = a.edge_begin();
      } else {
        ++it;
      }
    }

    // remove node from node_ids_
    size_type start_size = node_ids_.size();
    size_type remove_idx = a.index();

    node_ids_[remove_idx] = node_ids_[start_size-1];
    nodes_[node_ids_[remove_idx]].idx = remove_idx;
    nodes_[a.uid_].idx = -1;

    node_ids_.pop_back();
    assert(node_ids_.size() == start_size - 1);
    return 1;
  }

  /** Remove a node from the graph and return an iterator
   * @param[in] n_it iterator to node to be removed
   * @return a node_iterator pointing to node_begin() if removal
   *         was successful (if the node exists), else ++n_it
   * @post new num_nodes() == old num_nodes() - 1 if successful
   * @post all edges incident to @a a have been removed if successful
   *
   * Complexity: O(num_nodes + degree * num_edges) amortized operations 
   *             (note we make degree calls to remove_edge).
   */
  node_iterator remove_node(node_iterator n_it) {
    Node a = *n_it;
    size_type success = remove_node(a);
    if (success) {
      return this->node_begin();
    } else {
      return ++n_it;
    }
  }



  // Helpers to print all node and edge data
  void print_node_data() {
    for (node_iterator it = this->node_begin(); it != this->node_end(); ++it)
      std::cout << "Node id " << (*it).index() << ": " << (*it).position() <<(*it).value() << std::endl;
  }
  void print_edge_data() {
    for (edge_iterator it = this->edge_begin(); it != this->edge_end(); ++it)
      std::cout << "Edge id (" << (*it).node1().index() <<"," <<(*it).node2().index() << "): " <<(*it).value() << std::endl;
  }

 private:


  // Removing an edge takes O(1) + O(edges) + O(deg_lower) + O(deg_upper) = O(edges)
  void remove_edge(size_type lower, size_type upper) {

    std::tuple<size_type, size_type> edgeidx = std::make_tuple(lower, upper);
    
    // erase from edge map
    edges_.erase(edgeidx);

    // erase from edge_ids vector
    for (auto it = edge_ids_.begin(); it != edge_ids_.end(); ++it) {
      if (*it == edgeidx) {
        edge_ids_.erase(it);
        break;
      }
    }

    // erase edge from neighbor map
    for (auto it = neighbor_map_[lower].begin(); it != neighbor_map_[lower].end(); ++it) {
      if (*it == upper) {
        neighbor_map_[lower].erase(it);
        break;
      }
    }
    for (auto it = neighbor_map_[upper].begin(); it != neighbor_map_[upper].end(); ++it) {
      if (*it == lower) {
        neighbor_map_[upper].erase(it);
        break;
      }
    }
  }

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // node structure
  struct node_struct {
    Point point;
    size_type idx;
    node_value_type value;
    node_struct(Point point_, size_type idx_, node_value_type value_) 
      : point(point_), idx(idx_), value(value_) { }
  };

  // nodes
  std::vector<size_type> node_ids_; // maps index of current node to node id
  std::vector<node_struct> nodes_; // stores data for all nodes ever. indexed by node id

  // edge map
  std::map<std::tuple<size_type, size_type>, edge_value_type> edges_;
  std::vector<std::tuple<size_type, size_type>> edge_ids_;
  std::map<size_type, std::vector<size_type>> neighbor_map_;
  
  /** Helper function to update neighbor map when adding a new edge
   * @pre neighbor_map[lower] and neighbor_map[upper] have already 
   * been initalized */
  void add_neighbors(size_type lower, size_type upper) {
    // update maps for both lower and upper nodes
    neighbor_map_[lower].push_back(upper);
    neighbor_map_[upper].push_back(lower);
  }

};

#endif // CME212_GRAPH_HPP
