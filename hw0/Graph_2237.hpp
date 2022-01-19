#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// Debugging defines
// #define INFO_DEF
// #define WARNING_DEF
// #define ERROR_DEF

#if defined(INFO_DEF) || defined(WARNING_DEF) || defined(ERROR_DEF)
#define DEBUG_HEAD(level) std::cout << level << ": " << __FILE__ << ":" << __func__ << ":" << __LINE__ << ": "
#endif

// Info Macros
#ifdef INFO_DEF
#include <iostream>
#define INFO_EVAL(item) \
  do { \
    DEBUG_HEAD("INFO"); \
    std::cout << #item << " -> " << (item) << std::endl; \
  } while(0)
#define INFO(msg) \
  do { \
    DEBUG_HEAD("INFO"); \
    std::cout << msg << std::endl; \
  } while(0)
#else
#define INFO(msg)
#define INFO_EVAL(item)
#endif

// Warning Macros
#ifdef WARNING_DEF
#define WARNING_EVAL(item) \
  do { \
    DEBUG_HEAD("WARNING"); \
    std::cout << #item << " -> " << (item) << std::endl; \
  } while(0)
#define WARNING(msg) \
  do { \
    DEBUG_HEAD("WARNING"); \
    std::cout << msg << std::endl; \
  } while (0)
#else
#define WARNING_EVAL(item)
#define WARNING(msg)
#endif

// Error Macros (But does not terminate)
#ifdef ERROR_DEF
#define ERROR_EVAL(item) \
  do { \
    DEBUG_HEAD("ERROR"); \
    std::cout << #item << " -> " << (item) << std::endl; \
  } while(0)
#define ERROR(msg) \
  do { \
    DEBUG_HEAD("ERROR"); \
    std::cout << msg << std::endl; \
  } while (0)
#else
#define ERROR_EVAL(item)
#define ERROR(msg)
#endif

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
  private:
  /** Predeclare edge_internal struct */
  Graph* g_ptr = this; // Pointer to self object

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
  class Node {
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
      // Invalid node default constructor
      graph_ptr_ = nullptr;
      point_ptr_ = nullptr;
      idx_ = -1;
    }
    /** Return this node's position. */
    const Point& position() const {
      if (graph_ptr_ == nullptr) {throw std::runtime_error("Invalid Node");}
      return *point_ptr_;
      // return Point();
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (graph_ptr_ == nullptr) {throw std::runtime_error("Invalid Node");}
      return idx_;
      // return size_type(-1);
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
      // Graph pointer check
      bool gph_eq = this->graph_ptr_ == n.graph_ptr_;
      // Index check
      bool idx_eq = this->idx_ == n.idx_;
      // Point object check
      bool pnt_eq = this->point_ptr_ == n.point_ptr_;
      if (gph_eq && idx_eq && pnt_eq) {return true;}
      else {return false;}
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
      // Throw exception if not in same graph
      if (this->graph_ptr_  != n.graph_ptr_) {throw std::invalid_argument("Nodes not from same graph");}
      // Compare by index
      INFO_EVAL(this->idx_);
      INFO_EVAL(n.idx_);
      return this->idx_ < n.idx_;
      // (void) n;           // Quiet compiler warning
      // return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Private members
    Graph* graph_ptr_;   // Pointer to parent graph
    Point* point_ptr_;   // Pointer to corresponding point object
    int idx_;            // Index of node
    // Construction of valid node
    Node(Graph* g_ptr, Point position, unsigned int index) :
        graph_ptr_(g_ptr), point_ptr_(&position), idx_(index) {};
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    size_type num_nodes = node_vec_.size();
    return num_nodes;
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
  Node add_node(const Point& position) {
    // Construct an invalid node then update the attributes
    node_type* NewNode = new Node(this, position, this->size());
    // Update Node map
    INFO("Added new node index " << this->size() << " at " << NewNode);
    node_vec_.push_back(NewNode);
    return *NewNode;
    // (void) position;      // Quiet compiler warning
    // return Node();        // Invalid node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    bool is_inGraph = this == n.graph_ptr_;
    // Grab index, check if valid, and see if node pointer is same
    unsigned int index = n.idx_;
    if (index < this->size()) {
      Node* tmp_node = node_vec_[index];
      bool is_samePt = tmp_node->point_ptr_ == n.point_ptr_;
      bool is_sameIdx = tmp_node->idx_ == n.idx_;
      return (is_inGraph && is_samePt && is_sameIdx);
    } else {return false;} // Out of bounds error

    // (void) n;            // Quiet compiler warning
    // return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // Check out of bounds
    if (i >= this->size()) {return Node();}
    // Grab node object at position i
    return *node_vec_[i];
    // (void) i;             // Quiet compiler warning
    // return Node();        // Invalid node
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph_ptr_ = nullptr;
      node_1_ = nullptr;
      node_2_ = nullptr;
      idx_ = -1;
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (graph_ptr_ == nullptr) {return Node();} // Invalid Node
      return * node_1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (graph_ptr_ == nullptr) {return Node();} // Invalid Node
      return * node_2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // Check for same graph
      if (graph_ptr_ != e.graph_ptr_) {return false;}
      // Check for node equality permutations
      bool perm_1 = (this->node1() == e.node1()) && (this->node2() == e.node2());
      bool perm_2;
      if (graph_ptr_->is_directed_) {perm_2 = false;}
      else {perm_2 = (this->node2() == e.node1()) && (this->node1() == e.node2());}
      return perm_1 || perm_2;
      // (void) e;           // Quiet compiler warning
      // return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Throw exception if not in same graph
      if (this->graph_ptr_  != e.graph_ptr_) {throw std::invalid_argument("Nodes not from same graph");}
      // Compare by index
      INFO_EVAL(this->idx_);
      INFO_EVAL(e.idx_);
      return this->idx_ < e.idx_;
      // (void) e;           // Quiet compiler warning
      // return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Private members
    Graph* graph_ptr_;                  // Pointer to parent graph
    const Node* node_1_, * node_2_;     // Pointer to corresponding node objects
    int idx_;                           // Index of edge
    // Construction of valid edge
    Edge(Graph* g_ptr, const Node* n1,const  Node* n2, unsigned int index) :
        graph_ptr_(g_ptr), node_1_(n1), node_2_(n2), idx_(index) {};
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_vec_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    if (i >= edge_vec_.size()) {
      WARNING("Index out of bounds. Returning invalid edge.");
      return Edge();
    }
    // Loop through edges, return edge with index i
    for (const auto &e : edge_vec_) {
      if ((unsigned int) e->idx_ == i) {return *e;}
    }
    WARNING("EXCEPTION: No edge found.");
    return Edge();
    // (void) i;             // Quiet compiler warning
    // return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if Nodes are in graph
    if (!has_node(a) || !has_node(b)) {return false;}
    const Node* a_ptr = &a;
    const Node* b_ptr = &b;
    auto a_search = edge_map_.find(a_ptr);
    if (a_search != edge_map_.end()) {
      auto b_search = a_search->second.find(b_ptr);
      if (b_search != a_search->second.end()) {return true;}
    }
    // Search other permutation if graph is indirected
    if (!is_directed_) {
      auto b_search = edge_map_.find(b_ptr);
      if (b_search != edge_map_.end()) {
        auto a_search = b_search->second.find(a_ptr);
        if (a_search != b_search->second.end()) {return true;}
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
    // Check if nodes in graph
    if (!this->has_node(a) || !this->has_node(b)) {
      WARNING("Node(s) do not exist in graph. No edge will be added.");
      return Edge(); // Invalid Edge
    }
    // Check if edge already in graph
    if (has_edge(a, b)) {
      WARNING("Edge already exists. No edge will be added.");
      return Edge();
    }
    // Create new edge
    Edge* e = new Edge(this, &a, &b, edge_vec_.size());
    INFO("Added new edge index " << this->num_edges() << " at " << e);
    // Save to graph members
    edge_vec_.push_back(e);
    edge_map_[&a][&b] = e;
    return *e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // Loop through all nodes and edges, deleting the objects.
    for (const auto &n : node_vec_) {
      delete n;
    }
    for (const auto &e : edge_vec_) {
      delete e;
    }
    // Re-assign
    node_vec_ = std::vector<Node*>();
    edge_vec_ = std::vector<Edge*>();
    edge_map_ = std::unordered_map<const Node*, std::unordered_map<const Node*, Edge*>>();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator {
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
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
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

 private:
  /** Variables:
   *      Node vector: vector of node pointers
   *      Edge vector: vector of edge pointers
   *      Edge map   : Adjacency list (using maps to indicate what edge corresponds to the node pair) for searching
   */
  std::vector<Node*> node_vec_;
  std::vector<Edge*> edge_vec_;
  std::unordered_map<const Node*, std::unordered_map<const Node*, Edge*>> edge_map_;
  bool is_directed_ = false; // False: Graph is not directed.
};

#endif // CME212_GRAPH_HPP
