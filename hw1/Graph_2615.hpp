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
#define ERROR_DEF

#if defined(INFO_DEF) || defined(WARNING_DEF) || defined(ERROR_DEF)
#include <iostream>
#define DEBUG_HEAD(level) std::cout << level << ": " << __FILE__ << ":" << __func__ << ":" << __LINE__ << ": "
#endif

// Info Macros
#ifdef INFO_DEF
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
template <typename V = int>
class Graph {
  private:
  Graph* g_ptr = this; // Pointer to self object
  /** Predeclaration for node internal */
  struct node_internal;

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

  /** Assigns node value type from template */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // TODO: Unused variable: node_int_vec_
    node_vec_ = new std::vector<Node>;
    node_int_vec_ = new std::vector<node_internal>;
    edge_vec_ = new std::vector<Edge>;
    edge_map_ = new std::unordered_map<size_type, std::unordered_map<size_type, size_type>>;
  }

  /** Default destructor */
  ~Graph() {
    for (auto &n : *node_vec_) {
      delete n.intern_ptr_;
    }
    delete edge_map_;
    delete edge_vec_;
    delete node_int_vec_;
    delete node_vec_;
  };

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
      // Invalid node default constructor
      intern_ptr_ = nullptr;
    }
    /** Return this node's position. */
    const Point& position() const {
      if (intern_ptr_ == nullptr) {throw std::runtime_error("Invalid Node");}
      return intern_ptr_->point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (intern_ptr_ == nullptr) {throw std::runtime_error("Invalid Node");}
      return intern_ptr_->idx;
    }

    /** Returns the stored value as a mutable node_value_type */
    node_value_type& value() {return intern_ptr_->value;}

    /** Returns the stored value as a const node_value_type */
    const node_value_type& value() const {return intern_ptr_->value;};

    /** Returns the number of edges leaving this node */
    size_type degree() const {
      if (intern_ptr_->graph_ptr == nullptr) {throw std::runtime_error("Called degree() on an invalid node");}
      auto& e_map = *intern_ptr_->graph_ptr->edge_map_;
      auto n_search = e_map.find(index());
      if (n_search != e_map.end()) {
        INFO_EVAL(e_map[index()].size());
        return e_map[index()].size();
      }
      return 0;
    }

    /** Returns first iterator in the current Node entry of edge_map */
    incident_iterator edge_begin() const {
      auto& e_map = *intern_ptr_->graph_ptr->edge_map_;
      INFO("\t\t\tNum of adjacent nodes: " << e_map[index()].size());
      return IncidentIterator(this, e_map[index()].begin());
    }
    /** Returns iterator corresponding to end of incident edges */
    incident_iterator edge_end() const {
      auto& e_map = *intern_ptr_->graph_ptr->edge_map_;
      return IncidentIterator(this, e_map[index()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (intern_ptr_ == nullptr || n.intern_ptr_ == nullptr) {return false;} // One of the nodes is an invalid node
      // Graph pointer check
      bool gph_eq = (intern_ptr_->graph_ptr == n.intern_ptr_->graph_ptr);
      // Index check
      bool idx_eq = index() == n.index();

      return gph_eq && idx_eq;
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
      if (this->intern_ptr_->graph_ptr  != n.intern_ptr_->graph_ptr) {throw std::invalid_argument("Nodes not from same graph");}
      // Compare by index
      return index() < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Private members
    node_internal* intern_ptr_;   // Pointer to internal object
    // Construction of valid node
    Node(node_internal* ni) : intern_ptr_(ni) {};
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return node_vec_->size();
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
  Node& add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // TODO: Debug creation of node
    size_type i = size();
    node_internal* ni = new node_internal(val, g_ptr, position, i);
    Node n = Node(ni);
    node_vec_->push_back(n);
    // INFO("Created valid node " << i << " with internal object at " << n.intern_ptr_);
    // INFO("\t\tand graph ptr " << n.intern_ptr_->graph_ptr);
    return (*node_vec_)[i];
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Check if in same graph
    // INFO("Checking if graph has node ID at " << n.index());
    // INFO_EVAL(n.intern_ptr_->graph_ptr);
    // INFO_EVAL(g_ptr);
    if (this != n.intern_ptr_->graph_ptr) {return false;}
    // Grab index, check if valid, and see if node pointer is same
    unsigned int index = n.index();
    if (index < size()) {
      // INFO("Check if node is same as node stored");
      return ((*node_vec_)[index] == n);
    } else {
      ERROR("Index out of bounds");
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node& node(size_type i) const {
    // Check out of bounds
    if (i >= size()) {
      ERROR("Graph size: " << size() << "\tAccessing node out of bounds at index " << i );
      throw std::invalid_argument("Index out of bounds");
      // return Node();
      }
    // Grab node object at position i
    return (*node_vec_)[i];
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
      graph_ptr_ = nullptr;
      n1_idx_ = -1;
      n2_idx_ = -1;
      idx_ = -1;
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (graph_ptr_ == nullptr) {return Node();} // Invalid Node
      auto& n_vec = *graph_ptr_->node_vec_;
      return n_vec[n1_idx_];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (graph_ptr_ == nullptr) {return Node();} // Invalid Node
      auto& n_vec = *graph_ptr_->node_vec_;
      return n_vec[n2_idx_];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (graph_ptr_ == nullptr || e.graph_ptr_ == nullptr) {return false;} // Check for invalid edges
      // Check for same graph
      if (graph_ptr_ != e.graph_ptr_) {return false;}
      // Check for node equality permutations
      bool perm_1 = (node1() == e.node1()) && (node2() == e.node2());
      bool perm_2;
      if (graph_ptr_->is_directed_) {perm_2 = false;}
      else {perm_2 = (node2() == e.node1()) && (node1() == e.node2());}
      return perm_1 || perm_2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Throw exception if not in same graph
      if (graph_ptr_  != e.graph_ptr_) {throw std::invalid_argument("Nodes not from same graph");}
      // Compare by index
      INFO_EVAL(idx_);
      INFO_EVAL(e.idx_);
      return idx_ < e.idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Private members
    // TODO: Consider making Edge Internal obj
    Graph* graph_ptr_;                  // Pointer to parent graph
    unsigned int n1_idx_, n2_idx_;     // Pointer to corresponding node objects
    unsigned int idx_;                           // Index of edge
    // Construction of valid edge
    Edge(Graph* g_ptr, size_type n1, size_type n2, unsigned int index) :
        graph_ptr_(g_ptr), n1_idx_(n1), n2_idx_(n2), idx_(index) {};
  };

  /** Return the total number of edges in the graph.
   *
    * Complexity: O(1)
   */
  size_type num_edges() const {
    if (is_directed_) {return edge_vec_->size();}
    return (size_type) edge_vec_->size()/2; // No rounding should occur if size is undirected
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: O(1)
   */
  Edge edge(size_type i) const {
    if (i >= num_edges()) {
      WARNING("Index out of bounds. Returning invalid edge.");
      return Edge();
    }
    // Obtain correct index, and find in edge_vec
    if (!is_directed_) {i *= 2;}
    return (*edge_vec_)[i];
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
    // Search for node pair in map
    auto a_search = edge_map_->find(a.index());
    if (a_search != edge_map_->end()) {
      auto b_search = a_search->second.find(b.index());
      if (b_search != a_search->second.end()) {return true;}
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
    // INFO("Attempt to add edge between nodes " << a.index() << " and " << b.index());
    // Check if nodes in graph
    if (!has_node(a) || !has_node(b)) {
      WARNING("Node(s) do not exist in graph. No edge will be added.");
      return Edge(); // Invalid Edge
    }
    // Check if edge already in graph
    if (has_edge(a, b)) {
      WARNING("Edge already exists. No edge will be added.");
      size_type e_idx = (*edge_map_)[a.index()][b.index()];
      return (*edge_vec_)[e_idx];
    }
    // INFO("Create new edge between " << a.index() << " and " << b.index());
    // Grab node data from existing database
    Node n1 = (*node_vec_)[a.index()];
    Node n2 = (*node_vec_)[b.index()];
    // Check for node equality
    assert(n1 == a);
    assert(n2 == b);
    // Create new edge
    size_type idx =  num_edges();
    // Save to graph members
    (*edge_map_)[a.index()][b.index()] = edge_vec_->size();
    Edge e = Edge(this, a.index(), b.index(), idx);
    edge_vec_->push_back(e);
    // Save reverse map if not directed
    if (!is_directed_) {
      (*edge_map_)[b.index()][a.index()]  = edge_vec_->size();
      Edge e0 = Edge(this, b.index(), a.index(), idx);
      edge_vec_->push_back(e0);
    }
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // TODO: Unused variable node_int_vec_
    for (auto &n : *node_vec_) {
      delete n.intern_ptr_;
    }
    delete edge_map_;
    delete edge_vec_;
    delete node_int_vec_;
    delete node_vec_;
    // Re-assign
    node_vec_ = new std::vector<Node>;
    edge_vec_ = new std::vector<Edge>;
    node_int_vec_ = new std::vector<node_internal>;
    edge_map_ = new std::unordered_map<size_type, std::unordered_map<size_type, size_type>>;
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

    using vec_it_type = typename std::vector<Node>::iterator;   // Vector iterator

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** De-reference NodeIterator for the corresponding Node object */
    Node operator*() const {return *node_vec_it_;}

    /** Moves to next node */
    NodeIterator& operator++() {
      node_vec_it_++;
      return *this;
    }

    /** Checks for equality of NodeIterator */
    bool operator==(const NodeIterator& ni) const {return node_vec_it_ == ni.node_vec_it_;}

    bool operator!=(const NodeIterator& ni) const {return !(*this == ni);}

   private:
    friend class Graph;
    // Valid constructor for creating node iterator
    NodeIterator(const vec_it_type it) : node_vec_it_(it) {}
    // Store vector iterator
    vec_it_type node_vec_it_;
  };

  /** Creates a NodeIterator from the beginning of node_vec. */
  node_iterator node_begin() const {return NodeIterator(node_vec_->begin());}

  /** Creates a NodeIterator for one past the end of node_vec */
  node_iterator node_end() const {return NodeIterator(node_vec_->end());}

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

    using map_it_type = typename std::unordered_map<size_type, size_type>::iterator;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** De-references incident iterator into Edge object */
    Edge operator*() const {
      // Check if parent node is first node
      size_type e_idx = m_it_->second;
      auto& e_vec = *n_ptr_->intern_ptr_->graph_ptr->edge_vec_;
      Edge e0 = e_vec[e_idx];
      if (e0.n1_idx_ == n_ptr_->index()) {return e0;}
      INFO("Swapping nodes " << e0.n2_idx_ << " and " << e0.n1_idx_);
      return Edge(e0.graph_ptr_, e0.n2_idx_, e0.n1_idx_, e0.idx_);
    }
    /** Increment operator */
    IncidentIterator& operator++() {
      m_it_++;
      return *this;
    }
    /** Check equality operator */
    bool operator==(const IncidentIterator& ii) const {return m_it_ == ii.m_it_;}

    /** Check inequality operator */
    bool operator!=(const IncidentIterator& ii) const {return !(*this == ii);}

   private:
    friend class Graph;
    // Private members
    const Node* n_ptr_; // Pointer to parent node
    map_it_type m_it_; // Corresponding map iterator
    // Private constructor
    IncidentIterator(const Node* n, map_it_type m) : n_ptr_(n), m_it_(m) {}
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

    using vec_it_type = typename std::vector<Edge>::const_iterator;   // Vector iterator
    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereference to Edge object */
    Edge operator*() const {return *edge_vec_it_;}

    /** Increments to next edge */
    EdgeIterator& operator++() {
      edge_vec_it_++;
      return *this;
    }

    /** Check for equality of EdgeIterator */
    bool operator==(const EdgeIterator& ei) const {return edge_vec_it_ == ei.edge_vec_it_;}

    bool operator!=(const EdgeIterator& ei) const {return !(*this == ei);}

   private:
    friend class Graph;
    // Valid constructor
    EdgeIterator(const vec_it_type it) : edge_vec_it_(it) {}
    // Edge vector iterator
    vec_it_type edge_vec_it_;
  };

  /** Creates a EdgeIterator from beginning of edge_vec. */
  edge_iterator edge_begin() const {return EdgeIterator(edge_vec_->begin());}

  /** Creates a EdgeIterator from end of edge_vec */
  edge_iterator edge_end() const {return EdgeIterator(edge_vec_->end());}

 private:
  /** Variables:
   *      Node vector: vector of node objects
   *      Node internal vector: Vector of node internals
   *      Point vector: Vector of point objects
   *      Edge vector: vector of edge objects
   *      Edge map   : Adjacency list (using maps to indicate what edge corresponds to the node pair) for searching
   *                        The edge index returned corresponds to the Edge object stored in the Edge vector
   *  Struct:
   *      Node internal: Contains all of the attributes of Node
   */
  std::vector<Node>* node_vec_;
  std::vector<node_internal>* node_int_vec_;
  std::vector<Point>* point_vec_;
  std::vector<Edge>* edge_vec_;
  std::unordered_map<unsigned int, std::unordered_map<unsigned int, unsigned int>>* edge_map_;
  bool is_directed_ = false; // False: Graph is not directed.
  struct node_internal {
    node_value_type value;
    Graph* graph_ptr;
    Point point;
    size_type idx;
    node_internal(node_value_type val, Graph* g, Point p, size_type i) : value(val), graph_ptr(g), point(p), idx(i) {}
  };
};
#endif // CME212_GRAPH_HPP
