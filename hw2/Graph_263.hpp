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
#define INFO_DEF
#define WARNING_DEF
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

template <class T> class count_T;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
class Graph : private count_T<Graph<V, E>> {
  private:
  Graph* g_ptr = this; // Pointer to self object
  unsigned int gid_;    // Graph ID (for differentiating nodes and edges from different graphs)
  /** Predeclaration for internal objects*/
  struct node_internal;
  struct edge_internal;

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

  /** Assigns node and edge value type from template */
  using node_value_type = V;
  using edge_value_type = E;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    node_vec_ = new std::vector<Node>;
    node_int_vec_ = new std::vector<node_internal>;
    edge_vec_ = new std::vector<Edge>;
    edge_int_vec_ = new std::vector<edge_internal>;
    edge_map_ = new std::unordered_map<size_type, std::unordered_map<size_type, Edge>>;
    // Graph ID
    gid_ = count_T<Graph<V,E>>::get_count() - 1;
  }

  /** Default destructor */
  ~Graph() {
    delete edge_map_;
    delete edge_int_vec_;
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
      graph_ptr_ = nullptr;
    }
    ~Node() = default;
    /** Return this node's position. */
    const Point& position() const {
      if (!is_valid()) {throw std::runtime_error("Invalid Node, position.");}
      return internal_()->point;
    }

    /** Return mutable node position */
    Point& position() {
      if (!is_valid()) {throw std::runtime_error("Invalid Node, position (mutable).");}
      return internal_()->point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      if (!is_valid()) {throw std::runtime_error("Invalid Node, index.");}
      return internal_()->idx;
    }

    /** Returns the stored value as a mutable node_value_type */
    node_value_type& value() {
      assert(is_valid());
      return internal_()->value;
    }

    /** Returns the stored value as a const node_value_type */
    const node_value_type& value() const {
      assert(is_valid());
      return internal_()->value;
    };

    /** Returns the number of edges leaving this node */
    size_type degree() const {
      if (!is_valid()) {throw std::runtime_error("Called degree() on an invalid node");}
      auto& e_map = *graph_()->edge_map_;
      auto n_search = e_map.find(nint_idx_);
      if (n_search != e_map.end()) {
        // INFO_EVAL(e_map[index()].size());
        return e_map[nint_idx_].size();
      }
      return 0;
    }

    /** Returns FALSE if node is invalid, otherwise TRUE. */
    bool is_valid() const {
      // Check for valid pointer
      if (graph_ptr_ == nullptr) {return false;}
      // Check if node has been deleted
      if (internal_()->idx >= graph_ptr_->num_node_) {return false;}
      return true;
    }

    /** Returns first iterator in the current Node entry of edge_map */
    incident_iterator edge_begin() const {
      auto& e_map = *graph_()->edge_map_;
      // INFO("\t\t\tNum of adjacent nodes: " << e_map[index()].size());
      return IncidentIterator(this, e_map[nint_idx_].begin());
    }
    /** Returns iterator corresponding to end of incident edges */
    incident_iterator edge_end() const {
      auto& e_map = *graph_()->edge_map_;
      return IncidentIterator(this, e_map[nint_idx_].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph, same point, and the same index.
     */
    bool operator==(const Node& n) const {
      if (!is_valid() || !n.is_valid()) {return false;} // One of the nodes is an invalid node
      // Graph pointer check
      bool gph_eq = (graph_() == n.graph_());
      // Index check
      bool idx_eq = index() == n.index();
      // Points check
      bool pnt_eq = position() == n.position();

      return gph_eq && idx_eq && pnt_eq;
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
      // Compare by internal index and graph id
      return 1e8 * graph_ptr_->gid_ + nint_idx_ < 1e8 * n.graph_ptr_->gid_ + n.nint_idx_;
    }

    private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Returns memory address of associated graph
    Graph* graph_() const {
      return graph_ptr_;
    }
    // Returns pointer to internal object for ease-of-access
    node_internal* internal_() const {
      return &(*graph_ptr_->node_int_vec_)[nint_idx_];
    }
    // Private method to invalidate current node, assumes valid node
    void invalidate_() {
      internal_()->idx = graph_ptr_->num_node_;
      graph_ptr_ = nullptr;
    }
    // Private members
    Graph* graph_ptr_;
    size_type nint_idx_;   // Pointer to internal object
    // Construction of valid node
    Node(Graph* g, size_type i) : graph_ptr_(g), nint_idx_(i) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_node_;
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
    if (num_nodes() >= 1e8) {
      // Maximum nodes reached
      WARNING("Maximum node capacity reached, no node created.");
      return Node();
    }
    size_type i = num_node_++;
    try {
      node_vec_->push_back(Node(this, (*node_int_vec_).size()));
      node_int_vec_->push_back(node_internal(val, position, i));
    } catch(...) {
      std::cout << "Failed to create node. Unknown error.\n";
      throw;
    }
    // INFO("Created valid node " << i << " with internal object at " << n.nint_idx_);
    // INFO("\t\tand graph ptr " << n.graph_());
    return (*node_vec_)[i];
  }

  /** Removes node n from the graph
   * @pre @a n is a valid node
   *
   * @post Node @a n is invalidated (if it exists)
   * @post All edges incident to the node are removed
   * @post All incident iterators from @a n are invalidated
   * @return The number of nodes removed
   */
  size_type remove_node(Node& n) {
    if (!has_node(n)) {return 0;}
    // Removes all edges incident to the node
    for (auto it = n.edge_begin(); it != n.edge_end();) {
      auto n2 = (*it).node2();
      ++it; // Increment iterator before invalidating
      size_type i = remove_edge(n, n2);
      if (i == 0) {
        WARNING("\tRemoval unsuccessful");
      }
    }
    // Erase branch
    (*edge_map_).erase(n.nint_idx_);
    size_type idx = n.index();
    // Invalidate node
    n.invalidate_();
    // Swap with last valid node
    std::swap((*node_vec_)[idx], (*node_vec_)[--num_node_]);  // Also decrements node count
    // Make index match for swapped node
    size_type int_idx = (*node_vec_)[idx].nint_idx_;
    (*node_int_vec_)[int_idx].idx = idx;
    assert(!(*node_vec_)[num_node_].is_valid()); // Check node is invalid
    return 1;
  }

  /** Removes node by nodeiterator
   * @pre Node iterator is valid (not node_end())
   *
   * @post Node referred to by the iterator is invalidated
   * @post All edges incident to the node are removed
   * @post All incident iterators from the node are invalidated
   * @post The input node iterator is equivalent to the return iterator
   * @return The iterator to a valid node of the same index OR the end iterator
   */
  NodeIterator remove_node(NodeIterator n_it) {
    if (n_it == node_end()) {return node_end();} // No node removed
    // Grab node
    Node n = (*n_it);
    size_type idx = n.idx();
    assert(n.is_valid());
    // Remove node
    remove_node(n);
    // Return node iterator at same place
    return node_begin() + idx;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Check if in same graph
    if (this != n.graph_()) {return false;}
    assert(n.is_valid());
    // Grab index, check if valid, and see if node pointer is same
    unsigned int index = n.index();
    if (index < size()) {
      // INFO("Check if node is same as node stored");
      return ((*node_vec_)[index] == n);
    } else {
      ERROR("Index out of bounds.");
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
    }
    /** Default destructor */
    ~Edge() = default;

    /** Returns FALSE if edge is invalid, otherwise TRUE */
    bool is_valid() const {
      // Check for valid pointer
      if (graph_ptr_ == nullptr) {return false;}
      // Check if edge has been deleted
      if (internal_()->idx >= graph_ptr_->num_edge_) {return false;}
      return true;
    }

    /** Return a node of this Edge */
    Node node1() const {
      if (!is_valid()) {return Node();} // Invalid Node
      auto& n_vec = *graph_ptr_->node_vec_;
      auto& nint_vec = *graph_ptr_->node_int_vec_;
      size_type idx = nint_vec[n1int_idx_].idx;
      return n_vec[idx];
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (!is_valid()) {return Node();} // Invalid Node
      auto& n_vec = *graph_ptr_->node_vec_;
      auto& nint_vec = *graph_ptr_->node_int_vec_;
      size_type idx = nint_vec[n2int_idx_].idx;
      return n_vec[idx];
    }

    size_type index() const {
      if (!is_valid()) {throw std::runtime_error("Invalid edge index retrieval.");}
      return internal_()->idx;
    }

    /** Returns a reference to the current value */
    edge_value_type& value() {
      return internal_()->value;
    }

    /** Returns the current value  */
    const edge_value_type& value() const {
      return internal_()->value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (!is_valid() || !e.is_valid()) {
        INFO("Invalid edge");
        return false;
      } // Check for invalid edges
      // Check for same graph
      if (graph_ptr_ != e.graph_ptr_) {
        INFO("Not same graph");
        return false;
      }
      // Check for node equality permutations
      bool perm_1 = (node1() == e.node1()) && (node2() == e.node2());
      bool perm_2;
      if (graph_ptr_->is_directed_) {perm_2 = false;}
      else {perm_2 = (node2() == e.node1()) && (node1() == e.node2());}
      INFO_EVAL(perm_1);
      INFO_EVAL(perm_2);
      return perm_1 || perm_2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // Compare by function of internal index and graph id
      // 1e8 * gid + idx
      return 1e8 * graph_ptr_->gid_ + eint_idx_ < 1e8 * e.graph_ptr_->gid_ + e.eint_idx_;
    }

    private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Returns pointer to internal object for ease-of-access
    edge_internal* internal_() const {
      return &(*graph_ptr_->edge_int_vec_)[eint_idx_];
    }

    // Private method to invalidate current edge
    void invalidate_() {
      internal_()->idx = graph_ptr_->num_edge_;
      graph_ptr_ = nullptr;
    }
    // Private members
    Graph* graph_ptr_;
    size_type eint_idx_;                 // Index of edge internal
    size_type n1int_idx_, n2int_idx_;          // Index of node internals
    // Construction of valid edge
    Edge(Graph* g, unsigned int i, size_type n1, size_type n2) : graph_ptr_(g), eint_idx_(i), n1int_idx_(n1), n2int_idx_(n2) {}
  };

  /** Return the total number of edges in the graph.
   *
    * Complexity: O(1)
   */
  size_type num_edges() const {
    return num_edge_;
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
    auto a_search = edge_map_->find(a.nint_idx_);
    if (a_search != edge_map_->end()) {
      auto b_search = a_search->second.find(b.nint_idx_);
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& val = edge_value_type()) {
    // INFO("Attempt to add edge between nodes " << a.index() << " and " << b.index());
    // Check if nodes in graph
    if (!has_node(a) || !has_node(b)) {
      // WARNING("Node(s) do not exist in graph. No edge will be added.");
      return Edge(); // Invalid Edge
    }
    // Check if self loop
    if (a == b) {
      return Edge();
    }
    // Check if edge already in graph
    if (has_edge(a, b)) {
      // WARNING("Edge already exists. No edge will be added.");
      return (*edge_map_)[a.nint_idx_][b.nint_idx_];
    }
    if (num_edges() >= 1e8) {
      // Edge limit exceeded, no edge created
      return Edge();
    }
    // Grab node data from existing database
    Node n1 = (*node_vec_)[a.index()];
    Node n2 = (*node_vec_)[b.index()];
    // Check for node equality
    assert(n1 == a);
    assert(n2 == b);
    // Create new edge
    size_type idx =  num_edge_++;
    // Save to graph members
    edge_vec_->push_back(Edge(this, (*edge_int_vec_).size(), a.nint_idx_, b.nint_idx_));
    // If there were deleted edges, swap with last invalid element, such that index matches.
    if ((*edge_vec_).size() != idx) {
      std::swap((*edge_vec_)[idx], *((*edge_vec_).end() - 1));
    }
    // Save reverse map if not directed
    if (!is_directed_) {
      (*edge_map_)[b.nint_idx_][a.nint_idx_]  = Edge(this, (*edge_int_vec_).size(), b.nint_idx_, a.nint_idx_);
    }
    edge_int_vec_->push_back(edge_internal(val, idx));
    (*edge_map_)[a.nint_idx_][b.nint_idx_] = (*edge_vec_)[idx];
    return (*edge_vec_)[idx];
  }

  /** Remove an edge with two nodes from the graph
   * @pre @a n1 and @a n2 are valid nodes in the current graph
   *
   * @post The edge between @a n1 and @a n2 is removed (if it exists)
   * @post The incident iterator referring to this edge is invalidated
   * @return The number of edges removed, 0 if none are removed.
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // Check if nodes exists
    if (!has_node(n1) || !has_node(n2)) {
      INFO("Nodes do not exist");
      return 0;
    }
    // Check if edge exists
    if (!has_edge(n1, n2)) {
      INFO("Edge does not exist");
      return 0;
    }
    // Grab edge index
    size_type idx = (*edge_map_)[n1.nint_idx_][n2.nint_idx_].index();
    // Remove from map
    (*edge_map_)[n1.nint_idx_].erase(n2.nint_idx_);
    if (!is_directed_) {(*edge_map_)[n2.nint_idx_].erase(n1.nint_idx_);}
    // Swap with last valid element and invalidate
    (*edge_vec_)[idx].invalidate_();
    std::swap((*edge_vec_)[idx], (*edge_vec_)[--num_edge_]); // Also subtracts count for num_edge_
    // Make index match for swapped node
    size_type int_idx = (*edge_vec_)[idx].eint_idx_;
    (*edge_int_vec_)[int_idx].idx = idx;
    assert(!(*edge_vec_)[num_edge_].is_valid());  // Check that edge vector is invalidated
    return 1;
  }

  /** Remove an edge with an edge reference
   * @pre Edge is valid
   *
   * @post The edge is removed (if it exists)
   * @return The number of edges removed, 0 if none are removed.
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

    /** Removes an edge with a valid edge iterator
   * @pre Edge iterator is valid
   *
   * @post Edge referred by the edge iterator is removed and invalidated
   * @post The incident iterator referring to the edge is invalidated
   * @post The original iterator @a e_it is still valid (and is the same as the output iterator)
   * @return A valid edge iterator with the same index OR the new edge_end() iterator
   */
  EdgeIterator remove_edge(EdgeIterator e_it) {
    if (e_it == edge_end()) {return edge_end();} // No edge removed
    // Grab edge
    Edge e = (*e_it);
    size_type idx = e.index();
    assert(e.is_valid());
    // Remove edge
    remove_edge(e);
    // Return edge iterator at same place
    return edge_begin() + idx;
  }

  /** Removes an edge with from an edge iterator
   * @pre The edge iterator is valid
   *
   * @post The edge referred to by the edge iterator is removed
   * @return A valid edge iterator that can iterate through the remaining edges
   */

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    delete edge_map_;
    delete edge_int_vec_;
    delete edge_vec_;
    delete node_int_vec_;
    delete node_vec_;
    // Re-assign
    node_vec_ = new std::vector<Node>;
    node_int_vec_ = new std::vector<node_internal>;
    edge_vec_ = new std::vector<Edge>;
    edge_int_vec_ = new std::vector<edge_internal>;
    edge_map_ = new std::unordered_map<size_type, std::unordered_map<size_type, Edge>>;
    // Reset counters for node and edges
    num_node_ = 0;
    num_edge_ = 0;
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
    NodeIterator() {}

    /** De-reference NodeIterator for the corresponding Node object */
    Node operator*() const {return (*ptr_->node_vec_)[n_idx_];}

    /** Moves to next node */
    NodeIterator& operator++() {
      n_idx_++;
      return *this;
    }

    /** Post increment operator */
    NodeIterator& operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    /** Jumps to future node */
    NodeIterator& operator+(int i) {
      n_idx_ += i;
      return *this;
    }

    /** Jumps to previous node */
    NodeIterator& operator-(int i) {
      n_idx_ -= i;
      return *this;
    }

    /** Checks for equality of NodeIterator */
    bool operator==(const NodeIterator& ni) const {return n_idx_ == ni.n_idx_;}

    bool operator!=(const NodeIterator& ni) const {return !(*this == ni);}

   private:
    friend class Graph;
    // Valid constructor for creating node iterator
    NodeIterator(const size_type idx, const Graph* g) : n_idx_(idx), ptr_(g) {}
    // Store vector iterator
    size_type n_idx_;
    const Graph* ptr_;
  };

  /** Creates a NodeIterator from the beginning of node_vec. */
  node_iterator node_begin() const {return NodeIterator(0, this);}

  /** Creates a NodeIterator for one past the end of node_vec */
  node_iterator node_end() const {return NodeIterator(num_nodes(), this);}

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

    using map_it_type = typename std::unordered_map<size_type, Edge>::iterator;
    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** De-references incident iterator into Edge object */
    Edge operator*() const {
      return m_it_->second; // No need to check for swaps, the map always saves edge of correct permutation
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

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereference to Edge object
     * MUST be valid index entry to dereference
     */
    Edge operator*() const {return (*ptr_->edge_vec_)[edge_idx_];}

    /** Increments to next edge */
    EdgeIterator& operator++() {
      edge_idx_++;
      return *this;
    }

    /** Post increment operator */
    EdgeIterator& operator++(int) {
      auto tmp = *this;
      ++*this;
      return tmp;
    }

    /** Jumps to a edge further down */
    EdgeIterator& operator+(int i) {
      edge_idx_ += i;
      return *this;
    }

    /** Jumps to a previous edge */
    EdgeIterator& operator-(int i) {
      edge_idx_ -= i;
      return *this;
    }

    /** Check for equality of EdgeIterator */
    bool operator==(const EdgeIterator& ei) const {
      // Return equality of indices
      return edge_idx_ == ei.edge_idx_;}

    bool operator!=(const EdgeIterator& ei) const {return !(*this == ei);}

    private:
    friend class Graph;
    // Valid constructor
    EdgeIterator(const size_type idx, const Graph* g) : edge_idx_(idx), ptr_(g) {}
    // Private members
    size_type edge_idx_;  // Edge index
    const Graph* ptr_ = g_ptr;  // Pointer to original graph
  };

  /** Creates a EdgeIterator from beginning of edge_vec. */
  edge_iterator edge_begin() const {return EdgeIterator(0, this);}

  /** Creates a EdgeIterator from end of edge_vec */
  edge_iterator edge_end() const {return EdgeIterator(num_edges(), this);}

 private:
  /** Variables:
   *      Node vector: vector of node objects
   *      Node internal vector: vector of node_internal objects
   *      Point vector: Vector of point objects
   *      Edge vector: vector of edge objects
   *      Edge internal vector: vector of edge_internal objects
   *      Edge map   : Adjacency list (using maps to indicate what edge corresponds to the node pair) for searching
   *                        The edge index returned corresponds to the Edge object stored in the Edge vector
   *      num_node, num_edge: Number of nodes / edges in the graph. Can be different from corresponding vector.size()
   *  Struct:
   *      Node internal: Contains all of the attributes of Node
   */
  std::vector<Node>* node_vec_;             // Contains all nodes (up to the last clear())
  std::vector<node_internal>* node_int_vec_;  // Contains all node internal (up to the last clear())
  std::vector<Point>* point_vec_;         // TODO: REMOVE IF UNUSED
  std::vector<Edge>* edge_vec_;     // Contains all edges (up to last clear())
  std::vector<edge_internal>* edge_int_vec_;    // Contains all edge internal (up to last clear())
  std::unordered_map<size_type, std::unordered_map<size_type, Edge>>* edge_map_;  // Contains all valid edges (contains reverse pair if undirected)
  unsigned int num_node_ {},  num_edge_ {};
  bool is_directed_ = false; // False: Graph is not directed.
  struct node_internal {
    node_value_type value;
    Point point;
    size_type idx;
    node_internal(node_value_type val, Point p, size_type i) : value(val), point(p), idx(i) {}
  };
  struct edge_internal {
    edge_value_type value;
    unsigned int idx;
    edge_internal(edge_value_type val, unsigned int i) :
      value(val), idx(i) {}
  };
};

/** Counter object for curating graph ID
 * Inspired from https://stackoverflow.com/questions/7097679/simplest-way-to-count-instances-of-an-object
 */
template <class T>
class count_T {
  public:
  count_T() {++num_T_;} // Constructor count
  count_T(const count_T& obj) {if (this != &obj) {++num_T_;}}
  ~count_T() = default; // Leave default constructor
  static unsigned int get_count() {return num_T_;}
  private:
  static unsigned int num_T_;
};

template <class T>
unsigned int count_T<T>::num_T_ = 0;

#endif // CME212_GRAPH_HPP
