#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.

  /** Predeclaration of node internal struct. */
  struct internal_node;

  /** Predeclaration of edge helper structs. */
  struct NodePair;
  struct hash_func;

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
  using size_type = unsigned;  // evaluates to unsigned int

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {  // initialized using defaults of STL containers
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODE PROXY
  //

  /** @class Graph::Node
   * @brief Class representing the graph's nodes.
   *
   * Node objects are used to access information about the Graph's nodes.
   */
  class Node {
   private:
     // Allow Graph to access Node's private member data and functions.
     friend class Graph;

     // Use this space to declare private data members and methods for Node
     // that will not be visible to users, but may be useful within Graph.
     // i.e. Graph needs a way to construct valid Node objects

     // Pointer back to the Graph container
     Graph* graph_;

     // This node's unique identification number
     size_type uid_;

     /** Private Constructor */
     Node(const Graph* graph, size_type uid)
       : graph_(const_cast<Graph*>(graph)), uid_(uid) {
     }

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
      graph_ = nullptr;  // initialize node unattached to a graph
      uid_ = -1;         // initialize uid that is easily identified as invalid

    }

    /** Return this node's unique identifier. */
    size_type getid() const {
      return uid_;
    }

    /** Return this node's position. */
    const Point& position() const { // HW0: YOUR CODE HERE

      // extract uid from uid_to_idx hash map
      size_type idx = graph_->uid_to_idx[uid_];

      // extract point info from node internals
      Point& node_position = graph_->nodes_[idx].point_;

      return node_position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {  // HW0: YOUR CODE HERE

      // extract uid from uid_to_idx hash map
      size_type idx = graph_->uid_to_idx[uid_];

      // check if index valid
      size_type graph_size = graph_->nodes_.size();
      if (((int)idx >= 0) and (idx < graph_size)) {
        return idx;
      }

      return size_type(-1);
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
    bool operator==(const Node& n) const {  // HW0: YOUR CODE HERE
       (void) n;  // quiet compiler warning
       if (graph_->uid_to_idx[uid_] == n.index()) {  // check index equality
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
    bool operator<(const Node& n) const {  // HW0: YOUR CODE HERE

      // quiet compiler warning
      (void) n;

      // check less than condition, ordering of nodes
      if (uid_ < n.index()) {
         return true;
      }

      return false;
    }

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes_.size();  // vector container affords us size() method
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
    (void) position;  // quiet compiler warning

    // declare and initialize uid, idx
    size_type new_uid = next_uid_;
    size_type idx = nodes_.size();
    uid_to_idx[new_uid] = idx;

    // instantiate node internals and append to graph
    internal_node new_node(position, new_uid);
    nodes_.push_back(new_node);

    // increment next_uid_ tracker
    ++next_uid_;

    return Node(this, new_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;  // quiet compiler warning
    const size_type node_id = n.getid();  // access id from node class
    if (uid_to_idx.find(node_id) != uid_to_idx.end()) {  // search hash map
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
    (void) i;  // quiet compiler warning
    size_type node_id = nodes_[i].uid_;  // extract uid_

    return Node(this, node_id);  // instantiate proxy node
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
   private:
     // Allow Graph to access Edge's private member data and functions.
     friend class Graph;
     // HW0: YOUR CODE HERE
     // Use this space to declare private data members and methods for Edge
     // that will not be visible to users, but may be useful within Graph.
     // i.e. Graph needs a way to construct valid Edge objects

     // two nodes associated with edge
     Node a_, b_;

     /** Private Constructor */
     Edge(const Node &a, const Node &b)
       : a_(a), b_(b) {
     }


   public:
    /** Construct an invalid Edge. */
    Edge()
      : a_(), b_() {
    }

    /** Return a node of this Edge */
    Node node1() const {  // HW0: YOUR CODE HERE
      return a_;
    }

    /** Return the other node of this Edge */
    Node node2() const {  // HW0: YOUR CODE HERE
      return b_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {  // HW0: YOUR CODE HERE
      (void) e;   // quiet compiler warning

      // check equality of nodes in both directions
      bool case1 = (a_ == e.node1()) and (b_ == e.node2());
      bool case2 = (b_ == e.node1()) and (a_ == e.node2());

      if (case1 or case2) {
        return true;
      }

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {  // HW0: YOUR CODE HERE
      (void) e;           // quiet compiler warning

      // check index of first edge node
      if (a_ < e.node1()) {
        return true;
      }

      return false;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {  // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {  // HW0: YOUR CODE HERE
    (void) i;             // quiet compiler warning
    return Edge(edges_[i].node1(), edges_[i].node2());
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {  // HW0: YOUR CODE HERE
    (void) a; (void) b;   // quiet compiler warning

    // check both directions in nodes_to_eid hash map
    bool case1 = nodes_to_eid.find(NodePair(a, b)) != nodes_to_eid.end();
    bool case2 = nodes_to_eid.find(NodePair(b, a)) != nodes_to_eid.end();

    if (case1 or case2) {
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
  Edge add_edge(const Node& a, const Node& b) {  // HW0: YOUR CODE HERE
    (void) a, (void) b;   // quiet compiler warning

    // instantiate new edge
    Edge e(a, b);

    // check if preexisting
    if (!has_edge(a, b)) {

      // append new edge to edges_
      edges_.push_back(e);

      // declare and initialise eid, idx
      size_type new_eid = next_eid_;
      size_type idx = edges_.size();
      eid_to_idx[new_eid] = idx;
      nodes_to_eid[NodePair(a, b)] = new_eid;

      // increment next_eid
      ++next_eid_;

      return e;
    }

    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {  // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    uid_to_idx.clear();
    eid_to_idx.clear();
    nodes_to_eid.clear();
    next_uid_ = 0;
    next_eid_ = 0;
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

  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    Point point_;
    size_type uid_;

    internal_node(const Point& position, size_type uid) {
      point_ = position;
      uid_ = uid;

    }
  };

  struct NodePair {
    size_type min_uid;
    size_type max_uid;

    NodePair(const Node &a, const Node &b) {
      min_uid = std::min(a.uid_, b.uid_);
      max_uid = std::max(a.uid_, b.uid_);
    }

    // compare keys in hash collisiosn
    bool operator==(const NodePair &np) const {
        return min_uid == np.min_uid and max_uid == np.max_uid;
    }

    // less than operator for edges
    bool operator<(const NodePair &np) const {
      return min_uid < np.min_uid;  // compare smallest nodes
    }
  };

  struct hash_func {
    size_type operator()(const NodePair &np) const {
      // Hash function taken from boost library
      size_type seed = 0;
      std::hash<size_type> hasher;
      seed ^= hasher(np.min_uid) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= hasher(np.max_uid) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };

  /** Declare containers for nodes. */
  size_type next_uid_ = 0;
  std::vector<internal_node> nodes_;
  std::unordered_map<size_type, size_type> uid_to_idx;

  /** Declare containers for edges. */
  size_type next_eid_ = 0;
  std::vector<Edge> edges_;
  std::unordered_map<size_type, size_type> eid_to_idx;
  std::unordered_map<NodePair, size_type, hash_func> nodes_to_eid;

};

#endif // CME212_GRAPH_HPP
