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
template <typename V>
class Graph {
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

  /** Type of node value (mass, temperature, color, etc). */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator (following STL conventions). */
  using node_iterator = NodeIterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator (following STL conventions). */
  using edge_iterator = EdgeIterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator (following STL conventions). */
  using incident_iterator = IncidentIterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

 private:

  //
  // PRIVATE MEMBER DEFINITIONS
  //

  // Ordered map (indexed) containing key-value pair ID: Position
  std::map<size_type,Point> nodes_;

  // Ordered map (indexed) containing key-value pair ID: Value
  std::map<size_type,node_value_type> values_;

  // Ordered map (indexed) containing key-value pair ID: Incident Edges
  std::map<size_type,std::vector<edge_type>> incidents_;

  // Number of nodes currently on the graph
  size_type num_nodes_;

  // Number of nodes that have been added since the graph was created
  size_type num_nodes_gen_;
 
  // Number of edges currently on the graph
  size_type num_edges_;

  // Vector of each edge on the graph
  std::vector<Edge> edges_;

 public:

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : num_nodes_(0), num_nodes_gen_(0), num_edges_(0) {}

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
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph_->nodes_[id_];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return std::distance(graph_->nodes_.begin(),graph_->nodes_.find(id_));
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Return this node's value (with type defined by user). */
    node_value_type& value() {
      return graph_->values_[id_];
    }

    const node_value_type& value() const {
      return graph_->values_[id_];
    }

    /** Return the number of edges connected to this node. */
    size_type degree() const {
      return graph_->incidents_[id_].size();
    }

    /** Return incident_iterator pointing to first edge
        connected to this node.
     */
    IncidentIterator edge_begin() const {
      auto first = &graph_->incidents_[id_].front();
      return IncidentIterator(const_cast<edge_type*>(first));
    }

    /** Return incident_iterator pointing to one past
        the last edge connected to this node.
     */
    IncidentIterator edge_end() const {
      auto last = &graph_->incidents_[id_].back();
      return IncidentIterator(const_cast<edge_type*>(++last));
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return ((index() == n.index()) and (graph_ == n.graph()));
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
      if (graph_ == n.graph()) {return id_ < n.id();}
      else {return (graph_ < n.graph());}
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // Pointer to graph in which node exists
    graph_type* graph_;

    // Unique ID for node (different from current index on the graph)
    size_type id_;

    // Constructor for valid node
    Node(const graph_type* graph, size_type id)
      : graph_(const_cast<graph_type*>(graph)), id_(id) {
    }
    
    /** Return pointer to the graph for the node. */
    const graph_type* graph() const {
      return graph_;
    }

    /** Return unique ID of the node. */
    size_type id() const {
      return id_;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return num_nodes_;
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
  node_type add_node(const Point& position) {
    nodes_.insert(std::pair<size_type,Point>(num_nodes_gen_,position));
    values_.insert(std::pair<size_type,node_value_type>(num_nodes_gen_,0));
    num_nodes_ += 1;
    num_nodes_gen_ += 1;
    return Node(this, num_nodes_gen_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (this == n.graph());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    assert(i < size());
    auto it = nodes_.begin();
    std::advance(it,i);
    return Node(this, it->first);
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
    node_type node1() const {
      return node1_;
    }

    /** Return the other node of this Edge */
    node_type node2() const {
      return node2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (((node1_ == e.node1()) and (node2_ == e.node2())) or
        (node1_ == e.node2() and (node2_ == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1_ < e.node1()) {return true;}
      else if (e.node1() < node1_) {return false;}
      else {return node2_ < e.node2();}
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    
    // Edge stores two node objects
    node_type node1_;
    node_type node2_;

    // Constructor for valid edge
    Edge(const Node& a, const Node& b)
      : node1_(a), node2_(b) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  edge_type edge(size_type i) const {
    assert(i >= 0 and i < num_edges());
    return edges_[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (size_type i = 0; i < num_edges(); i++) {
      edge_type e = edges_[i];
      if ((e.node1() == a and e.node2() == b) or 
        (e.node2() == a and e.node1() == b)) {
        return true;}
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
  edge_type add_edge(const Node& a, const Node& b) {
    //assert(has_node(a) and has_node(b));
    //assert((a == b) == false); // cannot add edge between the same node
    edge_type e = Edge(a,b);
    if (has_edge(a,b) == false) {
      // only add the edge if it doesn't already exist
      num_edges_ += 1;
      edges_.push_back(e);
      // add this edge to incident container for Node a and Node b
      incidents_[a.id_].push_back(e);
      e = Edge(b,a);
      incidents_[b.id_].push_back(e);
    }
    return e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    num_nodes_ = 0;
    num_edges_ = 0;
    nodes_.clear();
    edges_.clear();
    values_.clear();
    incidents_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                                   // Element type
    using pointer           = std::pair<const size_type, Point>*;     // Pointers to elements
    using reference         = Node&;                                  // Reference to elements
    using difference_type   = std::ptrdiff_t;                         // Signed difference
    using iterator_category = std::forward_iterator_tag;                // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    /** Increment the node_iterator in place. */
    node_iterator& operator++() {
      std::pair<const size_type,Point> elem = *ptr_;
      auto it = graph_ptr_->nodes_.find(elem.first);
      ptr_ = const_cast<pointer>(&(*(++it)));
      return *this;
    }

    /** Test whether this iterator and @a iter are equal. */
    bool operator==(const node_iterator& iter) const {
      return ptr_ == iter.ptr_;
    }

    /** Dereference the node_iterator, returning a Node type. */
    value_type operator*() {
      std::pair<const size_type,Point> elem = *ptr_;
      return Node(graph_ptr_, elem.first);
    }

   private:
    friend class Graph;

    // Private data
    const graph_type* graph_ptr_;
    pointer ptr_;

    // Constructor
    NodeIterator(const graph_type* graph_ptr, pointer ptr) 
      : graph_ptr_{const_cast<graph_type*>(graph_ptr)}, ptr_{ptr} {}
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Return a node_iterator pointing to the first edge
      connected to this node.
   */
  node_iterator node_begin() const {
    auto first = &(*nodes_.begin());
    return NodeIterator(this, const_cast<std::pair<const size_type,Point>*>(first));
  }
  /** Return a node_iterator pointing to one past the last
      edge connected to this node.
   */
  node_iterator node_end() const {
    auto last = &(*nodes_.end());
    return NodeIterator(this, const_cast<std::pair<const size_type,Point>*>(last));
  }

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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Increment the iterator in place. */
    incident_iterator& operator++() {
      ++ptr_;
      return *this;
    }

    /** Test whether this iterator is equal to @a iter. */
    bool operator==(const incident_iterator& iter) const {
      return ptr_ == iter.ptr_;
    }

    /** Dereference this iterator, returning an Edge. */
    value_type operator*() {
      return *ptr_;
    }

   private:
    friend class Graph;
    
    pointer ptr_;

    IncidentIterator(pointer ptr)
      : ptr_{ptr} {}

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

    /** Increment the iterator in place. */
    edge_iterator& operator++() {
      ++ptr_;
      return *this;
    }

    /** Test whether this iterator is equal to @a iter. */
    bool operator==(const edge_iterator& iter) const {
      return ptr_ == iter.ptr_;
    }

    /** Dereference this iterator, returning an Edge. */
    value_type operator*() {
      return *ptr_;
    }

   private:
    friend class Graph;
    
    pointer ptr_;

    // Constructor
    EdgeIterator(pointer ptr)
      : ptr_{ptr} {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Return an edge_iterator pointing to the first edge
      of the graph.
   */
  edge_iterator edge_begin() const {
    auto first = &edges_.front();
    return EdgeIterator(const_cast<edge_type*>(first));
  }

  /** Return an edge_iterator pointing to one past
      the last edge of the graph.
   */
  edge_iterator edge_end() const {
    auto last = &edges_.back();
    return EdgeIterator(const_cast<edge_type*>(++last));
  }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
