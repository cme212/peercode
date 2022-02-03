#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <vector>
#include <map>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

#define make_pair_ordered(x, y) \
  (x < y ? std::make_pair(x, y) : std::make_pair(y, x))

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template<typename V>
class Graph {
 private:
  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Predeclare the internal types
  struct Node_;
  struct Edge_;

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
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : size_(0), num_edges_(0) {
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
  class Node: private totally_ordered<Node> {
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
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uindex_]->position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uindex_]->index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    node_value_type& value(){
      return graph_->nodes_[uindex_]->value;
    }

    // const node_value_type& value() const;
    const node_value_type& value() const{
      return graph_->nodes_[uindex_]->value;
    }

    // size_type degree() const;
    size_type degree() const{
      return graph_->nodes_[uindex_]->incident.size();
    }

    // incident_iterator edge_begin() const;
    incident_iterator edge_begin() const{
      return IncidentIterator(graph_, uindex_, graph_->nodes_[uindex_]->incident.begin());
    }

    // incident_iterator edge_end() const;
    incident_iterator edge_end() const{
      return IncidentIterator(graph_, uindex_,
                              graph_->nodes_[uindex_]->incident.end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_ == n.graph_ && uindex_ == n.uindex_;
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
      return graph_ < n.graph_ || uindex_ < n.uindex_;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type uindex_;

    Node(const graph_type* graph, size_type uindex)
        : graph_(const_cast<graph_type*>(graph)), uindex_(uindex) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
  }

  /** Synonym for size(). */
  size_type num_nodes() const { return size_; }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value =node_value_type()) {
    // HW0: YOUR CODE HERE
    size_type next_uindex = nodes_.size();
    nodes_.push_back(new Node_{size_, position, value, {}});
    index_to_uindex_.push_back(next_uindex);
    size_++;
    return Node(this, next_uindex);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n == node(n.index());
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < size_);
    return Node(this, index_to_uindex_[i]);
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
  class Edge: private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_uindex_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_uindex_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return make_pair_ordered(node1_uindex_, node2_uindex_) ==
             make_pair_ordered(e.node1_uindex_, e.node2_uindex_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE
      return graph_ < e.graph_ || node1_uindex_ < e.node1_uindex_ ||
             node2_uindex_ < e.node2_uindex_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    graph_type* graph_;
    size_type node1_uindex_;
    size_type node2_uindex_;

    Edge(const graph_type* graph, size_type node1_uindex,
         size_type node2_uindex)
        : graph_(const_cast<graph_type*>(graph)),
          node1_uindex_(node1_uindex),
          node2_uindex_(node2_uindex) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i >= 0 && i < num_edges_);
    size_type edge_uindex = edge_index_to_uindex_[i];
    return Edge(this, (edges_[edge_uindex])->node1_uindex,
                (edges_[edge_uindex])->node2_uindex);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    return has_edge_.find(make_pair_ordered(a.uindex_, b.uindex_)) !=
           has_edge_.end();
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
    assert(a != b);
    std::pair<size_type, size_type> this_edge =
        make_pair_ordered(a.uindex_, b.uindex_);
    std::map<std::pair<size_type, size_type>, size_type>::iterator tmp =
        has_edge_.find(this_edge);
    if (tmp != has_edge_.end()) return Edge(this, a.uindex_, b.uindex_);
    size_type next_edge_uindex = edges_.size();
    edges_.push_back(new Edge_{a.uindex_, b.uindex_});
    edge_index_to_uindex_.push_back(next_edge_uindex);
    has_edge_[this_edge] = next_edge_uindex;
    nodes_[a.uindex_]->incident[b.uindex_] = next_edge_uindex;
    nodes_[b.uindex_]->incident[a.uindex_] = next_edge_uindex;
    num_edges_++;
    return Edge(this, a.uindex_, b.uindex_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    num_edges_ = 0;
    nodes_.clear();
    index_to_uindex_.clear();
    edges_.clear();
    edge_index_to_uindex_.clear();
    has_edge_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private equality_comparable<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                            // Element type
    using pointer = Node*;                              // Pointers to elements
    using reference = Node&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    Node operator*() const { return graph_->node(index_); }

    // NodeIterator& operator++()
    NodeIterator& operator++() {
      index_++;
      return *this;
    }

    // bool operator==(const NodeIterator&) const
    bool operator==(const NodeIterator& niter) const {
      return graph_ == niter.graph_ && index_ == niter.index_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    graph_type* graph_;
    size_type index_;

    NodeIterator(const graph_type* graph, size_type index)
        : graph_(const_cast<graph_type*>(graph)), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  node_iterator node_begin() const { return NodeIterator(this, 0); }

  // node_iterator node_end() const
  node_iterator node_end() const { return NodeIterator(this, size_); }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private equality_comparable<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge*;                              // Pointers to elements
    using reference = Edge&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const{
      return Edge(graph_,node_uindex_,incident_iter_->first);
    }

    // IncidentIterator& operator++()
    IncidentIterator& operator++(){
      ++incident_iter_;
      return *this;
    }

    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& i) const{
      return graph_==i.graph_ && node_uindex_==i.node_uindex_ && incident_iter_==i.incident_iter_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph_;
    size_type node_uindex_;
    std::map<size_type,size_type>::iterator incident_iter_;

    IncidentIterator(const graph_type* graph, size_type node_uindex,
                     std::map<size_type, size_type>::iterator incident_iter)
        : graph_(const_cast<graph_type*>(graph)),
          node_uindex_(node_uindex),
          incident_iter_(incident_iter) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private equality_comparable<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                            // Element type
    using pointer = Edge*;                              // Pointers to elements
    using reference = Edge&;                            // Reference to elements
    using difference_type = std::ptrdiff_t;             // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    Edge operator*() const { return graph_->edge(edge_index_); }

    // EdgeIterator& operator++()
    EdgeIterator& operator++() {
      edge_index_++;
      return *this;
    }

    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& eiter) const {
      return graph_ == eiter.graph_ && edge_index_ == eiter.edge_index_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* graph_;
    size_type edge_index_;

    EdgeIterator(const graph_type* graph, size_type edge_index)
        : graph_(const_cast<graph_type*>(graph)), edge_index_(edge_index) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const { return EdgeIterator(this, 0); }

  // edge_iterator edge_end() const
  edge_iterator edge_end() const { return EdgeIterator(this, num_edges_); }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  size_type size_;  // Number of nodes

  size_type num_edges_;  // Number of edges

  struct Node_ {
    size_type index;
    Point position;
    node_value_type value;
    std::map<size_type,size_type> incident; //node2_uindex->edge_uindex
  };

  std::vector<Node_*> nodes_;  //
  std::vector<size_type> index_to_uindex_;

  struct Edge_ {
    size_type node1_uindex;
    size_type node2_uindex;
  };

  std::vector<Edge_*> edges_;
  std::vector<size_type> edge_index_to_uindex_;
  std::map<std::pair<size_type, size_type>, size_type> has_edge_; // (a,b)->edge_uindex, a<b

  // Disable copy and assignment (I don't know why but why not)
  Graph(const Graph&) = delete;
  Graph& operator=(const Graph&) = delete;
};

#endif  // CME212_GRAPH_HPP
