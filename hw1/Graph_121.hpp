#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph
{
private:
  struct internal_node;

public:
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V>;

  /** User specified node value. */
  using node_value_type = V;

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
  Graph() : internal_nodes_(), node_uid_to_index_(), next_node_uid_(0), edges_(), node_pair_to_index_(), incident_nodes_()
  {
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
  class Node : private totally_ordered<Node>
  {
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
    {
    }

    /** Return this node's position. */
    const Point &position() const
    {
      size_type point_index = graph_->node_uid_to_index_.at(uid_);
      return graph_->internal_nodes_.at(point_index).point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const
    {
      return graph_->node_uid_to_index_.at(uid_);
    }

    /**
     * @brief Get value of node by reference. Can be used to set the node value.
     * 
     * @return node_value_type& 
     */
    node_value_type &value()
    {
      size_type node_index = graph_->node_uid_to_index_.at(uid_);
      return graph_->internal_nodes_.at(node_index).node_value;
    }

    /** Return a constant reference to this node´s value. */
    const node_value_type &value() const
    {
      size_type node_index = graph_->node_uid_to_index_.at(uid_);
      return graph_->internal_nodes_.at(node_index).node_value;
    }

    /**
     * @brief Get number of incident nodes to the current node.
     * 
     * @return size_type 
     */
    size_type degree() const
    {
      if (graph_->incident_nodes_.find(uid_) == graph_->incident_nodes_.end())
      {
        return 0;
      }
      return graph_->incident_nodes_.at(uid_).size();
    }

    /**
     * @brief Get start of incident node iterator.
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_begin() const
    {
      return IncidentIterator(graph_, uid_, graph_->incident_nodes_.at(uid_).begin());
    }

    /**
     * @brief Get end of incident node iterator.
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_end() const
    {
      return IncidentIterator(graph_, uid_, graph_->incident_nodes_.at(uid_).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const
    {
      return uid_ == n.uid_ and graph_ == n.graph_;
    }

    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * The node ordering relation must obey trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     */
    bool operator<(const Node &n) const
    {
      return uid_ < n.uid_;
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    graph_type *graph_;
    size_type uid_;

    Node(const graph_type *graph, size_type uid)
        : graph_(const_cast<graph_type *>(graph)), uid_(uid)
    {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const
  {
    return internal_nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const
  {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type &node_value = node_value_type())
  {
    internal_node node = {position, next_node_uid_, node_value};
    internal_nodes_.push_back(node);
    node_uid_to_index_[next_node_uid_] = internal_nodes_.size() - 1;
    next_node_uid_++;
    return Node(this, next_node_uid_ - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const
  {
    return node_uid_to_index_.find(n.uid_) != node_uid_to_index_.end();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const
  {
    return Node(this, internal_nodes_.at(i).uid);
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
  class Edge : private totally_ordered<Edge>
  {
  public:
    /** Construct an invalid Edge. */
    Edge()
    {
    }

    /** Return a node of this Edge */
    Node node1() const
    {
      return a_;
    }

    /** Return the other node of this Edge */
    Node node2() const
    {
      return b_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const
    {
      return NodePair(a_, b_) == NodePair(e.a_, e.b_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const
    {
      return NodePair(a_, b_) < NodePair(e.a_, e.b_);
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Node a_;
    Node b_;

    Edge(const Node &a, const Node &b)
        : a_(a), b_(b)
    {
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const
  {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const
  {
    // HW0: YOUR CODE HERE
    return edges_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const
  {
    // HW0: YOUR CODE HERE
    return node_pair_to_index_.find(NodePair(a, b)) != node_pair_to_index_.end();
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
  Edge add_edge(const Node &a, const Node &b)
  {
    // HW0: YOUR CODE HERE
    NodePair node_pair = NodePair(a, b);

    if (node_pair_to_index_.find(node_pair) != node_pair_to_index_.end())
    {
      return Edge(a, b);
    }

    Edge edge = Edge(a, b);
    edges_.push_back(edge);
    node_pair_to_index_[node_pair] = edges_.size() - 1;

    incident_nodes_[a.uid_].insert(b.uid_);
    incident_nodes_[b.uid_].insert(a.uid_);

    return edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear()
  {
    internal_nodes_.clear();
    node_uid_to_index_.clear();
    edges_.clear();
    node_pair_to_index_.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                             // Element type
    using pointer = Node *;                              // Pointers to elements
    using reference = Node &;                            // Reference to elements
    using difference_type = std::ptrdiff_t;              // Signed difference
    using iterator_category = std::forward_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator()
    {
    }

    /**
     * @brief Dereference the node the iterator is currently pointing to.
     * 
     * @return Node 
     */
    Node operator*() const
    {
      return Node(graph_, internal_nodes_iterator_->uid);
    }

    /**
     * @brief Increment node iterator.
     * 
     * @return NodeIterator& 
     */
    NodeIterator &operator++()
    {
      ++internal_nodes_iterator_;
      return *this;
    }

    /**
     * @brief Check for equality to other node iterator. 
     * Checks that the iterators are iterating over the same container and at the same poisition.
     * 
     * @param node_iterator 
     * @return true 
     * @return false 
     */
    bool operator==(const NodeIterator &node_iterator) const
    {
      return internal_nodes_iterator_ == node_iterator.internal_nodes_iterator_;
    }

  private:
    friend class Graph;

    graph_type *graph_;
    typename std::vector<internal_node>::const_iterator internal_nodes_iterator_;

    NodeIterator(const graph_type *graph, typename std::vector<internal_node>::const_iterator internal_nodes_iterator)
        : graph_(const_cast<graph_type *>(graph)), internal_nodes_iterator_(internal_nodes_iterator)
    {
    }
  };

  /**
   * @brief Get iterator over nodes.
   * 
   * @return node_iterator 
   */
  node_iterator node_begin() const
  {
    return NodeIterator(this, internal_nodes_.begin());
  }

  /**
   * @brief Get end of node iterator.
   * 
   * @return node_iterator 
   */
  node_iterator node_end() const
  {
    return NodeIterator(this, internal_nodes_.end());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator()
    {
    }

    /**
     * @brief Get the edge between the node that was used to create the incident iterator and the incident node.
     * 
     * @return Edge 
     */
    Edge operator*() const
    {
      size_type other_uid = *incident_nodes_iterator_;
      return Edge(Node(graph_, uid_), Node(graph_, other_uid));
    }

    /**
     * @brief Increment the node that is incident to the node that was used to create the incident iterator.
     * 
     * @return IncidentIterator& 
     */
    IncidentIterator &operator++()
    {
      ++incident_nodes_iterator_;
      return *this;
    }

    /**
     * @brief Check for equality between two incident iterators.
     * Inherits equality from vector iterator.
     * 
     * @param incident_iterator 
     * @return true 
     * @return false 
     */
    bool operator==(const IncidentIterator &incident_iterator) const
    {
      return incident_nodes_iterator_ == incident_iterator.incident_nodes_iterator_;
    }

  private:
    friend class Graph;

    graph_type *graph_;
    size_type uid_;
    std::unordered_set<size_type>::iterator incident_nodes_iterator_;

    IncidentIterator(const graph_type *graph, size_type uid, std::unordered_set<size_type>::iterator incident_nodes_iterator)
        : graph_(const_cast<graph_type *>(graph)), uid_(uid), incident_nodes_iterator_(incident_nodes_iterator)
    {
    }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>
  {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator()
    {
    }

    /**
     * @brief Get the edge that the edge iterator is currently pointing to.
     * 
     * @return Edge 
     */
    Edge operator*() const
    {
      return *edges_iterator_;
    }

    /**
     * @brief Increment the edge iterator.
     * 
     * @return EdgeIterator& 
     */
    EdgeIterator &operator++()
    {
      ++edges_iterator_;
      return *this;
    }

    /**
     * @brief Equality between edge iterator.
     * Inherits from equality between two STL vector iterators.
     * 
     * @param edge_iterator 
     * @return true 
     * @return false 
     */
    bool operator==(const EdgeIterator &edge_iterator) const
    {
      return edges_iterator_ == edge_iterator.edges_iterator_;
    }

  private:
    friend class Graph;

    typename std::vector<Edge>::const_iterator edges_iterator_;

    EdgeIterator(typename std::vector<Edge>::const_iterator edges_iterator)
        : edges_iterator_(edges_iterator)
    {
    }
  };

  /**
   * @brief Get start of edge iterator.
   * 
   * @return edge_iterator 
   */
  edge_iterator edge_begin() const
  {
    return EdgeIterator(edges_.begin());
  }

  /**
   * @brief Get end of edge iterator.
   * 
   * @return edge_iterator 
   */
  edge_iterator edge_end() const
  {
    return EdgeIterator(edges_.end());
  }

private:
  // Internal node structure.
  struct internal_node
  {
    Point point;
    size_type uid;
    node_value_type node_value;
  };

  // Containers for nodes.
  std::vector<internal_node> internal_nodes_;
  std::unordered_map<size_type, size_type> node_uid_to_index_;
  size_type next_node_uid_;

  // Struct for containing an ordered pair of nodes.
  struct NodePair
  {
    Node a;
    Node b;

    NodePair(const Node &a, const Node &b) : a(a), b(b)
    {
    }

    size_type min_uid() const
    {
      return std::min(a.uid_, b.uid_);
    }
    size_type max_uid() const
    {
      return std::max(a.uid_, b.uid_);
    }

    // Used to compare keys in hash collisions.
    bool operator==(const NodePair &np) const
    {
      return (a == np.a and b == np.b) or (a == np.b and b == np.a);
    }

    // Used to compare edges.
    bool operator<(const NodePair &np) const
    {
      return true ? min_uid() < np.min_uid() : max_uid() < np.max_uid();
    }
  };

  // The specialized hash function for `unordered_map` keys
  struct hash_fn
  {
    size_type operator()(const NodePair &np) const
    {
      // Hash computation taken from boost library.
      size_type seed = 0;
      std::hash<size_type> hasher;
      seed ^= hasher(np.min_uid()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      seed ^= hasher(np.max_uid()) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
      return seed;
    }
  };

  // Containers for edges.
  std::vector<Edge> edges_;
  std::unordered_map<NodePair, size_type, hash_fn> node_pair_to_index_;
  std::unordered_map<size_type, std::unordered_set<size_type>> incident_nodes_;
};

#endif // CME212_GRAPH_HPP
