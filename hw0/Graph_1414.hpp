#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

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
  /** Predeclaration of internal_node and internal_edge types. */
  struct internal_node;
  struct internal_edge;

 public:
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
  Graph() : nodes_(), uid_to_idx(), uid_tracker(0) {}

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
  class Node { // Proxy Node actually
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
    Node() {}
    /** Define fetch() similarly to in proxy_example. Fetch helps us get the
     * internal node given a node.
     */
    internal_node& fetch() const {
      return parent_graph_->nodes_.at(parent_graph_->uid_to_idx.at(uid_));
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return parent_graph_->uid_to_idx[fetch().uid_];
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
      if (fetch().uid_ == n.uid_ &&
          this->parent_graph_ == n.parent_graph_) {
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
      size_type myuid = fetch().uid_;
      return myuid < n.uid_;
    }

   private:
    friend class Graph;
    size_type uid_;
    graph_type* parent_graph_;

    Node(const graph_type* graph, size_type uid)
      : uid_(uid), parent_graph_(const_cast<graph_type*>(graph)) {}
  };

  /** Return the number oft nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
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
  Node add_node(const Point& position) {
    size_type assign_uid = this->uid_tracker;
    uid_tracker += 1;
    size_type idx = nodes_.size();
    uid_to_idx[assign_uid] = idx;

    internal_node new_internal_node;
    new_internal_node.pos_ = position;
    new_internal_node.uid_ = assign_uid;
    this->nodes_.emplace_back(new_internal_node);

    return Node(this, assign_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type n_uid = n.uid_;
    if (uid_to_idx.find(n_uid) == uid_to_idx.end()) {
      return false;
    }

    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < num_nodes());
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
  class Edge {
   public:
    /** Construct an invalid Edge. */
    Edge() {}

    /** Help function to get internal_edge given an Edge object. Similar to
     *  fetch() in Node class, and fetch() in proxy_example.cpp.
     */
    internal_edge& edge_fetch() const {
      return edge_p_graph_->edges_.at(edge_p_graph_->eid_to_idx.at(eid_));
    }

    /** Return a node of this Edge */
    Node node1() const {
      return edge_fetch().n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return edge_fetch().n2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      size_type myeid = edge_fetch().eid_;
      if (myeid == e.eid_) {
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
      size_type myeid = edge_fetch().eid_;
      if (myeid < e.eid_) {
        return true;
      }

      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* edge_p_graph_;
    size_type eid_;

    // Initialized for Edge.
    Edge(const graph_type* graph, size_type eid)
      : edge_p_graph_(const_cast<graph_type*>(graph)), eid_(eid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    size_type myeid = edges_[i].eid_;
    return Edge(this, myeid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if nodes a and b exist
    assert(has_node(a));
    assert(has_node(b));

    std::set<size_type> edge_ab = {a.uid_, b.uid_};
    std::set<size_type> edge_ba = {b.uid_, a.uid_};

    // Since we add both edges to the graph, we need really only check one,
    // but this is more modular and robust!
    if (map_edges.find(edge_ab) == map_edges.end()
     && map_edges.find(edge_ba) == map_edges.end()) {
      return false;
    }

    return true;
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
    assert(a.uid_ != b.uid_); // a and b are distinct nodes
    assert(has_node(a)); // node a is valid
    assert(has_node(b)); // node b is valid

    std::set<size_type> edge_ab = {a.uid_, b.uid_};
    std::set<size_type> edge_ba = {b.uid_, a.uid_};

    // Check if edge already exists
    if (map_edges.find(edge_ab) != map_edges.end()) {
      size_type myidx = eid_to_idx.at(map_edges[edge_ab]);
      return Edge(this, myidx);
    } else if (map_edges.find(edge_ba) != map_edges.end()) {
      size_type myidx = eid_to_idx.at(map_edges[edge_ba]);
      return Edge(this, myidx);
    }

    internal_edge new_internal_edge;
    size_type assign_eid = this->eid_tracker;
    eid_tracker += 1;
    size_type idx = edges_.size();
    eid_to_idx[assign_eid] = idx;
    map_edges[edge_ab] = assign_eid; // undirected graph so a->b and b->a are
    map_edges[edge_ba] = assign_eid; // the same, so same eid_, too!

    new_internal_edge.n1_ = a;
    new_internal_edge.n2_ = b;
    new_internal_edge.eid_ = assign_eid;
    this->edges_.emplace_back(new_internal_edge);

    return Edge(this, assign_eid);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    delete this;
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
  // Nodes data
  std::vector<internal_node> nodes_;
  std::unordered_map<size_type, size_type> uid_to_idx;
  size_type uid_tracker;

  struct internal_node {
    internal_node() {}
    internal_node(size_type uid, const Point& pos) : uid_(uid), pos_(pos) {}
    size_type uid_;
    Point pos_;
  };

  // Edges data
  std::vector<internal_edge> edges_;
  size_type eid_tracker;
  std::unordered_map<size_type, size_type> eid_to_idx;
  std::map<std::set<size_type>, size_type> map_edges; // {(n1, n2): eid, etc}

  struct internal_edge {
    internal_edge() {}
    internal_edge(size_type eid, Node n1, Node n2)
      : eid_(eid), n1_(n1), n2_(n2) {}
    size_type eid_;
    Node n1_;
    Node n2_;
  };
};

#endif // CME212_GRAPH_HPP