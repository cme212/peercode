
#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

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

  /** User-specified attribute such as temperature, color, mass, etc. */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() : nodes_(), uid_to_idx(), uid_tracker(0), edges_(), eid_tracker(0),
            eid_to_idx(), map_edges() {}

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
  class Node : private totally_ordered<Node> { // Proxy Node actually
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

    // Return the degree of a node (i.e. number of edges it is incident to).
    size_type degree() const {
      return fetch().adjacent_nodes_.size();
    }

    /** Return an iterator pointing at the beginning of the adjacency list of
     * the given node.
     */
    incident_iterator edge_begin() const {
      return IncidentIterator(fetch().adjacent_nodes_.begin(), parent_graph_, *this);
    }

    /** Return an iterator pointing at the end of the adjacency list of the
     * given node.
     */
    incident_iterator edge_end() const {
      return IncidentIterator(fetch().adjacent_nodes_.end(), parent_graph_, *this);
    }

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

    /** Obtain the value from internal node while not requiring constantness. */
    node_value_type& value() {
      return fetch().value_;
    }

    /** Obtain the (const) value from internal node. */
    const node_value_type& value() const {
      return fetch().value_;
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
  Node add_node(const Point& position,
                const node_value_type& value = node_value_type()) {
    size_type assign_uid = this->uid_tracker;
    uid_tracker += 1;
    size_type idx = nodes_.size();
    uid_to_idx[assign_uid] = idx;

    internal_node new_internal_node;
    new_internal_node.pos_ = position;
    new_internal_node.uid_ = assign_uid;
    new_internal_node.value_ = value;
    new_internal_node.adjacent_nodes_ = {};
    this->nodes_.push_back(new_internal_node);

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
  class Edge : private totally_ordered<Edge> {
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

    std::pair<size_type, size_type> edge_ab = {a.uid_, b.uid_};
    std::pair<size_type, size_type> edge_ba = {b.uid_, a.uid_};

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

    std::pair<size_type, size_type> edge_ab = {a.uid_, b.uid_};
    std::pair<size_type, size_type> edge_ba = {b.uid_, a.uid_};

    /** Check if edge (a,b) already exists, and if it does: return it.
     * Otherwise, if (b,a) exists, delete (a,b) from map_edges and add (b,a)
     * with the same assign_eid.
     * Also update edges_ to swap n1_ and n2_ in the internal_edge.
     */
    if (map_edges.find(edge_ab) != map_edges.end()) {
      size_type myidx = eid_to_idx.at(map_edges[edge_ab]);
      return Edge(this, myidx);
    } else if (map_edges.find(edge_ba) != map_edges.end()) {
      // store eid and idx, so we can use it when adding (b,a)
      size_type eid_ab = map_edges[edge_ba];
      size_type myidx = eid_to_idx.at(eid_ab);
      // delete (a,b) from map_edges and add (b,a)
      map_edges.erase(edge_ba);
      map_edges[edge_ab] = eid_ab;
      // Change n1_ and n2_ in edges_ (we found (b,a), so add (a,b).)
      edges_[myidx].n1_ = a;
      edges_[myidx].n2_ = b;

      return Edge(this, myidx);
    }

    /** Add node b to adjacent_nodes of a is there is an edge a->b, and also
     * the other way around because the graph is undirected.
     */
    (a.fetch().adjacent_nodes_).push_back(b);
    (b.fetch().adjacent_nodes_).push_back(a);

    // Create new internal_edge and add do edges_ container.
    internal_edge new_internal_edge;
    size_type assign_eid = this->eid_tracker;
    eid_tracker += 1;
    size_type idx = edges_.size();
    eid_to_idx[assign_eid] = idx;
    map_edges[edge_ab] = assign_eid;

    new_internal_edge.n1_ = a;
    new_internal_edge.n2_ = b;
    new_internal_edge.eid_ = assign_eid;
    this->edges_.push_back(new_internal_edge);

    return Edge(this, idx);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects. Do this by clearing the
   * edges, nodes containers and the corresponding ids to indices maps. Also
   * reset the trackers :).
   */
  void clear() {
    edges_.clear();
    nodes_.clear();
    uid_to_idx.clear();
    eid_to_idx.clear();
    map_edges.clear();
    uid_tracker = 0;
    eid_tracker = 0;
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    // Dereference the iterator to obtain the node it is pointing at.
    Node operator*() const {
      size_type node_id = (*it_).uid_;
      return Node(ptr_graph_, node_id);
    }

    // Increment the iterator to the next node in the node container.
    NodeIterator& operator++() {
      it_++;
      return *this;
    }

    // Check if two node iterations are equal.
    bool operator==(const NodeIterator& node2) const {
      return it_ == node2.it_;
    }

   private:
    friend class Graph;
    /** A node iterator uses the underlying iterator of vectors, while also
     * keeping track of the parent graph, so we can refer to Node objects.
     */
    typename std::vector<internal_node>::iterator it_;
    graph_type* ptr_graph_;
    NodeIterator(typename std::vector<internal_node>::iterator it,
                 const graph_type* ptr_graph) :
        it_(it), ptr_graph_(const_cast<graph_type*>(ptr_graph)) {}
  };

  /** Given a node iterator, return an iterator to the beginning of the
   * node-container.
   */
  node_iterator node_begin() {
    return NodeIterator(nodes_.begin(), this);
  }

  /** Given a node iterator, return an iterator to the end of the
   * node-container.
   */
  node_iterator node_end() {
    return NodeIterator(nodes_.end(), this);
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
    IncidentIterator() {}

    // Dereference IncIt to obtain an incident edge.
    Edge operator*() const {
      Node snd = Node(inc_parent_graph_, (*inc_it_).uid_);
      return inc_parent_graph_->add_edge(orig_node, snd);
    }

    // Increment IncIt to point to next incident edge/the corresponding node.
    IncidentIterator& operator++() {
      inc_it_++;
      return *this;
    }

    // Check if two IncIt are equal.
    bool operator==(const IncidentIterator& it2) const {
      return inc_it_ == it2.inc_it_;
    }

   private:
    friend class Graph;

    typename std::vector<Node>::iterator inc_it_;
    graph_type* inc_parent_graph_;
    Node orig_node;
    IncidentIterator(typename std::vector<Node>::iterator inc_it,
                     const graph_type* ptr_graph, Node or_n) :
        inc_it_(inc_it), inc_parent_graph_(const_cast<graph_type*>(ptr_graph)),
        orig_node(or_n) {}
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
    EdgeIterator() {}

    // Dereference an edge iterator.
    Edge operator*() const {
      size_type edge_id = (*eit_).eid_;
      return Edge(eptr_graph_, edge_id);
    }

    // Increment an edge iterator to the next edge in the edge container.
    EdgeIterator& operator++() {
      eit_++;
      return *this;
    }

    // Check if two edge iterators are equal.
    bool operator==(const EdgeIterator& edge2) const {
      return eit_ == edge2.eit_;
    }

   private:
    friend class Graph;
    /** An edge iterator inherits its underlying iterator from the std::vector
     * iterator. It contains a pointer to its parent graph so that we can
     * dereference the iterator.
     */
    typename std::vector<internal_edge>::iterator eit_;
    graph_type* eptr_graph_;
    EdgeIterator(typename std::vector<internal_edge>::iterator it,
                 const graph_type* ptr_graph) :
          eit_(it), eptr_graph_(const_cast<graph_type*>(ptr_graph)) {}
  };

  // Obtain an iterator pointing to the start of the edges container.
  edge_iterator edge_begin() {
    return EdgeIterator(edges_.begin(), this);
  }

  // Obtain an iterator pointing to the end of the edges container.
  edge_iterator edge_end() {
    return EdgeIterator(edges_.end(), this);
  }

 private:
  // Nodes data
  std::vector<internal_node> nodes_;
  std::unordered_map<size_type, size_type> uid_to_idx;
  size_type uid_tracker;

  /** An internal node contains a position, a type value and an id (uid).
   * Furthermore, we keep track of all the adjacent nodes.
   */
  struct internal_node {
    node_value_type value_;
    size_type uid_;
    Point pos_;
    std::vector<Node> adjacent_nodes_;
  };

  // Edges data
  std::vector<internal_edge> edges_;
  size_type eid_tracker;
  std::unordered_map<size_type, size_type> eid_to_idx;
  std::map<std::pair<size_type, size_type>, size_type> map_edges; // {(n1, n2): eid, etc}

  // an internal edge contains two nodes, and an id (= eid).
  struct internal_edge {
    size_type eid_;
    Node n1_;
    Node n2_;
  };
};

#endif // CME212_GRAPH_HPP