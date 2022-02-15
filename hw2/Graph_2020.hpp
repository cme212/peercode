#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "CME212/Point.hpp"
#include "CME212/Util.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */

template <typename V, typename E> 
class Graph {

struct nodeinfo;
struct edgeinfo;

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
  using edge_value_type = E;


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph *graph_;
  Graph() {
    nodes_ = {};
    active_nodes_ = {};
    edges_ = {};
    active_edges_ = {};
    edgeSet_ = {};
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
    // Graph *graph_;
    Node() {
      uid_ = 0;
      graph_ = NULL;
    }
    

    /** Return this node's position (modifiable). */
    Point &position() { return graph_->nodes_[uid_].p_; }

    /** Return this node's position. */
    const Point &position() const { return graph_->nodes_[uid_].p_; }

    /** Return this node's initial position. */
    const Point &initial_position() const { 
      return graph_->nodes_[uid_].initial_p_; 
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const { return graph_->nodes_[uid_].idx_; }

    /** Return this node's value (modifiable). */
    node_value_type &value() { return graph_->nodes_[uid_].v_; }


    /** Return this node's value. */
    const node_value_type &value() const { return graph_->nodes_[uid_].v_; }

    /** Return this node's degree. */
    size_type degree() const {
      return (graph_->nodes_[uid_].incidents_).size();
    }

    /** begin() function for incident iterator. */
    incident_iterator edge_begin() const {
      return IncidentIterator(0, this, graph_);
    }

    /** end() function for incident iterator. */
    incident_iterator edge_end() const {
      return IncidentIterator(degree(), this, graph_);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node &n) const {
      if (graph_ == n.graph_ && index() == n.index()) {
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
    bool operator<(const Node &n) const { 
      if (index() == n.index()) {
        return (graph_ < n.graph_);
      } else {
       return index() < n.index();
      } 
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type uid_;
    Graph *graph_;

   /** Node constructor */
    Node(size_type n, const Graph *graph) {
      uid_ = n;
      graph_ = const_cast<Graph*> (graph);
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const { return active_nodes_.size(); }

  /** Synonym for size(). */
  size_type num_nodes() const { return size(); }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] v The new node's value (optional)
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position,
                const node_value_type& v = node_value_type()) {
    Node new_node = Node(nodes_.size(), this);
    assert(sizeof(new_node) <= 16);
    nodeinfo info;
    info.p_ = position;
    info.v_ = v;
    info.idx_ = size();
    info.incidents_ = {};
    info.active_ = true;
    info.initial_p_ = position;
    nodes_.push_back(info);
    active_nodes_.push_back(new_node.uid_);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node &n) const {
    if (this == n.graph_ && nodes_[n.uid_].active_) {
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
    assert(i < num_nodes());
    return Node(active_nodes_[i], this);
  }


  /** Delete a node to the graph and output 1.
   * If the node is not in the graph, output 0.
   * All edges that are incident to the node also
   * gets deleted.
   * @param[in] node The node to be deleted
   * @post new num_nodes() == old num_nodes() -1
   *       if the node was on the graph
   *
   * Complexity: O(num_edges)
   */

  size_type remove_node (const Node &n) {
    if (has_node(n)) {
      nodes_[n.uid_].active_ = false;
      //swap and pop back
      active_nodes_[n.index()] = active_nodes_[size() - 1];
      active_nodes_.pop_back();
      nodes_[active_nodes_[n.index()]].idx_ = n.index();
      // remove incident edges  
      incident_iterator cur_incident = n.edge_begin();
      incident_iterator end_incident = n.edge_end();
      while (cur_incident != end_incident) {
        remove_edge(*cur_incident);
        ++cur_incident;
      }
      return 1;
    } else {
      // grap has no node n
      return 0;
    }
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
      uid_ = 0;
      graph_ = NULL;
    }

    /** Return the value of this Edge (modifiable) */
    edge_value_type &value() { return graph_->edges_[uid_].v_; }

    /** Return the value of this Edge */
    const edge_value_type &value() const { return graph_->edges_[uid_].v_; }

    /** Return a node of this Edge */
    Node node1() const { return Node(graph_->edges_[uid_].n1_, graph_); }

    /** Return the other node of this Edge */
    Node node2() const { return Node(graph_->edges_[uid_].n2_, graph_); }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge &e) const {
      return (((node1() == e.node1()) && (node2() == e.node2())) ||
              ((node1() == e.node2()) && (node2() == e.node1())));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge &e) const {
      size_type min1 = std::min(node1().uid_, node2().uid_);
      size_type max1 = std::max(node1().uid_, node2().uid_);
      size_type min2 = std::min(e.node1().uid_, e.node2().uid_);
      size_type max2 = std::max(e.node1().uid_, e.node2().uid_);
      if (graph_ == e.graph_) {
        if (min1 < min2) {
          return true;
        } else if (min1 > min2) {
          return false;
        } else {
          return (max1 < max2);
        }
      } else {return (graph_ < e.graph_); }
    }
   
    /** Return the Euclidean length of this edge */
    double length() { return norm(node1().position() - node2().position()); }

    /** Return the other node of the edge */
    Node otherNode(Node given_node) {
      Node x = node1(), y = node2();

      return x == given_node ? y : x;
    }

  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    size_type uid_;
    Graph *graph_;

    Edge(size_type n, const Graph* graph) {
      uid_ = n;
      graph_ = const_cast<Graph*> (graph);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const { return active_edges_.size(); }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return Edge(active_edges_[i], this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node &a, const Node &b) const {
    size_type x = std::min(a.uid_, b.uid_), y = std::max(a.uid_, b.uid_);
    std::string edge_hash = std::to_string(x) + '-' + std::to_string(y);
    return (edgeSet_.find(edge_hash) != edgeSet_.end());
  }

  /** Add an edge to the graph, or return the current edge if it already
   * exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges()
   * + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node &a, const Node &b, const edge_value_type& v = edge_value_type()) {
    Edge new_edge(edges_.size(), this);
    assert(sizeof(new_edge) <= 32);
    if (!has_edge(a, b)) {
      size_type x = std::min(a.uid_, b.uid_), y = std::max(a.uid_, b.uid_);
      edgeinfo info;
      info.n1_ = a.uid_;
      info.n2_ = b.uid_;
      info.v_ = v;
      info.idx_ = active_edges_.size();
      edges_.push_back(info);
      active_edges_.push_back(new_edge.uid_);
      // create a string hash for edgeSet_
      edgeSet_[std::to_string(x) + '-' + std::to_string(y)] = edges_.size() - 1;
      nodes_[a.uid_].incidents_.push_back(b.uid_);
      if (a.uid_ != b.uid_) { // In case the nodes are the same
        nodes_[b.uid_].incidents_.push_back(a.uid_);
      }
    }
    return new_edge;
  }
 
  /** Delete an edge (a, b) to the graph and output 1.
   * If the edge is not in the graph, output 0.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @param[in] node a A node on the edge to be deleted
   * @param[in] node b The other node on the edge to be deleted
   * @post new num_edges() == old num_edg() -1
   *       if the edge was on the graph
   *
   * Complexity: O(num_edge)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a, b)) {
      size_type x = std::min(a.uid_, b.uid_), y = std::max(a.uid_, b.uid_);
      size_type cur_idx = edgeSet_[std::to_string(x) + '-' + std::to_string(y)];
      std::vector<size_type>& a_inc = nodes_[a.uid_].incidents_;
      std::vector<size_type>& b_inc = nodes_[b.uid_].incidents_;
      active_edges_[cur_idx] = active_edges_.back();
      active_edges_.pop_back();
      edges_[active_edges_[cur_idx]].idx_ = cur_idx;
      edgeSet_.erase(std::to_string(x) + '-' + std::to_string(y));
      a_inc.erase(std::remove(a_inc.begin(), a_inc.end(), cur_idx), a_inc.end());
      b_inc.erase(std::remove(b_inc.begin(), b_inc.end(), cur_idx), b_inc.end());
      return 1;
    } else {
      // grap has no edge (a, b)
      return 0;
    }
  }

  /** Delete an edge to the graph and output 1.
   * If the edge is not in the graph, output 0.
   * @param[in] edge e The edge to be deleted
   * @post new num_edges() == old num_edges() -1
   *       if the edge was on the graph
   *
   * Complexity: O(num_edge)
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes_ = {};
    active_nodes_ = {};
    edges_ = {};
    active_edges_ = {};
    edgeSet_ = {};
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Node;                           // Element type
    using pointer = Node *;                            // Pointers to elements
    using reference = Node &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {}

    /** Dereference */
    Node operator*() const { return Node(index_, graph_); };

    /** Increment */
    NodeIterator &operator++() {
    index_ ++;
      return *this;
    };

 
    /** Equality */
    bool operator==(const NodeIterator &nodeIter) const {
      return (graph_ == nodeIter.graph_) && (index_ == nodeIter.index_);
    };

  private:
    friend class Graph;
    size_type index_;
    Graph* graph_;
    //Private constructor that can be accessed by Graph 
    NodeIterator(size_type index, const Graph* graph) 
      : index_(index), graph_(const_cast<Graph*>(graph)) {
    }
  };

  /** begin() function for node iterator. */
  node_iterator node_begin() const { return NodeIterator(0, this); };

  /** end() function for node iterator. */
  node_iterator node_end() const { return NodeIterator(size(), this);};

  /** remove the node associated to the node iterator. */  
  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return NodeIterator(0, this);
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator.
   */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {}

    /** Dereference */
    Edge operator*() const {
      size_type a = node_->uid_;
      size_type b = graph_->nodes_[a].incidents_[index_];
      size_type x = std::min(a, b), y = std::max(a, b);
      size_type edge_id = graph_->edgeSet_[std::to_string(x) + '-' + std::to_string(y)];
      edgeinfo &info = graph_->edges_[edge_id];

      if (a != info.n1_) { 
        info.n2_ = info.n1_;
        info.n1_ = a;
        
      }
      return Edge(edge_id, graph_);

    };

    /** Increment */
    IncidentIterator &operator++() {
      index_ ++;
      return *this;
    };
 
    /** Equality */
    bool operator==(const IncidentIterator& IncidentIter) const {
      bool same_graph = (graph_ == IncidentIter.graph_);
      bool same_node = (node_ == IncidentIter.node_); 
      bool same_index = (index_ == IncidentIter.index_); 
      return (same_graph && same_node && same_index); 
    };

  private:
    friend class Graph;
    size_type index_;
    Node* node_;
    Graph* graph_;

    //Private constructor that can be accessed by the Node class.
    IncidentIterator(size_type index, const Node* node, const Graph* graph) 
        : index_(index), node_(const_cast<Node*>(node)), graph_(const_cast<Graph*>(graph)) {
    }

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  public:
    // These type definitions let us use STL's iterator_traits.
    using value_type = Edge;                           // Element type
    using pointer = Edge *;                            // Pointers to elements
    using reference = Edge &;                          // Reference to elements
    using difference_type = std::ptrdiff_t;            // Signed difference
    using iterator_category = std::input_iterator_tag; // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {}

    /** Dereference */
    Edge operator*() const { return Edge(index_, graph_); };

    /** Increment */
    EdgeIterator &operator++() {
      index_ ++;
      return *this;
    };
    
    /** Equality */
    bool operator==(const EdgeIterator& EdgeIter) const {
      return (index_ == EdgeIter.index_) && (graph_ == EdgeIter.graph_);
    };

  private:
    friend class Graph;
    size_type index_; 
    Graph* graph_;
    EdgeIterator(size_type index, const Graph* graph) 
      : index_(index), graph_(const_cast<Graph*>(graph)) {
    }
  };

  /** begin() function for edge iterator. */
  edge_iterator edge_begin() const {
    return EdgeIterator(0, this); 
  };

  /** end() function for edge iterator. */
  edge_iterator edge_end() const {
    return EdgeIterator(active_edges_.size(), this); 
  };

  /** remove the edge associated to the edge iterator. */  
  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return EdgeIterator(0, this);
  };

private:
  std::vector<nodeinfo> nodes_;
  std::vector<size_type> active_nodes_;
  std::vector<edgeinfo> edges_;
  std::vector<size_type> active_edges_;
  std::unordered_map<std::string, size_type> edgeSet_;

  struct nodeinfo {
    Point p_;
    node_value_type v_;
    size_type idx_;
    std::vector<size_type> incidents_;
    bool active_;
    Point initial_p_;
  };

  struct edgeinfo {
    size_type n1_;
    size_type n2_;
    edge_value_type v_;
    size_type idx_;
  };
};

#endif // CME212_GRAPH_HPP
