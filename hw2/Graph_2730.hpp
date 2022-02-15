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
template <typename V, typename E>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value */
  typedef V node_value_type;
  typedef E edge_value_type;

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
    // HW0: YOUR CODE HERE  
    num_nodes_ = 0;
    num_edges_ = 0;
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
    Node() {
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      id_ = 0;
    }

    /** Return this node's position. */
    Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[id_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[id_].ind;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return graph_->nodes_[id_].value;
    }

    const node_value_type& value() const {
      return graph_->nodes_[id_].value;
    }

    /** Return the number of incident edges. */
    size_type degree() const {
      return graph_->degrees[id_];
    }

    /** Start of the incident iterator */
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, id_, 0);
    }

    /** End of incident iterator */
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, id_, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph_ == n.graph_) && (id_ == n.id_);
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
      if (graph_ < n.graph_)
        return true;
      else if (graph_ > n.graph_)
        return false;
      else if (id_ < n.id_)
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type id_;
    Node(const Graph* graph, size_type id) 
      : graph_(const_cast<Graph*>(graph)), id_(id) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
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
  Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodeid_map.insert(std::make_pair(num_nodes_, nodes_.size()));
    nodes_.push_back(internal_node(position, num_nodes_, val));
    ++num_nodes_;
    adj_list.resize(num_nodes_);
    degrees.push_back(0);
    return Node(this, num_nodes_-1);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert((0 <= i) && (i < size()));
    return Node(this, nodeid_map.at(i));        
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
      graph_ = nullptr;
      node1_id_ = node2_id_ = 0;
      id_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node1_id_);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graph_, node2_id_);      
    }

    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    edge_value_type& value() {
      return graph_->edges_[id_].value;
    }

    const edge_value_type& value() const {
      return graph_->edges_[id_].value;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (graph_ == e.graph_ && id_ == e.id_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ < e.graph_)
        return true;
      else if (graph_ > e.graph_)
        return false;
      else if (id_ < e.id_)
        return true;
      else
        return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_;
    size_type node1_id_;
    size_type node2_id_;
    size_type id_;
    Edge(const Graph* graph, size_type node1_id, size_type node2_id, size_type(id)) 
      : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_id_(node2_id), id_(id) {
    }
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
    size_type eid = edgeid_map.at(i);
    return Edge(this, edges_[eid].node1_id_, edges_[eid].node2_id_, eid);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graph_ == this && b.graph_ == this);
    if (a.id_ == b.id_) 
      return false;
    for (size_type i = 0; i < a.degree(); i++) {
      size_type eid = adj_list[a.id_][i];
      if (b.id_ == edges_[eid].node1_id_ || b.id_ == edges_[eid].node2_id_) {
        return true;
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
    // HW0: YOUR CODE HERE
    assert((a.graph_ == this) && (b.graph_  == this) && (a.index() != b.index()));
    for (size_type i = 0; i < a.degree(); i++) {
      size_type eid = adj_list[a.id_][i];
      if (b.id_ == edges_[eid].node1_id_ || b.id_ == edges_[eid].node2_id_) {
        return Edge(this, a.id_, b.id_, eid);    
      }   
    }
    if (a.degree() == adj_list[a.id_].size()) 
      adj_list[a.id_].push_back(edges_.size());
    else
      adj_list[a.id_][a.degree()] = edges_.size();
    if (b.degree() == adj_list[b.id_].size()) 
      adj_list[b.id_].push_back(edges_.size());
    else
      adj_list[b.id_][b.degree()] = edges_.size();

    degrees[a.id_] += 1;
    degrees[b.id_] += 1;
    assert(a.id_ < nodes_.size() && b.id_ < nodes_.size());
    edges_.push_back(internal_edge(a.id_, b.id_, num_edges_, edge_value_type()));
    edgeid_map.insert(std::make_pair(num_edges_, edges_.size()-1));
    num_edges_++;
    return Edge(this, a.id_, b.id_, edges_.size()-1);       
  }

std::vector<size_type> update_adj_list(std::vector<std::vector<size_type>>& adj_list, 
                                                    size_type id1, size_type id2) {
    size_type eid;
    for (size_type i = 0; i < degrees[id1]; i++) {
      eid = adj_list[id1][i];
      assert (edges_[eid].node1_id_ < nodes_.size() && edges_[eid].node2_id_ < nodes_.size());
      if (id2 == edges_[eid].node1_id_ || id2 == edges_[eid].node2_id_) {
        size_type tmp = adj_list[id1][i];
        adj_list[id1][i] = adj_list[id1][degrees[id1]-1];
        adj_list[id1][degrees[id1]-1] = tmp;
        std::vector<size_type> res {1, eid};
        return res;
      }
    }
    std::vector<size_type> res {0, 0};
    return res;
  }

  /**
   * @brief Invalidate the edge connecting node a and node b. Function update_adj_list is called twice. Suppose
   * the max degree is K, then the calls take O(K) operations since each call takes O(K).
   * 
   * @param[in] Node a
   * @param[in] Node b
   * @return size_type Return 1 if there is an valid edge between a and b and we removed it; return 0 otherwise.
   * 
   * @pre Let d1 be old degree of node a, and d2 be old degree of node b. Then the first d1 elemenets of 
   *      adj_list[a.id_] are ids of valid incident edges on a, and the first d2 elements of adj_list[b.id_] 
   *      are ids of valid incident edges on b.
   * @pre Let ind be the index of the edge in the graph, edgeid_map[ind] = edge_id
   * @pre degrees[a.id_] = d1 and degrees[b.id_] = d2.
   * 
   * @post new @a a.degree() = old @a a.degree() - 1, new @a b.degree() = old @a b.degree() - 1
   * @post new @a num_edges_ = old @a num_edges_ - 1
   * @post @a edgeid_map.size() = new @a num_edges_, for all key in edgeid_map, 0 <= key < new num_edges_
   * @post @a adj_list[a.id] has the same size, but only the first new @a a.degree() elements are valid,
   *       and the id of the removed edge is not valid
   * @post @a adj_list[b.id] has the same size, but only the first new @a b.degree() elements are valid.
   * @post @a degrees[a.id_] = new @a a.degree(), @a degrees[b.id_] = new @a b.degree()
   */
  size_type remove_edge(const Node& a, const Node& b) {
    std::vector<size_type> res = update_adj_list(adj_list, a.id_, b.id_);
    if (res[0] == 0)
      return 0;

    size_type eid = res[1];
    update_adj_list(adj_list, b.id_, a.id_);
    size_type ind = edges_[eid].ind;
    if (ind < num_edges_ - 1) {
      size_type new_eid = edgeid_map[num_edges_-1];
      edgeid_map[ind] = new_eid;
      edges_[new_eid].ind = ind;
    }
    edgeid_map.erase(num_edges_);
    degrees[a.id_] -= 1;
    degrees[b.id_] -= 1;
    --num_edges_;
    return 1;
  }

  /**
   * @brief Invalidate the edge e. We call remove_edge(const Node& a, const Node& b) defined above, so 
   * it takes O(K) operations as well.
   * @param[in] Edge e
   * @return size_type Return 1 if e is valid and succesfully removed; return 0 otherwise .
   * 
   * @pre Let ind be the index of e in the graph, edgeid_map[ind] = e.id_
   * 
   * @post new @a n1.degree() = old @a n1.degree() - 1, new @a n2.degree() = old @a n2.degree() - 1
   * @post new @a num_edges_ = old @a num_edges_ - 1
   * @post @a edgeid_map.size() = new @a num_edges_, for all key in edgeid_map, 0 <= key < new num_edges_
   * @post @a adj_list[n1.id] has the same size, but only the first new @a n1.degree() elements are valid,
   *       and the id of the removed edge is not valid
   * @post @a adj_list[n2.id] has the same size, but only the first new @a n2.degree() elements are valid.
   * @post @a degrees[n1.id_] = new @a n1.degree(), @a degrees[n2.id_] = new @a n2.degree()
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(Node(this, e.node1_id_), Node(this, e.node2_id_));
  }

  /**
   * @brief Invalidate the edge iterator e_it. We call remove_edge(const Edge& e) defined above, so 
   * it takes O(K) operations as well.
   * @param edge_iterator e_it
   * @return edge_iterator Return edge_begin() if e is valid and succesfully removed; return e_it otherwise.
   * 
   * @pre Let ind be the index of e in the graph, edgeid_map[ind] = e.id_
   * 
   * @post new @a n1.degree() = old @a n1.degree() - 1, new @a n2.degree() = old @a n2.degree() - 1
   * @post new @a num_edges_ = old @a num_edges_ - 1
   * @post @a edgeid_map.size() = new @a num_edges_, for all key in edgeid_map, 0 <= key < new num_edges_
   * @post @a adj_list[n1.id] has the same size, but only the first new @a n1.degree() elements are valid,
   *       and the id of the removed edge is not valid
   * @post @a adj_list[n2.id] has the same size, but only the first new @a n2.degree() elements are valid.
   * @post @a degrees[n1.id_] = new @a n1.degree(), @a degrees[n2.id_] = new @a n2.degree()
   * @post @a e_it = edge_begin() if edge is successfully invalidated
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end()) {
      remove_edge(*e_it);
      e_it = edge_begin();
    } 
    return e_it;
  }

  /**
   * @brief Invalidate node n. Call remove_edge() and remove all the incident edges. Take O(K) as well.
   * 
   * @param[in] Node n 
   * @return size_type Return 1 if node is valid and succesfully removed; return 0 otherwise.
   * 
   * @post n.degree() = 0
   * @post new num_nodes_ = old num_nodes_ - 1
   * @post new nodeid_map.size() = old nodid_map.size() - 1
   * @post adj_list[n.id_] is empty
   */
  size_type remove_node(const Node& n) {
    if (!(has_node(n)))
      return 0;
    auto it = n.edge_begin();
    while (it != n.edge_end()) {
      remove_edge(*it);
      if (n.degree() == 0)
        break;
    }
    size_type ind = nodes_[n.id_].ind;
    if (ind < num_nodes_ - 1) {
      size_type new_id = nodeid_map[num_nodes_-1];
      nodeid_map[ind] = new_id;
      nodes_[new_id].ind = ind;
    }
    --num_nodes_;
    return 1;
  }

    /**
   * @brief Invalidate node_iterator n_it. Call remove_node() and remove all the incident edges. 
   * Take O(K) as well.
   * 
   * @param node_iterator n_it 
   * @return size_type Return node_begin() if n_it is valid and succesfully removed; return n_it otherwise.
   * 
   * @post n.degree() = 0
   * @post new num_nodes_ = old num_nodes_ - 1
   * @post new nodeid_map.size() = old nodid_map.size() - 1
   * @post adj_list[n.id_] is empty
   * @post n_it = node_begin() if removed
   */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it != node_end()) {
      remove_node(*n_it);
      n_it = node_begin();
    } 
    return n_it;
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
    nodeid_map.clear();
    edges_.clear(); 
    edgeid_map.clear();
    adj_list.clear();
    degrees.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(graph_, graph_->nodeid_map[ind]);
    }
    
    NodeIterator& operator++() {
      ++ind;
      return *this;
    }
    
    bool operator==(const NodeIterator& node_iter) const {
      return (graph_ == node_iter.graph_) && (ind == node_iter.ind);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type ind;
    NodeIterator(Graph* g, size_type i) : graph_(g), ind(i) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(const_cast<Graph*>(this), 0);
  }

  node_iterator node_end() const {
    return NodeIterator(const_cast<Graph*>(this),num_nodes_);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    Edge operator*() const {
      size_type eid = graph_->adj_list[node_id_][inc_id_];
      size_type node2_id_ = 0;
      if (graph_->edges_[eid].node2_id_ == node_id_) 
        node2_id_ = graph_->edges_[eid].node1_id_;
      else
        node2_id_ = graph_->edges_[eid].node2_id_;
      return Edge(graph_, node_id_, node2_id_, eid);
    }

    IncidentIterator& operator++() {
      inc_id_++;
      return *this;
    }

    bool operator==(const IncidentIterator& incident_iter) const {
      return (graph_ == incident_iter.graph_) && (node_id_ == incident_iter.node_id_) \
                                              && (inc_id_ == incident_iter.inc_id_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_;
    size_type inc_id_;
    IncidentIterator(Graph* g, size_type id1, size_type id2)
     : graph_(g), node_id_(id1), inc_id_(id2) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
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
    Edge operator*() const {
      return *inc_iter;
    }
    
    EdgeIterator& operator++() {
      ++inc_iter;
      size_type node_ind_ = node_iter.ind; 
      size_type inc_id_ = inc_iter.inc_id_;

      while (1==1) {
        size_type flag = 0;
        if (inc_id_ >= (*node_iter).degree()) {
          flag = 1;
          ++node_iter;
          inc_iter = (*node_iter).edge_begin();
        } else {
          size_type max_ind = std::max((*inc_iter).node1().index(), (*inc_iter).node2().index());
          if (node_ind_ >= max_ind) {
            ++inc_iter;
            flag = 1;
          }
        }
        node_ind_ = node_iter.ind; 
        inc_id_ = inc_iter.inc_id_;
        if (flag == 0 || node_ind_ == graph_->num_nodes())
          break;
      }
      return *this;
    }
    
    bool operator==(const EdgeIterator& edge_iter) const {
      return (graph_ == edge_iter.graph_) && (node_iter == edge_iter.node_iter)\
                                          && (inc_iter == edge_iter.inc_iter);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    node_iterator node_iter;
    incident_iterator inc_iter;
    EdgeIterator(Graph* g, node_iterator node, incident_iterator inc)
     : graph_(g), node_iter(node), inc_iter(inc) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    node_iterator n = node_begin();
    return EdgeIterator(const_cast<Graph*>(this), n, (*n).edge_begin());
  }

  edge_iterator edge_end() const {
    node_iterator n = node_end();
    return EdgeIterator(const_cast<Graph*>(this), n, (*n).edge_begin());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  struct internal_node {
    Point point;
    size_type ind;
    node_value_type value;
    internal_node(const Point& point, size_type id, node_value_type val)
      :  point(point), ind(id), value(val) {}
  };

  struct internal_edge {
    size_type node1_id_;
    size_type node2_id_;
    size_type ind;
    edge_value_type value;
    internal_edge(size_type id1, size_type id2, size_type id, edge_value_type val)
      :  node1_id_(id1), node2_id_(id2), ind(id), value(val) {}
  };

  size_type num_nodes_;
  std::vector<internal_node> nodes_;
  std::map<size_type, size_type> nodeid_map; // node index in the graph mapped to node id_

  size_type num_edges_;
  std::vector<internal_edge> edges_;
  std::map<size_type, size_type> edgeid_map; // edge index in the graph mapped to edge id_

  std::vector<std::vector<size_type> > adj_list;
  std::vector<size_type> degrees;
};

#endif // CME212_GRAPH_HPP
