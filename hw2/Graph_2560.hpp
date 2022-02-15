#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = double>
//template<typename V>
class Graph {
  //typedef V node_value_type;
  typedef E edge_value_type;
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct neighbor;
  struct internal_node;
  struct internal_edge;

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
  /** Declaration of Node Value Type */
  using node_value_type = V;
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
  Graph() 
    : nodes_(), nodes_i2u_(), edges_(), edges_i2u_(), all_edges_(), edges_ordered_(),
      nodes_size_(0), next_node_uid_(0), edges_size_(0), next_edge_uid_(0) {
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
      node_uid_ = 0;
    }

    /** Return this node's position. */
    Point& position() {
      return fetch().point_;
    }

    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().point_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      //return fetch().uid_;
      return fetch().idx_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return fetch().val_;
    }

    const node_value_type& value() const {
      return fetch().val_;
    }

    size_type degree() const {
      return fetch().neighbors_.size();
    }

    // std::vector<neighbor> neighbors() const {
    //   return fetch().neighbors_;
    // }

    // incident_iterator edge_begin() const;
    IncidentIterator edge_begin() {
      return IncidentIterator(graph_, 0, this);
    }
    // incident_iterator edge_end() const;
    IncidentIterator edge_end() {
      return IncidentIterator(graph_, degree(), this);
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (graph_ == n.graph_ && node_uid_ == n.node_uid_) {
        return true;
      } else {
        return false;
      }
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
      if (graph_ < n.graph_) {
        return true;
      } else if (graph_ == n.graph_) {
        if (node_uid_ < n.node_uid_) {
          return true;
        } else {
          return false;
        }
      } else {
        return false;
        
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Pointer back to the graph container
    Graph* graph_;
    // This node's unique identification number
    size_type node_uid_;
    /** private constructor */
    Node(const Graph* graph, size_type uid) 
      : graph_(const_cast<Graph*>(graph)), node_uid_(uid) {
      }
  
    
    /** Helper method to return the appropriate element. */
    internal_node& fetch() const {
      //size_type node_uid = graph_->nodes_i2u_[node_idx_];
      if (graph_->nodes_.find(node_uid_) == graph_->nodes_.end()) {
        assert(false);
      } else {
        return graph_->nodes_.at(node_uid_);
      }
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    //assert(nodes_.size() == nodes_size_);
    return nodes_i2u_.size();
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
    nodes_[next_node_uid_] = internal_node();
    nodes_[next_node_uid_].point_ = position;
    size_type next_idx = num_nodes();
    nodes_[next_node_uid_].idx_ = next_idx;
    nodes_[next_node_uid_].val_ = val;
    nodes_i2u_.push_back(next_node_uid_);
    ++nodes_size_;
    ++next_node_uid_;
    return Node(this, next_node_uid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if (this == n.graph_ && n.index() < num_nodes()) {
      return true;
    } else {
      return false;
    }
    // if (nodes_.find(n.node_uid_) == nodes_.end()) {
    //   return false;
    // } else {
    //   return true;
    // }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < num_nodes()) {
      size_type uid = nodes_i2u_[i];
      return Node(this, uid);
    } else {
      assert(false);
    }
    
    // if (nodes_.find(i) == nodes_.end()) {
    //   assert(false);
    // } else {
    //   return Node(this, i);
    // }
  }


  size_type remove_node(const Node& cn) {
    if (has_node(cn)) {
      // remove all incident edges
      Node& n = const_cast<Node&>(cn);

      while (n.degree() != 0) {
        auto it = n.edge_begin();
        remove_edge(*it);
      }

      size_type cur_idx = n.index();
      size_type cur_uid = nodes_i2u_[cur_idx];

      if (cur_idx == nodes_i2u_.size()-1) {
        nodes_i2u_.pop_back();
      } else {
        size_type back_uid = nodes_i2u_.back();
        nodes_i2u_.back() = cur_uid;
        nodes_i2u_.at(cur_idx) = back_uid;
        // update the index for the element being swapped
        nodes_[back_uid].idx_ = cur_idx;
        nodes_[cur_uid].idx_ = nodes_i2u_.size()-1;
        // remove the element from active nodes
        nodes_i2u_.pop_back();
      }
      return 1;
    } else {
      return 0;
    }
  }

  node_iterator remove_node(node_iterator n_it) {
    remove_node(*n_it);
    return node_begin();
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
    Edge() {
      // HW0: YOUR CODE HERE
      graph_ = nullptr;
      edge_uid_ = 0;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return fetch().node1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return fetch().node2_;
    }

    edge_value_type& value() {
      return fetch().val_;
    }

    const edge_value_type& value() const {
      return fetch().val_;
    }

    double length() const {
      Point p1 = fetch().node1_.position();
      Point p2 = fetch().node2_.position();
      double length = norm_2((p1 - p2));
      return length;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      internal_edge this_edge = fetch();
      internal_edge other_edge = e.fetch();
      if (graph_ == e.graph_
          && this_edge.node1_.index() == other_edge.node1_.index() 
          && this_edge.node2_.index() == other_edge.node2_.index()) {
        return true;
      } else if (graph_ == e.graph_
                 && this_edge.node1_.index() == other_edge.node2_.index() 
                 && this_edge.node2_.index() == other_edge.node1_.index()) {
        return true;
      } else {
        return false;
      }
    }

    bool operator!=(const Edge& e) const {
      return !(*this == e);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (graph_ < e.graph_) {
        return true;
      } else if (graph_ == e.graph_) {
        if (edge_uid_ < e.edge_uid_) {
          return true;
        } else {
          return false;
        }
      } else {
        return false;
        
      }

    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    Graph* graph_;
    // This edge's user-facing index
    size_type edge_uid_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type uid)
      : graph_(const_cast<Graph*>(graph)), edge_uid_(uid) {
    }

    /** Helper method to return the appropriate element. */
    internal_edge& fetch() const {
      //size_type edge_uid = graph_->edges_i2u_[edge_idx_];
      if (graph_->edges_.find(edge_uid_) == graph_->edges_.end()) {
        assert(false);
      } else {
        return graph_->edges_[edge_uid_];
      }
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    //assert(edges_.size() == edges_size_);
    //assert(all_edges_.size() == edges_size_);
    return edges_i2u_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    size_type uid = edges_i2u_[i];
    if (edges_.find(uid) == edges_.end()) {
       assert(false);
     } else {
       return Edge(this, uid);
     }
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    size_type a_uid = nodes_i2u_.at(a.index());
    size_type b_uid = nodes_i2u_.at(b.index());
    std::set<size_type> cur = {a_uid, b_uid};
    if (all_edges_.find(cur) == all_edges_.end()) {
      return false;
    } else {
      return true;
    }
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
    //Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    //std::set<size_type> cur = {a.index(), b.index()};

    size_type a_uid = nodes_i2u_.at(a.index());
    size_type b_uid = nodes_i2u_.at(b.index());
    std::set<size_type> cur = {a_uid, b_uid};
    
    if (has_edge(a, b)) {
      size_type edge_uid = all_edges_[cur];
      // update node 1 and node 2
      edges_[edge_uid].node1_ = a;
      edges_[edge_uid].node2_ = b;
      //size_type edge_idx = edges_[edge_uid].idx_;
      return Edge(this, edge_uid);
    } else {
      edges_[next_edge_uid_] = internal_edge();
      size_type next_idx = num_edges();
      edges_[next_edge_uid_].idx_ = next_idx;
      edges_[next_edge_uid_].node1_ = a;
      edges_[next_edge_uid_].node2_ = b;
      //edges_[next_edge_uid_].length_ = norm_2((a.position() - b.position()));
      edges_[next_edge_uid_].val_ = val;
      all_edges_[cur] = next_edge_uid_;
      edges_ordered_.push_back(next_edge_uid_);
      // nodes_[a.node_uid_].neighbors_.push_back(neighbor(b.node_uid_, next_edge_uid_));
      // nodes_[b.node_uid_].neighbors_.push_back(neighbor(a.node_uid_, next_edge_uid_));
      nodes_[a_uid].neighbors_.push_back(neighbor(b, Edge(this, next_edge_uid_)));
      nodes_[b_uid].neighbors_.push_back(neighbor(a, Edge(this, next_edge_uid_)));
      edges_i2u_.push_back(next_edge_uid_);
      assert(edges_i2u_.size() == next_idx + 1);
      ++edges_size_;
      ++next_edge_uid_;
      return Edge(this, next_edge_uid_-1);
    }
    
  }

  size_type remove_edge(const Node& a, const Node& b) {
    size_type a_uid = nodes_i2u_.at(a.index());
    size_type b_uid = nodes_i2u_.at(b.index());
    std::set<size_type> cur = {a_uid, b_uid};
    if (has_edge(a, b)) {
      size_type edge_uid = all_edges_[cur];
      internal_edge e = edges_[edge_uid];
      size_type cur_idx = e.idx_;
      assert(edge_uid == edges_i2u_.at(cur_idx));
      if (cur_idx == edges_i2u_.size()-1) {
        edges_i2u_.pop_back();
      } else {
        size_type back_uid = edges_i2u_.back();
        edges_i2u_.back() = edge_uid;
        edges_i2u_.at(e.idx_) = back_uid;
        // update the index for the element being swapped
        edges_[back_uid].idx_ = cur_idx;
        edges_[edge_uid].idx_ = edges_i2u_.size()-1;
        // remove from active edges
        edges_i2u_.pop_back();
      }
      // remove edge from all_edges
      all_edges_.erase(cur);

      // remove b from a's neighbors
      //a_node.neighbors.erase(a_node.neighbors.find())
      internal_node& a_node = nodes_[a_uid];
      for (auto it = a_node.neighbors_.begin(); it != a_node.neighbors_.end();) {
        Node cur_nb = (*it).nb_;
        if (cur_nb == b) {
          it = a_node.neighbors_.erase(it);
          break;
        } else {
          ++it;
        }
      }
      // remove a from b's neighbors
      internal_node& b_node = nodes_[b_uid];
      for (auto it = b_node.neighbors_.begin(); it != b_node.neighbors_.end();) {
        Node cur_nb = (*it).nb_;
        if (cur_nb == a) {
          it = b_node.neighbors_.erase(it);
          break;
        } else {
          ++it;
        }
      }
      return 1;

    } else {
      return 0;
    }

  }

  size_type remove_edge(const Edge& e) {
    size_type res = remove_edge(e.node1(), e.node2());
    return res;
  }

  edge_iterator remove_edge(edge_iterator e_it) {
    remove_edge(*e_it);
    return edge_begin();
  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    nodes_i2u_.clear();
    edges_.clear();
    edges_i2u_.clear();
    all_edges_.clear();
    nodes_size_ = 0;
    next_node_uid_ = 0;
    edges_size_ = 0;
    next_edge_uid_ = 0;
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }
    

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    Node operator*() const {
      return graph_->node(idx_);
    }

    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
      if (graph_ == node_iter.graph_ && idx_ == node_iter.idx_) {
        return true;
      } else {
        return false;
      }
    }

    bool operator!=(const NodeIterator& node_iter) const {
      return !(*this == node_iter);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    const Graph* graph_;
    size_type idx_;

    /** Private Constructor */
    NodeIterator(const Graph* graph, size_type idx) 
      : graph_{graph}, idx_{idx} {}

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const;
  // node_iterator node_end() const;
  NodeIterator node_begin() const {
    return NodeIterator(this, 0);
  }

  NodeIterator node_end() const {
    size_type n = num_nodes();
    return NodeIterator(this, n);
  }

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
    //Edge operator*() const
    Edge operator*() const {
      size_type node_uid = graph_->nodes_i2u_[(*node_).index()];
      //size_type edge_uid = graph_->nodes_.at(node_uid_).neighbors_[idx_].edge_uid_;
      //Edge e = Edge(graph_, edge_uid);
      Edge e = graph_->nodes_.at(node_uid).neighbors_[idx_].edge_;
      if (!(e.node1() == *node_)) {
        (*graph_).add_edge(e.node2(), e.node1());
      }
      //return Edge(graph_, edge_uid);
      return e;
    }

    // IncidentIterator& operator++()
    IncidentIterator& operator++() {
      ++idx_;
      return *this;
    }

    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& inc_iter) const {
      if (graph_ == inc_iter.graph_ && idx_ == inc_iter.idx_ && node_ == inc_iter.node_) {
        return true;
      } else {
        return false;
      }
    }

    bool operator!=(const IncidentIterator& inc_iter) const {
      return !(*this == inc_iter);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type idx_;
    Node* node_;
    IncidentIterator(Graph* graph, size_type idx, Node* node) 
      : graph_{graph}, idx_{idx}, node_{node} {}

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
    //EdgeIterator() {
    //}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const {
      //size_type edge_uid = graph_->edges_ordered_.at(idx_);
      size_type edge_uid = graph_->edges_i2u_.at(idx_);
      return Edge(graph_, edge_uid);
    }
    
    EdgeIterator& operator++() {
      ++idx_;
      return *this;
    }
    
    bool operator==(const EdgeIterator& edge_iter_) const {
      if (graph_ == edge_iter_.graph_ && idx_ == edge_iter_.idx_) {
        return true;
      } else {
        return false;
      }
    }

    bool operator!=(const EdgeIterator& edge_iter_) const {
      return !(*this == edge_iter_);
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_;
    size_type idx_;

    /** Private Constructor */
    EdgeIterator(const Graph* graph, size_type idx) 
      : graph_{graph}, idx_{idx} {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, 0);
  };
  
  EdgeIterator edge_end() const {
    return EdgeIterator(this, edges_i2u_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct neighbor{
    //size_type neighbor_uid_;
    //size_type edge_uid_;
    Node nb_;
    Edge edge_;
    // neighbor(size_type neighbor_uid, size_type edge_uid) 
    //   : neighbor_uid_(neighbor_uid), edge_uid_(edge_uid) {}
    neighbor(Node nb, Edge edge): nb_(nb), edge_(edge) {}
  };

  struct internal_node {
    Point point_; // The point held by a node
    size_type idx_; 
    std::vector<neighbor> neighbors_;
    node_value_type val_;
  };

  struct internal_edge {
    size_type idx_;
    Node node1_;
    Node node2_;
    edge_value_type val_;
  };

  std::unordered_map<size_type, internal_node> nodes_;
  std::vector<size_type> nodes_i2u_;
  std::unordered_map<size_type, internal_edge> edges_;
  std::vector<size_type> edges_i2u_;
  std::map<std::set<size_type>, size_type> all_edges_;
  std::vector<size_type> edges_ordered_;
  size_type nodes_size_ = 0;
  size_type next_node_uid_ = 0;
  size_type edges_size_ = 0;
  size_type next_edge_uid_ = 0;
};

#endif // CME212_GRAPH_HPP
