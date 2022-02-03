#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <iostream>
#include <map>
#include <set>
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
template <typename V = int>
class Graph {
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;
  
  using node_value_type = V;

  struct node_info {
    Point p;
    node_value_type value;
    size_type degree;
    std::vector<size_type> edges;
  };

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


  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
      // Uses unitialized pointer and integer, invalid unless called by
      // graph.add_node().
    }

    /** Return this node's position. 
     */
    const Point& position() const {
      return graph_->node_info_.at(uid_).p;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // Return private property UID
      return uid_;
    }

    // Return actual storage of value for that node
    node_value_type& value() {
      return graph_->node_info_[uid_].value;
    }

    const node_value_type& value() const {
      return graph_->node_info_.at(uid_).value;
    }

    size_type degree() const {
      return graph_->node_info_.at(uid_).degree;
    }
 
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (this->uid_ != n.uid_) {
        return false;
      }
      return true;
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
      if (this->uid_ > n.uid_) {
        return false;
      }
      if (this->uid_ == n.uid_) {
        return false;
      }
      return true;
    }
    
    incident_iterator edge_begin() const {
      incident_iterator result;
      result.node_id_ = uid_;
      result.graph_ = graph_;
      result.edge_ = graph_->node_info_.at(uid_).edges.begin();
      return result;
    }

    incident_iterator edge_end() const {
      incident_iterator result;
      result.node_id_ = uid_;
      result.graph_ = graph_;
      result.edge_ = graph_->node_info_.at(uid_).edges.end();
      return result;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Address of graph that contains this node.
    Graph* graph_;
    // Index of node in graph.
    size_type uid_;
    // Validity indicator
    bool valid_ = false;
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    if (node_storage_.empty()) {
      return 0;
    }
    else {
      return node_storage_.size();
    }
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
                const node_value_type& v = node_value_type()) {
    Node newNode = Node();
    node_info newInfo = node_info();
    Point point(position.x, position.y, position.z);
    newNode.graph_ = this;
    newNode.uid_ = this->size();
    newNode.valid_ = true;
    newInfo.p = point;
    newInfo.value = v;
    newInfo.degree = 0;
    node_storage_[newNode.uid_] = newNode;
    node_info_[newNode.uid_] = newInfo;
    return newNode;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // Access map and check index
    if (!n.valid_) {
      return false;
    }
    if (!(n.uid_ < size())) {
      return false;
    }
    if (n.graph_ != this) {
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
  const Node node(size_type i) const {
    return node_storage_.at(i);
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
      // Uses unitialized pointer and integer, invalid unless called by
      // graph.add_edge().
    }

    /** Return a node of this Edge */
    Node node1() const {
      return graph_->node_storage_.at(node1_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return graph_->node_storage_.at(node2_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (this->graph_ != e.graph_) {
        return false;
      }
      if (this->uid_ != e.uid_) {
        return false;
      }
      return true;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (this->graph_ != e.graph_) {
        return false;
      }
      if (this->uid_ > e.uid_) {
        return false;
      }
      return true;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Pointer to graph object
    Graph* graph_;
    // Unique identifier
    size_type uid_;
    // Node indentifiers
    size_type node1_;
    size_type node2_;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_storage_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edge_storage_.at(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::set<size_type> node_ids {a.index(), b.index()};
    return edges_.count(node_ids);
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
    std::set<size_type> nodes {a.uid_, b.uid_};
    if (this->has_edge(a, b)) {
      Edge newEdge = Edge();
      newEdge.graph_ = this;
      newEdge.node1_ = a.uid_;
      newEdge.node2_ = b.uid_;
      return newEdge;
      // Waiting for ed response
      //return edge_storage_[edges_[nodes]];
    }
    Edge newEdge = Edge();
    newEdge.uid_ = edge_storage_.size();
    newEdge.node1_ = a.uid_;
    newEdge.node2_ = b.uid_;
    newEdge.graph_ = this;
    edge_storage_[newEdge.uid_] = newEdge;
    edges_[nodes] = newEdge.uid_;
    node_info_[a.uid_].degree += 1;
    node_info_[b.uid_].degree += 1;
    node_info_[a.uid_].edges.push_back(newEdge.uid_);
    node_info_[b.uid_].edges.push_back(newEdge.uid_);
    return newEdge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    node_storage_.clear();
    node_info_.clear();
    edge_storage_.clear();
    edges_.clear();
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
    NodeIterator() {
    }

   // Node operator*()
   Node operator*() const {
     return graph_->node_storage_.at(uid_);
   }

   // NodeIterator& operator++()
   NodeIterator& operator++() {
     uid_++;
     return *this;
   }

   // bool operator==(const NodeIterator&) const
   bool operator==(const NodeIterator& nI) const {
     if (this->graph_ != nI.graph_) {
       return false;
     }
     if (this->uid_ != nI.uid_) {
       return false;
     }
     return true;
   }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type uid_;
    const Graph* graph_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    node_iterator result;
    result.graph_ = this;
    result.uid_ = 0;
    return result;
  }

  node_iterator node_end() const {
    node_iterator result;
    result.graph_ = this;
    result.uid_ = this->size();
    return result;
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

    // Edge operator*() const
    Edge operator*() const {
      Edge newEdge;
      newEdge.graph_ = this->graph_;
      newEdge.node1_ = node_id_;
      if (newEdge.node1_ == graph_->edge_storage_[*edge_].node1().index()) {
        newEdge.node2_ = graph_->edge_storage_[*edge_].node2().index();
      }
      else {
        newEdge.node2_ = graph_->edge_storage_[*edge_].node1().index();
      }
      return newEdge;
    }

    // IncidentIterator& operator++()
    IncidentIterator& operator++() {
      edge_++;
      return *this;
    }

    // bool operator==(const IncidentIterator&) const
    bool operator==(const IncidentIterator& incI) const {
      if (this->graph_ != incI.graph_) {
        return false;
      }
      if (this->node_id_ != incI.node_id_) {
        return false;
      }
      if (this->edge_ != incI.edge_) {
        return false;
      }
      return true;
    }

   private:
    friend class Node;
    friend class Graph;
    Graph* graph_;
    size_type node_id_;
    std::vector<size_type>::iterator edge_;
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

    // Edge operator*() const
    Edge operator*() const {
      return graph_->edge_storage_.at(uid_);
    }
    // EdgeIterator& operator++()
    EdgeIterator& operator++() {
      uid_++;
      return *this;
    }
    // bool operator==(const EdgeIterator&) const
    bool operator==(const EdgeIterator& eI) const {
      if (this->graph_ != eI.graph_) {
        return false;
      }
      if (this->uid_ != eI.uid_) {
        return false;
      }
      return true;
    }

    void print_vals() {
      std::cout << "Edge ID: " << uid_ << std::endl;
    }

   private:
    friend class Graph;
    // Current edge id
    size_type uid_;
    // Graph object pointer
    const Graph* graph_;
  };

  // edge_iterator edge_begin() const
  edge_iterator edge_begin() const {
    edge_iterator result;
    result.graph_ = this;
    result.uid_ = 0;
    return result;
  }
  // edge_iterator edge_end() const 
  edge_iterator edge_end() const {
    edge_iterator result;
    result.graph_ = this;
    result.uid_ = this->num_edges();
    return result;
  }

 private:
  
  std::map<size_type, Node> node_storage_;
  std::map<size_type, node_info> node_info_;
  std::map<size_type, Edge> edge_storage_;
  std::map<std::set<size_type>, size_type> edges_;
};

#endif // CME212_GRAPH_HPP
