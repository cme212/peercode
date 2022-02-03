#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>

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

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
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

  using node_value_type = V; // added

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
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
  class Node : private totally_ordered<Node>{
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
     * 
     */
    
    // added for HW1
    node_value_type& value () {
      return graph_->nodes_[uid_].value;
    }

    // added for HW1
    const node_value_type& value () const {
      return graph_->nodes_[uid_].value;
    }

    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph_->nodes_[uid_].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(0 <= uid_ && uid_ < graph_->nodes_.size());
      return graph_->nodes_[uid_].uid;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    size_type degree() const {
      return graph_->edges_[uid_].size();
    }
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, uid_, 0);
    }
    incident_iterator edge_end() const {
      return IncidentIterator(graph_, uid_, this->degree()); // degree -1??
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return graph_==n.graph_ and uid_==n.uid_;

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
      return graph_==n.graph_ && uid_ < n.uid_;

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* graph_;
    size_type uid_; // uid_ belongs to node class, uid belongs to node_element

    /** Private Constructor */
    Node(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodes_.size();

  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  // /** Add a node to the graph, returning the added node.
  //  * @param[in] position The new node's position
  //  * @post new num_nodes() == old num_nodes() + 1
  //  * @post result_node.index() == old num_nodes()
  //  *
  //  * Complexity: O(1) amortized operations.
  //  */
  // Node add_node(const Point& position) {
  //   // HW0: YOUR CODE HERE
  //   size_type size_ = size();
  //   // Set the text and uid for the new element
  //   node_element new_node;
  //   new_node.point = position;
  //   new_node.uid = size_;
  //   // Delete the old elements and reassign its value
  //   nodes_.push_back(new_node);
  //   // Returns a SimpleElement that points to the new element
  //   return Node(this, new_node.uid);
  //   // (void) position;      // Quiet compiler warning
  //   // return Node();        // Invalid node
  // }

  // added for HWK1
  Node add_node(const Point& position, const node_value_type& node_val = node_value_type ()) { // what is this doing ??????
    size_type size_ = size();
    // Set the text and uid for the new element
    node_element new_node;
    new_node.point = position;
    new_node.uid = size_;
    new_node.value = node_val;
    // Delete the old elements and reassign its value
    nodes_.push_back(new_node);
    // add new vector to edges_
    std::vector<size_type> new_vec;
    edges_.push_back(new_vec);
    // Returns a SimpleElement that points to the new element
    return Node(this, new_node.uid);
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.graph_ == this and n.uid_ < size();

  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
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
  class Edge : private totally_ordered<Edge>{
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
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

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (e.graph_ == graph_) and ((e.node1_id_ == node1_id_ and e.node2_id_ == node2_id_) or (e.node1_id_ == node2_id_ and e.node2_id_ == node1_id_));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (node1_id_ == e.node1_id_) return node2_id_ < e.node2_id_;
      return node1_id_ < e.node1_id_;
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

    /** Private Constructor */
    Edge(const Graph* graph, size_type node1_id, size_type node2_id)
        : graph_(const_cast<Graph*>(graph)), node1_id_(node1_id), node2_id_(node2_id) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return id2edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return Edge(this, id2edges_.at(i)[0], id2edges_.at(i)[1]);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::vector<size_type> pair1;
    pair1.push_back(a.uid_);
    pair1.push_back(b.uid_);
    std::vector<size_type> pair2;
    pair2.push_back(b.uid_);
    pair2.push_back(a.uid_);
    return edges2id_.count(pair1) == 1 or edges2id_.count(pair2) == 1;

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
    std::vector<size_type> pair;
    pair.push_back(a.uid_);
    pair.push_back(b.uid_);
    std::vector<size_type> pair2;
    pair2.push_back(b.uid_);
    pair2.push_back(a.uid_);


    if(edges2id_.count(pair) == 0 and edges2id_.count(pair2) == 0){
      id2edges_[curr_edge_id_] = pair;
      edges2id_[pair] = curr_edge_id_;
      curr_edge_id_ += 1;

      edges_[a.uid_].push_back(b.uid_);
      edges_[b.uid_].push_back(a.uid_);
    }
    return Edge(this, a.uid_, b.uid_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges2id_.clear(); // maps ids of two nodes in an edge to unique id of edge
    id2edges_.clear(); // maps unique id of edge to ids of two nodes in an edge
    curr_edge_id_ = 0;
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

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return Node(graph_, graph_->nodes_[uid_].uid);
    }

    NodeIterator& operator++() {
      uid_ ++;
      return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
      return graph_ == node_iter.graph_ && uid_ == node_iter.uid_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type uid_; 

    // Constructor
    NodeIterator(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator(this, this->nodes_.size());
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
      return Edge(graph_, node1_, graph_->edges_[node1_][curr_idx_]);
    }
    IncidentIterator& operator++() {
      curr_idx_ ++;
      return *this;
    }
    bool operator==(const IncidentIterator& node_iter) const {
      return graph_ == node_iter.graph_ && curr_idx_ == node_iter.curr_idx_ && node1_ == node_iter.node1_;
    }


   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type node1_; 
    size_type curr_idx_;

    // Constructor
    IncidentIterator(const Graph* graph, size_type uid, size_type cur_idx)
        : graph_(const_cast<Graph*>(graph)), node1_(uid), curr_idx_(cur_idx) {         
    }

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
      Edge e = Edge(graph_, node1_, graph_->edges_[node1_][node2_idx_]);
      return e;
    }

    EdgeIterator& operator++() {
      if (node2_idx_ < graph_->edges_[node1_].size()-1) {
          node2_idx_ ++;
      } else {
        node1_++;
        node2_idx_ = 0;
      }
      while (node1_ > graph_->edges_[node1_][node2_idx_] && node1_ < graph_->edges_.size() && node2_idx_ < graph_->edges_[node1_].size()) {
        if (node2_idx_ < graph_->edges_[node1_].size()-1) {
          node2_idx_ ++;
        } else {
          node1_++;
          node2_idx_ = 0;
        }
      }
      return *this;
    }

    bool operator==(const EdgeIterator& edge_iter) const {
      return graph_ == edge_iter.graph_ && node1_ == edge_iter.node1_ && node2_idx_ == edge_iter.node2_idx_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type node1_; 
    size_type node2_idx_;

    // Constructor
    EdgeIterator(const Graph* graph, size_type node1, size_type node2_idx)
        : graph_(const_cast<Graph*>(graph)), node1_(node1), node2_idx_(node2_idx) {        
    }

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->edges_.size(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct node_element {
    Point point;
    size_type uid;      // The unique identifcation for an element
    node_value_type value;
  };

  std::vector<node_element> nodes_;
  std::map<std::vector<size_type>, size_type> edges2id_; // maps ids of two nodes in an edge to unique id of edge
  std::map<size_type, std::vector<size_type>> id2edges_; // maps unique id of edge to ids of two nodes in an edge
  size_type curr_edge_id_ = 0;
  std::vector<std::vector<size_type>> edges_;// each idx represent node's uid (node1), the vector at each idx represents all nodes ids that are adjacent to node1

};

#endif // CME212_GRAPH_HPP
