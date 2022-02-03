#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <set>

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
  struct internal_edge;
  struct internal_node;
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
  Graph(): node_size_(0), next_node_id_(0),
    edge_size_(0), next_edge_id_(0){
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
      index_ = 0;
      parentG_ = nullptr;
    }
      

    
    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      if (parentG_->node_.count(index_) > 0 ){
          return (parentG_->node_[index_]).position;
      }
      assert(false);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
//     node_value_type& value();
//     const node_value_type& value() const;
//     size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
      node_value_type& value(){
          if (parentG_->node_.count(index_) > 0 ){
              return (parentG_->node_[index_]).value;
          }
          assert(false);
      }
      const node_value_type& value() const{
          if (parentG_->node_.count(index_) > 0 ){
              return (parentG_->node_[index_]).value;
          }
          assert(false);
      }
      
      size_type degree() const{
          return (parentG_->node_edges_[index_]).size();
      }
      incident_iterator edge_begin() const{
          return IncidentIterator(parentG_, index_, 0);
      }
      incident_iterator edge_end() const{
          size_type deg = (parentG_->node_edges_[index_]).size();
          return IncidentIterator(parentG_, index_, deg);
      }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if (parentG_ == n.parentG_ and index_ == n.index_){
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
      // HW0: YOUR CODE HERE
        if ( parentG_ < n.parentG_){
            return true;
        }else if(parentG_ == n.parentG_){
            if(index_ < n.index_){
                return true;
            }else{
                return false;
            }
        }else{
            return false;
        }
    }

   private:
    Graph* parentG_;
    size_type index_;
    
    /**Private Constructor*/
    Node(const Graph* graph, size_type index):
      parentG_(const_cast<Graph*>(graph)), index_(index){}
    

    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
      
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_size_;
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
      //create a new node
      internal_node new_node;
      new_node.position = position;
      new_node.index = next_node_id_;
      new_node.value = val;
      //update the graph attributes
      node_[next_node_id_] = new_node;
      node_size_ += 1;
      next_node_id_ += 1;
      return Node(this, next_node_id_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
      return node_.count(n.index_) > 0;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
      assert(i < size());
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
      index_ = 0;
      parentG_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (parentG_->edge_.count(index_) > 0){
        return (parentG_->edge_[index_]).node1;
      }
      assert(false);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (parentG_->edge_.count(index_) > 0){
        return (parentG_->edge_[index_]).node2;
      }
      assert(false);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //since each edge is unique(with unique pair of nodes and index)
      //edges connecting the same pair of nodes must have the same index
      return index_ == e.index_;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //edges are in the same graph
      return index_ < e.index_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    Graph* parentG_;
    size_type index_;
    Edge(const Graph* graph, size_type index)
      : parentG_(const_cast<Graph*>(graph)), index_(index){}
      
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    //check a, b are valid nodes of this graph
    assert(has_node(a));
    assert(has_node(b));
    //make a pair of node index
    std::set<size_type> pair = {a.index(), b.index()};
    //check if the pair already exists in the map key
    return edge_id_.count(pair) > 0;
    
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
      assert(!(a == b));
      size_type a_idx = a.index();
      size_type b_idx = b.index();
      //if edge already exists
      std::set<size_type> pair = {a_idx, b_idx};
      if (has_edge(a, b)){
          size_type cur_edge = edge_id_[pair];
          edge_[cur_edge].node1 = a;
          edge_[cur_edge].node2 = b;
          return Edge(this, cur_edge);
      }
      //create a new edge
      internal_edge new_edge;
      new_edge.node1 = a;
      new_edge.node2 = b;
      new_edge.index = next_edge_id_;
      //update graph attributes
      node_edges_[a_idx].push_back(next_edge_id_);
      node_edges_[b_idx].push_back(next_edge_id_);
      edge_[next_edge_id_] = new_edge;
      edge_id_[pair] = next_edge_id_;
      edge_size_ += 1;
      next_edge_id_ += 1;
      return Edge(this, next_edge_id_-1);
      
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_.clear();
    edge_.clear();
    node_edges_.clear();
    edge_id_.clear();
    node_size_ = 0;
    next_node_id_ = 0;
    edge_size_ = 0;
    next_edge_id_ = 0;
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
          gp_ = nullptr;
          node_id_ = 0;
    }
    

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{
        return gp_->node(node_id_);
    }
    NodeIterator& operator++(){
      node_id_ += 1;
      return *this;
    }
      
      bool operator==(const NodeIterator& nit) const{
          return gp_ == nit.gp_ && node_id_ == nit.node_id_;
      }

   private:
    friend class Graph;
    Graph* gp_;
    size_type node_id_;
    
    NodeIterator(const Graph* graph, size_type node_id): gp_(const_cast<Graph*>(graph)), node_id_(node_id) {
    }
      
      
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
    NodeIterator node_begin() const {
        return NodeIterator(this, 0);
    }
    NodeIterator node_end() const{
        return NodeIterator(this, node_size_);
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
        gp_ = nullptr;
        nidx_ = 0;
        ie_idx_ = 0;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
      Edge operator*() const{
          Edge this_eg = gp_->edge(gp_->node_edges_[nidx_][ie_idx_]);
          if (this_eg.node1().index() == nidx_){
              return this_eg;
          }else{
              //change edge orientation
              Node n2 = gp_->node(this_eg.node1().index());
              Node n1 = gp_->node(nidx_);
              return gp_->add_edge(n1, n2);
          }
      }
      IncidentIterator& operator++(){
          ie_idx_ += 1;
          return *this;
      }
      bool operator==(const IncidentIterator& ie_it) const {
          return gp_ == ie_it.gp_ && nidx_ == ie_it.nidx_ && ie_idx_ == ie_it.ie_idx_;
      }
   private:
    friend class Graph;
    friend class Node;
    Graph* gp_;
    size_type nidx_;
    size_type ie_idx_;
      
    IncidentIterator(const Graph* graph, size_type nidx, size_type ie_idx):gp_(const_cast<Graph*>(graph)), nidx_(nidx), ie_idx_(ie_idx){}
    // HW1 #3: YOUR CODE HERE
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
        gp_ = nullptr;
        eidx_ = 0;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
      Edge operator*() const {
          return gp_->edge(eidx_);
      }
      EdgeIterator& operator++(){
          eidx_ += 1;
          return *this;
      }
      bool operator==(const EdgeIterator& e_it) const{
          return gp_ == e_it.gp_ && eidx_ == e_it.eidx_;
      }
   private:
    friend class Graph;
    Graph* gp_;
    size_type eidx_;
      EdgeIterator(const Graph* graph, size_type eidx): gp_(const_cast<Graph*>(graph)), eidx_(eidx) {
      }
    
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
    edge_iterator edge_begin() const{
        return EdgeIterator(this, 0);
    }
    edge_iterator edge_end() const{
        return EdgeIterator(this, edge_size_);
    }
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
    struct internal_node {
        Point position;
        size_type index;
        node_value_type value;
    };
    struct internal_edge {
        Node node1;
        Node node2;
        size_type index;
    };
    std::unordered_map<size_type, internal_node> node_;
    std::unordered_map<size_type, internal_edge> edge_;
    std::map<std::set<size_type>, size_type> edge_id_;
    std::unordered_map<size_type, std::vector<size_type>> node_edges_;
    std::unordered_set<size_type> reachable;
    size_type node_size_;
    size_type next_node_id_;
    size_type edge_size_;
    size_type next_edge_id_;
    
};

#endif // CME212_GRAPH_HPP
