#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
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
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of node value. */
  // using node_value_type = V;
  typedef V node_value_type;
  /** Type of edge value. */
  // using edge_value_type = E;
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
  Graph() 
  // HW0: YOUR CODE HERE
  : maxnode_(0), maxedge_(0), nodeset_size_(0), edgeset_size_(0){}
    

  /** Default Destructor */
  ~Graph() {
    // // delete all internal nodes
    // for (auto it = nodes_.begin(); it != nodes_.end(); it++) {
    //   delete it->second;
    // }
    // // delete all internal edges
    // for (auto it = edges_.begin(); it != edges_.end(); it++) {
    //   delete it->second;
    // }
  };

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
      graphptr_ = nullptr;
      idx_ = size_type(-1);
    }
    
    /** Return this node's position. Non-const case */
    Point& position() {
      // HW0: YOUR CODE HERE
      return fetch().pos_;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().pos_;
    }

    /** Return this node's index, a number in the range [0, nodeset_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return idx_;
    }

    // HW1: YOUR CODE HERE
    /** Return this node's value */
    node_value_type& value () {
      return fetch().val_;
    }

    /** Deal with const graph and return this node's value */
    const node_value_type& value () const {
      return fetch().val_;
    }

    /** Return this node's degree */
    size_type degree() const {
      size_type node_uid = graphptr_->node_idx2uid_.at(idx_);
      size_type deg{};
      if (graphptr_->node2edge_.count(node_uid) == 0) {
        deg = 0;
      } else {
        deg = graphptr_->node2edge_.at(node_uid).size();
      }
      return deg;
    }

    /**
     * @brief Begin function for IncidentIterator
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_begin() const {
      // if no edge with this node
      size_type deg = degree();
      if (deg== 0) {
        return IncidentIterator();
      }
      size_type node_uid = graphptr_->node_idx2uid_.at(idx_);
      return IncidentIterator(graphptr_, *this, 
        graphptr_->node2edge_.at(node_uid).begin());
    }

    /**
     * @brief End function for IncidentIterator
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_end() const {
      // if no edge with this node
      if (degree() == 0) {
        return IncidentIterator();
      }
      size_type node_uid = graphptr_->node_idx2uid_.at(idx_);
      return IncidentIterator(graphptr_, *this, 
        graphptr_->node2edge_.at(node_uid).end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE

      return (graphptr_ == n.graphptr_ && idx_ == n.idx_);
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
      return graphptr_< n.graphptr_ || idx_ < n.idx_;  // Comparing by index
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    // Pointer back to the Graph container
    Graph* graphptr_;
    // This element's unique identification number
    size_type idx_;
    /** Private Constructor */
    Node(const Graph* graphptr, size_type idx)
    : graphptr_(const_cast<Graph*>(graphptr)), idx_(idx) {}
    
    internal_node& fetch() const {
      assert(idx_<graphptr_->num_nodes());
      size_type my_uid = graphptr_->node_idx2uid_.at(idx_);
      return *(graphptr_->nodes_.at(my_uid));
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return nodeset_size_;
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
  Node add_node(const Point& position, const node_value_type& node_value = node_value_type ()){
    // HW0: YOUR CODE HERE
    std::unique_ptr<internal_node> new_node_ptr = make_unique<internal_node>
      (position, maxnode_, nodeset_size_, node_value);
    nodes_[maxnode_] = std::move(new_node_ptr);
    node_idx2uid_[nodeset_size_] = maxnode_;
    maxnode_++;
    nodeset_size_++;
    return Node(this, nodeset_size_-1);
  }

  /** Remove a node from the graph, if it exists.
   * @return 1 if a node has been removed; 0 if not
   * @pre node @n is a valid node
   * @post has_node(const Node& @n) == false
   * @post If old has_node(const Node& @n), new num_nodes() == old num_nodes() - 1.
   *       Else,                        new num_nodes() == old num_nodes().
   *
   * The NodeIterator will skip removed node.
   * 
   * Complexity: O(n.degree()), approximately O(1) when we consider the 
   * sparsity of the graph.
   */
  size_type remove_node(const Node& n) {
    // check if has node
    if (!has_node(n)) {
      return 0;
    }

    // get idx and uid for the node
    size_type n_idx = n.index();
    assert(n_idx < num_nodes());
    size_type n_uid = node_idx2uid_.at(n_idx);
    
    // modify node2edge_ (remove all adjacent edges)
    std::vector<Edge> edges_remove{};
    for (auto e = n.edge_begin(); e != n.edge_end(); ++e) { 
      edges_remove.push_back(*e);
    }
    for (auto e:edges_remove) {
      const Edge edge = e;
      remove_edge(edge);
    }


    // modify node_idx2uid_ (swap its idx with the last and pop)
    size_type last_idx = nodeset_size_-1;
    size_type last_uid = node_idx2uid_[last_idx];
    node_idx2uid_[n_idx] = last_uid;
    node_idx2uid_.erase(last_idx);
    // change node_idx
    (*nodes_[last_uid]).idx_ = n_idx;  
    // change adjacent edges' node1_idx
    if (node2edge_.count(last_uid) != 0) {
      std::unordered_map<size_type, size_type> adj_nodes_e = node2edge_.at(last_uid);

      for (auto ei = adj_nodes_e.begin(); ei != adj_nodes_e.end(); ++ei) { 
        size_type edge_uid = ei->second;
        size_type n1_idx = edges_.at(edge_uid)->node1_idx_;
        size_type n2_idx = edges_.at(edge_uid)->node2_idx_;
        assert (n1_idx == last_idx || n2_idx == last_idx);
        if (n1_idx == last_idx) {
          edges_.at(edge_uid)->node1_idx_ = n_idx;
        }else {
          edges_.at(edge_uid)->node2_idx_ = n_idx;
        }
      }
    }

    // modify nodes_
    nodes_.erase(n_uid);

    // modify nodeset_size_
    nodeset_size_--;

    return 1;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
   // (void) n;            // Quiet compiler warning
    if (this == n.graphptr_ && n.idx_ < this->nodeset_size_) return true;  
    return false;
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
    Node result_node = Node(this, i);
    assert (result_node.index() == i);
    return result_node;     
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      graphptr_ = nullptr;
      idx_ = size_type(-1);
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(graphptr_, node1_idx_); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(graphptr_, node2_idx_); 
    }
    /** Return the length of this Edge */
    double length() const {
      return norm(node1().position() - node2().position());
    }

    /** Return this edge's idx */
    size_type index() {
      return fetch().idx_;
    }

        /** Return this edge's uid */
    size_type uid() {
      return fetch().uid_;
    }


    /** Return this edge's value */
    edge_value_type& value() {
      return fetch().val_;
    }

    /** Deal with const graph and return this edge's value */
    const edge_value_type& value () const {
      return fetch().val_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return (graphptr_ == e.graphptr_ && idx_ == e.idx_);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //(void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      return graphptr_< e.graphptr_ || idx_ < e.idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* graphptr_;         // The graph the edge belongs to
    size_type idx_;        // The position of the edge in the map edges_
    size_type node1_idx_;   // The idx of Node1()
    size_type node2_idx_;   // The idx of Node2()

    /** Private Constructor */
    Edge(const Graph* graphptr, size_type idx, size_type node1_idx, 
      size_type node2_idx): graphptr_(const_cast<Graph*>(graphptr)), idx_(idx), 
      node1_idx_(node1_idx), node2_idx_(node2_idx) {}

    internal_edge& fetch() const {
      assert(idx_<graphptr_->num_edges());
      size_type my_uid = graphptr_->edge_idx2uid_.at(idx_);
      return *(graphptr_->edges_.at(my_uid));
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edgeset_size_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    Edge result_edge = Edge(this, i, edges_.at(edge_idx2uid_.at(i))->node1_idx_, 
      edges_.at(edge_idx2uid_.at(i))->node2_idx_);
    return result_edge;   
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    assert(a.graphptr_ != nullptr && b.graphptr_ != nullptr);  
    assert(this->has_node(a) && this->has_node(b));

    size_type a_uid = this->node_idx2uid_.at(a.idx_);
    size_type b_uid = this->node_idx2uid_.at(b.idx_);

    if (this->node2edge_.count(a_uid)>0) {
      if (this->node2edge_.at(a_uid).count(b_uid)>0) {
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& edge_value 
    = edge_value_type ()) {
    // HW0: YOUR CODE HERE
    size_type a_idx = a.index();
    size_type b_idx = b.index();
    size_type a_uid = this->node_idx2uid_.at(a_idx);
    size_type b_uid = this->node_idx2uid_.at(b_idx);

    // check if has_edge 
    if (this->has_edge(a, b)) {
      size_type edge_uid = this->node2edge_.at(a_uid).at(b_uid);
      size_type edge_idx = (this->edges_.at(edge_uid))->idx_;
      return Edge(this, edge_idx, a_idx, b_idx);
    }
    
    assert(a != b);
    
    // internal_edge* new_edge_ptr = new internal_edge(maxedge_, edgeset_size_, a_idx, 
    //   b_idx, edge_value, true);
    std::unique_ptr<internal_edge> new_edge_ptr = make_unique<internal_edge>
      (maxedge_, edgeset_size_, a_idx, b_idx, edge_value);
    
    // modify edges_
    edges_[maxedge_] = std::move(new_edge_ptr);
    
    // modify node2edge_
    if (node2edge_.count(a_uid)>0) {
      node2edge_[a_uid][b_uid] = maxedge_;
    } else {
      std::unordered_map<size_type, size_type> partial_edge ({{b_uid, maxedge_}});
      node2edge_[a_uid] = partial_edge;
    }
    if (node2edge_.count(b_uid)>0) {
      node2edge_[b_uid][a_uid] = maxedge_;
    } else {
      std::unordered_map<size_type, size_type> partial_edge ({{a_uid, maxedge_}});
      node2edge_[b_uid] = partial_edge;
    }

    // modify edge_idx2uid_
    edge_idx2uid_[edgeset_size_] = maxedge_;

    // mdodify maxedge_, edgeset_size_
    maxedge_++;
    edgeset_size_++;

    return Edge(this, edgeset_size_-1, a_idx, b_idx);
  }

  /** Remove an edge from the graph, if it exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return 1 if an edge has been removed; 0 if not
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   *
   * Edges will be reindexed -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   *
   * The EdgeIterator will skip removed edge.
   * The IncidentIterator will skip removed edge.
   * 
   * Complexity: O(1)
   */
  size_type remove_edge(const Node& a, const Node& b) {
    size_type a_idx = a.index();
    size_type b_idx = b.index();
    size_type a_uid = node_idx2uid_.at(a_idx);
    size_type b_uid = node_idx2uid_.at(b_idx);

    // check if has_edge 
    if (!has_edge(a, b)) {
      return 0;
    }

    assert(a != b);

    // get uid and idx for the edge
    size_type r_uid = node2edge_[a_uid][b_uid];
    size_type r_idx = (*edges_[r_uid]).idx_;

    // modify edge_idx2uid_ (swap its idx with the last and pop)
    size_type last_uid = edge_idx2uid_[edgeset_size_-1];
    edge_idx2uid_[r_idx] = last_uid;  
    edges_[last_uid]->idx_ = r_idx;  // modify idx of the swapped       
    edge_idx2uid_.erase(edgeset_size_-1);
    
    // modify edges_ 
    edges_.erase(r_uid);

    // modify node2edge_
    node2edge_[a_uid].erase(b_uid);
    if (node2edge_[a_uid].size() == 0) node2edge_.erase(a_uid);

    node2edge_[b_uid].erase(a_uid);
    if (node2edge_[b_uid].size() == 0) node2edge_.erase(b_uid);
    
    // modify edgeset_size_
    edgeset_size_--;

    return 1;
  }

    /** Remove an edge from the graph, if it exists.
   * @pre @a e is an edge object
   * @return 1 if an edge has been removed; 0 if not
   * @post has_edge(@a a, @a b) == false
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges() - 1.
   *       Else,                        new num_edges() == old num_edges().
   * 
   * Edges will be reindexed -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   *
   * The EdgeIterator will skip removed edge.
   * The IncidentIterator will skip removed edge.
   * 
   * Complexity: O(1)
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
    // HW0: YOUR CODE HERE
    maxnode_ = 0;
    nodeset_size_ = 0;
    nodes_.clear();
    node_idx2uid_.clear();

    maxedge_ = 0;
    edgeset_size_ = 0;
    edges_.clear();
    node2edge_.clear();
    edge_idx2uid_.clear();

    std::cout<<"clear all pointers"<<std::endl;
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered <NodeIterator>{
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
    
    /**
     * @brief Dereference operator for NodeIterator
     * 
     * @pre The node being dereferenced is valid
     * @return Node 
     */
    Node operator*() const {
      //"Invalid node to derefence."
      assert(graphptr_ != nullptr && iteridx_ < graphptr_->nodeset_size_);
      return Node(graphptr_, iteridx_);
    }

    /**
     * @brief Advance operator for NodeIterator
     * 
     * @return NodeIterator& 
     */
    NodeIterator& operator++() {
      iteridx_ += 1;
      return *this;
    }

    /**
     * @brief Equal to operator for NodeIterator
     * 
     * @param n another NodeIterator
     * @return true if it is equal to @a n
     * @return false if it is not equal to @a n
     */
    bool operator==(const NodeIterator& n) const {
      return ((graphptr_ == n.graphptr_) && (iteridx_ == n.iteridx_));
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graphptr_;
    size_type iteridx_;

    NodeIterator(const Graph* graphptr, size_type iteridx) {
      graphptr_ = const_cast<Graph*>(graphptr);
      iteridx_ = iteridx;
    }
  };

  // HW1 #2: YOUR CODE HERE
  /**
   * @brief Begin of this NodeIterator
   * 
   * @return node_iterator 
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /**
   * @brief End of this NodeIterator
   * 
   * @return node_iterator 
   */
  node_iterator node_end() const {
    return NodeIterator(this, nodeset_size_);
  }
  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered <IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      graphptr_ = nullptr;
      node_ = Node(nullptr, size_type(-1));
    }

    // HW1 #3: YOUR CODE HERE
    /**
     * @brief Dereference operator for IncidentIterator
     * 
     * @pre the Edge it is dereferencing into is valid 
     * @return Edge 
     */
    Edge operator*() const {
      //"Invalid IncIter."
      assert(graphptr_ != nullptr && graphptr_->has_node(node_));

      size_type center_idx = node_.index();
      size_type adjacent_uid = (*mapiter_).first;
      size_type adjacent_idx = (graphptr_->nodes_.at(adjacent_uid))->idx_;
      size_type edge_uid = (*mapiter_).second;
      size_type edge_idx = (graphptr_->edges_.at(edge_uid))->idx_;

      return Edge(graphptr_, edge_idx, center_idx, adjacent_idx);
     }

    /**
     * @brief Advance operator for IncidentIterator
     * 
     * @return IncidentIterator& 
     */
    IncidentIterator& operator++() {
      ++mapiter_;
      // if edge not valid find the next valid edge
      size_type node_idx = node_.index();
      size_type node_uid = graphptr_->node_idx2uid_.at(node_idx);
      std::unordered_map<size_type, size_type> adj_nodes = graphptr_->node2edge_.at(node_uid);

      return *this;
    }

    /**
     * @brief Equal to operator for IncidentIterator
     * 
     * @param inc Another IncidentIterator
     * @return true if it is equal to @a inc
     * @return false if it is not equal to @a inc
     */
    bool operator==(const IncidentIterator& inc) const {
      return (graphptr_==inc.graphptr_ && node_==inc.node_ && mapiter_==inc.mapiter_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graphptr_;
    node_type node_;
    std::unordered_map<size_type, size_type>::const_iterator mapiter_;

    // private constructor
    IncidentIterator(const Graph* graphptr, const node_type &node, 
      std::unordered_map<size_type, size_type>::const_iterator mapiter)
      :graphptr_(const_cast<Graph*>(graphptr)), 
      node_(const_cast<node_type&>(node)), mapiter_(mapiter){
      // find the first valid edge
      // size_type node_uid = graphptr_->node_idx2uid_.at(node.index());
      // while(mapiter_ != graphptr_->node2edge_.at(node_uid).end() &&
      //   !(graphptr_->edges_.at((*mapiter_).second)->valid_) ) {
      //     ++mapiter_;
      //   }
    }
      
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered <EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      graphptr_ = nullptr;
      idx_ = size_type(-1);
    }

    // HW1 #5: YOUR CODE HERE
    /**
     * @brief Dereference operator for EdgeIter
     * 
     * @pre @a idx_ is in [0, number of edges]
     * @return Edge 
     */
    Edge operator*() const {
      // "Invalid EdgeIter invalid edge"
      // if (graphptr_ == nullptr || idx_ > graphptr_->num_edges() ||
      // !(graphptr_->edges_.at(graphptr_->edge_idx2uid_.at(idx_))->valid_) ) {
      //   std::cerr << "Error: Dereferencing to an invalid edge." << std::endl;
      //   exit(0);
      // }

      size_type edge_uid = graphptr_->edge_idx2uid_.at(idx_);
      size_type node1_idx = (graphptr_->edges_.at(edge_uid))->node1_idx_;
      size_type node2_idx = (graphptr_->edges_.at(edge_uid))->node2_idx_;

      return Edge(graphptr_, idx_, node1_idx, node2_idx);
    }

    /**
     * @brief Advance operator for EdgeIterator
     * 
     * @return EdgeIterator& 
     */
    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }

    /**
     * @brief Equal to operator for EdgeIterator
     * 
     * @param e another EdgeIterator
     * @return true if it is equal to @a e
     * @return false if it is not equal to @a e
     */
    bool operator==(const EdgeIterator& e) const {
      return (graphptr_==e.graphptr_ && idx_==e.idx_);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graphptr_;
    size_type idx_;

    // private constructor
    EdgeIterator(const Graph* graphptr, size_type idx)
      :graphptr_(const_cast<Graph*>(graphptr)), idx_(idx){}

  };

  // HW1 #5: YOUR CODE HERE
  /**
   * @brief Begin function of EdgeIterator
   * 
   * @return edge_iterator 
   */
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  /**
   * @brief End function of EdgeIterator 
   * 
   * @return edge_iterator 
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, edgeset_size_);
  }

  template<typename T, typename... Args>
std::unique_ptr<T> make_unique(Args&&... args) {
    return std::unique_ptr<T>(new T(std::forward<Args>(args)...));
}

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal type for set nodes
  struct internal_node {
    Point pos_;            // The position of a node
    size_type uid_;        // The unique identifcation for a node
    size_type idx_;        // The position in map
    node_value_type val_;  // The template value for node 

    internal_node(Point pos, size_type uid, size_type idx, node_value_type val)
    :pos_(pos),uid_(uid), idx_(idx), val_(val){}
  };

  // Internal type for set edges
  struct internal_edge {
    size_type uid_;        // The unique identifcation for an edge
    size_type idx_;        // The position in map
    size_type node1_idx_;  // The index of node1
    size_type node2_idx_;  // The index of node2
    edge_value_type val_;  // The template value for edge

    internal_edge(size_type uid, size_type idx, size_type node1_idx, 
      size_type node2_idx, edge_value_type val)
    :uid_(uid), idx_(idx), node1_idx_(node1_idx), node2_idx_(node2_idx),
      val_(val){}
  };

  // number of nodes added so far
  size_type maxnode_;     
  // number of edges added so far
  size_type maxedge_;
  // size of nodeset
  size_type nodeset_size_;
  // size of edgeset
  size_type edgeset_size_;
  // the map from node uid to internal node
  std::unordered_map<size_type, std::unique_ptr<internal_node>> nodes_;
  // the map from edge uid to internal edge
  std::unordered_map<size_type, std::unique_ptr<internal_edge>> edges_;
  // the map from node1_uid -> node2_uid -> edge_uid
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> node2edge_;
  // Node: the map from idx to uid
  std::unordered_map<size_type, size_type> node_idx2uid_;
  // Node: the map from idx to uid
  std::unordered_map<size_type, size_type> edge_idx2uid_;
  
};

#endif // CME212_GRAPH_HPP
