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

  /** Declare for class template */
  typedef V node_value_type;
  typedef E edge_value_type;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), 
        active_nodes_(),
        size_(0), 
        tol_num_nodes_(0),
        edges_(), 
        active_edges_(), 
        num_edges_(0), 
        tol_num_edges_(0) {
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
      this->graph_ = nullptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return ((graph_->nodes_).at(uid_)).pt_;
    }

    /** Return this node's position but not constant */
    Point& position() {
      return ((graph_->nodes_).at(uid_)).pt_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      size_t idx = ((graph_->nodes_).at(uid_)).idx_;
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return ((graph_->nodes_).at(uid_)).attri_;
    };

    const node_value_type& value() const {
      const internal_node& node = ((graph_->nodes_).at(uid_));
      return node.attri_;
    };
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ( (this->graph_ == n.graph_) && (this->uid_ == n.uid_) );
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
      return (this->uid_ < n.uid_);
    }

    // return the number of incident edges.
    size_type degree () const {
      return (
        (graph_->nodes_).at(uid_).outgoing_edges_
      ).size();
    }

    // Start of the incident iterator.
    incident_iterator edge_begin () const {
      size_t idx = graph_->nodes_.at(uid_).idx_;
      return IncidentIterator(graph_, idx, 0);
    }

    // End of incident iterator.
    incident_iterator edge_end () const {
      size_t deg = degree();
      size_t idx = graph_->nodes_.at(uid_).idx_;
      return IncidentIterator(graph_, idx, deg);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer to the graph container
    Graph* graph_;

    // This node's index in the graph
    size_t uid_;

    /* Private constructor */
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
    assert(size_==active_nodes_.size());
    return this->size_;
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
  Node add_node(
    const Point& position, const node_value_type& attri = node_value_type()) {
    // HW0: YOUR CODE HERE
        
    // add this point object to the internal representation of nodes
    std::unordered_set <size_t> connected;
    std::vector <size_t> outgoing_edges;
    internal_node node = internal_node
        {size(), position, attri, connected, outgoing_edges};
    nodes_.insert(std::pair<size_t, internal_node> (tol_num_nodes_, node));
    active_nodes_.push_back(tol_num_nodes_);

    // Increment size
    size_++;
    tol_num_nodes_++;

    return Node(this, tol_num_nodes_-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ( (n.graph_ == this) && 
             (n.graph_->nodes_.at(n.uid_).idx_ < size() ));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<size());
    size_t uid = active_nodes_.at(i);
    return Node(this, uid);        // Invalid node
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
      this->graph_ = nullptr;
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return graph_->edges_.at(uid_).node1_; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return graph_->edges_.at(uid_).node2_; 
    }

    /** Return value associated with the edge */
    edge_value_type& value () {
      return (graph_->edges_.at(uid_).attri_);
    };

    /** Return value associated with the edge */
    const edge_value_type& value () const {
      const internal_edge& edge = graph_->edges_.at(uid_);
      return edge.attri_;
    };

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return ( (graph_ == e.graph_) 
            && (node1()==e.node1()) 
            && (node2()==e.node2()) );
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (&graph_ != &e.graph_) {
        return (&graph_ < &e.graph_);
      };

      return ( 
        uid_ < e.uid_
      );
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer to the graph container
    Graph* graph_;

    // This node's index in the graph
    size_t uid_;

    // Private Constructor for Edge
    Edge(const Graph* graph, size_type uid)
        : graph_(const_cast<Graph*>(graph)), uid_(uid) {
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    assert(num_edges_ == active_edges_.size());
    return num_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<num_edges());
    size_t uid = active_edges_.at(i);
    return Edge(this, uid);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    // get id of both nodes
    size_t a_uid = active_nodes_.at(a.index());
    size_t b_uid = active_nodes_.at(b.index());
    
    // return is a and b has been connected
    if (nodes_.at(a_uid).connected_nodes_.find(b_uid) 
      != nodes_.at(a_uid).connected_nodes_.end()) {
      return true;
    } else {
      return false;
    };

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
    
    bool dont_add = has_edge(a,b);

    if (dont_add) {

      // if a and b are connected, return the edge

      size_t i = 0;
      bool found = false;
      std::vector <size_t>& outgoing_edges 
        = nodes_.at(active_nodes_.at(a.index())).outgoing_edges_;
      size_t max_iter = outgoing_edges.size();

      size_t edge_idx;

      // loop through the edges
      while (!found && i < max_iter) {
        internal_edge edge = edges_.at(outgoing_edges[i]);
        edge_idx = edges_.at(outgoing_edges[i]).idx_;
        if (edge.node1_ == a && edge.node2_ == b) {
          found = true;
        } else if (edge.node1_ == b && edge.node2_ == a) {
          found = true;
          edges_.at(outgoing_edges[i]) = internal_edge{edge_idx,a,b,E()};
        } else {
          i++;
        }
      }; 

      return Edge(this, outgoing_edges[i]);

    } else {

      // create internal edge and add it to the graph
      internal_edge edge = internal_edge{num_edges(), a, b, E()};
      edges_.insert(std::pair<size_t, internal_edge> (tol_num_edges_, edge));
      active_edges_.push_back(tol_num_edges_);


      // get uid of both nodes
      size_t a_uid = active_nodes_.at(a.index());
      size_t b_uid = active_nodes_.at(b.index());

      // indicate two nodes are connected in internal representation
      nodes_.at(a_uid).connected_nodes_.insert(b_uid);
      nodes_.at(b_uid).connected_nodes_.insert(a_uid);

      // increment outgoing edges of each node
      nodes_.at(a_uid).outgoing_edges_.push_back(tol_num_edges_);
      nodes_.at(b_uid).outgoing_edges_.push_back(tol_num_edges_);

      // Increment size
      num_edges_++;
      tol_num_edges_++;

      return Edge(this, tol_num_edges_-1);   

    };
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    num_edges_ = 0;
    tol_num_edges_ = 0;
    tol_num_nodes_ = 0;
    nodes_.clear();
    edges_.clear();
    active_edges_.clear();
    active_nodes_.clear();
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
    NodeIterator(const Graph* graph, size_t idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const {
      return (*graph_).node(idx_);
    }
    
    NodeIterator& operator++() {
      idx_++;
      return *this;
    }
    
    bool operator==(const NodeIterator& node_iter) const {
      return (graph_ == node_iter.graph_ && idx_ == node_iter.idx_);
    }

    bool operator!=(const NodeIterator& node_iter) const {
      return !(*this == node_iter);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_t idx_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  
  node_iterator node_end() const {
    return NodeIterator(this, size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator(const Graph* graph, size_t node_idx, size_t idx)
      : graph_(const_cast<Graph*>(graph)), node_idx_(node_idx), idx_(idx) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {

      // retrieve the edge id
      size_t node_uid = graph_->active_nodes_.at(node_idx_);
      internal_node& i_e = (graph_->nodes_).at(node_uid);
      size_t edge_uid = (i_e.outgoing_edges_).at(idx_);
      
      // swap the two ends
      if (node_idx_!=graph_->edges_.at(edge_uid).node1_.index()){
        Node temp = graph_->edges_.at(edge_uid).node1_;
        graph_->edges_.at(edge_uid).node1_ = graph_->edges_.at(edge_uid).node2_;
        graph_->edges_.at(edge_uid).node2_ = temp;
      };

      return Edge(graph_, edge_uid);
    }

    IncidentIterator& operator++() {
      ++idx_;
      return *this;
    }

    bool operator==(const IncidentIterator& iter) const {
      return (graph_==iter.graph_ 
          && node_idx_==iter.node_idx_ 
          && idx_==iter.idx_);
    }

    // bool operator!=(const IncidentIterator& iter) const {
    //   return !(*this == iter);
    // }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    
    Graph* graph_;
    size_t node_idx_;
    size_t idx_;    

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
    EdgeIterator(const Graph* graph, size_t idx)
      : graph_(const_cast<Graph*>(graph)), idx_(idx) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return (*graph_).edge(idx_);
    }

    EdgeIterator& operator++() {
      idx_++;
      return *this;
    }
    
    bool operator==(const EdgeIterator& iter) const {
      return (graph_==iter.graph_ && idx_==iter.idx_);
    }

    bool operator!=(const EdgeIterator& iter) const {
      return !(*this==iter);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_t idx_;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

  /** Helper function to check if a edge is valid to be removed 
   *
   *  A valid edge points to this graph g, has an 0 =< idx_ < n(active 
   *  edges)an uid_ less than n(edges ever added to graph), and is 
   *  stored in sync in edges_ with matching idx_ */
  bool valid_edge(const Edge& e) {

    return (
      (e.graph_) == this &&
      e.uid_ < tol_num_edges_ &&
      edges_.at(e.uid_).idx_ >= 0 &&
      edges_.at(e.uid_).idx_ < num_edges() &&
      active_edges_.at(edges_.at(e.uid_).idx_) == e.uid_
    );

  };


  /** Remove an edge from the graph, return 1 if successful and 0 otherwise.
   *
   * @param[in, out] e is a valid edge object with two valid ends a and b 
   * @return 1 if the input edge has been successfully removed, 0 otherwise
   *
   * @post If @ e is a valid object of the graph:
   *           new num_edges() == old num_edges() - 1
   *           new active_edges_.size() == old active_edges_.size() - 1
   *           new tol_num_edges_ == old tol_num_edges_
   *           (b's uid_) in a's connected_nodes_ == false
   *           (a's uid_) in b's connected_nodes_ == false
   *           (@ e's uid_) in a's outgoing_edges == false
   *           (@ e's uid_) in b's outgoing_edges == false
   *           
   * @post If @ e is NOT a valid object of the graph:
   *           no internal structure is modified
   *           new num_edges() == old num_edges()
   *           new active_edges_.size() == old active_edges_.size()
   *           new tol_num_edges_ == old tol_num_edges_
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   *
   * Can invalidate edge iterator -- while iterating through graph edges, 
   * remove edge may cause old (*it) != new (*it)  
   *
   * Complexity: expect O(# connected nodes and outgoing edges of @ e's ends),
   *             at most O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Edge& e) {

    if (this->valid_edge(e)) {

      size_t edge_uid = e.uid_;
      size_t edge_idx = edges_.at(edge_uid).idx_;
      size_t back_uid = active_edges_.at(this->num_edges()-1);
      size_t back_idx = this->num_edges()-1;

      size_t n1_uid = edges_.at(edge_uid).node1_.uid_;
      size_t n2_uid = edges_.at(edge_uid).node2_.uid_;

      // disconnect the two ends of e
      nodes_.at(n1_uid).connected_nodes_.erase(
        nodes_.at(n1_uid).connected_nodes_.find(n2_uid)
      );

      nodes_.at(n2_uid).connected_nodes_.erase(
        nodes_.at(n2_uid).connected_nodes_.find(n1_uid)
      );

      // remove e as an outgoing edge of its two ends
      nodes_.at(n1_uid).outgoing_edges_.erase(
        std::find(
          nodes_.at(n1_uid).outgoing_edges_.begin(),
          nodes_.at(n1_uid).outgoing_edges_.end(),
          edge_uid
        )
      );

      nodes_.at(n2_uid).outgoing_edges_.erase(
        std::find(
          nodes_.at(n2_uid).outgoing_edges_.begin(),
          nodes_.at(n2_uid).outgoing_edges_.end(),
          edge_uid
        )
      );
      
      // swap idx in internal edges_ to prepare for removal
      edges_.at(back_uid).idx_ = edge_idx;
      edges_.at(edge_uid).idx_ = back_idx;

      // swap position in active_edges_ to prepare for removal
      active_edges_.at(edge_idx) = active_edges_.back();

      // remove edge from active_edges_
      active_edges_.pop_back();
      num_edges_--;

      return 1;

    } else {

      return 0;

    }

  };

  // returns the edge object that connects n1 and n2
  Edge find_edge(const Node& n1, const Node& n2) {

    size_t i = 0;
    bool found = false;
    std::vector <size_t>& outgoing_edges 
      = nodes_.at(active_nodes_.at(n1.index())).outgoing_edges_;
    size_t max_iter = outgoing_edges.size();

    // loop through the edges
    while (!found && i < max_iter) {
      internal_edge edge = edges_.at(outgoing_edges[i]);
      if (edge.node1_ == n1 && edge.node2_ == n2) {
        found = true;
      } else if (edge.node1_ == n2 && edge.node2_ == n1) {
        found = true;
      } else {
        i++;
      }
    }; 

    return Edge(this, outgoing_edges[i]);

  };


  /** Remove an edge from the graph using its two end points
   *
   * @param[in, out] a, b are two valid node objects 
   * @return 1 if the edge connecting the input nodes has been 
   *         successfully removed, 0 otherwise
   *
   * @post If has_edge(a,b):
   *           new num_edges() == old num_edges() - 1
   *           new active_edges_.size() == old active_edges_.size() - 1
   *           new tol_num_edges_ == old tol_num_edges_
   *           (b's uid_) in a's connected_nodes_ == false
   *           (a's uid_) in b's connected_nodes_ == false
   *           (@ edge's uid_) in a's outgoing_edges == false
   *           (@ edge's uid_) in b's outgoing_edges == false
   *           
   * @post Else:
   *           no internal structure is modified
   *           new num_edges() == old num_edges()
   *           new active_edges_.size() == old active_edges_.size()
   *           new tol_num_edges_ == old tol_num_edges_
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   *
   * Can invalidate edge iterator -- while iterating through graph edges, 
   * remove edge may cause old (*it) != new (*it)  
   *
   * Complexity: expect O(# connected nodes and outgoing edges of @ e's ends),
   *             at most O(num_nodes() + num_edges())
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    
    bool exist = has_edge(n1, n2);
    
    if (exist) {
      Edge e = find_edge(n1, n2);
      return remove_edge(e);
    } else {
      return 0;
    }

  };


  /** Remove an edge from the graph using its iterator
   *
   * @param[in, out] e_it is a edge iterator 
   * @return edge_begin() if the edge referenced by @ e_it has been 
   *         successfully removed, edge_end() otherwise
   *
   * @post If e_it != edge_end():
   *           new num_edges() == old num_edges() - 1
   *           new active_edges_.size() == old active_edges_.size() - 1
   *           new tol_num_edges_ == old tol_num_edges_
   *           (b's uid_) in a's connected_nodes_ == false
   *           (a's uid_) in b's connected_nodes_ == false
   *           (@ edge's uid_) in a's outgoing_edges == false
   *           (@ edge's uid_) in b's outgoing_edges == false
   *           
   * @post Else:
   *           no internal structure is modified
   *           new num_edges() == old num_edges()
   *           new active_edges_.size() == old active_edges_.size()
   *           new tol_num_edges_ == old tol_num_edges_
   *
   * Can invalidate edge indexes -- in other words, old edge(@a i) might not
   * equal new edge(@a i). 
   *
   * Can invalidate edge iterator -- while iterating through graph edges, 
   * remove edge may cause old (*it) != new (*it)  
   *
   * Complexity: expect O(# connected nodes and outgoing edges of @ e's ends),
   *             at most O(num_nodes() + num_edges())
   */
  edge_iterator remove_edge(edge_iterator e_it) {

    if (e_it == this->edge_end()) {

      return e_it;

    } else {

      Edge e = *e_it;
      size_t temp = remove_edge(e);
      (void) temp;

      return this->edge_begin();
    }; 

  }


  /** Helper function to check if a edge is node to be removed 
   *
   *  A valid node belongs to this graph, has an 0 =< idx_ < n(active 
   *  nodes), an uid_ less than n(nodes ever added to graph), 
   *  and is stored in sync in nodes_ with matching idx_ */
  bool valid_node(const Node& n) {

    return (
      (n.graph_) == this &&
      n.uid_ < tol_num_nodes_ &&
      nodes_.at(n.uid_).idx_ >= 0 &&
      nodes_.at(n.uid_).idx_ < size_ &&
      active_nodes_.at(nodes_.at(n.uid_).idx_) == n.uid_
    );

  };

  /** Remove a node and all its outgoing edges from the graph
   *
   * @param[in, out] n is a valid node object
   * @return 1 if the input node has been successfully removed, 0 otherwise
   *
   * @post If @ n is a valid object of the graph:
   *           new num_nodes() == old num_edges() - 1
   *           new active_nodes_.size() == old active_nodes_.size() - 1
   *           new tol_num_nodes_ == old tol_num_nodes_
   *           new num_edges() == old num_edges() - n(outgoing edges from n)
   *           new active_edges_.size() 
   *                  == old active_edges_.size() - n(outgoing edges from n)
   *           new tol_num_edges_ == old tol_num_edges_
   *           
   * @post If @ e is NOT a valid object of the graph:
   *           no internal structure is modified
   *           new num_nodes() == old num_edges() 
   *           new active_nodes_.size() == old active_nodes_.size() 
   *           new tol_num_nodes_ == old tol_num_nodes_
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   *
   * Can invalidate node iterator -- while iterating through graph nodes, 
   * remove node may cause old (*it) != new (*it)  
   *
   * Complexity: expect n(outgoing edges) * O(remove_edge()),
   *             at most O(num_nodes())
   */
  size_type remove_node(const Node& n) {

    if (this->valid_node(n)) {

      size_t node_uid = n.uid_;
      size_t node_idx = nodes_.at(node_uid).idx_;

      // iterate to remove all incident edges
      for (auto ei = n.edge_begin(); ei != n.edge_end();) {
        this->remove_edge(*ei);
        ei = n.edge_begin();
      };

      size_t back_uid = active_nodes_.at(this->size()-1);
      size_t back_idx = this->size()-1;

      // swap in internal nodes_ to prepare for removal
      nodes_.at(back_uid).idx_ = node_idx;
      nodes_.at(node_uid).idx_ = back_idx;
      
      // swap in active_nodes to prepare for removal
      active_nodes_.at(node_idx) = back_uid;

      // remove
      active_nodes_.pop_back();
      size_--;

      return 1;

    } else {

      return 0;

    };

  };


  /** Remove a node from the graph through its iterator
   *
   * @param[in] n_it is a node iterator
   * @return node_begin() if (*n_it) has been successfully removed, 
   *         node_end() otherwise
   *
   * @post If n_it != node_end():
   *           new num_nodes() == old num_edges() - 1
   *           new active_nodes_.size() == old active_nodes_.size() - 1
   *           new tol_num_nodes_ == old tol_num_nodes_
   *           new num_edges() == old num_edges() - n(outgoing edges from n)
   *           new active_edges_.size() 
   *                  == old active_edges_.size() - n(outgoing edges from n)
   *           new tol_num_edges_ == old tol_num_edges_
   *           
   * @post Else:
   *           no internal structure is modified
   *           new num_nodes() == old num_edges() 
   *           new active_nodes_.size() == old active_nodes_.size() 
   *           new tol_num_nodes_ == old tol_num_nodes_
   *
   * Can invalidate node indexes -- in other words, old node(@a i) might not
   * equal new node(@a i). 
   *
   * Can invalidate node iterator -- while iterating through graph nodes, 
   * remove node may cause old (*it) != new (*it)  
   *
   * Complexity: expect n(outgoing edges) * O(remove_edge()),
   *             at most O(num_nodes())
   */
  node_iterator remove_node(node_iterator n_it) {

    if (n_it == node_end()) {

      return n_it;

    } else {

      Node n = *n_it;
      size_t temp = remove_node(n);
      (void) temp;

      return this->node_begin();

    }

  };


 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    size_t idx_; // used to index active_nodes_
    Point pt_;
    V attri_;
    std::unordered_set <size_t> connected_nodes_;
    std::vector <size_t> outgoing_edges_;
  };

  struct internal_edge {
    size_t idx_; // used to index active_edges_
    Node node1_;
    Node node2_;
    E attri_;
  };

  // information of all nodes
  std::unordered_map<size_t, internal_node> nodes_;

  // currently active nodes
  std::vector <size_t> active_nodes_;

  // number of active nodes
  size_t size_;

  // total number of nodes added so far
  size_t tol_num_nodes_; 

  // information of all edges
  std::unordered_map<size_t, internal_edge> edges_;

  // currently active edges
  std::vector <size_t> active_edges_;

  // number of active edges
  size_t num_edges_;

  // total number of edges added so far
  size_t tol_num_edges_;


};

#endif // CME212_GRAPH_HPP
