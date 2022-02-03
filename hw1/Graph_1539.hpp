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
template <typename V>
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
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
      : nodes_(), edges_(), size_(0), num_edges_(0) {
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
      return ((graph_->nodes_).at(id_)).pt_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return ((graph_->nodes_).at(id_).attri_);
    };

    const node_value_type& value() const {
      const internal_node& node = (graph_->nodes_).at(id_);
      return node.attri_;
    };
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ( (this->graph_ == n.graph_) && (this->id_ == n.id_) );
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
      return (this->id_ < n.id_);
    }

    // return the number of incident edges.
    size_type degree () const {
      return (((graph_->nodes_).at(id_)).outgoing_edges_).size();
    }

    // Start of the incident iterator.
    incident_iterator edge_begin () const {
      return IncidentIterator(graph_, id_, 0);
    }

    // End of incident iterator.
    incident_iterator edge_end () const {
      size_t deg = degree();
      return IncidentIterator(graph_, id_, deg);
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
    size_t id_;

    /* Private constructor */
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
        {position, attri, connected, outgoing_edges};
    nodes_.insert(std::pair<size_t, internal_node> (size_, node));

    // Increment size
    size_++;

    return Node(this, size_-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ( (n.graph_ == this) && (n.id_ < size_ ));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i<size_);
    return Node(this, i);        // Invalid node
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
      return (graph_->edges_).at(id_).node1_;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return (graph_->edges_).at(id_).node2_;      // Invalid Node
    }

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
      return (id_ < e.id_);
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
    size_t id_;

    // Private Constructor for Edge
    Edge(const Graph* graph, size_type id)
        : graph_(const_cast<Graph*>(graph)), id_(id) {
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
    assert(i<num_edges_);
    return Edge(this, i);        // Invalid Edge
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
    size_t a_idx = a.index();
    size_t b_idx = b.index();

    // get all the nodes a is connected to
    //const std::unordered_set <size_t>& a_nodes_ 
        //= nodes_.at(a_idx).connected_nodes_;
    
    // return is a and b has been connected
    if (nodes_.at(a_idx).connected_nodes_.find(b_idx) 
      != nodes_.at(a_idx).connected_nodes_.end()) {
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
          = nodes_.at(a.index()).outgoing_edges_;
      size_t max_iter = outgoing_edges.size();

      // loop through the edges
      while (!found && i < max_iter) {
        internal_edge edge = edges_.at(outgoing_edges[i]);
        if (edge.node1_ == a && edge.node2_ == b) {
          found = true;
        } else if (edge.node1_ == b && edge.node2_ == a) {
          found = true;
          edges_.at(outgoing_edges[i]) = internal_edge{a,b};
        } else {
          i++;
        }
      }; 

      return Edge(this, outgoing_edges[i]);

    } else {

      // create internal edge and add it to the graph
      internal_edge edge = internal_edge{a,b};
      edges_.insert(std::pair<size_t, internal_edge> (num_edges_, edge));

      // indicate two nodes are connected in internal representation
      nodes_.at(a.index()).connected_nodes_.insert(b.index());
      nodes_.at(b.index()).connected_nodes_.insert(a.index());

      // increment outgoing edges of each node
      nodes_.at(a.index()).outgoing_edges_.push_back(num_edges_);
      nodes_.at(b.index()).outgoing_edges_.push_back(num_edges_);

      // Increment size
      num_edges_++;

      return Edge(this, num_edges_-1);   

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
    nodes_.clear();
    edges_.clear();
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
    return NodeIterator(this, size_);
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
      size_t edge_id = 
          (((graph_->nodes_).at(node_idx_)).outgoing_edges_).at(idx_);

      // swap the two ends
      if (node_idx_!=graph_->edges_.at(edge_id).node1_.index()){
        Node temp = graph_->edges_.at(edge_id).node1_;
        graph_->edges_.at(edge_id).node1_ = graph_->edges_.at(edge_id).node2_;
        graph_->edges_.at(edge_id).node2_ = temp;
      };

      return Edge(graph_, edge_id);
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
    return EdgeIterator(this, num_edges_);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node {
    Point pt_;
    node_value_type attri_;
    std::unordered_set <size_t> connected_nodes_;
    std::vector <size_t> outgoing_edges_;
  };

  struct internal_edge {
    Node node1_;
    Node node2_;
  };

  std::unordered_map<size_t, internal_node> nodes_;
  std::unordered_map<size_t, internal_edge> edges_;
  size_t size_;
  size_t num_edges_;


};

#endif // CME212_GRAPH_HPP
