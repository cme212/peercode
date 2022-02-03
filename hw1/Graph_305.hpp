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

//Template to generalize nodes
template <typename V>

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

  // Declare struct used in Node class
  struct internal_node;

  // Declare struct used in Edge class
  struct internal_pair;
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Predeclaration of generalized Node type*/
  using node_value_type = V;

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
  Graph() : size_nodes_(0), size_edges_(0) {
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
    Node() : graph_(NULL) {
    
    }
    
  
    /**
     * @brief Returns number of edges incident to node
     * 
     */    
    size_type degree() const {
      return graph_->adjacencies_[idx_].size();
    }
   
    /**
     * @brief Returns iterator to first of edges incident to node.
     */ 
    incident_iterator edge_begin() const {
      return IncidentIterator(this, 0);
    }
   
    /**
     * @brief Returns end iterator of edges incident to node. 
     */ 
    incident_iterator edge_end() const {
      return IncidentIterator(this, this.degree());
    }


    /**
     * @brief Returns value of node's attribute
     * 
     * Default initializes if no value is provided.
     */ 

    node_value_type& value() {
      return graph_->nodes_[idx_].node_value;
    };

    /**
     * @brief Returns value of node's attribute
     * 
     * Default initializes if no value is provided.
     */ 
    const node_value_type& value() const {
      return graph_->nodes_[idx_].node_value;
    };

    /** Return this node's position. */
    const Point& position() const {
      assert(graph_ != NULL);   // Invalid Node
      return graph_->nodes_[idx_].node_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      
      return idx_;
    }

   
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (graph_ == n.graph_ && idx_ == n.index()) return true;
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
      assert(graph_ != NULL && n.graph_ != NULL);   // Invalid Node
      
      return idx_ < n.index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
       
    // Links node to access graph class private data
    Graph* graph_;

    // Identifies node with index
    size_type idx_;

    // Private constructor reserved for internal initialization with .add_node()
    Node(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), idx_(index) {
    }

   
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return size_nodes_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    internal_node new_node;
    new_node.node_point = position;
    new_node.node_value = value;

    nodes_.push_back(new_node);
    ++size_nodes_;
    return Node(this, size_nodes_-1);       
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {    
    return n.index() < size_nodes_ && n.graph_ == this;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
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
    Edge() : graph_(NULL) {
    }

    /** Return a node of this Edge */
    Node node1() const {
     
      internal_pair node_pair = graph_->edges_[idx_];
      return node_pair.node_a;     
    }

    /** Return the other node of this Edge */
    Node node2() const {
     
      internal_pair node_pair = graph_->edges_[idx_];
      return node_pair.node_b; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
     
      internal_pair node_pair = graph_->edges_[idx_]; 
      if (node_pair.node_a == e.node1() && node_pair.node_b == e.node2()) return true;
      if (node_pair.node_a == e.node2() && node_pair.node_b == e.node1()) return true;

      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      assert(graph_ != NULL);     // Invalid Edge
      return idx_ < e.idx_;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
  
    // Links edge to access graph class private data 
    Graph* graph_;

    // Identifies edge with index
    size_type idx_;


    Edge(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), idx_(index) {
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return size_edges_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return Edge(this, i);       
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.graph_ != NULL && b.graph_ != NULL);   // Invalid nodes
    for (size_type i = 0; i < size_edges_; i++) {
      internal_pair edge = edges_[i];
      if (edge.node_a == a && edge.node_b == b) return true;
      if (edge.node_a == b && edge.node_b == a) return true;
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
   
    assert(a.graph_ != NULL && b.graph_ != NULL);
    
    //Implement feature to avoid adding edges if already existing
    for (size_type i = 0; i < size_edges_; i++) {
      internal_pair edge = edges_[i];
      if (edge.node_a == a && edge.node_b == b) return Edge(this, i);
      if (edge.node_a == b && edge.node_b == a) {
        edge.node_a = b; edge.node_b = a;
        return Edge(this,i);
      }
    } 
      //Add edges with nodes in correct order
    internal_pair edge_nodes;
    edge_nodes.node_a = a;
    edge_nodes.node_b = b;
    edges_.push_back(edge_nodes);

    // Update adjacencies
    adjacencies_[a.index()].push_back(size_edges_);
    adjacencies_[b.index()].push_back(size_edges_);

    ++size_edges_;

    return Edge(this, size_edges_-1);       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
       
    nodes_.clear();
    edges_.clear();

    adjacencies_.clear();
    
    size_nodes_ = 0;
    size_edges_ = 0;
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
    NodeIterator() : graph_(NULL) {
    }

    /*
     * @brief Defines dereference operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */
    Node operator*() const {
      return Node(graph_, idx_);
    }

    /*
     * @brief Defines dereference operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    Node operator*() {
      return Node(graph_, idx_);
    }

    /*
     * @brief Defines prefix increment  operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */
    NodeIterator& operator++() {
      ++idx_;
      return *this;
    }

   /*
     * @brief Defines equality operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */
    bool operator==(const NodeIterator& it) const {
      return Node(graph_, idx_) == *it;
    }


   private:
    friend class Graph;

    Graph* graph_;

    size_type idx_;

    NodeIterator(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), idx_(index) {}

    
  };

  /*
   * @brief Returns iterator to first node in graph
   *
   * @pre Assumes that graph is defined and nonempty
   */
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }

  /*
   * @brief Returns end iterator
   *
   * @pre Assumes that graph is defined and nonempty
   */
  node_iterator node_end() const{
   i return NodeIterator(this, size_nodes_);
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
    IncidentIterator() : graph_(NULL) {
    }
  
    /*
     * @brief Defines dereference operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    Edge operator*() const {
      size_type edge_idx = graph_->adjacencies_[node_.index()][idx_];
      if (graph_->edges_[edge_idx].node_a != node_) {
        graph_->edges_[edge_idx].node_b = graph_->edges_[edge_idx].node_a;
        graph_->edges[edge_idx].node_a = node_; 
      }
      return Edge(graph_, edge_idx);
    }
   
    /*
     * @brief Defines prefix increment  operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    IncidentIterator& operator++() {
      ++idx_;
      return *this;
    }
   
   /*
     * @brief Defines equality operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    bool operator==(const IncidentIterator& it) const {
      return graph_ == it.graph_ && node_ == it.node_ && idx_ == it.idx_;
    }

   private:
    friend class Graph;
    
    Graph* graph_;

    Node node_;

    size_type idx_;

    IncidentIterator(const Graph* graph, Node node, size_type index) : graph_(const_cast<Graph*>(graph)), node_(node), idx_(index) {}
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
    EdgeIterator() : graph_(NULL) {
    }
  
    /*
     * @brief Defines dereference operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    Edge operator*() const {
      return Edge(graph_, idx_);
    }
    
    /*
     * @brief Defines prefix increment  operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    EdgeIterator& operator++() {
      ++idx_;
      return *this;
    }
    

   /*
     * @brief Defines equality operator.
     *
     * @pre Assumes that graph is defined and node is valid
     */

    bool operator==(const EdgeIterator& it) const {
      return graph_ == it.graph_ && idx_ == it.idx_;
    }


   private:
    friend class Graph;
    
    Graph* graph_;

    size_type idx_; 

    EdgeIterator(const Graph* graph, size_type index) : graph_(const_cast<Graph*>(graph)), idx_(index) {}


  };

  /*
   * @brief Returns iterator to first edge in graph
   *
   * @pre Assumes that graph is defined and nonempty
   */

  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  
  /*
   * @brief Returns end iterator
   *
   * @pre Assumes that graph is defined and nonempty
   */
  edge_iterator edge_end() const {
    return EdgeIterator(this, size_edges_);
  }

 private:

  // Stores nodes' points, attributes, and number of nodes
  struct internal_node {
    Point node_point;
    node_value_type node_value;
  };

  std::vector<internal_node> nodes_;
  size_type size_nodes_;

  // Stores edges' nodes and number of edges

  size_type size_edges_; 

  // Vector of pairs of nodes to represent each edge
  struct internal_pair {
    node_type node_a;
    node_type node_b;
  };

  // Stores indices of edges adjacent to key node
  std::map<size_type, std::vector<size_type>> adjacencies_;

  std::vector<internal_pair> edges_; 

};

#endif // CME212_GRAPH_HPP
