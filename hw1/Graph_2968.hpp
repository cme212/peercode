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

using namespace std;


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
  
  // keep nodes in a vector of Point pointer
  // use node id to access points. The index is the node_id
  vector<pair<Point*, V>> nodes_;
  vector<vector<unsigned>> edges_;
  

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
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
  Graph() {
    // HW0: YOUR CODE HERE
    nodes_ = vector<pair<Point*, V>>(0);
    edges_ = vector<vector<size_type>>(0);
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
      // nodes_ = vector(Point*);

    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      assert(graph_ != NULL);
      // use nodes index to return Point
      return *((graph_ -> nodes_[node_id_]).first);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(this -> graph_ != NULL);
      assert(graph_ -> size() != 0 and graph_ ->size() > node_id_);
      return node_id_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    //node_value_type& value () {
    //  assert(this->graph_ != NULL);
    //  return (*graph_).nodes_[node_id_].second;
    //}

    const node_value_type& value () const {
      assert(this->graph_ != NULL);
      return (*graph_).nodes_[node_id_].second;
    }

    node_value_type& value ()  {
      assert(this->graph_ != NULL);
      return graph_ -> nodes_[node_id_].second;
    }

    size_type degree() const {
       return graph_->edges_[node_id_].size();  

    }

    incident_iterator edge_begin() const {
      return IncidentIterator(graph_, node_id_, 0);
    }

    incident_iterator edge_end() const {
      return IncidentIterator(graph_, node_id_, degree());
    }



    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      
      if (graph_ == n.graph_ and node_id_ == n.node_id_) {
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
      if (graph_ == n.graph_ and node_id_ < n.node_id_) {
        return true;
      } else if (graph_ < n.graph_) {
        return true;
      }
      return false;
      
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    graph_type* graph_;
    size_type node_id_; // which node
    // used in add_node
    Node(const graph_type* graph, size_type node_id) :
      graph_(const_cast<Graph*>(graph)), node_id_(node_id){
    };
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const{
    return nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return nodes_.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    // since position is const, need to create non-const one
    Point* new_Point = new Point;
    *new_Point = position;
    pair<Point*, node_value_type> new_pair(new_Point, node_value_type());
    nodes_.push_back(new_pair);

    // edges
    vector<size_type> tmp;
    edges_.push_back(tmp);
    return Node(this, size()-1);

  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    
    assert(n.graph_ != NULL);
    if (n.graph_ == this and n.node_id_ < size()){
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
    
    assert( i < nodes_.size());
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      
      return Node(graph_, node_id1_);      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      
      return Node(graph_, graph_->edges_[node_id1_][vec_id2_]); 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if (graph_ == e.graph_) {
        if (node1() == e.node1() and node2() == e.node2()) {
          return true;
        } else if (node1() == e.node2() and node2() == e.node1()) {
          return true;
        }
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (graph_ == e.graph_ and edge_id_ < e.edge_id_) {
        return true;
      } else if (graph_ < e.graph_) {
        return true;
      }
      return false;
      
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    const graph_type* graph_;
    
    size_type edge_id_;
    size_type node_id1_;
    size_type vec_id2_;
    // used in add_edge
    Edge(const Graph* graph, size_type node_id1, size_type vec_id2)
        : graph_(const_cast<Graph*>(graph)), node_id1_(node_id1), vec_id2_(vec_id2){
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    size_type numedges = 0;
    for(auto i : edges_){
      numedges += i.size();
    }
    return numedges/2;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < num_edges());
    return *std::next(edge_begin(), i);     
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(this == a.graph_ and this == b.graph_);
    size_type a_id = a.node_id_;
    size_type b_id = b.node_id_;
    assert(a_id < size() and b_id < size());
    size_type deg = a.degree();

    for(size_type i = 0; i < deg; ++i){
      if(edges_[a_id][i] == b_id) return true;
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
    if (has_edge(a,b)) {
      size_type a_id = a.node_id_;
      size_type b_id = b.node_id_;
      size_type deg = a.degree();
      for(size_type i = 0; i < deg; ++i){
        if(edges_[a_id][i] == b_id) {
          return Edge(this, a_id, i);
        }
      }
    }
  
    size_type a_id = a.node_id_;
    size_type b_id = b.node_id_;
    
    //size_type deg = a.degree();

    edges_[a_id].push_back(b_id);
    edges_[b_id].push_back(a_id);
    return Edge(this, a_id, edges_[a_id].size()-1);
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
  }
  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator: private totally_ordered<NodeIterator> {
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
    Node operator*() const  {
      return Node(graph_, node_id_);
    }

    NodeIterator& operator++() {
      ++node_id_;
      return *this;
    }

    bool operator==(const NodeIterator& new_node) const {
      return graph_ == new_node.graph_ and node_id_ == new_node.node_id_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type node_id_;
    NodeIterator(const Graph* graph, size_type node_id)
        : graph_(const_cast<Graph*>(graph)), node_id_(node_id) {
    }
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
  class IncidentIterator : private totally_ordered<IncidentIterator>{
   private:
    friend class Graph;
    graph_type* graph_;
    size_type node_id1_;
    size_type vec_id2_;
    IncidentIterator(const Graph* graph, size_type node_id1, size_type vec_id2)
        : graph_(const_cast<Graph*>(graph)), node_id1_(node_id1), vec_id2_(vec_id2) {
    }

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
      return Edge(graph_, node_id1_, vec_id2_);

    }
    
    IncidentIterator& operator++() {
      ++vec_id2_;
      return *this;
    }
    bool operator==(const IncidentIterator& ii) const {
      if (graph_ == ii.graph_ and node_id1_ == ii.node_id1_ and vec_id2_ == ii.vec_id2_) {
        return true;
      }
      return false;
    }

   
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
  private:
    friend class Graph;
    graph_type* graph_;
    size_type node_id1_;
    size_type vec_id2_;
    EdgeIterator(const Graph* graph, size_type node_id1, size_type vec_id2)
        : graph_(const_cast<Graph*>(graph)), node_id1_(node_id1), vec_id2_(vec_id2) {
    }
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
      return Edge(graph_, node_id1_, vec_id2_);
    }
    EdgeIterator& operator++() {
       ++vec_id2_;
       size_type deg = graph_->edges_[node_id1_].size();
       if (vec_id2_ == deg) {
        ++node_id1_;
        vec_id2_ = 0;
       }
        return *this;
    }
    bool operator==(const EdgeIterator& ii) const {
      if (graph_ == ii.graph_ and node_id1_ == ii.node_id1_ and vec_id2_ == ii.vec_id2_) {
        return true;
      }
      return false;
    }

   
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, size(), 0);
  }


};

#endif // CME212_GRAPH_HPP
