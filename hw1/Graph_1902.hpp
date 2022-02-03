#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <tuple>
#include <unordered_map>
#include <vector>

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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  /* Design principles
  * Nodes stored in a map of indexes to a tuple containing Points,
  * adjacent node indexes, and the value() attribute.
  * (node indexes will not change with insertions or removals)
  * Edges stored in two ways:
  * 1. edge_map is keyed by unordered node pairs for fast lookup of existence
  * and indices
  * 2. reverse_map is keyed by index and returns ordered node pairs for use
  * in Edge methods
  */
  std::unordered_map<unsigned,
  std::tuple<Point, std::vector<unsigned>, V>> node_map;
  std::unordered_map<std::string, std::vector<unsigned>> edge_map;
  std::unordered_map<unsigned, std::tuple<unsigned, unsigned>> reverse_map;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;
  using node_map_type = std::unordered_map<unsigned,
  std::tuple<Point, std::vector<unsigned>, V>>;

  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;
  using edge_map_type = std::unordered_map<std::string, std::vector<unsigned>>;

  /** Type of node iterators, which iterate over all graph nodes. */
  class NodeIterator;
  /** Synonym for NodeIterator */
  using node_iterator = NodeIterator;
  using node_map_iterator = typename node_map_type::const_iterator;

  /** Type of edge iterators, which iterate over all graph edges. */
  class EdgeIterator;
  /** Synonym for EdgeIterator */
  using edge_iterator = EdgeIterator;
  using edge_map_iterator = edge_map_type::const_iterator;

  /** Type of incident iterators, which iterate incident edges to a node. */
  class IncidentIterator;
  /** Synonym for IncidentIterator */
  using incident_iterator = IncidentIterator;
  using node_vector_type = std::vector<unsigned int>;
  using node_vector_iterator = node_vector_type::const_iterator;

  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;

  /** Type of node contents */
  using node_value_type = V;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph using default constructors for maps. */
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
    Node() : graph_(nullptr), node_index_(0){
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return std::get<0>(graph_->node_map.at(node_index_));
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_index_;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    /** Return this node's type to be modified. */
    node_value_type& value() {
      return std::get<2>(graph_->node_map.at(node_index_));
    };

    /** Return this node's type to be accessed, not modified. */
    const node_value_type& value() const {
      return std::get<2>(graph_->node_map.at(node_index_));
    };

    /** Return number of adjacent nodes */
    size_type degree() const
    {
      if(graph_){
        return std::get<1>(graph_->node_map.at(node_index_)).size();
      }else{
        return 0;
      }
    };

    /** Return iterator pointing at beginning of adjacent node vector */
    incident_iterator edge_begin() const
    {
      unsigned node_degree = std::get<1>(graph_->
        node_map.at(node_index_)).size();
      node_vector_iterator it = std::get<1>(graph_->
        node_map.at(node_index_)).begin();
      return IncidentIterator(graph_, it, node_degree);
    };

    /** Return iterator pointing at end of adjacent node vector */
    incident_iterator edge_end() const
    {
      unsigned node_degree = std::get<1>(graph_->
        node_map.at(node_index_)).size();
      node_vector_iterator it = std::get<1>(graph_->
        node_map.at(node_index_)).end();
      return IncidentIterator(graph_, it, node_degree);
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if(graph_){
        return (node_index_ == n.node_index_ && graph_ == n.graph_);
      }else if(!graph_ && !n.graph_){ //two invalid nodes
        return true;
      }else{ //one valid, one invalid
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
      (void) n;           // Quiet compiler warning

      /* first order by graph pointers (ordered by declaration order)
         then within the same graph, compare node indices
      */
      if(graph_ < n.graph_){
        return true;
      }else if(graph_ == n.graph_){
        return node_index_ >= n.node_index_ ? false : true;
      }else{
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

    // Pointer back to the Graph container
    Graph* graph_;
    // This node's index
    size_type node_index_;
    /** Private Constructor */
    Node(const Graph* graph, size_type node_index)
        : graph_(const_cast<Graph*>(graph)), node_index_(node_index) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_map.size();
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
    const node_value_type& val = node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    size_type new_index = node_map.size();
    std::get<0>(node_map[new_index]) = position;
    std::get<2>(node_map[new_index]) = val;
    return Node(this, new_index);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    (void) n;            // Quiet compiler warning
    return n.graph_ == this && n.node_index_ < this->size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
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
    Edge() : graph_(nullptr), edge_index_(0){
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if(graph_){
        return Node(graph_, std::get<0>(graph_->reverse_map[edge_index_]));
      }else{
        return Node();
      }
    }


    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if(graph_){
        return Node(graph_, std::get<1>(graph_->reverse_map[edge_index_]));
      }else{
        return Node();
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      if(graph_ && e.graph_){ //two valid edges
        return (this->node1() == e.node1() && this->node2() == e.node2()) ||
        (this->node1() == e.node2() && this->node2() == e.node1()) ;
      }else if(!graph_ && !e.graph_){ //two invalid Edges
        return true;
      }else{
        return false; //one valid, one not valid
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //Get min and max node indices for edges
      size_type this_min = this->get_min_index();
      size_type this_max = this->get_max_index();
      size_type e_min = e.get_min_index();
      size_type e_max = e.get_max_index();

      //First order based on graph memory address
      if(graph_ < e.graph_){
        return true;

      }else if(graph_ == e.graph_){

        //If edges are part of the same graph, compare based on min node index
        if(this_min < e_min){
          return true;
        }else if(this_min == e_min){
          //If same min node index, compare based on max node index
          return this_max >= e_max ? false : true;
        }else{
          return false;
        }

      }else{
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
    // This edge's node pair
    size_type edge_index_;
    /** Private Constructor */
    Edge(const Graph* graph, size_type edge_index)
        : graph_(const_cast<Graph*>(graph)), edge_index_(edge_index) {
    }

    /**
    * @brief Helper function to get minimum of node indices associated
    * with an edge
    * @return min(node1.index(), node2.index())
    */
    size_type get_min_index() const {
        return this->node1().node_index_ <= this->node2().node_index_ ?
        this->node1().node_index_ : this->node2().node_index_;
    }
    /**
    * @brief Helper function to get minimum of node indices associated
    * with an edge
    * @return max(node1.index(), node2.index())
    */
    size_type get_max_index() const {
        return this->node1().node_index_ >= this->node2().node_index_ ?
        this->node1().node_index_ : this->node2().node_index_;
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
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
    (void) a; (void) b;   // Quiet compiler warning
    //get unordered edge id key for quick lookups
    std::string edge_id = this->get_edge_id(a.node_index_, b.node_index_);
    return edge_map.count(edge_id);
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
    (void) a, (void) b;   // Quiet compiler warning
    //check if edge exists
    std::string edge_id = this->get_edge_id(a.node_index_, b.node_index_);
    bool has_edge = this->has_edge(a, b);

    //return existing object if edge exists
    if(has_edge){
      //find correct ordered edge in reverse map
      size_type node1, node2;
      size_type ordered_index = 0;
      for(auto d : edge_map[edge_id]){
        node1 = std::get<0>(reverse_map[d]);
        node2 = std::get<1>(reverse_map[d]);
        if(node1 == a.node_index_ && node2 == b.node_index_){
          ordered_index = d;
        }
      }
      return Edge(this, ordered_index);

    }else{
      //otherwise add unique edge to edge map
      size_type edge_index = reverse_map.size();

      //first, add edge(a, b)
      edge_map[edge_id].push_back(edge_index);
      reverse_map[edge_index] = std::tuple<unsigned, unsigned>{a.node_index_,
        b.node_index_};
      std::get<1>(node_map.at(a.node_index_)).push_back(edge_index);

      //next, add edge(b, a)
      edge_map[edge_id].push_back(edge_index + 1);
      reverse_map[edge_index + 1] = std::tuple<unsigned,
      unsigned>{b.node_index_, a.node_index_};
      std::get<1>(node_map.at(b.node_index_)).push_back(edge_index + 1);

      return Edge(this, edge_index);
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    node_map.erase(node_map.begin(), node_map.end());
    edge_map.erase(edge_map.begin(), edge_map.end());
    reverse_map.erase(reverse_map.begin(), reverse_map.end());
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), it_(nullptr), size_(0) {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const
    {
      if(size_ == 0){
        return Node();
      }else{
        return Node(graph_, (*it_).first);
      }
    }
    NodeIterator& operator++()
    {
      ++it_;
      return *this;
    }
    bool operator==(const NodeIterator& node_iter) const
    {
      return it_ == node_iter.it_;
    }
    bool operator!=(const NodeIterator& node_iter) const
    {
      return it_ != node_iter.it_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    // Pointer back to the Graph container
    Graph* graph_;
    node_map_iterator it_;
    unsigned size_;

    //Private constructor that can be accessed by Graph class
    NodeIterator(const Graph* graph, node_map_iterator it, unsigned size) :
    graph_(const_cast<Graph*>(graph)), it_(it), size_(size) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const
  {
    return NodeIterator(this, node_map.begin(), this->num_nodes());
  }
  node_iterator node_end() const
  {
    return NodeIterator(this, node_map.end(), this->num_nodes());
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
    IncidentIterator() : graph_(nullptr), it_(nullptr), size_(0) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const
    {
      if(size_ == 0){
        return Edge();
      }else{
        return Edge(graph_, (*it_));
      }
    }
    IncidentIterator& operator++()
    {
      ++it_;
      return *this;
    }
    bool operator==(const IncidentIterator& inc_iter) const
    {
      return it_ == inc_iter.it_;
    }
    bool operator!=(const IncidentIterator& inc_iter) const
    {
      return it_ != inc_iter.it_;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    node_vector_iterator it_;
    unsigned size_;

    //Private constructor that can be accessed by Graph class
    IncidentIterator(const Graph* graph, node_vector_iterator it,
    unsigned size) :
    graph_(const_cast<Graph*>(graph)), it_(it), size_(size) {}

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
    EdgeIterator() : graph_(nullptr), it_(nullptr), size_(0) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const
    {
      if(size_ == 0){
        return Edge();
      }else{
        return Edge(graph_, (*it_).second.at(0));
      }
    }
    EdgeIterator& operator++()
    {
      ++it_;
      return *this;
    }
    bool operator==(const EdgeIterator& edge_iter) const
    {
      return it_ == edge_iter.it_;
    }
    bool operator!=(const EdgeIterator& edge_iter) const
    {
      return it_ != edge_iter.it_;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // Pointer back to the Graph container
    Graph* graph_;
    edge_map_iterator it_;
    unsigned size_;

    //Private constructor that can be accessed by Graph class
    EdgeIterator(const Graph* graph, edge_map_iterator it,
    unsigned size) :
    graph_(const_cast<Graph*>(graph)), it_(it), size_(size) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const
  {
    return EdgeIterator(this, edge_map.begin(), this->num_edges());
  }
  edge_iterator edge_end() const
  {
    return EdgeIterator(this, edge_map.end(), this->num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /**
  * @brief Helper function to create a unique edge id for edge a-b and b-a
  * @param[in] a unsigned int, node1.index()
  * @param[in] b unsigned int, node2.index()
  * @return string formatted as "node[_first_].index()-node[_second_].index()",
  * where _first_ = min(a,b) and _second_ = max(a,b)
  *
  * @pre a and b are indices of valid nodes
  */
  std::string get_edge_id(const size_type a, const size_type b) const{
    std::vector<size_type> id_vec = {a, b};
    std::sort(id_vec.begin(), id_vec.end());
    return std::to_string(id_vec[0]) + "-" + std::to_string(id_vec[1]);
  }

};

#endif // CME212_GRAPH_HPP
