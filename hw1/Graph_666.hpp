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

public:

  
  /** Predeclaration of Node type. */
  class Node;
  
  using node_value_type = V;
  node_value_type& value ();
  const node_value_type& value () const;

  
  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

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

  /** Type of unordered map of node-edge pairs */
  using Pair = std::unordered_map<size_type, size_type>;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    Edges = {};
    Nodes = {};
    node_edge_pairs = {};
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
      this->graph_node = nullptr;
      this->node_num = -1;
    }

    Node(const Graph* graph_node_, size_type node_num_) {
      this->graph_node = graph_node_;
      this->node_num = node_num_;
    }
    

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      const Point* point = graph_node->Nodes[index()];
      return *point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_num;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      if ((index() == n.index())&&(graph_node == n.graph_node)) {return true;}
      else {return false;}
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
      if (n.index() < index()) {return true;}
      else {return false;}
    }

  private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    size_type node_num;
    const Graph* graph_node;
  }; // End node

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */

  size_type size() const {
    return Nodes.size();
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
  Node add_node(const Point& position,				\
		const node_value_type& = node_value_type ()) {
    // HW0: YOUR CODE HERE
    size_type node_total = this->num_nodes();
    Node* new_node = new Node(this, node_total);
    node_edge_pairs[node_total] = {};
    Point* position_ptr = new Point(position);
    Nodes.push_back(position_ptr);
    return *new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if ((n.graph_node == this) && (n.index() < num_nodes())) {return true;}
    else {return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    Node* node = new Node(this, i);
    return *node;
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
  class Edge : totally_ordered<Edge> {
  public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      node_1 = Node();
      node_2 = Node();
    }

    /** Return a node of this Edge */
    Node node1() const {
      return node_1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {

      if ((e.node1() == node1() && e.node2() == node2())
	  || (e.node1() == node2() && e.node2() == node1()))
	{return true;}
      else {return false;}
    }
    

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    bool operator<(const Edge& e) const {
      if (e.node1() < node1() || \
	  (e.node1() == e.node1() && e.node2() < e.node2()))
	{return true;}
      else {return false;}
    }


    
  private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Node node_1;
    Node node_2;

    void set_nodes(const Node a, const Node b) {
      node_1 = a;
      node_2 = b;
    }
    
  }; // End edge

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return Edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return *Edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    Pair node_edge_pair = get_node_edges(a.node_num); 
    return node_edge_pair.find(b.node_num) != node_edge_pair.end();
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
      Pair a_edge = get_node_edges(a.node_num);
      size_type node_num_ = a_edge.at(b.node_num);
      return edge(node_num_);
    }
    else {
      Edge* new_edge = new Edge();
      new_edge->set_nodes(a,b);
      add_node_edge(a.node_num, b.node_num, num_edges());
      Edges.push_back(new_edge);
      return *new_edge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */

  
  void clear() {
    // HW0: YOUR CODE HERE
    Nodes.clear();
    Edges.clear();
    node_edge_pairs.clear();
    
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    NodeIterator(Graph graph_, size_type node_num_) {
      graph_ptr = &graph_;
      node_num = node_num_;
    }
    
    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /**
     * @brief Return the Node that the NodeIterator currently represents
     * @return Node of the current graph that has a node_num of node_num
     * @pre The Nodes vector of graph_ptr has a Node whose node_num is node_num
     */

    Node operator*() const {
      return Node(graph_ptr, node_num);
    }
    
    /**
     * @brief Iterate the NodeIterator class
     * @return NodeIterator++ Returns NodeIterator whose node_num has been 
     *  incremented by 1
     * @post The returned NodeIterator has a node_num that is greater than the 
     *  argument NodeIterator's node_num by 1
     */    
    
    NodeIterator operator++() {
      size_type node_num_plus = node_num+1;
      Graph graph_ = *graph_ptr;
      return NodeIterator(graph_, node_num_plus);
    }

    /**
     * @brief Check if argument NodeIterator and current NodeIterator are 
     *  identical
     * @return Return a boolean that is true if the NodeIterator's parameters 
     *  are equal, and false otherwise
     * @post The argument NodeIterator and current NodeIterator have been 
     *  compared
     */

    bool operator==(const NodeIterator& node_iter) const {
      return ((node_iter.node_num == node_num)	\
	      && (node_iter.graph_ptr == graph_ptr));
    }

  private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_ptr;
    size_type node_num;
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /**
   * @brief Create node iterator that represents beginning of list of nodes
   * @return NodeIterator at the beginning
   * @pre Nodes are defined in the graph
   */
  
  NodeIterator node_begin() const {
    const Graph graph_ = this;
    size_type node_num_ = 0;
    return NodeIterator(graph_, node_num_);
  }

  /**
   * @brief Create node iterator that represents end of list of nodes
   * @return NodeIterator at the end
   * @pre Nodes are defined in the graph
   */
  
  NodeIterator node_end() const {
    const Graph graph_ = this;
    size_type node_num_ = num_nodes();
    return NodeIterator(graph_, node_num_);
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

    IncidentIterator(Graph graph_, size_type node_num_, \
		     Pair::iterator map_iter_) {
      graph_ptr = &graph_;
      node_num = node_num_;
      map_iter = map_iter_;
      edge_num = map_iter->second;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /**
     * @brief Check current edge and return
     * @return Edge incident with current node
     * @pre Node has an edge connected to it
     */

    Edge operator*() const {
      return graph_ptr->edge(edge_num);
    }

    /**
     * @brief Increment IncidentIterator
     * @return By incrementing the map iterator, return the IncidentIterator
     *  with incremented map iterator
     * @pre map_iter should not be at the end of the map
     * @post map_iter iterates one past its current value
     */

    IncidentIterator& operator++() {
      map_iter++;
      return IncidentIterator(*graph_ptr, node_num, map_iter);
    }

    /**
     * @brief Check if current IncidentIterator is equal to argument 
     * IncidentIterator
     * @return Return true if equal, false if not
     */
    
    bool operator==(const IncidentIterator& II) const {
      return ((graph_ptr==II.graph_ptr) && (node_num==II.node_num) && \
	      (edge_num==II.edge_num));
    }

  private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_ptr;
    size_type node_num;
    Pair::iterator map_iter;
    size_type edge_num;
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
    EdgeIterator() {
    }

    EdgeIterator(Graph graph_, size_type edge_num_) {
      graph_ptr = &graph_;
      edge_num = edge_num_;
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /**
     * @brief Return the Edge that the EdgeIterator currently represents
     * @return Edge of the current graph that has a edge_num of edge_num
     * @pre The Edges vector of graph_ptr has a edge whose edge_num is edge_num
     */

    Edge operator*() const {
      return edge(edge_num);
    }
    
    /**
     * @brief Iterate the EdgeIterator class
     * @return EdgeIterator++ returns EdgeIterator whose edge_num has been 
     *  incremented by 1
     * @post The returned EdgeIterator has a edge_num that is greater than the 
     *  argument EdgeIterator's edge_num by 1
     */    
    
    EdgeIterator operator++() {
      size_type edge_num_plus = edge_num+1;
      Graph graph_ = *graph_ptr;
      return EdgeIterator(graph_, edge_num_plus);
    }

    /**
     * @brief Check if argument EdgeIterator and current EdgeIterator are 
     *  identical
     * @return Return a boolean that is true if the EdgeIterator's parameters 
     *  are equal, and false otherwise
     * @post The argument EdgeIterator and current EdgeIterator have been 
     *  compared
     */

    bool operator==(const EdgeIterator& edge_iter) const {
      return ((edge_iter.edge_num == edge_num)	\
	      && (edge_iter.graph_ptr == graph_ptr));
    }

  private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_ptr;
    size_type edge_num;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /**
   * @brief Create edge iterator that represents beginning of list of edges
   * @return EdgeIterator at the beginning
   * @pre Edges are defined in the graph
   */

  EdgeIterator edge_begin() const {
    const Graph graph_ = *this;
    size_type edge_num_ = 0;
    return EdgeIterator(graph_, edge_num_);
  }

  /**
   * @brief Create edge iterator that represents end of list of edges
   * @return EdgeIterator at the end
   * @pre Edges are defined in the graph
   */
  
  EdgeIterator edge_end() const {
    const Graph graph_ = *this;
    size_type edge_num_ = num_edges();
    return EdgeIterator(graph_, edge_num_);
  }

private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  std::vector<Point*> Nodes;
  std::vector<Edge*> Edges;

  // Stores nodes and edges connected to each node
  std::unordered_map<size_type, Pair> node_edge_pairs;

  Pair get_node_edges(size_type node_num_) const {
    return node_edge_pairs.at(node_num_);
  }

  void add_node_edge(size_type node_num_1, size_type node_num_2, \
		     size_type edge_num_) {
    node_edge_pairs[node_num_1].insert({node_num_2, edge_num_});
    node_edge_pairs[node_num_2].insert({node_num_1, edge_num_});
  }
  
  
};

#endif // CME212_GRAPH_HPP
