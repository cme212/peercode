#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <list>
#include <unordered_map>
#include <unordered_set> 
#include <tuple>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = short>
class Graph {
 private:
   /** Forward declation of InternalNode type**/
  struct InternalNode; 

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
  using node_value_type = V; 

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
  Graph() {}

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
    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // verify that the node is valid (i.e. belongs to an existing graph)
      assert(graph != nullptr and *graph != nullptr);

      return fetch().position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // if a node is invalid (i.e. does not belong to an active graph object),
      // it does not have a notion of an index
      assert(graph != nullptr and *graph != nullptr);

      // if a node is valid, its id can be interpreted as an index by which it
      // can be accessed from the associated graph, as graph.node(index)
      return id; 
    }

    /** Return an l-value reference to the node's associated value. */
    node_value_type& value() { 
      assert(graph != nullptr and *graph != nullptr);

      return fetch().value;
    }

    /** Returns an r-value reference to the node's associated value */
    const node_value_type& value() const { 
      assert(graph != nullptr and *graph != nullptr);

      return fetch().value;
    }

    /** Return the degree of the node, defined as the number of nodes adjacent
    to the current one, i.e. connected to the current node through an edge. */
    size_type degree() const {
      assert(graph != nullptr and *graph != nullptr); 

      const auto& adjacent_nodes = (**graph).connected_nodes.at(id);
      return adjacent_nodes.size(); 
    }

    /** Return an iterator to the first element of the collection of nodes
    adjacent to the current node. */
    incident_iterator edge_begin() const { 
      const auto& adjacent_nodes = (**graph).connected_nodes.at(id);

      return IncidentIterator(*this, adjacent_nodes.begin()); 
    }

    /** Return an iterator to the element one past the end of the collection of
    nodes adjacent to the current node. Attempting to dereference this iterator
    results in undefined behavior. */
    incident_iterator edge_end() const {
      const auto& adjacent_nodes = (**graph).connected_nodes.at(id);

      return IncidentIterator(*this, adjacent_nodes.end()); 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
     
      // if either or both nodes are invalid (i.e. do not belong to an existing
      // graph object, the nodes compare false
      if (this->graph == nullptr or n.graph == nullptr)
          return false; 
      if (*graph == nullptr or *(n.graph) == nullptr)
          return false; 
      
      // valid nodes must belong to the same graph and be found at the same
      // index in order to compare true
      if (*graph != *(n.graph))
        return false; 

      return (this->id == n.id);  
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
      // if both nodes are valid and belong to the same graph, the node with
      // the smaller index compares lower. If one or both nodes are invalid, 
      // or belong to different graphs, then the comparison is still made 
      // based on node id, but this id does not have any significance w.r.t. 
      // any graph object (i.e. cannot be interpreted as an index)      

      return id < n.id; 

   }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    friend class NodeIterator; 
    friend class IncidentIterator; 

    graph_type** graph = nullptr; 
    size_type id; 

    /** Private member function that returns the InternalNode associated with
     *  the current node. 
     * 
     * @pre the node is valid (i.e. graph != nullptr and *graph != nullptr)
     * @return an l-value reference to the InternalNode object containing the
     *         data associated with the current node
     */
    InternalNode& fetch() const {
      auto& nodes = (**graph).nodes; 
      return nodes[id]; 
    } 

    /** Private constructor for use by the Graph class and associated 
     *  iterators.
    */
    Node(graph_type** graph_, size_type id): graph(graph_), id(id) {}
   
 };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodes.size();
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
                const node_value_type& value = node_value_type()) {

    nodes.emplace_back(position, next_id, value);
   
    // add a new entry to the adjacency map to store any nodes to be connected
    // to the newly added node in the future
    connected_nodes[next_id]; 
 
    return Node(const_cast<graph_type**>(&self_ptrs.back()), next_id++);
 }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph == nullptr) {
      return false; 
    }

    return (this == *(n.graph));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < nodes.size()); 
    
    return Node(const_cast<graph_type**>(&self_ptrs.back()), nodes[i].id);
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
    Edge() {}

    /** Return a node of this Edge */
    Node node1() const {
      return node_a; 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node_b; 
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if ((this->node_a == e.node_a) and (this->node_b == e.node_b)) 
        return true; 

      return ((this->node_a == e.node_b) and (this->node_b == e.node_a)); 
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {

      // by trichotomy, if two edges compare equal, then e1 < e2 is false
      if (*this == e)
        return false;
 
      // if two edges are not equal, then at least one pair of nodes, with
      // one node from e1 and one node from e2, will not compare equal.               // Define comparison between edges as a comparison of the first 
      // encountered pair of unequal nodes. 
      if (node_a != e.node_a) 
        return (node_a < e.node_a); 

      return (node_b < e.node_b); 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator; 

    Node node_a; 
    Node node_b;

    /** Private constructor for use by the graph class and associated 
        iterators.*/
    Edge(const node_type& node_a, const node_type& node_b) : node_a(node_a), 
                                                             node_b(node_b) {}

 };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < edges.size()); 

    return edges[i]; 
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1)
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) and has_node(b)); 

    return connected_nodes.at(a.id).count(b.id); 
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
   * Complexity: average O(1), worse case O(num_nodes())
   */
  Edge add_edge(const Node& a, const Node& b) {
    // verify preconditions
    assert(has_node(a) and has_node(b)); 
    assert(a != b); 

    Edge edge(a, b); 

    if (not has_edge(a, b)) {
      edges.push_back(edge); 

      connected_nodes[a.id].insert(b.id); 
      connected_nodes[b.id].insert(a.id); 
    } 

    return edge; 
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges.clear(); 
    connected_nodes.clear(); 

    // invalidate the graph pointer associated with all former Node objects
    // released by this Graph instance
    self_ptrs.back() = nullptr; 

    // add a new valid graph pointer to be released to future, valid Nodes
    // constructed by the graph
    self_ptrs.push_back(this); 
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
    NodeIterator() {}

    /** Return the underlying Node object pointed to by the NodeIterator
     * @pre iterator != graph.node_end()  
     */
    Node operator*() const {
      return Node(graph, internal_iter->id);        
    }

    /** Increment the NodeIterator to point to the next element in the 
     *  container of graph Nodes
     * @pre old iterator != graph.node_end() 
     * @post either new iterator  == graph.node_end() or new iterator can be
     *       safely dereferenced to access the next Node in the graph
     * @return l-value reference to the updated NodeIterator (enables chaining
     *         (of expressions) 
     */
    NodeIterator& operator++() { 
      internal_iter++;
      return *this; 
    }

    /** Check whether two NodeIterator objects compare equal
     * @return true if the two iterators point to the same element in the graph
     *         Node container; otherwise, false
     */
    bool operator==(const NodeIterator& n) const {
      return (internal_iter == n.internal_iter); 
    }

    /** Check whether two NodeIterator objects are unequal 
     * @return !(*this == _n_)  
     */
    bool operator!=(const NodeIterator& n) const {
      return !(*this == n); 
    }

   private:
    friend class Graph;

    using internal_iterator = 
          typename std::vector<InternalNode>::const_iterator;
    
    // an internal iterator to traverse the graph's vector of InternalNodes
    internal_iterator internal_iter; 

    // the graph associated with the NodeIterator
    graph_type** graph; 

    /** Private constructor for use by the Graph class */
    NodeIterator(const internal_iterator& internal_iter, graph_type** graph) :
     internal_iter(internal_iter), graph(graph) {} 
  };

  /** Return an iterator to the first element in the container of Node objects
      associated with the present graph */
  node_iterator node_begin() const {
    return NodeIterator(nodes.begin(), 
                        const_cast<graph_type**>(&self_ptrs.back())); 
  }

  /** Return an iterator to one-past the last element in the container of Node
      objects associated with the present graph. Attempting to derference this
      iterator results in undefined behavior. */
  node_iterator node_end() const {
    return NodeIterator(nodes.end(), 
                        const_cast<graph_type**>(&self_ptrs.back())); 
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

    /** Return the underlying Edge object pointed to by the IncidentIterator,
     * @pre iterator != node.edge_end()  
     * @post _return_.node1() is equal to the Node that spawned the 
     *       IncidentIterator
     */
    Edge operator*() const {
      Node node_2 = Node(node.graph, *adj_node_iter); 
      return Edge(node, node_2); 
    }

    /** Increment the IncidentIterator to point to the next incident Edge
     * @pre old iterator != node.edge_end() 
     * @post either new iterator  == node.edge_end() or new iterator can be
     *       safely dereferenced to return the next incident Edge
     * @return l-value reference to the updated IncidentIterator (enables
     *         (chaining of expressions) 
     */
    IncidentIterator& operator++() {
      adj_node_iter++; 
      return *this; 
    }

    /** Checks whether the two IncidentIterator objects are equal
     * @return true if both iterators point to the same incident edge of the
     *         same node; otherwise, false 
     */
    bool operator==(const IncidentIterator& it) const {
      return (adj_node_iter == it.adj_node_iter); 
    }

    /** Check whether two InicidentIterator objects are unequal
     * @return !(*this == _it_)  
     */
    bool operator!=(const IncidentIterator& it) const {
      return !(*this == it);
    }

   private:
    friend class Node;

    using int_iterator = std::unordered_set<size_type>::const_iterator; 

    // the node that spawned the IncidentIterator 
    Node node;  
  
    // an iterator to traverse the container of nodes adjacent to the spawning 
    // node
    int_iterator adj_node_iter; 
 
    /** Private constructor for use by the Graph class*/
    IncidentIterator(const Node& node, int_iterator iter) : 
     node(node), adj_node_iter(iter) {}
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

    /** Return the underlying Edge object pointed to by the EdgeIterator
     * @pre iterator != graph.edge_end()  
     */
    Edge operator*() const {
      return *internal_iterator; 
    }

    /** Increment the EdgeIterator to point to the next element in the 
     *  container of graph edges
     * @pre old iterator != graph.edge_end() 
     * @post either new iterator  == graph.edge_end() or new iterator can be
     *       safely dereferenced to access the next Edge in the graph
     * @return l-value reference to the updated EdgeIterator (enables chaining
     *         (of expressions) 
     */
 
    EdgeIterator& operator++() {
      internal_iterator++; 
      return *this; 
    }

    /** Check whether two EdgeIterator objects are equal
     * @return true if the two iterators point to the same element in the Edge
     *         container; otherwise, false
     */
    bool operator==(const EdgeIterator& ei) const {
      return (internal_iterator == ei.internal_iterator); 
    }

    /** Check whether two EdgeIterator objects are unequal 
     * @return !(*this == _ei_)  
     */ 
    bool operator!=(const EdgeIterator& ei) const {
      return !(*this == ei); 
    }

   private:
    friend class Graph;

    using int_iter = typename std::vector<Edge>::const_iterator; 

    // iterator to traverse the associated graph's vector of Edge objects
    int_iter internal_iterator; 

    /** Private constructor for use by the Graph class */
    EdgeIterator(const int_iter& iter) : internal_iterator(iter) {}
  };


  /** Return an iterator to the first element in the graph's collection of Edge
      objects */
  edge_iterator edge_begin() const {
    return EdgeIterator(edges.begin()); 
  }

  /** Return an iterator to one past the last element in the graph's collection
      of Edge objects */
  edge_iterator edge_end() const {
    return EdgeIterator(edges.end()); 
  }

 private:
  // internal attributes of the Graph class

  /** graph data: vectors of internal nodes and edges, as well as
      an adjacency map to facilitate O(1) lookups of existing edges */
  std::vector<InternalNode> nodes;  
  std::vector<Edge> edges;
  std::unordered_map<size_type, std::unordered_set<size_type>> connected_nodes;
 
  // index to assign to the next node added to the graph 
  size_type next_id = 0; 

  // pointers to the current graph to be released with Node objects
  std::list<Graph*> self_ptrs {this}; 

  /** InternalNode objects hold the position and value data associated with a
      graph Node. Private type for sole use by the Graph class. 
   */
  struct InternalNode {
    Point position; 
    size_type id; 
    node_value_type value{}; 

    /** InternalNode constructor */
    InternalNode(const Point& position, size_type id, 
                 const node_value_type& value) : position(position), 
                                                 id(id),
                                                 value(value) {} 
 };

};

#endif // CME212_GRAPH_HPP
