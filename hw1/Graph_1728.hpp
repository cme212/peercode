#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
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
  unsigned int size_;
  unsigned int edge_size;
  unsigned int edge_counter;
  struct nodeProxy;
  struct edgeProxy;
  std::vector<nodeProxy> Nodes;
  std::vector<edgeProxy> Edges;
  std::map<std::vector<unsigned int>,unsigned int> Etoi; 
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
  Graph() {
    // HW0: YOUR CODE HERE
    size_ = 0;
    edge_size = 0;
    edge_counter = 0;
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
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return graph->Nodes[ind].pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      return ind;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    size_type degree() const{
      return graph->Nodes[ind].neighbors.size();
    }

    IncidentIterator edge_begin() const{
      return IncidentIterator(ind,graph,0);
    }

    IncidentIterator edge_end() const{
      return IncidentIterator(ind,graph,this->degree());
    }

    const node_value_type& value() const{
      return graph->Nodes[ind].val;
    }

    node_value_type& value(){
      return graph->Nodes[ind].val;
    }
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return this->ind == n.ind && n.graph == this->graph;
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
      return this->ind < n.ind;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    Graph* graph; 
    size_type ind;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Node(const Graph* g, size_type i): graph{const_cast<Graph*>(g)}, ind(i){}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return size_;
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
  Node add_node(const Point& position, const node_value_type& v = node_value_type()) {
    // HW0: YOUR CODE HERE
    std::vector<size_type> neighbor = std::vector<size_type>();
    size_++;
    nodeProxy new_node = nodeProxy(position,size_-1,v,neighbor);
    Nodes.push_back(new_node);
    return Node(this,size_-1); 
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return n.ind < size_;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i < size_ ){
      return Node(this,i);
    }
    else{return Node();}     // Invalid node
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
  class Edge:private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      size_type nodeId = graph->Edges[id].node1;
      return Node(graph,nodeId);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      size_type nodeId = graph->Edges[id].node2;
      return Node(graph,nodeId);     
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return id == e.id;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return id < e.id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type id;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Edge(const Graph* g,size_type i): graph(const_cast<Graph*>(g)), id(i) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_size;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
  
    return Edge(this,2*i);        
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    std::vector<unsigned int> key1;
    key1.push_back(a.ind);
    key1.push_back(b.ind);
    if (Etoi.find(key1) == Etoi.end()){
      return false;
    }
    else{
      return true;
    }
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
    std::vector<unsigned int> key1,key2;
    key1.push_back(a.ind);
    key1.push_back(b.ind);
    key2.push_back(b.ind);
    key2.push_back(a.ind);
    if (has_edge(a,b)){
      return Edge(this,Etoi.find(key1)->second);
    }
    else{
      edge_size++;
      edge_counter += 2;
      edgeProxy new_edge = edgeProxy(a.ind,b.ind,edge_counter-2);
      edgeProxy new_edge_reverse = edgeProxy(b.ind,a.ind,edge_counter-1);
      Nodes[a.ind].neighbors.push_back(b.ind);
      Nodes[b.ind].neighbors.push_back(a.ind);
      Edges.push_back(new_edge);
      Edges.push_back(new_edge_reverse);
      Etoi[key1] = edge_counter-2;
      Etoi[key2] = edge_counter-1; 
    }
    return Edge(this,edge_counter-2);        
    
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    Edges.clear();
    Nodes.clear();
    Etoi.clear();
    size_ = 0;
    edge_size = 0;
    edge_counter = 0;
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const
    NodeIterator& operator++(){
      if (n_id  == this->graph->size()) {
        return *this;
      }
      else{
        n_id++;
        return *this;
      }
    }

    Node operator*() const{
      return this->graph->node(n_id);
    }

    bool operator==(const NodeIterator& n_iter) const{
      return n_id == n_iter.n_id;
    }

    bool operator!=(const NodeIterator& n_iter) const{
      return n_id != n_iter.n_id;
    }
    
    bool operator<=(const NodeIterator& n_iter) const {
      return n_id <= n_iter.n_id;
}

   private:
    friend class Graph;

    size_type n_id;
    Graph* graph;
    NodeIterator(size_type i,const Graph* g): n_id(i), graph(const_cast<Graph*>(g)) {};
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  NodeIterator node_begin() const{
    return NodeIterator(0,this);
  }

  NodeIterator node_end() const{
    return NodeIterator(this->Nodes.size(),this);
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

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const
    Edge operator*() const{
      std::vector<unsigned int> key;
      key.push_back(a_id);
      key.push_back(graph->Nodes[a_id].neighbors[current]);
      return Edge(graph,graph->Etoi.find(key)->second);
    }

    IncidentIterator& operator++(){
      current++;
      return *this;
    }

    bool operator==(const IncidentIterator& ii) const{
      return this->current == ii.current && a_id == ii.a_id;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    size_type a_id;
    Graph* graph;
    size_type current;
    IncidentIterator(size_type a, const Graph* g, size_type c): a_id(a),graph(const_cast<Graph*>(g)), current(c) {}
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

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    Edge operator*() const{
      size_type e_id = graph->Edges[current].ind;
      return Edge(graph,e_id);
    }

    EdgeIterator& operator++(){
      if (current > graph->Edges.size()){
	return *this;
      }
      else{
        current += 2;
        return *this;
      }
    }

    bool operator==(const EdgeIterator& eI) const{
      return current == eI.current;
    }

    bool operator!=(const EdgeIterator& eI) const{
      return current != eI.current;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    size_type current;
    Graph* graph;

    EdgeIterator(size_type cid, const Graph* g): current(cid), graph(const_cast<Graph*>(g)) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  EdgeIterator edge_begin() const{
    return EdgeIterator(0,this);
  }

  EdgeIterator edge_end() const {
    return EdgeIterator(this->Edges.size(),this);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  struct nodeProxy{
    Point pt;
    size_type ind; 
    node_value_type val;
    std::vector<size_type> neighbors;
    nodeProxy(Point p,size_type i, node_value_type v, std::vector<size_type> nei){
      pt = p;
      ind = i;
      val = v;
      neighbors = nei;
    }
  };

  struct edgeProxy{
    size_type node1;
    size_type node2;
    size_type ind;
    edgeProxy(size_type n1, size_type n2, size_type id){
      node1 = n1;
      node2 = n2;
      ind = id; 
    }
  };
  
};

#endif // CME212_GRAPH_HPP
