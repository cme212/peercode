#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <iterator>
#include <map>

#include <unordered_map>

#include <tuple>
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

  struct internal_element {
    Point pt;
    V v;

    internal_element(Point point, V value){
      this->pt = Point(point.x, point.y, point.z);
      this->v = value;
    }

    internal_element(Point point){
      this->pt = Point(point.x, point.y, point.z);
      this->v = 0;
    }
  };




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

  // HW1 added
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

  
  // additional data structure
  std::vector<internal_element*> elements; // store the nodes
  std::vector<Node*> node_list; // store the node
  std::vector<Edge*> edge_list; // store the edge =
  std::map<std::pair<unsigned, unsigned>, unsigned> edges_map;// store edges
  std::vector<std::pair<unsigned, unsigned>> edges_vector; // store the edges
  unsigned graph_size; // number of nodes
  unsigned numEdges; // number of edges
  
  // store the edges connected to a node
  std::map<unsigned, std::vector<Edge*>> node_edge_map; 

  
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() 
  :elements(), node_list(),edge_list(), edges_map(), edges_vector(), 
  node_edge_map(){
    graph_size = 0;
    numEdges = 0;
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
  class Node : private totally_ordered<Node>{
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

    Node(const Graph* new_graph, size_type new_idx) 
    : graph(const_cast<Graph*>(new_graph)), idx(new_idx) {
    }

    /** Return this node's position. */
    const Point& position() const{
      assert(idx < graph->graph_size);
      return graph->elements[idx]->pt;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;


    // HW1 Implementation
    node_value_type& value() {
      assert(idx < graph->graph_size);
      return graph->elements[idx]->v; 
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      assert(idx < graph->graph_size);
      return const_cast<node_value_type&>(graph->elements[idx]->v);
    }

    /** Return the number of edges connected this node */
    size_type degree() const{
      std::vector<Edge*> incident = graph->node_edge_map[this->idx];
      return incident.size(); // to be changed later
    }

    /** Return the begin of IncidentIterator of a given node. */
    IncidentIterator edge_begin() const{
      return IncidentIterator(graph, this->idx,  0);

    }

    /** Return the end of IncidentIterator of a given node. */
    IncidentIterator edge_end() const{
      std::vector<Edge*> incident = graph->node_edge_map[this->idx];
      return IncidentIterator(graph, this->idx,  degree());

    }
   

    /** Return the graph of the given node */
    Graph node_graph() const{
      return *(this->graph);
    }



    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // check if index is the same
      if (n.idx != this->idx)
        return false;
      // check if in same graph
      if (n.graph != this->graph)
        return false;
      return true;

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
      if(this->idx < n.idx)
        return true;
      return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
    Graph* graph;
    size_type idx;
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return this->graph_size;
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
  Node add_node(const Point& position) {
    internal_element* add_element = new internal_element(position);
    elements.push_back(add_element); // add to elements
    graph_size++; // increment graph_size by 1
    Node* added = new Node(this, size()-1);
    node_list.push_back(added);
    return *added; 
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& value){
    internal_element* add_element = new internal_element(position, value);
    elements.push_back(add_element); // add to elements
    graph_size++; // increment graph_size by 1
    Node* added = new Node(this, size()-1);
    node_list.push_back(added);
    return *added; 
  }


  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph != this)
      return false;
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < this->size());
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
    }


   
    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, node1_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, node2_idx);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      // check the same direction
      if ((e.node1() == this->node1()) && (e.node2() == this->node2()))
        return true;

      // check the opposite direction
      if ((e.node1() == this->node2()) && (e.node2() == this->node1()))
        return true;
      return false;
    }

    size_type id() const{
      return this->uid;
    }
   
    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return (!((*this) == e) && (this->uid < e.uid));
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    size_type node1_idx;
    size_type node2_idx;
    size_type uid;
    Graph* graph;

    Edge(size_type n1, size_type n2, size_type uid, const Graph* new_graph){
      this->node1_idx = n1;
      this->node2_idx = n2;
      this->uid = uid;
      this->graph = const_cast<Graph*>(new_graph);
    }





  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return numEdges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < this->num_edges());
    return Edge(edges_vector[i].first, edges_vector[i].second, i, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // modified version after reviewing the peer code
    auto tmp_1 = edges_map.find(std::make_pair(a.index(), b.index()));
    auto tmp_2 = edges_map.find(std::make_pair(b.index(), a.index())); 

    if ((tmp_1 != edges_map.end()) || (tmp_2 != edges_map.end())) {
      return true;
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

     // modified version after reviewing the peer code
    std::map<std::pair<unsigned, unsigned>, unsigned>::iterator tmp_1 = edges_map.find(std::make_pair(a.index(), b.index()));
    std::map<std::pair<unsigned, unsigned>, unsigned>::iterator tmp_2 = edges_map.find(std::make_pair(b.index(), a.index()));

    if (tmp_1 != edges_map.end())
      return Edge(a.index(), b.index(), edges_map[std::make_pair(a.index(), b.index())] , this);

    if (tmp_2 != edges_map.end())
      return Edge(b.index(), a.index(), edges_map[std::make_pair(b.index(), a.index())] , this);

    // add the new edge to the graph
    Edge* add_edge = new Edge(a.index(), b.index(), numEdges, this);
    Edge* add_edge_2 = new Edge(b.index(), a.index(), numEdges, this);
    std::pair<size_type, size_type> pair_idx = std::make_pair(a.idx, b.idx);
    edges_map[pair_idx] = numEdges;
    edges_vector.push_back(pair_idx);


    // add the new edge to the node_edge_map
    node_edge_map[a.index()].push_back(add_edge);
    node_edge_map[b.index()].push_back(add_edge_2);

    // add the new edge to the edge_list
    edge_list.push_back(add_edge);

    this->numEdges++;

    return *add_edge;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    elements.clear();
    edges_map.clear();
    edges_vector.clear();
    node_list.clear();
    node_edge_map.clear();
    edge_list.clear();
    numEdges = 0;
    graph_size = 0;
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

    /** Construct an valid NodeIterator. */
    NodeIterator(const Graph* graph, unsigned idx) 
    : m_graph(const_cast<Graph*>(graph)), m_idx(idx) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    
    /** Construct the dereferencing. */
    Node operator*() const {
      //return *(m_graph->node_list[m_idx]);
      return m_graph->node(m_idx);
    }

    /** Add the increment function*/
    NodeIterator& operator++(){
      m_idx++;
      return *this;
    }

    /** Add the comparator. */
    bool operator==(const NodeIterator& itr) const {
      return ((itr.m_idx == this->m_idx) && (itr.m_graph == this->m_graph));
    }

    /** Add the comparator. */
    bool operator!=(const NodeIterator& itr) const {
      return ((itr.m_idx != this->m_idx) || (itr.m_graph != this->m_graph));
    }

   private:
    friend class Graph;
    Graph* m_graph;
    unsigned m_idx;

    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Begin of the node iterator. */
  NodeIterator node_begin() const{
    return NodeIterator(this, 0);
  }

  /** End of the node iterator */
  NodeIterator node_end() const {
    return NodeIterator(this, this->size());
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
    using iterator_category = std::forward_iterator_tag;// Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    /** Construct the IncidentIterator. */
    IncidentIterator(const graph_type* _graph, size_type _node_idx,
     size_type _edge_idx) 
    : graph(const_cast<Graph*>(_graph)), node_idx(_node_idx),
     edge_idx(_edge_idx) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // bool operator==(const IncidentIterator
    // IncidentIterator& operator++()&) const

    /** Dereferencing operator. */
    Edge operator*() const {
      return *(graph->node_edge_map[node_idx][edge_idx]);
    } 

    /** Increment operator */
    IncidentIterator& operator++() {
      edge_idx++;
      return *this;
    }

    /** Comparator operator */
    bool operator==(const IncidentIterator& itr2) const{
      return ((this->graph == itr2.graph) && 
        (this->node_idx == itr2.node_idx) && 
        (this->edge_idx == itr2.edge_idx));
    }

    /** Comparator operator */
    bool operator!=(const IncidentIterator& itr2) const{
      return ((this->graph != itr2.graph) || 
        (this->node_idx != itr2.node_idx) || 
        (this->edge_idx != itr2.edge_idx));    
    }
    

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    graph_type* graph;
    size_type node_idx; // node index
    size_type edge_idx; // edge index of the edges connected to the given node
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

    /** Construct an valid EdgeIterator. */
    EdgeIterator(const graph_type* graph, unsigned idx): 
    m_graph(const_cast<Graph*>(graph)), m_idx(idx) {
    }
    
     /** Add the dereferencing operator. */
    Edge operator*() const{
      return *(this->m_graph->edge_list[m_idx]);
    }

    /**  Add the increment operator. */
    EdgeIterator& operator++(){
      m_idx++;
      return *this;
    }

    /** Add the comparator operator. */
    bool operator==(const EdgeIterator& edge_itr) const {
      return ((this->m_idx == edge_itr.m_idx) && 
        (this->m_graph == edge_itr.m_graph));
    }

    /** Add the comparator operator. */
    bool operator!=(const EdgeIterator& edge_itr) const {
      return  ((this->m_idx != edge_itr.m_idx) || 
        (this->m_graph != edge_itr.m_graph));
    }




    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    graph_type* m_graph;
    unsigned m_idx;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Begin of the EdgeIterator */
  EdgeIterator edge_begin() const{
    return EdgeIterator(this, 0);
  }
  
  /** End of the EdgeIterator. */
  EdgeIterator edge_end() const{
    return EdgeIterator(this, num_edges()); 
  }

 private:
  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
