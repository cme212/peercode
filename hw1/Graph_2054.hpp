#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

// Newly added include file.
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @ class Graph
 * @ brief A template for 3D undirected graphs.
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
  /** For HW0, this part is empty. */


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
    /** Do nothing here to construct a empty graph. */
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
      // Doesn't do anything for the invalid constructor.
    }
    

    /** Return this node's position. */
    const Point& position() const {
      // The node needs to be valid.
      if(id >= graph->size()){
        std::cout<<this->id<<" versus "<<graph->size();
      }
      // assert(id < graph->size());
      return graph -> point_list[id];
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id;
    }


    // HW1: YOUR CODE HERE
    /** Return this node's value.*/
    node_value_type &value(){
      assert(id < graph->size());
      return graph -> value_list[id];
    }

    /** Return this node's value. */
    const node_value_type& value() const{
      assert(id < graph->size());
      return graph -> value_list[id];
    }


    /**Set the value
     */
    void setValue(node_value_type new_value){
      assert(id < graph->size());
      graph -> value_list[id] = new_value;
    }
    

    // Supply definitions AND SPECIFICATIONS for:
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    /**@brief Degree of node is defined as the number of
     * edges that connects to the current node.
     *@return The number of edges that incident to the current node.
     */
    size_type degree() const{
      return graph -> node_edge[id].size();
    }


    /**@brief Get the begin iterator for the incident iterator.
     * In my data strucuture, each node maps a list for the edges.
     * The edge list is traversed and the first element is the starting
     * element and the last one is the end element.
     *@return An incident iterator that refers to the begining
     * value. 
     */
    incident_iterator edge_begin() const{
      return IncidentIterator(graph, id, 0);
    }


    /**@brief get the end iterator for the incident iterator.
     *@return An incident iterator that refers to the end value.
     */
    incident_iterator edge_end() const{
      return IncidentIterator(graph, id, this->degree());
    }
    


    /** Test whether this node and _n_ are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->id == n.id && this->graph == n.graph;
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
      return this->id < n.id;
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    Graph *graph;
    size_type id;

    // Constructor only available in graph class.
    Node(const Graph* graph_input, const size_type id_input)
         : graph(const_cast<Graph*>(graph_input)), id(id_input){}
  };


  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return point_list.size();
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
  Node add_node(const Point& position, const node_value_type & v = node_value_type()) {
    size_type num = size();
    Node new_node = Node(this, num);
    point_list.push_back(position);
    value_list.push_back(v);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.id < size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < size());
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
      /* Doesn't do anything for an invalid edge. */
    }

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph, node1_id); 
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph, node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return id == e.id;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return id < e.id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type node1_id;
    size_type node2_id;

    // We define a id for edge like the id for node.
    size_type id;

    // Constructor only available for graph class. 
    Edge(const Graph* graph_input, const size_type node1_id_input, 
      const size_type node2_id_input, const size_type id_input):
        graph(const_cast<Graph*>(graph_input)),
        node1_id(node1_id_input), 
        node2_id(node2_id_input), 
        id(id_input){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_list.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    return edge_list[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(has_node(a) && has_node(b) && ! (a==b));
    // find_edge is another function defined below.
    // It is used to find the id of edge that connects two nodes.
    // The complexity is O(num_edges()).
    size_type i = find_edge(a, b);
    if(i < num_edges())
      return true;
    return false;
  }



  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre _a_ and _b_ are distinct valid nodes of this graph
   * @return an Edge object e with _e_._node1_() == _a_ and _e_._node2_() == _b_
   * @post has_edge(_a_, _b_) == true
   * @post If old has_edge(_a_, _b_), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indexes -- in other words, old edge(_i_) might not
   * equal to new edge(_i_). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    assert(has_node(a) && has_node(b) && ! (a==b));
    size_type num = num_edges();
    size_type id_a = a.id;
    size_type id_b = b.id;
    size_type i = find_edge(a, b);
    if(i == num_edges()){
      Edge new_edge = Edge(this, id_a, id_b, num);
      edge_list.push_back(new_edge);
      node_edge[id_a].push_back(i);
      node_edge[id_b].push_back(i);
      return new_edge;
    }else{
      return edge_list[i];
    }

  }



  /** Given two node and find the id of the edge that connects them.
   * This is a helpful function and
   * @pre _node_a_ and  _node_b_ are different valid nodes
   * @return the id of edge if it exists. Otherwise, return the number of edges.
   * 
   * Complexity: O(num_edges()).
   */
  size_type find_edge(const Node& node_a, const Node& node_b) const{
    // We need to ensure that node_a/b are different and valid.
    assert(has_node(node_a) && has_node(node_b) && !(node_a == node_b));
    if(node_edge.find(node_a.id) != node_edge.end()){
      for(size_type i:node_edge.at(node_a.id)){
        if(edge_list[i].node1().id == node_b.id)
          return i;
      }
    }
    return num_edges();
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    edge_list.clear();
    point_list.clear();
    node_edge.clear();
    value_list.clear();
  }
  
  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator:private totally_ordered<NodeIterator>{

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
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** @brief Get the value that the current iterator refers to.
     * @return A node object. 
     * Complexity: O(1)
     */
    Node operator*() const{
      return Node(graph, id_node);
    }
    

    /**@brief Increment the current iterator.
     *@return A iterator with different value to reference.
     *complexity: O(1)
     */
    NodeIterator &operator++(){
      id_node = id_node + 1;
      return (*this);
    }
 

    /**@brief Check if two iterator are the same.
     *@param[in] The iterator (*this) and it. 
     *@return return true if the two are the same, otherwise, 
     * return false.
     *Complexity: O(1)
     */
    bool operator==(const NodeIterator & it) const{
      return this->graph == it.graph && this->id_node == it.id_node;
    }

   private:
    friend class Graph;
    Graph *graph;
    size_type id_node;

    // Constructor for the graph and id.
    NodeIterator(const Graph *graph_input, const size_type id_input):
                 graph(const_cast<Graph*>(graph_input)), id_node(id_input){}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
  /**@brief Define the start iterator.
   *@return return a node that serves as the start iterator.
   *complexity:O(1)
   */
  node_iterator node_begin() const{
    return NodeIterator(this, 0);
  }
  


  /**@brief Define the end iterator.
   *@return return a node that serves as the end iterator.
   *complexity:O(1)
   */
  node_iterator node_end() const{
    return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator> {
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

    /** @brief Get the value that the current iterator refers to.
     * @return An edge object. 
     * @post The node1 of the returned edge must be equal to the 
     * id_node_ below.
     * Complexity: O(1)
     */
    Edge operator*() const{
      size_type id_curr_edge = graph_ -> node_edge[id_node_][id_edge_];
      Edge e = graph_->edge_list[id_curr_edge];
      if(e.node1_id == id_node_) return e;
      else return Edge(graph_, id_node_, e.node1_id, id_curr_edge);
    }

    /**@brief Increment the iterator.
     *@return An incident iterator object with incremented iterator.
     * complexity: O(1)
     */
    IncidentIterator &operator++(){
      id_edge_++;
      return (*this);
    }

    /**@brief Check if the two incident iterators are the same.
     *@return True if they are the same or false.
     */
    bool operator==(const IncidentIterator& it) const{
      return this->graph_ == it.graph_ && this->id_node_ == it.id_node_
             && this->id_edge_ == it.id_edge_;
    }


   private:
    friend class Graph;
    Graph *graph_;
    // Current node that edges incident to.
    size_type id_node_;
    // Note that each node maps to a list of edges with their ids.
    // id_edge_ refers to the index of current edge in the list.
    // It is not the id of edge.
    size_type id_edge_;
    IncidentIterator(const Graph* graph, size_type id_node, size_type id_edge):
                    graph_(const_cast<Graph*>(graph)), 
                    id_node_(id_node), 
                    id_edge_(id_edge){}
  };


  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
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
    /**@brief Get the current edge
     *@return Get the edge that the current edge iterator refers to
     * Complexity: O(1)
     */
    Edge operator*() const{
      return graph_->edge_list[id_edge_];
    }


    /**@brief Increment the current edge iterator.
     *@return Get the edge iterator that refers to the next value.
     * Complexity: O(1)
     */
    EdgeIterator &operator ++(){
      id_edge_++;
      return (*this);
    }
    

    /**@brief Check if two Edge iterators are the same.
     *@return True if the same or false.
     */
    bool operator == (const EdgeIterator & it) const{
      return this -> graph_ == it.graph_ 
             && this->id_edge_ == it.id_edge_;
    }


   private:
    friend class Graph;
    Graph *graph_;
    size_type id_edge_;
    // Constructor
    EdgeIterator(const Graph* graph, size_type id_edge)
                : graph_(const_cast<Graph*>(graph)), 
                  id_edge_(id_edge){}
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  /**@brief In the data structure, edges are stored with unique id
   * in a list. The first element in the list is used as the starting
   * element and the last one is the end element.
   * @return Starting edge iterator.
   * Complexity: O(1)
   */
  edge_iterator edge_begin() const{
    return EdgeIterator(this, 0);
  }


  /**brief Return the end edge iterator.
   *@return End edge iterator.
   *Complexity: O(1)
   */
  edge_iterator edge_end() const{
    return EdgeIterator(this, this->edge_list.size());
  }

 private:
  /** Here, we define the member variables of graph.
   * There are multiple nodes and edges in one graph. 
   * The container here should faciliate the functionality
   * to add/remove node/edges. Two vectors are used to store the
   * points and edges. Note that the edge has attributes of the two
   * nodes information. To relate node with the edge. We define another
   * map to associate the node with the edges that use the node as one end.
   * This map is useful for other functions, like finding the edge 
   * between two nodes, removing/adding node/edge.
   * For the point_list and edge_list, the ith elements represent the 
   * positions of ith node and edge, respectively.
   */
  std::vector <Point> point_list;
  std::vector <edge_type> edge_list;

  /** Added for homework 1.
   * Note that value_list stores the information of value at each node.
   * Also, ith element represents the value of ith node in value_list.
   */
  std::vector <node_value_type> value_list; 

  /** The map node_edge has the key as the index of node. The corresponding 
   * item of the key is the index of edges that connect the node. This helps 
   * to lower the time complexity of many operations later on.
   */
  std::unordered_map<size_type, std::vector<size_type>> node_edge;
};

#endif // CME212_GRAPH_HPP
