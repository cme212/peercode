#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <list>
#include <map>
#include <unordered_map>
#include <unordered_set> 
#include <utility>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = short, typename E = short>
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
  using edge_value_type = E; 

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

    /** Return this node's position as an l-value reference. */
    Point& position() {
      // verify that the node is valid (i.e. belongs to an existing graph)
      assert(graph != nullptr and graph->active_nodes.count(id)); 

      return graph->node_positions.at(id);
    }

    /** Return this node's position as an r-value reference. */
    const Point& position() const {
      assert(graph != nullptr and graph->active_nodes.count(id)); 

      return graph->node_positions.at(id);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {

      // if a node is invalid (i.e. does not belong to an active graph object),
      // it does not have a notion of an index
      assert(graph != nullptr and graph->active_nodes.count(id)); 
      
      return graph->node_indices.at(id); 

    }

    /** Return an l-value reference to the node's associated value. */
    node_value_type& value() { 
      assert(graph != nullptr and graph->active_nodes.count(id)); 

      return graph->node_values.at(id);
    }

    /** Returns an r-value reference to the node's associated value */
    const node_value_type& value() const { 
      assert(graph != nullptr and graph->active_nodes.count(id)); 

      return graph->node_values.at(id);
    }

    /** Return the degree of the node, defined as the number of nodes adjacent
    to the current one, i.e. connected to the current node through an edge. */
    size_type degree() const {
      assert(graph != nullptr and graph->active_nodes.count(id)); 

      const auto& adjacent_nodes = (*graph).connected_nodes.at(id);
      return adjacent_nodes.size(); 
    }

    /** Return an iterator to the first element of the collection of nodes
    adjacent to the current node. */
    incident_iterator edge_begin() const { 
      const auto& adjacent_nodes = (*graph).connected_nodes.at(id);

      return IncidentIterator(*this, adjacent_nodes.begin()); 
    }

    /** Return an iterator to the element one past the end of the collection of
    nodes adjacent to the current node. Attempting to dereference this iterator
    results in undefined behavior. */
    incident_iterator edge_end() const {
      const auto& adjacent_nodes = (*graph).connected_nodes.at(id);

      return IncidentIterator(*this, adjacent_nodes.end()); 
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
     
      // if either or both nodes are invalid (i.e. do not belong to an existing
      // graph object), the nodes compare unequal
      if (graph == nullptr or n.graph == nullptr)
          return false; 

      // valid nodes must belong to the same graph to compare equal
      if (graph != n.graph)
        return false; 

      // valid nodes must have the same id to compare equal
      if (this->id != n.id) 
        return false; 

      // if all comparisons above are passed, return true if the node is still
      // active in the graph 
      return (graph->active_nodes).count(id); 
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
      /** if both nodes are valid and belong to the same graph, the node with
       *  the smaller id compares lower. If one or both nodes are invalid, 
       *  or belong to different graphs, then the comparison is still made 
       *  based on node id, but this id does not have any significance w.r.t. 
       *  any graph object. If two nodes belong to different graphs, graph    
       *  pointers are compared. 
       */

      if (graph == nullptr or n.graph == nullptr) { 
        return id < n.id; 
      }

      if (graph == n.graph) {
        return id < n.id; 
      }
      
      return graph < n.graph; 
   }

   private:
    /** Allow Graph, Edge, and relevant iterators to access Node's private 
    member data and functions. */ 
    friend class Graph;
    friend class Edge; 
    friend class NodeIterator; 
    friend class IncidentIterator; 

    graph_type* graph; 
    size_type id; 

   /** Private constructor for use by the Graph class and associated 
     *  iterators.
    */
    Node(graph_type* graph, size_type id): graph(graph), id(id) {}
   
 };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return active_nodes.size();
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

    active_nodes.insert(next_id);

    // add node to graph's node data containers
    node_positions[next_id] = position;
    node_values[next_id] = value;  
    node_indices[next_id] = nodes.size(); 
    nodes.push_back(next_id); 

    // add a new entry to the adjacency map to store any nodes to be connected
    // to the newly added node in the future
    connected_nodes[next_id]; 
 
    return Node(const_cast<graph_type*>(this), next_id++);
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

    return (n.graph == this and active_nodes.count(n.id));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(i < nodes.size()); 
    
    return Node(const_cast<graph_type*>(this), nodes[i]); 
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

    /** Return an l-value reference to the value of this Edge */
    edge_value_type& value() {
      return fetch_value();  
    }

    /** Return an r-value reference to the value of this Edge */
    const edge_value_type& value() const {
      return fetch_value(); 
    }

    /** Return the vector representation of this Edge, from node1 to node2 
      * Note that the chosen orientation of the vector has no physical meaning,
      * since the graph is undirected. i.e. vector() and -1*vector() are 
      * equivalent representations. 
      */
    Point vector() const { 
      return node_b.position() - node_a.position(); 
    }

    /** Return the length of this Edge */ 
    double length() const {
      return norm(vector()); 
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
      if (*this == e) 
        return false; 

      if (node_a != e.node_a) 
        return node_a < e.node_a; 

      return node_b < e.node_b; 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    friend class IncidentIterator; 

    Node node_a; 
    Node node_b;

    /** Private constructor for use by the graph class and associated 
        iterators.*/
    Edge(const node_type& node_a, const node_type& node_b) :
      node_a(node_a), node_b(node_b)  {}

    /** Retrieve an l-value reference to the value associated with the Edge 
        object */
    edge_value_type& fetch_value() const {
      auto& graph = node_a.graph;
      
      assert(graph->has_edge(node_a, node_b)); 

      auto& edge_values = graph->edge_values; 
      
      return edge_values.at(std::min(node_a.id, node_b.id)).at(
                            std::max(node_a.id, node_b.id)); 
    }
    
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
    if (not has_node(a) or not has_node(b)) 
      return false; 
    
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

    Edge edge(a,b); 

    if (not has_edge(a, b)) {

      // insert each node into the other's adjacency list
      connected_nodes[a.id].insert(b.id); 
      connected_nodes[b.id].insert(a.id); 

      size_type min_id = std::min(a.id, b.id); 
      size_type max_id = std::max(a.id, b.id); 

      // insert the new edge into the Graph's edge data containers
      edge_values[min_id][max_id] = edge_value_type{}; 
      edge_indices[min_id][max_id] = edges.size(); 
      edges.push_back(edge); 
    } 

    return edge; 
  }

  /** Remove an edge given an EdgeIterator 
   *  Time complexity: O(1) 
   * 
   *  @pre @a e_it is a valid iterator to an edge in the graph 
   *  @post the edge pointed to by old @a e_it is removed from the graph
   *  @post graph.has_edge(*(old @a e_it)) returns false
   *  @post Iterators to the last graph Edge (i.e. old graph.edge_end()) are
   *        invalidated. Generally, @a e_it remains valid but points to a 
   *        different graph  edge. The Edge object pointed to by old @a e_it is
   *         rendered invalid, as are any equivalent Edges (Edges for which 
   *         *(old e_it) == e returns true).  
   *  @return A valid EdgeIterator to the start of the Graph's edge collection 
   *          (i.e. graph.edge_begin() )
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    Edge edge = *e_it; 
    Node a = edge.node_a; 
    Node b = edge.node_b; 

    // erase each of the edge's nodes from the other's adjacency list
    // O(1) 
    connected_nodes.at(a.id).erase(b.id); 
    connected_nodes.at(b.id).erase(a.id); 

    // erase the data pertaining to the edge object (i.e. value and index) 
    // O(1)
    edge_values.at(std::min(a.id, b.id)).erase(std::max(a.id, b.id)); 
    edge_indices.at(std::min(a.id, b.id)).erase(std::max(a.id, b.id)); 

    // locate the index of the edge in the Graph's ordered collection     
    unsigned index = e_it.internal_iterator - edges.begin(); 

    // O(1) swap-and-pop operation from the vectorized collection of Edges
    std::swap(edges.back(), edges[index]); 
    edges.pop_back(); 
    
    return edge_begin();      
  }

  /** Remove an edge from the graph given two Node objects 
   *  Time complexity: O(1) 
   * 
   *  @post If @a a and @b b are valid and connected nodes of the Graph object,
   *        the graph edge connecting the two nodes is removed 
   *  @post graph.has_edge(a, b) returns false 
   *  @post EdgeIterators to the last graph Edge (i.e. old graph.edge_end()) 
   *        are invalidated. IncidentIterators spawned by node @a a or node @b 
   *        b pointing to the edge connecting @a a to @b b are invalidated.     
   *  @post Any Edge object connecting @a a to @b b is invalidated. I.e. any 
   *        Edge for which edge.node1() == a and edge.node2() == b or vice 
   *        versa is rendered invalid.  
   *  @return If an edge was removed from the graph, return 1; otherwise, 0.  
   */
  size_type remove_edge(const Node& a, const Node& b) {

    if (not has_edge(a, b)) 
      return 0;

    // remove each node from the other's adjacency list: O(1) 
    connected_nodes.at(a.id).erase(b.id); 
    connected_nodes.at(b.id).erase(a.id); 

    // erase corresponding edge value: O(1)
    edge_values.at(std::min(a.id, b.id)).erase(std::max(a.id, b.id)); 

    // identify index of edge to be removed and then erase the index information
    // O(1)
    size_type edge_index = edge_indices[std::min(a.id, b.id)]
                                       [std::max(a.id, b.id)]; 
    edge_indices.at(std::min(a.id, b.id)).erase(std::max(a.id, b.id)); 


    // assign the index of the last element in the Edge container equal to the
    // index of the edge to be removed: O(1)
    Edge& back = edges.back(); 
    size_type a2_id = back.node_a.id, b2_id = back.node_b.id; 
    edge_indices[std::min(a2_id, b2_id)][std::max(a2_id, b2_id)] = edge_index;

    // swap-and-pop operation to remove the given Edge: O(1)
    std::swap(edges[edge_index], back);
    edges.pop_back(); 

    return 1; 
  }

  /** Remove an edge from the graph given the corresponding Edge object 
   *  Time complexity: O(1) 
   * 
   *  @post If edge.node1() and edge.node2() are valid and connected nodes of 
   *        the Graph object, the graph edge is removed 
   *  @post graph.has_edge(edge.node1(), edge.node2()) returns false 
   *  @post EdgeIterators to the last graph edge (i.e. old graph.edge_end()) 
   *        are invalidated. IncidentIterators spawned by old edge.node(1) and
   *        old edge.node(2) pointing to old @a edge are invalidated.
   *  @post All Edge objects that compare equal to @a edge are invalidated. 
   *  @return If an edge was removed from the graph, return 1; otherwise, 0.  
   */
  size_type remove_edge(const Edge& edge) {
    return remove_edge(edge.node_a, edge.node_b); 
  }

  /** Remove a node from the graph given a Node object. 
   *  Time complexity: O(node.degree()) 
   *  
   *  @post if old graph.has_node(@a node) returns true, the corresponding 
   *        node is removed from the graph
   *  @post new graph.has_node(@a node) returns false  
   *  @post Any NodeIterator equal to old graph.node_end() is invalidated. 
   *        EdgeIterators pointing to the last old node.degree() graph edges are
   *        invalidated. All IncidentIterators spawned by @a node are 
   *        invalidated. All IncidentIterators pointing to an edge which 
   *        connects the removed node to another graph node are invalidated.  
   *  @post all Edge objects for which edge.node1() == old @a node or 
   *        edge.node2() == old @a node are invalidated 
   *  @post @a node and any Node objects which compare equal to it are 
   *        invalidated. 
   *  @return if a node was removed from the graph, return 1; otherwise, 0
   */
  size_type remove_node(const Node& node) { 

    if (not has_node(node)) { 
      return 0;
    }

    // remove all incident edges: O(node.degree()) 
    while (node.edge_begin() != node.edge_end()) {
      remove_edge(*node.edge_begin()); 
    }
    
    // erase all information pertaining to the node object: O(1)
    connected_nodes.erase(node.id);     
    active_nodes.erase(node.id);
    node_positions.erase(node.id); 
    node_values.erase(node.id); 

    // replace the index of the last node in the vectorized collection with the 
    // index of the node to be removed: O(1)
    node_indices[nodes.back()] = node_indices[node.id];
 
    // swap-and-pop operation: O(1)
    std::swap(nodes[node_indices[node.id]], nodes.back());
    nodes.pop_back();  

    node_indices.erase(node.id); 

    return 1;    
  }

  /** Remove a node from the graph given a NodeIterator. 
   *  Time complexity: O(node.degree()) 
   *  
   *  @pre @a n_it is a valid, dereferenceable NodeIterator 
   *  @post the node pointed to by old @n_it is removed from the graph
   *  @post new graph.has_node(*(old @a n_it)) returns false  
   *  @post Any NodeIterator equal to old graph.node_end() is invalidated. 
   *        EdgeIterators pointing to the last old node.degree() graph edges are
   *        invalidated. All IncidentIterators spawned by the node are 
   *        invalidated. All IncidentIterators pointing to an edge which 
   *        connects the removed node to another graph node are invalidated. 
   *  @post all Edge objects for which edge.node1() == *(old @a n_it) or 
   *        edge.node2() == *(old @a n_it) are invalidated 
   *  @post *(old @a n_it) and any Node objects which compare equal to it are 
   *        invalidated. 
   */
  node_iterator remove_node(node_iterator n_it) {
    Node node = *n_it; 
    remove_node(node); 
    return node_begin(); 
  }

  /** Remove all nodes and edges from this graph.
   * Time complexity: O(1) 
   * 
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects, as well as all 
   * NodeIterators, EdgeIterators, and IncidentIterators. 
   */
  void clear() {
    // clear all node data
    active_nodes.clear(); 
    node_positions.clear();
    node_values.clear(); 
    node_indices.clear();
    nodes.clear(); 
 
    // clear all edge and adjacency data
    edges.clear(); 
    connected_nodes.clear(); 
    edge_values.clear();
    edge_indices.clear(); 
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
      return Node(graph, *internal_iterator);        
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
      internal_iterator++;
      return *this; 
    }

    /** Check whether two NodeIterator objects compare equal
     * @return true if the two iterators point to the same element in the graph
     *         Node container; otherwise, false
     */
    bool operator==(const NodeIterator& n) const {
      return (internal_iterator == n.internal_iterator); 
    }

    /** Check whether two NodeIterator objects are unequal 
     * @return !(*this == _n_)  
     */
    bool operator!=(const NodeIterator& n) const {
      return !(*this == n); 
    }

   private:
    friend class Graph;

    using iter_type  = typename std::vector<size_type>::const_iterator;
    
    // an internal iterator to traverse the graph's vector of InternalNodes
    iter_type internal_iterator; 

    // the graph associated with the NodeIterator
    graph_type* graph; 

    /** Private constructor for use by the Graph class */
    NodeIterator(const iter_type& internal_iterator, graph_type* graph) :
      internal_iterator(internal_iterator), graph(graph) {} 
  };

  /** Return an iterator to the first element in the container of Node objects
      associated with the present graph */
  node_iterator node_begin() const {
    return NodeIterator(nodes.begin(), const_cast<graph_type*>(this)); 
    std::cout << "nodes.end() - nodes.begin() " << nodes.end() - nodes.begin() 
              << std::endl;  
  }

  /** Return an iterator to one-past the last element in the container of Node
      objects associated with the present graph. Attempting to derference this
      iterator results in undefined behavior. */
  node_iterator node_end() const {
    return NodeIterator(nodes.end(), const_cast<graph_type*>(this)); 
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

  /** graph data: vectors of nodes and edges, containers mapping node and edge 
   *  IDs to respective node and edge data attributes, a set of active graph 
   *  nodes, and an adjacency map to facilitate O(1) lookups of existing edges 
   */

  std::unordered_set<size_type> active_nodes; 
  std::unordered_map<size_type, Point> node_positions; 
  std::unordered_map<size_type, node_value_type> node_values; 
  std::unordered_map<size_type, size_type> node_indices; 
  std::vector<size_type> nodes;  
  
  std::vector<Edge> edges;
  std::unordered_map<size_type, std::unordered_set<size_type>> connected_nodes;
  std::unordered_map<size_type, 
                     std::unordered_map<size_type, edge_value_type>> 
                                                               edge_values; 
  std::unordered_map<size_type, 
                     std::unordered_map<size_type, size_type>> edge_indices; 
  

 // id to assign to the next node added to the graph 
  size_type next_id = 0; 

};

#endif // CME212_GRAPH_HPP
