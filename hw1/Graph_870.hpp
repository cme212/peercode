#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <set>

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
  Graph() : nodelst(), edgelst(), edgesets() {}

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
      _graph_obj = nullptr;
      _node_id = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return _graph_obj -> nodelst.at(_node_id)._position;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return _graph_obj -> nodelst.bucket(_node_id);
    }
  
    /** Return this node's value. */
    node_value_type& value(){
      return _graph_obj-> nodelst.at(_node_id)._value;
    }

    /** Return this node's value. */
    const node_value_type& value() const {
      return _graph_obj-> nodelst.at(_node_id)._value;
    }

    /** Return this node's degree. */
    size_type degree() const{
      assert(*this != Node()); //asserting that the node is valid
      return _graph_obj-> nodelst.at(_node_id).incident_edge_ids.size();
    }

    /** Return the begin iterator for the incident iterator */
    incident_iterator edge_begin() const{
      return IncidentIterator(_graph_obj, _node_id, 0);
    }

   /** Return the end iterator for the incident iterator */
    incident_iterator edge_end() const{
      return IncidentIterator(_graph_obj, _node_id, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return (n._node_id == _node_id) && (_graph_obj == n._graph_obj);
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
      return _node_id < n._node_id;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Pointer back to the Graph container
    Graph* _graph_obj;
    // This nodes's unique index number
    size_type _node_id;

    /** Private Constructor */
    Node(const Graph* _graph, size_type node_id)
        : _graph_obj(const_cast<Graph*>(_graph)), _node_id(node_id) {}
    
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return nodelst.size();
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
  Node add_node(const Point &position, const node_value_type& val = node_value_type()) {
    //inserting new node_id & internal_node pair into node map
    size_type new_id = num_nodes();
    internal_node new_node = internal_node(position, new_id, val);
    nodelst.insert(std::pair<size_type,internal_node>(new_id, new_node));
    return Node(this, new_id);        
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return (n._graph_obj == this);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    if (i == size()){
      return Node();
    }
    assert(i < num_nodes() && i >= 0);
    auto it = nodelst.begin(i);
    return Node(this, it->first); 
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
      _graph_obj = nullptr;
      _edge_id = 0;
    }

    /** Return a node of this Edge */
    Node node1() const{
      if (_node1_idx == -1) {
        return Node(_graph_obj, _graph_obj -> edgelst.at(_edge_id).node1_id);
      }
      else{
        return Node(_graph_obj, _node1_idx);
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      if (_node2_idx == -1) {
        return Node(_graph_obj, _graph_obj -> edgelst.at(_edge_id).node2_id);
      }
      else{
        return Node(_graph_obj, _node2_idx);
      }
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      if (_graph_obj != e._graph_obj) {return false;}
      return ((node1() == e.node1() && node2() == e.node2()) ||
      (node1() == e.node2() && node2() == e.node1()));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return _edge_id < e._edge_id;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    Graph* _graph_obj;
    // This edge's unique index number
    size_type _edge_id;
    //private constructor for edge
    Edge(const Graph* _graph, size_type edge_id)
        : _graph_obj(const_cast<Graph*>(_graph)), _edge_id(edge_id) {}

    int _node1_idx = -1;
    int _node2_idx = -1;

    Edge(const Graph* _graph, size_type edge_id, size_type node1_idx, \
      size_type node2_idx)
        : _graph_obj(const_cast<Graph*>(_graph)), _edge_id(edge_id), 
        _node1_idx(node1_idx), _node2_idx(node2_idx) {}  
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edgelst.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges() && i >= 0);
    auto it = edgelst.begin(i);
    return Edge(this, it->first);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    //using overload operator below for comparisons
    assert(a._graph_obj == this && b._graph_obj == this);
    assert(a._node_id != b._node_id);
    //compares the (a, b) set with the existing map of sets
    std::set<size_type> compare_set;
    compare_set.insert(a._node_id);
    compare_set.insert(b._node_id);
    if (edgesets.find(compare_set) != edgesets.end()){
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
    assert(a._graph_obj == b._graph_obj);
    assert(a._node_id != b._node_id);
    size_type new_id = num_edges();
    std::set<size_type> new_set; 
    new_set.insert(a._node_id); 
    new_set.insert(b._node_id);

    if (!has_edge(a, b)) { 
      //inserting new edge_id & internal_edge pair into edge map 
      //inserting a new (a, b) pair into the edgesets map
      internal_edge new_edge{a._node_id, b._node_id, new_id};
      edgelst.insert(std::pair<size_type,internal_edge>(new_id, new_edge));
      edgesets.insert({new_set, new_id});
      //adding values to the incident edge_lst
      nodelst.at(a._node_id).incident_edge_ids.push_back(new_id);
      nodelst.at(b._node_id).incident_edge_ids.push_back(new_id);
    }
    else{
      size_type existing_edge_id = edgesets.at(new_set);
      internal_edge* existing_edge = &edgelst.at(existing_edge_id);

      if (existing_edge->node1_id == a.index() && existing_edge->node2_id == b.index()) {
        return Edge(this, edgesets.at(new_set));
      }
      else {
        existing_edge->node1_id = a.index();
        existing_edge->node2_id = b.index();
        return Edge(this, edgesets.at(new_set));
      }
    }
  return Edge(this, new_id);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodelst.clear();
    edgelst.clear();
    edgesets.clear();
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
      _graph_obj = nullptr;
      _iterator_id = 0;
    }

    //dereference the node iterator
    Node operator*() const{
      return _graph_obj->node(_iterator_id);
    }

    //increment to the next node_iterator id
    NodeIterator& operator++() {
      _iterator_id += 1;
      return *this;
    }

    //equality operator within the node iterator class 
    bool operator==(const NodeIterator& iterator_2) const {
      return (_graph_obj == iterator_2._graph_obj && \
        _iterator_id == iterator_2._iterator_id);
    }

   private:
    friend class Graph;

    Graph* _graph_obj;
    size_type _iterator_id;

    /*private constructor for the NodeIterator class*/
    NodeIterator(const Graph* _graph, size_type it_val)
        : _graph_obj(const_cast<Graph*>(_graph)), _iterator_id(it_val) {}
  };

  /*returns the beginning iterator for NodeIterator*/
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  } 

  /*returns the end iterator for NodeIterator*/
  node_iterator node_end() const {
    return NodeIterator(this, nodelst.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
   class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
      _graph_obj = nullptr;
      _iterator_id = 0;
      _node_id = 0;
    }

    //dereference operator for the incident iterator class
    Edge operator*() const{
        Edge e = Edge(_graph_obj, _graph_obj->nodelst.at(_node_id).incident_edge_ids.at(_iterator_id));
        if (e.node1() != Node(_graph_obj, _node_id)){
          return Edge(_graph_obj, _graph_obj->nodelst.at(_node_id).incident_edge_ids.at(_iterator_id), \
          e.node2().index(), e.node1().index());
        }
        else{
          return e;
        }
    }

    //increment to the next iterator_id
    IncidentIterator& operator++() {
      _iterator_id += 1;
      return *this;
    }

    //equality operator within the incident_iterator class 
    bool operator==(const IncidentIterator& iterator_2) const {
      return (_graph_obj == iterator_2._graph_obj && \
      _iterator_id == iterator_2._iterator_id && \
      _node_id == iterator_2._node_id);
    }

   private:
    friend class Graph;

    Graph* _graph_obj;
    size_type _node_id;
    size_type _iterator_id;

    /*private constructor for the IncidentIterator class*/
    IncidentIterator(const Graph* _graph, size_type node_id, size_type it_val) \
      : _graph_obj(const_cast<Graph*>(_graph)), _node_id(node_id), _iterator_id(it_val) {}
  };
  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
      _graph_obj = nullptr;
      _iterator_id = 0;
    }

    //dereference operator for the edge iterator
    Edge operator*() const{
      return _graph_obj->edge(_iterator_id);
    }

    //increment to the next iterator_id
    EdgeIterator& operator++() {
      _iterator_id += 1;
      return *this;
    }

    //equality operator for the edge iterator
    bool operator==(const EdgeIterator& iterator_2) const {
      return (this-> _graph_obj == iterator_2._graph_obj && \
      this -> _iterator_id == iterator_2._iterator_id);
    }

   private:
    friend class Graph;

    Graph* _graph_obj;
    size_type _iterator_id;

    /*private constructor for the EdgeIterator class*/
    EdgeIterator(const Graph* _graph, size_type it_val)
        : _graph_obj(const_cast<Graph*>(_graph)), _iterator_id(it_val) {}
  };

  /*returns the beginning iterator for EdgeIterator*/
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  } 

  /*returns the end iterator for EdgeIterator*/
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }
 private:

  //internal node definition 
  struct internal_node {

  Point _position;
  size_type _idx;
  node_value_type _value;

  std::vector<size_type> incident_edge_ids;

  internal_node(const Point& position, \
    const size_type idx, node_value_type val) : _position(position), 
    _idx(idx), _value(val) {};
  };

  //internal edge definition
  struct internal_edge{
  size_type node1_id;
  size_type node2_id;
  size_type _edge_id;
  
  internal_edge(size_type node1, size_type node2, \
    size_type edge_id) : node1_id(node1), \
    node2_id(node2), _edge_id(edge_id) {};
  };

  std::unordered_map<size_type, internal_node> nodelst; //id, internal pair
  std::unordered_map<size_type, internal_edge> edgelst; //id, internal pair
  std::map<std::set<size_type>, size_type> edgesets; //set, id pair
};

#endif // CME212_GRAPH_HPP