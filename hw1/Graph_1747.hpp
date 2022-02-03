#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <cassert>

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

  //
  // PUBLIC TYPE DEFINITIONS
  //
  
  /** Type of a node : user specified. */
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
    v2e = {};
    edges = {};
    vertices = {};
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
      this -> g = nullptr;
      this -> idx = -1;
    }

    Node(const Graph* g, size_type i) {
      this -> g = g;
      this -> idx = i;
    }

    /** Return this node's position. */
    const Point& position() const {
      const Point* vertex = g -> vertices[index()].first;
      return *vertex;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    /**
     * @brief gets a reference to a node's value.
     * @pre node is valid
     * @return node_value_type& 
     */
    node_value_type& value() {
      node_value_type* value = g -> vertices[index()].second;
      return *value;
    }

    /**
     * @brief gets a const reference to a node's value
     * @pre node is valid
     * @return const node_value_type& 
     */
    const node_value_type& value() const {
      const node_value_type*  value = g -> vertices[index()].second;
      return *value;
    }
    /**
     * @brief gets degree of a node
     * i.e the number of edges incident
     * 
     * @pre node is valid
     * 
     * @return size_type : number
     */
    size_type degree() const {
      std::unordered_map<size_type, size_type> incident_vertices = g -> getIncidentEdgesUnsafe(index());
      return incident_vertices.size();
    }
    
    /**
     * @brief gets the first edge as an iterator
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_begin() const {
      std::unordered_map<size_type, size_type> incident_vertices = \
        g -> getIncidentEdgesUnsafe(index());

      
      incident_iterator iter = IncidentIterator(incident_vertices, 0, this -> g, this);
      return iter;
    }

    /**
     * @brief gets the last edge as an iterator.
     * 
     * @return incident_iterator 
     */
    incident_iterator edge_end() const {
      std::unordered_map<size_type, size_type> incident_vertices = \
        g -> getIncidentEdgesUnsafe(index());
      
      incident_iterator iter = IncidentIterator(incident_vertices, degree(), this -> g, this);
      return iter;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (index() == n.index()) && (g == n.g);
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
      return n.index() < index();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;

    size_type idx;

    // tells node what graph they're apart of
    const Graph* g;


    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    
  }; // end Node Class Definition

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    return vertices.size();
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
  Node add_node(const Point& position, const node_value_type& = node_value_type()) {
    // initialize node
    // add it to vertices.

    size_type idx = this -> num_nodes();
    
    // create an entry for it in the v2e data structure so later we can just add to it.
    Node n = Node(this, idx);
    v2e[idx] = {};
    Point* point = new Point(position);
    node_value_type* val = new node_value_type();
    auto pair = std::make_pair(point, val);
    vertices.push_back(pair);
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return n.g == this && n.index() < num_nodes();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    Node n = Node(this, i);
    return n;
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
    Edge(const Node& start, const Node& end) : start(start), end(end) {

    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return start;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return end;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (e.node1() == node1() && e.node2() == node2())
       || (e.node1() == node2() && e.node2() == node1());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return e.node1() < node1() || \
       (e.node1() == e.node1() && e.node2() < e.node2());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // start and end vertices
    Node start;
    Node end;



    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    return *edges[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {

    // guaranteed to exist - it just may be empty
    std::unordered_map<size_type, size_type> incident_vertices = getIncidentEdgesUnsafe(a.idx); 
    return incident_vertices.find(b.idx) != incident_vertices.end();
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
    if (has_edge(a, b)) {

      std::unordered_map<size_type, size_type> a_edges = getIncidentEdgesUnsafe(a.idx);
      size_type idx = a_edges.at(b.idx);

      // swaps order of vertices: FIX as a result of hw0.
      // Order doesn't really matter except for returning what we 
      // expect to the user.

      // TODO Probably a data leak.
      edges[idx] = new Edge(a, b);
      return edge(idx);
    } else {
      Edge* e = new Edge(a, b);
      addVertexPairAndEdgeToV2E(a.idx, b.idx, num_edges());
      edges.push_back(e);
      return *e;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE

    // clears nodes
    vertices.clear();

    edges.clear();

    v2e.clear();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator  : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
      g = nullptr;
      current_vertex_index = -1;
    }

    /**
     * @brief Overload * operator
     * @pre 0 <= current_vertex_index < vertices.size()
     * @return Node 
     */
    Node operator*() const {
      return g -> node(current_vertex_index);
    }


    /**
     * @brief Overload ++ operator
     * @post current_vertex_index = current_vertex_index + 1
     * @return NodeIterator& 
     */
    NodeIterator& operator++() {
        current_vertex_index++;
        return *this;
    }

    /**
     * @brief 
     * 
     * @param other 
     * @return true if graphs and indices are same
     * @return false if not
     */
    bool operator==(const NodeIterator& other) const {
      return current_vertex_index == other.current_vertex_index && g == other.g;
    }

   private:
    friend class Graph;

    size_type current_vertex_index;

    const Graph* g;

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  /**
   * @brief gets node_iterator for first node.
   * @post current_vertex_index = 0
   * 
   * @return node_iterator 
   */
  node_iterator node_begin() const {
    node_iterator* it = new NodeIterator();
    it -> g = this;
    it -> current_vertex_index = 0;
    return *it;
  }

  /**
   * @brief gets node_iterator for end iterator
   * @post current_vertex_index = num_nodes()
   * @return node_iterator 
   */
  node_iterator node_end() const {
    node_iterator* it = new NodeIterator();
    it -> g =  this;
    it -> current_vertex_index = num_nodes();
    return *it;
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
    IncidentIterator() {

    }


    /**
     * @brief Construct a new Incident Iterator object
     * 
     * @param map : vertex : edge map
     * @param idx : iterated index
     * @param g : graph reference
     * @param n : reference to node it came from.
     */
    IncidentIterator(const std::unordered_map<size_type, size_type> &map, size_type idx, const Graph* g, const Node* n) : g(g), map(map) {
      

      // Set up order that we're going to iterate through incident edges.
      for (auto kv : map) {
        keys.push_back(kv.first);
      }

      // Set node pointer.
      this -> n = n;


      // Set index for iterator.
      index = idx;
    } 



    /**
     * @brief Overload * Operator
     * 
     * @return edge s.t edge.node1 == n
     * 
     * @return Edge 
     */
    Edge operator*() const {

      auto key = keys[index];
      auto edge_idx = map.at(key);
      auto edge = g -> edge(edge_idx);
      if (edge.node1() != *n) {
        auto new_edge = new Edge(edge.node2(), edge.node1());
        return *new_edge;
      } else {
        return edge;
      }
    }


    /**
     * @brief Overload ++ operator
     * @post index = index + 1;
     * 
     * @return IncidentIterator& 
     */
    IncidentIterator& operator++() {
      index++;
      return *this;
    }
    

    /**
     * @brief Equality Operator
     * 
     * @param other 
     * @return true if indices match
     * @return false else
     */
    bool operator==(const IncidentIterator& other) const {
      return (index == other.index);
    }

   private:
    friend class Graph;

    const Graph* g;
    const std::unordered_map<size_type, size_type> map;

    std::vector<size_type> keys;

    size_type index = -1;

    const Node* n;

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
    
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph* g, size_type index) : g(g) {
      this -> index = index;
    }


    /**
     * @brief Overloaded * Operator
     * @pre 0 <= index < g -> edges.size();
     * @return Edge at correct index 
     */
    Edge operator*() const {
      return g -> edge(index);
    }

    /**
     * @brief Overloaded ++ operator
     * @post index = index + 1
     * @return EdgeIterator& 
     */
    EdgeIterator& operator++() {
      index++;
      return *this;
    }

    /**
     * @brief Equality Operator
     * 
     * @param other 
     * @return true if indices match
     * @return false else
     */
    bool operator==(const EdgeIterator& other) const {
      return other.index == index;
    }

   private:
    friend class Graph;

    const Graph* g;


    // 
    size_type index = -1;

  };


  /**
   * @brief Get edge begin
   * 
   * 
   * 
   * @return edge_iterator, s.t edge_iterator.index = 0
   */
  edge_iterator edge_begin() const {
    return edge_iterator(this, 0);
  }

  /**
   * @brief Get edge end
   * 
   * 
   * 
   * @return edge_iterator, s.t edge_iterator.index = num_edges()
   */
  edge_iterator edge_end() const {
    return edge_iterator(this, num_edges());
  }

 private:
  std::vector<std::pair<Point*, node_value_type*>> vertices;
  std::vector<Edge*> edges;

  // map v1 : {v2 : v1-v2 edge index, v3: v1-v3 edge index}
  std::unordered_map<size_type, std::unordered_map<size_type, size_type>> v2e;

  /* Function adds a vertex pair of indices and edge index into the v2e data structure
  / @return void 
  / that is guaranteed to exist 
  */
  void addVertexPairAndEdgeToV2E(size_type a_idx, size_type b_idx, size_type e_idx) {
      v2e[a_idx].insert({b_idx, e_idx});
      v2e[b_idx].insert({a_idx, e_idx});
  }

  /** Get an unordered map of incident edges around a vertex
   * @pre a_idx exists in v2e.
   * @return unordered map of incident {vertex_idx : edge_idx} around a vertex index a_idx
   */
  std::unordered_map<size_type, size_type> getIncidentEdgesUnsafe(size_type a_idx) const {
    return v2e.at(a_idx);
  }

};

#endif // CME212_GRAPH_HPP
