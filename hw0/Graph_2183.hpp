#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 private:

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
    idx_node_map = {};
    idx_edge_map = {};
    edge_idx_map = {};
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
  class Node {
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

    // Node(Graph* this_graph, size_type this_idx) {
    //   // instantiate new node with pointer to graph and index
    //   my_graph = this_graph;
    //   idx = this_idx;
    // }

    Node(const Graph* this_graph, size_type this_idx) 
    : my_graph(const_cast<Graph*>(this_graph)), idx(this_idx) {
    }

    //     SimpleElement(const Ex_Graph* set, size_type uid)
    //     : set_(const_cast<Ex_Graph*>(set)), uid_(uid) {
    // }

    Node() {}

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // get internal node
      return this->my_graph->idx_node_map.at(this->idx).my_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->idx;
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
      bool graph_equal = (this->my_graph == n.my_graph);
      bool idx_equal = (this->idx == n.idx);       
      return (graph_equal & idx_equal);
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
      // Compare indices, see if equal
      return this->idx < n.idx;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // the node has a pointer to graph 
    Graph* my_graph;
    // the node has an index
    size_type idx;
  };



  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    size_type my_size = this->idx_node_map.size();
    return my_size;
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
  // called on the graph object, self
  Node add_node(const Point& position) {
      // HW0: YOUR CODE HERE
      // Create internal node, add to set 
      size_type idx = num_nodes();
      internal_node my_int_node;
      // (in order to not copy over pointer)
      my_int_node.my_point = Point(position.x, position.y, position.z);
      // add to map
      idx_node_map[idx] = my_int_node;
      // Create node and return
      Node my_node = Node(this, idx);
      return my_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    size_type my_idx = n.index();
    if ((my_idx < num_nodes()) & (this == n.my_graph)){
        if (idx_node_map.at(my_idx).my_point == n.position()){
            return true;
        }
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
    // HW0: YOUR CODE HERE
    assert(i < num_nodes());
    return Node(this, i);        
  }
















  //
  // EDGES  legooooo
  //

  /** @class Graph::Edge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge {
   public:

    /** Construct a valid edge to return **/
    Edge(const Graph* this_graph, size_type this_idx) 
    : my_graph(const_cast<Graph*>(this_graph)), idx(this_idx) {
    }


    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return this->my_graph->idx_edge_map.at(this->idx).node_1;      // Invalid Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return this->my_graph->idx_edge_map.at(this->idx).node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE      
      bool same_graph = (this->my_graph == e.my_graph);
      // case 1: n1=n1, n2=n2
      bool node_same = ((this->node1() == e.node1()) & \
      (this->node2() == e.node2()));
      bool node_flipped = ((this->node1() == e.node2()) & \
      (this->node2() == e.node1()));
      
      return (same_graph & (node_flipped | node_same));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // compare indices of the nodes
      return this->idx < e.idx;;
    }

    /** Return the index of this edge 
     * 
     */
     size_type index() const {
         return this->idx;
     } 

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* my_graph;
    size_type idx;

  };



  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->idx_edge_map.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // assert we are asking for valid edge
    assert(i < num_edges());
    // return new node
    return Edge(this, i);        
  }

  

 /** Returns the index of the edge connecting two valid nodes.
  * @pre assumes valid distinct nodes
  * @return index if valid, UINT_MAX if not valid **/
  size_type index(const Node& a, const Node& b) const {
      // HW0: my code
          // case 1: idx_a less than idx_b
    if (a < b) {
        // if a -> b in graph
        if (a.my_graph->edge_idx_map.count(a.idx) == 1) {
            if (a.my_graph->edge_idx_map.at(a.idx).count(b.idx) == 1) {
                return a.my_graph->edge_idx_map.at(a.idx).at(b.idx);
            }
        }
        // case 2: idx_b less than idx_a
    } else if (b < a) {
        // if b -> a in graph
        if (b.my_graph->edge_idx_map.count(b.idx) == 1) {
            if (b.my_graph->edge_idx_map.at(b.idx).count(a.idx) == 1) {
                return b.my_graph->edge_idx_map.at(b.idx).at(a.idx);
            }
        }
    }
    return UINT_MAX;
  }

/** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // see if index is invalid
    if (index(a,b) == UINT_MAX) {
        return false;
    }
    return true;
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
    // check if edge already exists 
    if (index(a, b) != UINT_MAX) {
        // return the current edge
        return Edge(this,index(a,b));
    };
    // decleare index and new internal edge
    size_type idx = num_edges();
    internal_edge my_int_edge;
    // add to idx_edge_map
    // add to edge_idx_map
    if (a < b) {
        // if there's not already an unordered map, make one
        if (a.my_graph->edge_idx_map.count(a.idx) == 0) {
            a.my_graph->edge_idx_map[a.idx] = {};
        }
        a.my_graph->edge_idx_map.at(a.idx)[b.idx] = idx;
        my_int_edge.node_1 = a;
        my_int_edge.node_2 = b;
        a.my_graph->idx_edge_map[idx] = my_int_edge;
    } else if (b < a) {
        // if there's not already an unordered map, make one
        if (b.my_graph->edge_idx_map.count(b.idx) == 0) {
            b.my_graph->edge_idx_map[b.idx] = {};// 
        }
        b.my_graph->edge_idx_map.at(b.idx)[a.idx] = idx;
        my_int_edge.node_1 = a;
        my_int_edge.node_2 = b;
        a.my_graph->idx_edge_map[idx] = my_int_edge;
    } else {
        assert(!(a == b)); // can't add edge between identical nodes
    };
    // return the new edge
    return Edge(this, idx);        
  }



  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // clear all datastructures
    this->idx_node_map = {};
    this->idx_edge_map = {};
    this->edge_idx_map = {};
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

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

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
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

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const













 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.
  
  // internal_node DATASTRUCTURE
  struct internal_node {
    Point my_point;
  };

    // idx_node_map
  std::unordered_map<size_type, internal_node> idx_node_map;

// internal_edge DATASTRUCTURE
  struct internal_edge {
    Node node_1;
    Node node_2;
};

// idx_edge_map
std::unordered_map<size_type,internal_edge> idx_edge_map;

// edge_idx_map
std::unordered_map<size_type, \
std::unordered_map<size_type, size_type>> edge_idx_map;

};

#endif // CME212_GRAPH_HPP