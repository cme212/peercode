#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

// Newly added include file.
#include <map>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @ class Graph
 * @ brief A template for 3D undirected graphs.
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
    Node() {
      // Doesn't do anything for the invalid constructor.
    }

    /** Return this node's position. */
    const Point& position() const {
      // The node needs to be valid.
      assert(id < graph->size());
      return graph -> point_list[id];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return id;
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
  Node add_node(const Point& position) {
    size_type num = size();
    Node new_node = Node(this, num);
    point_list.push_back(position);
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
  class Edge {
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
        node1_id(node1_id_input), node2_id(node2_id_input), 
        graph(const_cast<Graph*>(graph_input)), id(id_input){}
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
    assert(has_node(a) && has_node(b) && ! (a==b));
    size_type num = num_edges();

    // It is meaningful to mention that I alway let the node1
    // of the edge has lower id. This facilitates the traversal
    // process of edges. The complexity is O(num_edges()).
    size_type id_min = std::min(a.id, b.id);
    size_type id_max = std::max(a.id, b.id);
    size_type i = find_edge(a, b);
    if(i == num_edges()){
      Edge new_edge = Edge(this, id_min, id_max, num);
      edge_list.push_back(new_edge);
      node_edge[id_min].push_back(i);
      node_edge[id_max].push_back(i);
      return new_edge;
    }else{
      return edge_list[i];
    }
  }



  /** Given two node and find the id of the edge that connects them.
   * This is a helpful function and
   * @pre @a node_a and @a node_b are different valid nodes
   * @return the id of edge if it exists. Otherwise, return the number of edges.
   * 
   * Complexity: O(num_edges()).
   */
  size_type find_edge(const Node& node_a, const Node& node_b) const{
    // We need to ensure that node_a/b are different and valid.
    assert(has_node(node_a) && has_node(node_b) && !(node_a == node_b));
    size_type id_max = std::max(node_a.id, node_b.id);
    size_type id_min = node_a.id + node_b.id - id_max;
    if(node_edge.find(id_max) != node_edge.end()){
      for(size_type i:node_edge.at(id_max)){
        if(edge_list[i].node1().id == id_min)
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

  /** The map node_edge has the key as the index of node. The corresponding 
   * item of the key is the index of edges that connect the node. This helps 
   * to lower the complexity of many operations later on.
   */
  std::map<size_type, std::vector<size_type>> node_edge;
};

#endif // CME212_GRAPH_HPP
