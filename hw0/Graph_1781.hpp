#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#define MAX_NODES_LONG_INT 92682 // if size_t is to be able to index the edges, this is the largest number of nodes we can contain (for a complete graph)
#define MAX_NODES_LONG_LONG_INT 6074001000

#include <algorithm>
#include <vector>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <list>
#include <string>
#include <functional>
#include <limits.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
using std::vector;
using std::pair;
using std::unordered_map;
using std::unordered_set;
using std::list;

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
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
  using size_type = unsigned long int;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
  }

  /** Default destructor */
  ~Graph() = default;

  //
  // NODES
  //

  const Point& get_point(size_type index) const {
    return points[index].first;
  }

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
    Node () {
      idx = ULONG_MAX;
    }

    /** Return this node's position. */
    const Point& position() const {
      return parent->get_point(idx);
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

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return parent == n.parent && idx==n.idx;
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
      return idx < n.idx;
    }

    const Graph* get_parent() const {
      return parent;
    }
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Declaring private data members and methods for Node
    const Graph *parent; // 8 bytes (x64)
    size_type idx; // 4 bytes
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
  Node add_node(const Point& position) {
    if (nodes.size() == MAX_NODES_LONG_INT)
      return Node();
    Node new_node;
    new_node.idx = nodes.size();
    new_node.parent = this;
    incident Edges;
    nodes[new_node] = Edges;
    pair<Point, const Node*> point_node_pair (position, &nodes.find(new_node)->first);
    points.push_back(point_node_pair);
    return nodes.find(new_node)->first; // equivalent to `new_node`
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    return nodes.find(n) != nodes.end();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    return *points[i].second;
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
    }

    /** Return a node of this Edge */
    Node node1() const {
      return n1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return n2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return e.node1()==node1() && e.node2()==node2();
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      return node1()==e.node1() ? node2()<e.node2() : node1()<e.node1();
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Declaring private data members and methods for Edge
    Node n1; // 12
    Node n2; // 12
    size_type idx;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edge_index.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    return edge_index[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    return find_edge(a,b) != NULL;
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
    if (has_edge(a,b)) {
      return *find_edge(a,b);
    }
    else if (a.parent==this && b.parent==this) { // we want the edge to connect only vertices from the same parent graph
      Edge new_edge;
      new_edge.idx = edge_index.size();
      new_edge.n1 = a<b ? a : b;
      new_edge.n2 = a<b ? b : a;
      edge_index.push_back(new_edge);
      nodes[a].emplace(&edge_index[edge_index.size()-1]);
      nodes[b].emplace(&edge_index[edge_index.size()-1]);
      return new_edge;
    }
    return Edge();
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    points.clear();
    nodes.clear();
    edge_index.clear();
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

 void print_edges() {
   for (unsigned long int k = 0; k < edge_index.size() ; k++) {
      std::cout<< "Edge index: " << edge(k).idx<< "  n1 index: " <<edge(k).node1().index() << std::endl;
   }
 }
 private:
  // Graph class's private attributes/data members and helper functions
  using incident = unordered_set<const Edge*>;
// 
  const Edge* find_edge(const Node& a, const Node& b) const {
    // if (!has_node(a) && !has_node(b)) {// nodes not in graph, i.e., invalid
    // if (nodes.find(a) != nodes.end() || nodes.find(b) != nodes.end()) {// nodes not in graph, i.e., invalid
    //   return NULL;
    // }
    if (nodes.at(a).size() < nodes.at(b).size()) {
      for (auto incident_edge : nodes.at(a)) { // iterate over the vertex with a smaller degree
        if (nodes.at(b).find(incident_edge) != nodes.at(b).end()) {
          return incident_edge;
        }
      }
    }
    else {
      for (auto incident_edge : nodes.at(b)) {// iterate over the vertex with a smaller degree
        if (nodes.at(a).find(incident_edge) != nodes.at(a).end()) {
          return incident_edge;
        }
      }
    }
   return NULL;
  }
  
  // code for hash functions referenced from https://en.cppreference.com/w/cpp/utility/hash
  struct NodeHash
  {
      std::size_t operator()(Node const& node) const noexcept
      {
          return std::hash<size_type>{}(node.idx);
      }
  };

  vector<pair<Point, const Node*>> points;
  unordered_map<Node,unordered_set<const Edge*>,NodeHash> nodes;
  vector<Edge> edge_index;

};

#endif // CME212_GRAPH_HPP