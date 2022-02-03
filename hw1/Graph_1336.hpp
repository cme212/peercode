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
  /* Forward Declarations of Itnernal Data Structures */
  struct internal_node;
  struct internal_incident;

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
  using node_value_type = int;

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
      graph = nullptr;
      idx = 0;
    }

    /** Return this node's position. */
    const Point& position() const {
      return graph->nodes[idx].pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      assert(idx < graph->size());
      return idx;
    }

    /** Supply definitions AND SPECIFICATIONS for:
     *  node_value_type& value();
     *  const node_value_type& value() const;
     *  size_type degree() const;
     *  incident_iterator edge_begin() const;
     *  incident_iterator edge_end() const;
     */

    /** Check the size of the incident map for the corresponding Node */
    size_type degree() const {
      return graph->incident_map[*this].incidents.size();
    }

    /** Get the beginning of the incident map for the corresponding Node*/
    incident_iterator edge_begin() const {
      incident_iterator begin = IncidentIterator(graph, this, graph->incident_map[*this].incidents.begin());
      return begin;
    }

    /** Get the end of the incident map for the corresponding Node*/
    incident_iterator edge_end() const {
      incident_iterator end = IncidentIterator(graph, this, graph->incident_map[*this].incidents.end());
      return end;
    }

    /** Modify the Node's value data member */
    node_value_type& value(node_value_type val) {
      graph->nodes[idx].value = val;
      return graph->nodes[idx].value;
    }

    /** Get the Node's value data member */
    const node_value_type& value() const {
      return graph->nodes[idx].value;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      (void) n;          // Quiet compiler warning
      if (graph == n.graph && idx == n.idx) {
        return true;
      } else {
        return false;
      }
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
      (void) n;          // Quiet compiler warning
      if (this->index() < n.index()) {
        return true;
      } else {
        return false;
      }
    }

   private:
    /** Allow Graph to access Node's private member data and functions. */
    friend class Graph;

    /** Private data members needed for a valid Node */
    Graph* graph;
    size_type idx;

    /** Private Constructor for a valid Node */
    Node(const Graph* graph_, size_type idx_)
      : graph(const_cast<Graph*>(graph_)), idx(idx_) {}
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

  Node add_node(const Point& position, const node_value_type& val = 0) {
    (void) position;      // Quiet compiler warning

    /* Create new Node object with an index and assign it to an
     * internal Node with its corresponding position */
    Node new_node = Node(this, nodes.size());
    internal_node new_element;
    new_element.pos = position;
    new_element.node = new_node;

    /** Add internal Node to the graph's Node vector, then assign the Node a
     *  value, then return the Node */
    nodes.push_back(new_element);
    new_node.value(val);
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    (void) n;            // Quiet compiler warning
    if (n.graph == this && n.index() < nodes.size()) {
      return true;
    } else {
      return false;
    }
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    (void) i;             // Quiet compiler warning
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
  class Edge : private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      graph = nullptr;
      idx = 0;
      node1 = Node();
      node2 = Node();
    }

    /** Return a node of this Edge */
    Node node1() const {
      return node_1;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return node_2;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (this->node1() == e.node1() && this->node2() == e.node2()) {
        return true;
      } else if (this->node1() == e.node2() && this->node2() == e.node1()) {
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      if (this->node1() < e.node1()) {
        return true;
      } else {
        return false;
      }
    }

   private:
    /** Allow Graph to access Edge's private member data and functions. */
    friend class Graph;

    /** Private Data Members needed for a valid Edge */
    Graph* graph;
    size_type idx;
    const Node& node_1;
    const Node& node_2;

    /** Private Constructor of a valid Edge */
    Edge(const Graph* graph_, size_type idx_, const Node& node_1_, const Node& node_2_)
      : graph(const_cast<Graph*>(graph_)), idx(idx_), node_1(node_1_), node_2(node_2_) {
      }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return edges_vec.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes()), hopefully less
   */
  Edge edge(size_type i) const {
    (void) i;             // Quiet compiler warning
    return edges_vec[i];
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(log(num_nodes) + log(node_degree))
   */
  bool has_edge(const Node& a, const Node& b) {
    (void) a; (void) b;   // Quiet compiler warning
    if (incident_map.count(a) && (incident_map[a]).incidents.count(b)) {
      return true;
    } else if (incident_map.count(b) && (incident_map[b]).incidents.count(a)) {
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
   * Complexity: No more than O(log(num_nodes) + log(node_degree))
   */
  Edge add_edge(const Node& a, const Node& b) {
    (void) a, (void) b;   // Quiet compiler warning
    if (this->has_edge(a, b)) {
      /** Find and return existing edge */
      size_type idx = incident_map[a].incidents[b];
      return edges_vec[idx];
    } else {
      /** Create new valid edge */
      size_type idx = num_edges();
      Edge new_edge = Edge(this, idx, a, b);

      /** Add edge index to incident map for lookup by either node */
      if (incident_map.count(a)) {
        incident_map[a].incidents[b] = idx;
      } else {
        internal_incident new_element;
        new_element.incidents[b] = idx;
        incident_map[a] = new_element;
      }

      if (incident_map.count(b)) {
        incident_map[b].incidents[a] = idx;
      } else {
        internal_incident new_element;
        new_element.incidents[a] = idx;
        incident_map[b] = new_element;
      }

      /** Add edge to edge vector, return the new Edge */
      edges_vec.push_back(new_edge);
      return new_edge;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    edges_vec.clear();
    incident_map.clear();
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

    /** Incorporate Iterator from graph's node vector */
    using node_iter         = typename std::vector<internal_node>::const_iterator;

    /** Dereference to current corresponding node */
    Node operator*() const {
      return (*idx).node;
    }

    /** Increment to next node in node vector */
    NodeIterator& operator++() {
      idx++;
      return *this;
    }

    /** Compare graph's node vector iterators */
    bool operator==(const NodeIterator& it) const {
      if (graph == it.graph && idx == it.idx) {
        return true;
      }
      return false;
    }

    /** Compare graph's node vector iterators */
    bool operator!=(const NodeIterator& it) const {
      if (graph != it.graph || idx != it.idx) {
        return true;
      }
      return false;
    }

   private:
    /** Allow Graph to access NodeIterator's private member data and functions. */
    friend class Graph;

    /** Private data members to link to graph data */
    Graph* graph;
    node_iter idx;

    /** Construct an invalid NodeIterator. */
    NodeIterator(const Graph* graph_, node_iter idx_)
      : graph(const_cast<Graph*>(graph_)), idx(idx_) {}

  };

  /** Beginning of graph's node vector */
  NodeIterator node_begin() const {
    NodeIterator begin = NodeIterator(this, nodes.begin());
    return begin;
  }

  /** End of graph's node vector */
  NodeIterator node_end() const {
    NodeIterator end = NodeIterator(this, nodes.end());
    return end;
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator {
   public:
    /** These type definitions let us use STL's iterator_traits. */
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Incorporate Iterator from graph's incident map*/
    using incident_iter     = typename std::map<Node, size_type>::iterator;

    /** Construct an Invalid IncidentIterator */
    IncidentIterator(const Graph* graph_) {
      graph = const_cast<Graph*>(graph_);
      node = nullptr;
      idx = (graph->incident_map[(graph->nodes[0].node)]).incidents.begin();
    }

    /** Dereference corresponding Edge with correct internal node order*/
    Edge operator*() const {
      // Rebuild valid edge with root node as edge->node1()
      const Graph* graph_ = graph;
      size_type idx_ = (*idx).second;
      const Node& node_1_ = *node;
      const Node& node_2_ = (*idx).first;
      // Return rebuilt edge
      return Edge(graph_, idx_, node_1_, node_2_);
    }

    /** Increment Iterator to next incident node and corresponding Edge*/
    IncidentIterator& operator++(){
      idx++;
      return *this;
    }

    /** Compare graph's incident map iterators */
    bool operator==(const IncidentIterator& it) const {
      if (graph == it.graph && idx == it.idx) {
        return true;
      }
      return false;
    }

    /** Compare graph's incident map iterators */
    bool operator!=(const IncidentIterator& it) const {
      if (graph != it.graph || idx != it.idx) {
        return true;
      }
      return false;
    }

   private:
    /** Allow Graph to access IncidentIterator's private member data and functions. */
    friend class Graph;

    /** Private data members to link to Graph data */
    Graph* graph;
    const Node* node;
    incident_iter idx;

    /** Construct an valid IncidentIterator. */
    IncidentIterator(const Graph* graph_, const Node* node_, incident_iter idx_)
      : graph(const_cast<Graph*>(graph_)), node(const_cast<Node*>(node_)), idx(idx_) {}

  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator {
   public:
    /** These type definitions let us use STL's iterator_traits. */
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    /** Dereference internal incident_iterator to current Edge */
    Edge operator*() const {
      return (*incident_iter);
    }

    /** Increment internal incident and node iterators */
    EdgeIterator& operator++() {
      /** Loop until arrive at untraveled Edge */
      bool repeat = true;
      while(repeat) {
        /** Increment IncidentIterator First */
        ++incident_iter;
        if (incident_iter == (*node_iter).edge_end()) {
          /** Increment NodeIterator after all Incidents traversed */
          ++node_iter;
          if (node_iter != graph->node_end()) {
            /** Reset Incident Iterator to beginning of next Node */
            incident_iter = (*node_iter).edge_begin();
          } else {
            /** EdgeIterator end condition */
            incident_iter = IncidentIterator(graph);
            break;
          }
        }

        /* Check if current Edge is untouched */
        if (pass.count((*incident_iter).idx)) {
          continue;
        } else {
          /* Add current Edge to traversed Edges and end loop condition */
          pass.insert((*incident_iter).idx);
          repeat = false;
        }
      }
      return *this;
    }

    /** Check for equality by incident and node iterator equality */
    bool operator==(const EdgeIterator& it) const {
      if (incident_iter == it.incident_iter && node_iter == it.node_iter) {
        return true;
      }
      return false;
    }

    /** Check for inequality by incident and node iterator inequality */
    bool operator!=(const EdgeIterator& it) {
      if (incident_iter != it.incident_iter || node_iter != it.node_iter) {
        return true;
      }
      return false;
    }

   private:
    /** Allow Graph to access EdgeIterator's private member data and functions. */
    friend class Graph;

    /** Private data members to link to corresponding node and incident iterators */
    Graph* graph;
    NodeIterator node_iter;
    IncidentIterator incident_iter;
    std::set<size_type> pass;

    /** Construct a valid EdgeIterator. */
    EdgeIterator(const Graph* graph_, NodeIterator node_iter_, IncidentIterator incident_iter_)
      : graph(const_cast<Graph*>(graph_)), node_iter(node_iter_), incident_iter(incident_iter_) {}
  };

  /** Find beginning node and coressponding beginning incident edge */
  edge_iterator edge_begin() const {
    NodeIterator node_begin = this->node_begin();
    edge_iterator begin = EdgeIterator(this, node_begin, nodes[0].node.edge_begin());
    return begin;
  }

  /** Create ending edge condition */
  edge_iterator edge_end() const {
    NodeIterator node_end = this->node_end();
    IncidentIterator edge_end = IncidentIterator(this);
    edge_iterator end = EdgeIterator(this, node_end, edge_end);
    return end;
  }

 private:
  /** Internal Node structure to track geometric positions of each node */
  struct internal_node {
    Point pos;
    Node node;
    node_value_type value;
  };

  /** Internal Incident structure to track node pairings and their corresponding
   *  Edge index (bi-directional) */
  struct internal_incident {
    std::map<Node, size_type> incidents;
  };

  /** Private data members to track and link Edges, Nodes, and their
   *  geometric positions */
  std::vector<internal_node> nodes;
  std::vector<Edge> edges_vec;
  std::map<Node, internal_incident> incident_map;
};

#endif // CME212_GRAPH_HPP
