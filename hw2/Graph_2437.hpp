#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 *  @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

// To complete HW1 I combined several ideas that were given by
// students in whiteboard collaboration in office hours, as well as code 
// directly suggested by TAs in 1-on-1s and again whiteboard collaboration.
// Examples of are: BFS implementation, adjusting ids for edge_iterator, 
//                  filter_iterator set-up and ++ adjustment.
// The office hours I attended were Axel and Matteo both first and sec week
// of the assignment. Some of my code is also inspired by posted peer code,
// i.e. codes 976, 1038.

// To complete HW2 I combined several ideas that were given by
// students in whiteboard collaboration in office hours, as well as code 
// directly suggested by TAs in 1-on-1s and again whiteboard collaboration.
// Examples of are: polymorphic design for combined forces and constraints
//                  constraint implementation and simulate it
//                  remove_node (diagrams and pseucode from office hours)
// The office hours I attended were Matteo and Aasavari on the second week of the assignment.


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
class Graph {

 private:

  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  struct int_node; // internal representation of a node (proxy design pattern)
 
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
  // HW0: YOUR CODE HERE
  
  Graph(): nodes_vec(), edges_adj_vec(), edges_val_vec(), corresp(), edges_count(0) {}

  /** Default destructor */
  ~Graph() = default;

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_vec.clear();
    edges_adj_vec.clear();
    edges_val_vec.clear();
    corresp.clear();
    edges_count = 0;
  }

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
    // HW0: YOUR CODE HERE
    Node(): graph_ptr(nullptr), node_id(0) {}

    /** Return this node's position. */
    const Point& position() const {
      return graph_ptr->nodes_vec[node_id].point;
    }

    /** Return this node's position. Modifiable version. */
    Point& position() {
      return graph_ptr->nodes_vec[node_id].point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return graph_ptr->nodes_vec[node_id].node_id;
    }

    /** return a reference to the new value for the node so it can be assigned */
    node_value_type& value() {
      return const_cast<node_value_type&>(static_cast<const Node*>(this)->value());
    }

    /** directly return the value stored in the node (field!) */
    const node_value_type& value() const {
      return graph_ptr->nodes_vec[node_id].value;
    }

    size_type degree() const {
      if (node_id >= graph_ptr->edges_adj_vec.size()) {
        return 0;
      }
      return graph_ptr->edges_adj_vec[node_id].size();
    }

    incident_iterator edge_begin() const {
      return incident_iterator(graph_ptr, node_id, 0);
    }

    incident_iterator edge_end() const {
      // end at the last adjcent
      return incident_iterator(graph_ptr, node_id, degree());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // (void) n;          // Quiet compiler warning
      return (node_id == n.node_id) && (graph_ptr == n.graph_ptr);
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
      if (node_id != n.node_id) {
        return node_id < n.node_id;
      } else {
        return graph_ptr < n.graph_ptr;
      }
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    Graph* graph_ptr;
    size_type node_id;

    Node(const Graph* ptr, size_type id) : graph_ptr(const_cast<Graph*>(ptr)), node_id(id) {}

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return corresp.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return corresp.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point &position, const node_value_type& v = node_value_type()) {
    // emplace_back is more efficient as it doesn't create then destroy temp objs
    nodes_vec.emplace_back(position, v, corresp.size());
    corresp.push_back(nodes_vec.size() - 1);
    edges_adj_vec.resize(nodes_vec.size() + 1);
    edges_val_vec.resize(nodes_vec.size() + 1);
    return Node(this, nodes_vec.size() - 1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return this == n.graph_ptr && n.node_id < nodes_vec.size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(const_cast<Graph*>(this), corresp[i]);
  }

  /**
   * Remove passed-in node from the graph and edges that contain it
   * 
   * @complexity: O(num_nodes)
   * 
   * @pre @a n is a valid node
   * @pre @a g.node(i).index() == i
   * 
   * @post g.node(i).index() == i
   * @post num_nodes() is decremented by 1
   * @post edges_count is decremented by degree of @a n
   * @post edges that contained @a n are removed
   * @post has_node(n) will return False
   * @post associated NodeIterators, EdgeIterators, Incident iterators are invalidated
   * 
   */
  size_type remove_node(const Node& n) {
    if (!has_node(n)) {
      return (size_type) 0;
    }
    size_type removed_pos;
    for (unsigned i = 0; i < corresp.size(); i++) {
      if (corresp[i] == n.node_id) {
        removed_pos = i;
        break;
      }
    }
    corresp.erase(corresp.begin() + removed_pos);
    while (edges_adj_vec[n.node_id].size() > 0) {
      remove_edge(n, Node(n.graph_ptr, edges_adj_vec[n.node_id][0]));
    }
    for (unsigned i = removed_pos; i < corresp.size(); i++) {
      nodes_vec.at(corresp[i]).node_id = i;
    }
    return (size_type) 1;
  }

  /**
   * Remove the node referenced by the passed-in iterator
   * 
   * @Complexity: O(num_nodes)
   * 
   * @pre @a n is a valid node, where n = (*node_iter)
   * @pre @a g.node(i).index() == i
   * 
   * @post g.node(i).index() == i
   * @post num_nodes() is decremented by 1
   * @post edges_count is decremented by degree of @a n
   * @post edges that contained @a n are removed
   * @post has_node(n) will return False
   * @post associated NodeIterators, EdgeIterators, Incident iterators are invalidated
   */ 
  node_iterator remove_node(node_iterator node_iter) {
    return NodeIterator(this, remove_node(*node_iter));
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
    Edge() : graph_ptr(nullptr), node1_id(0), node2_id(0) {}

    /** Return a node of this Edge */
    Node node1() const {
      return Node(graph_ptr, node1_id);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(graph_ptr, node2_id);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool cond1 = graph_ptr == e.graph_ptr;
      bool cond2 = node1_id == e.node1_id;
      bool cond3 = node2_id == e.node2_id;
      bool cond4 = node1_id == e.node2_id;
      bool cond5 = node2_id == e.node1_id;
      return cond1 && ((cond2 && cond3) || ((cond4 && cond5)));
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      if (node1_id != e.node1_id) {
        return node1_id < e.node1_id;
      } else { // same first node id
        if (node2_id != e.node2_id) {
          return (node2_id < e.node2_id);
        }
        return true;
      }
    }

    float edge_len() {
      return norm(node2().position() - node1().position());
    }

    /**
     * Returns modifiable reference to edge value
     */
    edge_value_type& value() {
      Edge ed = Edge(graph_ptr, std::min(node1().node_id, node2().node_id),
          std::max(node1().node_id, node2().node_id));
      IncidentIterator it = ed.node1().edge_begin();
      while (it != ed.node1().edge_end() && (*it).node2_id != ed.node2_id) {
        ++it;
      }
      return graph_ptr->edges_val_vec[ed.node1_id][it.adj_id];
    }

    /**
     * Returns unmodifiable reference to edge value
     */
    const edge_value_type& value() const {
      Edge ed = Edge(graph_ptr, std::min(node1_id, node2_id), std::max(node1_id, node2_id));
      IncidentIterator it = ed.node1().edge_begin();
      while (it != ed.node1().edge_end() && (*it).node2_id != ed.node2_id) {
        ++it;
      }
      return graph_ptr->edges_val_vec[ed.node1_id][it.adj_id];
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* graph_ptr;
    size_type node1_id;
    size_type node2_id;

    Edge(Graph* ptr, size_type id_1, size_type id_2) : graph_ptr(ptr),
      node1_id(id_1), node2_id(id_2) {}

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_count;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Now we can use an iterator to do this faster!
    EdgeIterator iterator = edge_begin();
    while (i > 0) {
      ++iterator;
      i--;
    }
    // now we de-reference the iterator at the correct position
    Edge e = *iterator;
    return Edge(const_cast<Graph*>(this), e.node1_id, e.node2_id);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    for (unsigned int i = 0; i < edges_adj_vec[a.node_id].size(); i++) {
      if (edges_adj_vec[a.node_id][i] == b.node_id) {
        return true;
      }
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
    size_type i = a.node_id;
    size_type j = b.node_id;
    if (not has_edge(a, b)) {
      if (i >= edges_adj_vec.size()) {
        edges_adj_vec.resize(i+1);
      }
      if (j >= edges_adj_vec.size()) {
        edges_adj_vec.resize(j+1);
      }
      edge_type empty_val;
      edges_adj_vec[i].push_back(j);
      edges_val_vec[i].push_back(edge_value_type());
      edges_adj_vec[j].push_back(i);
      edges_val_vec[j].push_back(edge_value_type());
      edges_count++;
    }
    return Edge(this, a.node_id, b.node_id);
  }

  /**
   * Remove the edges associated with the passed-in nodes from the graph
   * 
   * @Complexity: O(degree(a) + degree(b))
   * 
   * @pre @a and @b are valid nodes
   * 
   * Assuming the edge existed:
   *     @post edges_count is decremented by 1
   *     @post degree of @a a is decremented by 1
   *     @post degree of @a b is decremented by 1
   *     @post associated EdgeIterator is invalidated
   * 
   */ 
  size_type remove_edge(const Node& a, const Node& b) {
      size_type removed = (size_type) 0;
      if (has_edge(a, b) && has_node(a) && has_node(b)) {
        IncidentIterator it = a.edge_begin();
        while (it != a.edge_end()) {
          if ((*it).node2_id == b.node_id) {
            if (edges_adj_vec[a.node_id].size() == 1) {
              edges_adj_vec[a.node_id].clear();
              edges_val_vec[a.node_id].clear();
            } else {
              edges_adj_vec[a.node_id].erase(edges_adj_vec[a.node_id].begin() + it.adj_id);
              edges_val_vec[a.node_id].erase(edges_val_vec[a.node_id].begin() + it.adj_id);
            }
            edges_count--;
            removed = (size_type) 1;
            break;
          }
          ++it;
        }
        // reset iterator for removing edge b->a
        it = b.edge_begin();
        while (it != b.edge_end()) {
          if ((*it).node2_id == a.node_id) {
            if (edges_adj_vec[b.node_id].size() == 1) {
              edges_adj_vec[b.node_id].clear();
              edges_val_vec[b.node_id].clear();
            } else {
              edges_adj_vec[b.node_id].erase(edges_adj_vec[b.node_id].begin() + it.adj_id);
              edges_val_vec[b.node_id].erase(edges_val_vec[b.node_id].begin() + it.adj_id);
            }
            break;
          }
          ++it;
        }
      }
      return removed;

  }


  /**
   * Remove the passed-in edge from the graph
   * 
   * @Complexity: O(degree(a) + degree(b)) where e = <AB>
   * 
   * @pre @a e is a valid edge
   * Assuming the edge existed:
      @post edges_count is decremented by 1
      @post degree of @a a is decremented by 1
      @post degree of @a b is decremented by 1
      @post associated EdgeIterator is invalidated
   */ 
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  }

  /**
   * Remove edge pointed by the passed-in iterator from the graph
   * 
   * @Complexity: O(degree(a) + degree(b)), where (*edge_iter) = <AB>
   * 
   * @pre @a edge_iter is a valid edge iterator
   * Assuming the edge existed:
        @post edges_count is decremented by 1
        @post degree of @a a is decremented by 1
        @post degree of @a b is decremented by 1
        @post associated EdgeIterator(s) and IncidentIterator(s) is invalidated 
   */ 
  edge_iterator remove_edge(edge_iterator edge_iter) {
    return EdgeIterator(this, remove_edge((*edge_iter).node1(), (*edge_iter).node2()));
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_ptr(nullptr), node_id(0) {}

    Node operator*() const {
      // return a node pointing to the graph and with its node_id
      return Node(graph_ptr, graph_ptr->corresp[node_id]);
    }

    NodeIterator& operator++() {
      node_id++;
      return *this;
    }

    bool operator==(const NodeIterator& other) const {
      return (graph_ptr == other.graph_ptr) && (node_id == other.node_id);
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    Graph* graph_ptr;
    size_type node_id;

    // Construct Valid Iterator:
    NodeIterator(const Graph* ptr, size_type id):
      graph_ptr(const_cast<Graph*>(ptr)), node_id(id) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return node_iterator(this, 0);
  }

  node_iterator node_end() const {
    return node_iterator(this, num_nodes());
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
    IncidentIterator() : graph_ptr(nullptr), node_id(0), adj_id(0) {}

    // De-reference the iterator
    Edge operator*() const {
      return Edge(graph_ptr, node_id, graph_ptr->edges_adj_vec[node_id][adj_id]);
    }

    // Step towards the next adjcent node to the original
    IncidentIterator& operator++() {
      adj_id++;
      return *this;
    }

    bool operator==(const IncidentIterator& other) const {
      if (graph_ptr == other.graph_ptr) {
        return (node_id == other.node_id) && (adj_id == other.adj_id);
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    Graph* graph_ptr;
    size_type node_id;
    size_type adj_id; // id of the adj_node in respect to the first node

    // The valid constructor
    IncidentIterator(const Graph* ptr, size_type id1, size_type id2):
        graph_ptr(const_cast<Graph*>(ptr)), node_id(id1), adj_id(id2) {}

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
    EdgeIterator() : graph_ptr(nullptr), node_id(0), adj_id(0) {}

    Edge operator*() const {
      return Edge(graph_ptr, node_id, graph_ptr->edges_adj_vec[node_id][adj_id]);
    }

    // Pre-Increment Operator
    EdgeIterator& operator++() {
      re_adjust_ids();
      return *this; // silence compiler warning
    }

    // Post-Increment operator
    EdgeIterator& operator++(int) {
      re_adjust_ids();
      return *this; // silence compiler warning
    }

    bool operator==(const EdgeIterator& other) const {
      if (graph_ptr == other.graph_ptr) {
        return (node_id == other.node_id) && (adj_id == other.adj_id);
      } else {
        return false;
      }
    }

   private:
    friend class Graph;

    // Same as described for the IncidentIterator we keep a graph pointer, 
    // an id to the fist node of the edge, and an id to the adjecant node in 
    // terms of the first node. 

    Graph* graph_ptr;
    size_type node_id;
    size_type adj_id;

    // Ensure that the ids pointed by the iterator are left at valid positions
    void re_adjust_ids(bool step = true) {
      if (step) {
        adj_id++;
      }
      while (node_id < graph_ptr->edges_adj_vec.size()) {
        while (adj_id < graph_ptr->edges_adj_vec[node_id].size()) {
          if (graph_ptr->edges_adj_vec[node_id][adj_id] > node_id) {
            return;
          }
          ++adj_id;
        }
        adj_id = 0;
        ++node_id;
      }
    }

    // Construct Valid EdgeIterator
    EdgeIterator(const Graph* ptr, size_type id1, size_type id2):
        graph_ptr(const_cast<Graph*>(ptr)), node_id(id1), adj_id(id2) {
          re_adjust_ids(false);
        }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0, 0);
  }
  
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_adj_vec.size(), 0);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::vector<int_node> nodes_vec;
  std::vector<std::vector<size_type>> edges_adj_vec;
  std::vector<std::vector<edge_value_type>> edges_val_vec;  
  std::vector<size_type> corresp;
  size_type edges_count;

  struct int_node {
    Point point;
    node_value_type value;
    size_type node_id;

    int_node(): point(Point()), value(node_value_type()), node_id(0) {}
    int_node(size_type id): point(Point()), value(node_value_type()), node_id(id) {}
    int_node(Point p, size_type id): point(p), value(node_value_type()), node_id(id) {}
    int_node(Point p, node_value_type v, size_type id) : point(p), value(v), node_id(id) {}

    bool operator==(const int_node& other) const {
      return (node_id == other.node_id) && (point == other.point);
    }

  };

};

#endif // CME212_GRAPH_HPP
