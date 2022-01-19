#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"
#define assert_GRAPH_PTR_EXIST_TX (assert(this->graph_ptr != nullptr))
#define assert_GRAPH_PTR_EXIST_TX_general(x) (assert((x).graph_ptr != nullptr))

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
class Graph {
 public:
    /** Predeclaration of Node and edge type. */
    class Node;
    class Edge;
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // TX: Add sets of nodes and edges
  std::vector <Point*> node_collection; //TODO: see if this fix works
  std::vector<std::array<size_type, 2> > edge_collection;


 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph;

  /** Predeclaration of Node type. */
  // class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  /** Predeclaration of Edge type. */
  // class Edge;
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



  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    // TX: Here I will intialize with zero vector
    this->node_collection = std::vector<Point*>(); //TODO: check correctness
    this->edge_collection = std::vector<std::array<size_type, 2> >();

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
      // Predeclaration
  private:
   // Allow Graph to access Node's private member data and functions.
   friend class Graph;
   Graph* graph_ptr;
   size_type node_idx;


   // HW0: YOUR CODE HERE
   // Use this space to declare private data members and methods for Node
   // that will not be visible to users, but may be useful within Graph.
   // i.e. Graph needs a way to construct valid Node objects

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

      //TX: assume it has nullptr
      this->graph_ptr = nullptr;

    }

    Node(const Graph* ptr_tmp, size_type idx_tmp) {
          // HW0: YOUR CODE HERE

          //TX: Done
          this->graph_ptr = const_cast<Graph*>(ptr_tmp);
          this->node_idx = idx_tmp;

      }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      // TX: Add very elementary read of a point position.
        assert_GRAPH_PTR_EXIST_TX;
        return *((*(this->graph_ptr)).node_collection[this->node_idx]);
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
        assert_GRAPH_PTR_EXIST_TX;
      // HW0: YOUR CODE HERE
      // TX: We assume this is a maintained node_idx type
      return this->node_idx;
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

      // TX: added one more macro to deal with comparing with nullptr
      assert_GRAPH_PTR_EXIST_TX;
      assert_GRAPH_PTR_EXIST_TX_general(n);
      return this->node_idx == n.node_idx;
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

      //TX: copied from == part
        assert_GRAPH_PTR_EXIST_TX;
        assert_GRAPH_PTR_EXIST_TX_general(n);
        return this->node_idx < n.node_idx;
    }


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE

    // TX: add a very easy function

    return this->node_collection.size();
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
    // HW0: YOUR CODE HERE


    // TX: Constructs node, add to node_collection. Then return
    Node node_tmp;
    node_tmp.graph_ptr = this;
    node_tmp.node_idx  = this->size(); // Only works because of zero indexing
    Point* tmp_ptr =  const_cast<Point*>(&position);
    (this->node_collection).push_back(tmp_ptr); // Add the point
    // TODO: is the above line legal? Do I encounter illegal change to the Point instance?
    return node_tmp;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    assert_GRAPH_PTR_EXIST_TX_general(n);

    // TX: Determine if node is a valid index
    return n.graph_ptr == this && n.node_idx < this->size();
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert (i < this->size());
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
  private:
      // TX: add graph pointer and two node indices. Same logic here.
      // Allow Graph to access Edge's private member data and functions.
      friend class Graph;
      Graph* graph_ptr;
      size_type node_idx_1;
      size_type node_idx_2;
      size_type edge_idx; // Order in the graph.edge_collection attribute
      // HW0: YOUR CODE HERE
      // Use this space to declare private data members and methods for Edge
      // that will not be visible to users, but may be useful within Graph.
      // i.e. Graph needs a way to construct valid Edge objects
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      this->graph_ptr = nullptr;
    }

    Edge(const Graph* ptr_tmp, size_type idx_tmp) {
      // HW0: YOUR CODE HERE

      //TX: Done
      this->graph_ptr = const_cast<Graph*>(ptr_tmp);
      this->node_idx_1 = (*ptr_tmp).edge_collection[idx_tmp][0];
      this->node_idx_2 = (*ptr_tmp).edge_collection[idx_tmp][1];
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // TX: added access
        assert_GRAPH_PTR_EXIST_TX;

        return (*this->graph_ptr).node(this->node_idx_1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // TX: added access
        assert_GRAPH_PTR_EXIST_TX;

        return (*this->graph_ptr).node(this->node_idx_2);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //TX: added my code here
        assert_GRAPH_PTR_EXIST_TX;
        assert_GRAPH_PTR_EXIST_TX_general(e);
        bool cond1 = (this->node_idx_1 == e.node_idx_1) && (this->node_idx_2 == e.node_idx_2);
        bool cond2 = (this->node_idx_1 == e.node_idx_2) && (this->node_idx_2 == e.node_idx_1);

      return cond1 || cond2;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      //TX: added my code
        assert_GRAPH_PTR_EXIST_TX;
        assert_GRAPH_PTR_EXIST_TX_general(e);
      return (this->edge_idx < e.edge_idx);
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // TX: added my code
    return (this->edge_collection).size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    //TX: Done by defining constructor above
    return Edge(this, i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    // TX: added code by simple enumeration
      assert_GRAPH_PTR_EXIST_TX_general(a);
      assert_GRAPH_PTR_EXIST_TX_general(b);
      size_type node_idx_a = a.node_idx;
      size_type node_idx_b = b.node_idx;

    for (unsigned int idx = 0; idx < this->num_edges(); idx++){
        size_type node_1 = this->edge_collection[idx][0];
        size_type node_2 = this->edge_collection[idx][1];

        bool cond1 = (node_1 == node_idx_a) && (node_2 == node_idx_b);
        bool cond2 = (node_1 == node_idx_b) && (node_2 == node_idx_a);
        if (cond1 || cond2) {
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
    // HW0: YOUR CODE HERE
      // TX: this block of code is from has edge
      assert_GRAPH_PTR_EXIST_TX_general(a);
      assert_GRAPH_PTR_EXIST_TX_general(b);
      size_type node_idx_a = a.node_idx;
      size_type node_idx_b = b.node_idx;

      for (unsigned int idx = 0; idx < this->num_edges(); idx++){
          size_type node_1 = this->edge_collection[idx][0];
          size_type node_2 = this->edge_collection[idx][1];

          bool cond1 = (node_1 == node_idx_a) && (node_2 == node_idx_b);
          bool cond2 = (node_1 == node_idx_b) && (node_2 == node_idx_a);
          if (cond1 || cond2) {
              return Edge(this , idx);
          }
      }
      //TX: if runs to this line, then no match found
      // TX: Constructs node, add to node_collection. Then return
      Edge edge_tmp;
      edge_tmp.graph_ptr = this;
      edge_tmp.edge_idx  = this->num_edges(); // Only works because of zero indexing
      edge_tmp.node_idx_1  = node_idx_a; // Only works because of zero indexing
      edge_tmp.node_idx_2  = node_idx_b; // Only works because of zero indexing


      (this->edge_collection).push_back({node_idx_a, node_idx_a}); // Add the point
      return edge_tmp;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
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

};

#endif // CME212_GRAPH_HPP
