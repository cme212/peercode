#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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
  // TODO: maybe declare <size_type> instead of <unsigned>?
  
  // Vectors to hold the index and the Point of each node in the Graph
  // The immutable node number inside the Node class informs the vector index
  std::vector<unsigned> node_index_;
  std::vector<Point*> node_Point_;

  /* Vector to hold edges; an edge is defined as an unordered set since graph
  is undirected and has no self-loops. Set here will always have size 2, but
  could be easily generalized to hypergraphs with sets of bigger size. */
  std::vector<std::unordered_set<unsigned>> edges_;
  /* Unordered map to unordered sets for O(1) complexity for has_edge and
  add_edge. The map key is always the smallest node number in the edge. */
  std::unordered_map<unsigned, std::unordered_set<unsigned>> pairs_;

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
  Graph() {}
    // HW0: YOUR CODE HERE
    // TODO
    // No need of anything here actually

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
      // HW0: YOUR CODE HERE
      // TODO
      // No need to do anything here for now
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // TODO
      // Get pointer to Point stored in the Graph class
      const Point* point_ptr = (this->graph_ptr_->node_Point_)[this->
                               node_number_];
      return *point_ptr;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // TODO
      // Get index from Graph class and assure that is is within the range
      size_type index = (this->graph_ptr_->node_index_)[this->node_number_];
      assert(index < graph_ptr_->size());
      return index;
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
      // TODO
      // Equal only if both graph and node index are equal
      bool same_graph = (this->graph_ptr_ == n.graph_ptr_);
      bool same_node = ((this->graph_ptr_->node_index_)[this->node_number_]
                         == (n.graph_ptr_->node_index_)[n.node_number_]);
      return (same_graph and same_node);
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
      // TODO
      // First compare Graph pointers
      if (this->graph_ptr_ < n.graph_ptr_) { return true; }
      else if (this->graph_ptr_ > n.graph_ptr_) { return false; }
      // If graph is the same, just compare node index
      return ((this->graph_ptr_->node_index_)[this->node_number_]
               < (n.graph_ptr_->node_index_)[n.node_number_]);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // TODO
    // Pointer to Graph and unique node_number_ (index for vectors in Graph)
    size_type node_number_;
    const Graph* graph_ptr_;

    /* Private constructor for Graph class to create valid Nodes */
    Node(size_type node_number, const Graph* graph_ptr):
      node_number_(node_number), graph_ptr_(graph_ptr) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // TODO
    // Both vectors must be of equal size and equal the number of nodes
    assert(this->node_index_.size() == this->node_Point_.size());
    return this->node_index_.size();
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
    // TODO
    //std::cout << "\nAdding Point:"; // TODO
    //std::cout << "\nx = " << position.elem[0];
    //std::cout << "\ny = " << position.elem[1];
    //std::cout << "\nz = " << position.elem[2] << "\n\n";
    // Create new point to make address will not get fluxed
    Point* new_point = new Point;
    *new_point = Point(position.elem[0], position.elem[1], position.elem[2]);
    // Add new node with index equal to current size of the graph
    size_type old_graph_size = this->size();
    this->node_index_.push_back(old_graph_size);
    this->node_Point_.push_back(new_point);
    // Test whether the addition worked as expected
    assert(this->size() == old_graph_size+1);
    //std::cout << "\nCurrent list size: " << this->node_Point_.size(); //TODO
    //for (unsigned i = 0; i < this->node_Point_.size(); i++) {
        //std::cout << "\ni = " << i; // TODO
        //std::cout << "\n" << this->node_Point_[i];
        //std::cout << "\nx = " << (*(this->node_Point_[i])).elem[0];
        //std::cout << "\ny = " << (*(this->node_Point_[i])).elem[1];
        //std::cout << "\nz = " << (*(this->node_Point_[i])).elem[2] <<"\n\n";
    //}
    return Node(old_graph_size, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // TODO
    // Graph must be the same and node number below graph size
    bool same_graph = (this == n.graph_ptr_);
    bool contain_node = (this->size() > n.node_number_);
    return (same_graph and contain_node);
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // TODO
    // Test if node does exist and then return it
    assert(i < this->num_nodes());
    return Node(i, this);
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
      // HW0: YOUR CODE HERE
      // TODO
      // No need to implement anything here for now
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // TODO
      return Node(this->node_num1_, this->graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // TODO
      return Node(this->node_num2_, this->graph_ptr_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // TODO
      // If not in same graph, cannot be equal
      bool same_graph = (this->graph_ptr_ == e.graph_ptr_);
      if (not same_graph) { return false; }
      // Graph is undirected so we need to test 4 equalities
      bool node1_eq1 = (this->node_num1_ == e.node_num1_);
      bool node1_eq2 = (this->node_num1_ == e.node_num2_);
      bool node2_eq1 = (this->node_num2_ == e.node_num1_);
      bool node2_eq2 = (this->node_num2_ == e.node_num2_);
      // Both cases in which edge is the same (undirected graph)
      bool node_match = (node1_eq1 and node2_eq2) or (node1_eq2 and node2_eq1);
      return node_match;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      // TODO
      /* Definition of "<" here: the bigger edge is the one with the "bigger"
      Graph pointer. In case they are equal, the bigger is the one with the
      biggest lowest node number. In case the smallest node numbers are the
      same, the bigger is the one with biggest highest node number. If all
      these 3 are equal, then the edges are equal. */
      // First compare Graph pointers
      if (this->graph_ptr_ < e.graph_ptr_) { return true; }
      else if (this->graph_ptr_ > e.graph_ptr_) { return false; }
      // If equal graph pointers (same Graph), compare lowest node
      size_type this_min = std::min(this->node_num1_, this->node_num2_);
      size_type e_min = std::min(e.node_num1_, e.node_num2_);
      if (this_min < e_min) { return true; }
      else if (this_min > e_min) { return false; }
      // If lowest node is equal, compare highest node
      size_type this_max = std::max(this->node_num1_, this->node_num2_);
      size_type e_max = std::max(e.node_num1_, e.node_num2_);
      if (this_max < e_max) { return true; }
      // If did not return true so far, it must be false
      return false;
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // TODO
    // Pointer to Graph and the node numbers of both nodes that form the edge
    const Graph* graph_ptr_;
    size_type node_num1_;
    size_type node_num2_;

    /* Private constructor for Graph to create valid edges */
    Edge(const Graph* graph_ptr, size_type node_num1, size_type node_num2):
        graph_ptr_(graph_ptr), node_num1_(node_num1), node_num2_(node_num2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // TODO
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // TODO
    // Make sure that the i-th edge exists and get the associated set
    assert(i < this->num_edges());
    std::unordered_set<unsigned> set_nodes = this->edges_[i];
    // Make sure edge is well-defined (2 nodes) and get the node numbers
    assert(set_nodes.size() == 2);
    unsigned node_num1 = *(set_nodes.begin());
    unsigned node_num2 = *(++set_nodes.begin());
    return Edge(this, node_num1, node_num2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // TODO
    // Get min/max node number due to organization of 'pairs_'
    size_type min_node = std::min(a.node_number_, b.node_number_);
    size_type max_node = std::max(a.node_number_, b.node_number_);
    // If no edge from min_node to a node of higher number, edge not existent
    if (pairs_.find(min_node) == pairs_.end()) { return false; }
    // If there are edges from min_node, check if one of the is with max_node
    std::unordered_set<unsigned> set_nodes = pairs_.find(min_node)->second;
    if (set_nodes.find(max_node) == set_nodes.end()) { return false; }
    // If both if statements were false, the edge exists, so we return true
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
    // TODO
    size_type old_edge_count = this->num_edges();
    if (not this->has_edge(a,b)) {
        // If edge does not exist, we have to add it to the Graph
        // First we add it to 'edges_'
        std::unordered_set<unsigned> new_edge ({a.node_number_,
                                                b.node_number_});
        this->edges_.push_back(new_edge);
        // Now we add it to 'pairs_'
        size_type min_node = std::min(a.node_number_, b.node_number_);
        size_type max_node = std::max(a.node_number_, b.node_number_);
        this->pairs_[min_node].insert(max_node);
        // Make sure everything went well
        assert(this->num_edges() == old_edge_count + 1);
    }
    return Edge(this, a.node_number_, b.node_number_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // TODO
    // Allocated memory for points, so we free the memory now
    for (unsigned i = 0; i < this->node_Point_.size(); i++) {
        delete this->node_Point_[i];
    }
    // Simply clear all data to remove all nodes and edges
    this->node_index_.clear();
    this->node_Point_.clear();
    this->edges_.clear();
    this->pairs_.clear();
    // Assure it worked as expected
    assert(num_nodes() == 0 && num_edges() == 0);
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
  // TODO

};

#endif // CME212_GRAPH_HPP
