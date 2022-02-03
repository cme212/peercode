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
template <typename V = int>
class Graph {
 public:
  // HW1 addition // TODO
  using node_value_type = V;
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  // Make type declaration less verbose
  using my_set = std::unordered_set<size_type>;
  using my_map = std::unordered_map<size_type, my_set>;
 
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
  // Vectors to hold the value and the Point of each node in the Graph
  // The immutable node number inside the Node class informs the vector index
  std::vector<node_value_type> node_value_;
  std::vector<Point*> node_Point_;

  /* Vector to hold edges; an edge is defined as an unordered set since graph
  is undirected and has no self-loops. Set here will always have size 2, but
  could be easily generalized to hypergraphs with sets of bigger size. */
  std::vector<my_set> edges_;
  /* Unordered map to unordered sets for O(1) complexity for has_edge and
  add_edge. The map key is always the smallest node number in the edge. */
  std::unordered_map<size_type, my_set> pairs_;

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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {}
    // HW0: YOUR CODE HERE
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
  class Node: private totally_ordered<Node> {
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
      // No need to do anything here for now
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // Get pointer to Point stored in the Graph class
      const Point* point_ptr = (this->graph_ptr_->node_Point_)[this->
                               node_index_];
      return *point_ptr;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // Get index from Graph class and assure that is is within the range
      size_type index = this->node_index_;
      assert(index < graph_ptr_->size());
      return index;
    }

    // TODO: not necessary anymore; for debugging
    //my_set all_neighbours() {
    //  if (graph_ptr_->pairs_.find(node_index_) != graph_ptr_->pairs_.end()) {
    //    return graph_ptr_->pairs_[node_index_];
    //  } else { my_set empty; return empty;}
    //}

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // TODO
    node_value_type& value() {
      return (*this->graph_ptr_).node_value_[this->node_index_];
    }
    const node_value_type& value() const {
      return (*this->graph_ptr_).node_value_[this->node_index_];
    }
    size_type degree() const {
      my_map& pairs = (*(this->graph_ptr_)).pairs_;
      if (pairs.find(this->node_index_) == pairs.end()) { return 0; }
      else { return pairs[this->node_index_].size(); }
    };
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_ptr_, node_index_, 0);
    };
    incident_iterator edge_end() const {
      return IncidentIterator(graph_ptr_, node_index_, this->degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Equal only if both graph and node index are equal
      bool same_graph = (this->graph_ptr_ == n.graph_ptr_);
      bool same_node = (this->node_index_ == n.node_index_);
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
      // First compare Graph pointers
      if (this->graph_ptr_ < n.graph_ptr_) { return true; }
      else if (this->graph_ptr_ > n.graph_ptr_) { return false; }
      // If graph is the same, just compare node index
      return (this->node_index_ < n.node_index_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Pointer to Graph and unique node_index_ (index for vectors in Graph)
    size_type node_index_;
    Graph* graph_ptr_;

    /* Private constructor for Graph class to create valid Nodes */
    Node(size_type node_index, Graph* graph_ptr):
      node_index_(node_index), graph_ptr_(graph_ptr) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // Both vectors must be of equal size and equal the number of nodes
    assert(this->node_value_.size() == this->node_Point_.size());
    return this->node_value_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post if user passed new_index, result_node.index() == new_index
   * Otherwise, result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, 
    const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Create new point to make address will not get fluxed
    // TODO: HW1 changes
    Point* new_point = new Point;
    *new_point = Point(position.elem[0], position.elem[1], position.elem[2]);
    // Get current size of graph
    size_type old_graph_size = this->size();
    // Push back given point and value for new node
    this->node_Point_.push_back(new_point);
    this->node_value_.push_back(value);
    // Test whether the addition worked as expected
    assert(this->size() == old_graph_size+1);
    return Node(old_graph_size, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // Graph must be the same and node number below graph size
    bool same_graph = (this == n.graph_ptr_);
    bool contain_node = (this->size() > n.node_index_);
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
    // Test if node does exist and then return it
    assert(i < this->num_nodes());
    return Node(i, const_cast<Graph*>(this));
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
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      // No need to implement anything here for now
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(this->node_num1_, this->graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(this->node_num2_, this->graph_ptr_);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
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
    // Pointer to Graph and the node numbers of both nodes that form the edge
    Graph* graph_ptr_;
    size_type node_num1_;
    size_type node_num2_;

    /* Private constructor for Graph to create valid edges */
    Edge(Graph* graph_ptr, size_type node_num1, size_type node_num2):
        graph_ptr_(graph_ptr), node_num1_(node_num1), node_num2_(node_num2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Make sure that the i-th edge exists and get the associated set
    assert(i < this->num_edges());
    my_set set_nodes = this->edges_[i];
    // Make sure edge is well-defined (2 nodes) and get the node numbers
    assert(set_nodes.size() == 2);
    size_type node_num1 = *(set_nodes.begin());
    size_type node_num2 = *(++set_nodes.begin());
    return Edge(const_cast<Graph*>(this), node_num1, node_num2);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // If no edge from a/b to some other node, edge not existent
    if (pairs_.find(a.node_index_) == pairs_.end()) { return false; }
    if (pairs_.find(b.node_index_) == pairs_.end()) { return false; }
    // Check if such edges are between a and b
    my_set edges_a = pairs_.find(a.node_index_)->second;
    my_set edges_b = pairs_.find(b.node_index_)->second;
    if (edges_a.find(b.node_index_) == edges_a.end()) { return false; }
    if (edges_b.find(a.node_index_) == edges_b.end()) { return false; }
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
    size_type old_edge_count = this->num_edges();
    if (not this->has_edge(a,b)) {
      // If edge does not exist, we have to add it to the Graph
      // First we add it to 'edges_'
      my_set new_edge ({a.node_index_, b.node_index_});
      this->edges_.push_back(new_edge);
      // Now we add it to 'pairs_' in duplicate (this helps eval degree)
      this->pairs_[a.node_index_].insert(b.node_index_);
      this->pairs_[b.node_index_].insert(a.node_index_);
      // Make sure everything went well
      assert(this->num_edges() == old_edge_count + 1);
    }
    return Edge(this, a.node_index_, b.node_index_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Allocated memory for points, so we free the memory now
    for (unsigned i = 0; i < this->node_Point_.size(); i++) {
        delete this->node_Point_[i];
    }
    // Simply clear all data to remove all nodes and edges
    this->node_value_.clear();
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
  class NodeIterator: private totally_ordered<NodeIterator> {
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
    // TODO
    Node operator*() const {
      return graph_ptr_->node(index_);
    }
    NodeIterator& operator++() {
      index_++;
      return *this;
    }
    bool operator==(const NodeIterator& iter) const {
      bool eq_graph = (this->graph_ptr_ == iter.graph_ptr_);
      bool eq_index = (this->index_ == iter.index_);
      return (eq_graph and eq_index);
    }

   private:
    const Graph* graph_ptr_;
    size_type index_;
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /** Construct a valid NodeIterator. */
    // TODO
    NodeIterator(const Graph* graph_ptr, size_type index):
      graph_ptr_(graph_ptr), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // TODO
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator(this, this->size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered <IncidentIterator> {
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
    // TODO
    Edge operator*() const {
      if (edge_number_ < neighbours_.size()) {
        return Edge(const_cast<Graph*>(graph_ptr_),
                    node_index_, neighbours_[edge_number_]);
      } else { return Edge(); }
    // if ( iter_ == empty_.end() ) { return Edge(); }
    // return Edge(const_cast<Graph*>(graph_ptr_), node_index_, *iter_);
    }
    IncidentIterator& operator++() {
      edge_number_++;
      return *this;
    }
    bool operator==(const IncidentIterator& i_iter) const {
        bool eq_graph = (this->graph_ptr_ == i_iter.graph_ptr_);
        bool eq_node = (this->node_index_ == i_iter.node_index_);
        bool eq_edge = (this->edge_number_ == i_iter.edge_number_);
        return (eq_graph and eq_node and eq_edge);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    // TODO
    Graph* graph_ptr_;
    size_type node_index_;
    size_type edge_number_;
    std::vector<size_type> neighbours_;
    /** Construct a valid IncidentIterator. */
    IncidentIterator(Graph* graph_ptr, size_type node_index,
      size_type edge_number): graph_ptr_(graph_ptr),
      node_index_(node_index), edge_number_(edge_number) {
        // If node has edges, initialize a vector of its neighbours
        if (graph_ptr_->pairs_.find(node_index) != graph_ptr_->pairs_.end()) {
          my_set set_neighbours = graph_ptr_->pairs_[node_index];
          // Store unordered data as vector for easy index access
          neighbours_.insert(neighbours_.end(), set_neighbours.begin(),
                             set_neighbours.end());
        }
      }
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered <EdgeIterator> {
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
    // TODO
    Edge operator*() const {
      my_set nodes_set = graph_ptr_->edges_[edge_number_];
      std::vector<size_type> nodes_vec(nodes_set.begin(), nodes_set.end());
      //nodes_vec.insert(nodes_vec.end(), nodes_set.begin(), nodes_set.end());
      assert(nodes_vec.size() == 2);
      return Edge(const_cast<Graph*>(graph_ptr_), nodes_vec[0], nodes_vec[1]);
    }
    EdgeIterator& operator++() {
      edge_number_++;
      return *this;
    }
    bool operator==(const EdgeIterator& iter) const {
        bool eq_graph = (this->graph_ptr_ == iter.graph_ptr_);
        bool eq_edge = (this->edge_number_ == iter.edge_number_);
        return (eq_graph and eq_edge);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    // TODO
    const Graph* graph_ptr_;
    size_type edge_number_;
    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph* graph_ptr, size_type edge_number):
      graph_ptr_(graph_ptr), edge_number_(edge_number) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // TODO
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, this->edges_.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
