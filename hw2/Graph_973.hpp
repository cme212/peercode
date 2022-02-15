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
template <typename V = int, typename E = double>
class Graph {
 public:
  // TODO: HW1-2 additions
  typedef V node_value_type;
  typedef E edge_value_type;
  /** Type of indexes and sizes.
      Return type of Graph::Node::index(), Graph::num_nodes(),
      Graph::num_edges(), and argument type of Graph::node(size_type) */
  using size_type = unsigned;
  // Make type declaration less verbose
  using my_set = std::unordered_set<size_type>;
  using pairs_map_type = std::unordered_map<size_type, my_set>;
  using value_map_type = std::unordered_map<size_type, std::unordered_map<
                           size_type, edge_value_type>>;

 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // Vectors to hold the value and the Point of each node in the Graph
  std::vector<node_value_type> node_values_;
  std::vector<Point*> node_Points_;

  // TODO: HW2 addition
  // Map to map the immutable node ID inside the Node class to vector index
  std::unordered_map<size_type, int> id_to_index_;
  // Map to map the vector index to the immutable node ID inside the Node class
  std::unordered_map<size_type, size_type> index_to_id_;

  /* Vector to hold edges; an edge is defined as an unordered set since graph
  is undirected and has no self-loops. Set here will always have size 2, but
  could be easily generalized to hypergraphs with sets of bigger size. */
  std::vector<my_set> edges_;

  /* Unordered map to unordered sets for O(1) complexity for has_edge and
  add_edge. The map key is either node ID number in the edge. */
  pairs_map_type pairs_;

  /* Unordered map to unordered map to edge_value_type for O(1) edge value
  addition, deletion and lookup. Keys are respectively min and max node ID */
  value_map_type edge_values_;

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

    // TODO: HW2 addition
    /** Return this node's position. */
    Point& position() {
      // HW0: YOUR CODE HERE
      // Get pointer to Point stored in the Graph class
      Point* point_ptr = (graph_ptr_->node_Points_)[index()];
      return *point_ptr;
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // Get pointer to Point stored in the Graph class
      const Point* point_ptr = (graph_ptr_->node_Points_)[index()];
      return *point_ptr;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // TODO: HW change
      // HW0: YOUR CODE HERE
      // Get index from Graph class and assure that is is within the range
      size_type index = graph_ptr_->id_to_index_[n_id_];
      assert((0 <= index) and (index < graph_ptr_->size()));
      return index;
    }

    // TODO: not necessary anymore; for debugging
    //my_set all_neighbours() {
    //  if (graph_ptr_->pairs_.find(n_id_) != graph_ptr_->pairs_.end()) {
    //    return graph_ptr_->pairs_[n_id_];
    //  } else { my_set empty; return empty;}
    //}

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value() {
      return (*graph_ptr_).node_values_[index()];
    }
    const node_value_type& value() const {
      return (*graph_ptr_).node_values_[index()];
    }
    size_type degree() const {
      pairs_map_type& pairs = (*graph_ptr_).pairs_;
      if (pairs.find(n_id_) == pairs.end()) { return 0; }
      else { return pairs[n_id_].size(); }
    };
    incident_iterator edge_begin() const {
      return IncidentIterator(graph_ptr_, n_id_, 0);
    };
    incident_iterator edge_end() const {
      return IncidentIterator(graph_ptr_, n_id_, degree());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Equal only if both graph and node index are equal
      bool same_graph = (graph_ptr_ == n.graph_ptr_);
      bool same_node = (n_id_ == n.n_id_);
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
      if (graph_ptr_ < n.graph_ptr_) { return true; }
      else if (graph_ptr_ > n.graph_ptr_) { return false; }
      // If graph is the same, just compare node index
      return (n_id_ < n.n_id_);
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    // Pointer to Graph and unique, immutable n_id_
    size_type n_id_;
    Graph* graph_ptr_;

    /* Private constructor for Graph class to create valid Nodes */
    Node(size_type n_id, Graph* graph_ptr):
      n_id_(n_id), graph_ptr_(graph_ptr) {}
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // Both vectors must be of equal size and equal the number of nodes
    assert(node_values_.size() == node_Points_.size());
    return node_values_.size();
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
    // TODO: HW2 changes
    // HW0: YOUR CODE HERE
    // Create new point to make address will not get fluxed
    Point* new_point = new Point;
    *new_point = Point(position.elem[0], position.elem[1], position.elem[2]);
    // Get current size of graph and current number of node IDs
    size_type old_graph_size = size();
    size_type old_number_ids = id_to_index_.size();
    // Push back given point and value for new node
    node_Points_.push_back(new_point);
    node_values_.push_back(value);
    // Associate new node ID with vector index
    id_to_index_[old_number_ids] = old_graph_size;
    index_to_id_[old_graph_size] = old_number_ids;
    // Test whether the addition worked as expected
    assert(size() == old_graph_size+1);
    assert(index_to_id_.size() == old_graph_size+1);
    assert(id_to_index_.size() == old_number_ids+1);
    return Node(old_graph_size, this);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // TODO: HW2 changes
    // HW0: YOUR CODE HERE
    // Graph must be the same and node number below graph size
    int index = n.index();
    bool same_graph = (this == n.graph_ptr_);
    // We set index as -1 when node has been removed
    bool contain_node = (index != -1);
    assert(index < int(size()));
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
    assert(i < num_nodes());
    return Node(index_to_id_.at(i), const_cast<Graph*>(this));
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
      return Node(node_id1_, graph_ptr_);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(node_id2_, graph_ptr_);
    }

    // TODO: HW2 additions
    double length() {
        return norm(node2().position() - node1().position());
    }
    edge_value_type& value() {
      size_type min_id = std::min(node_id1_, node_id2_);
      size_type max_id = std::min(node_id1_, node_id2_);
      return graph_ptr_->edge_values_[min_id][max_id];
    }
    const edge_value_type& value() const {
      size_type min_id = std::min(node_id1_, node_id2_);
      size_type max_id = std::min(node_id1_, node_id2_);
      return graph_ptr_->edge_values_[min_id][max_id];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // If not in same graph, cannot be equal
      bool same_graph = (graph_ptr_ == e.graph_ptr_);
      if (not same_graph) { return false; }
      // Graph is undirected so we need to test 4 equalities
      bool node1_eq1 = (node_id1_ == e.node_id1_);
      bool node1_eq2 = (node_id1_ == e.node_id2_);
      bool node2_eq1 = (node_id2_ == e.node_id1_);
      bool node2_eq2 = (node_id2_ == e.node_id2_);
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
      if (graph_ptr_ < e.graph_ptr_) { return true; }
      else if (graph_ptr_ > e.graph_ptr_) { return false; }
      // If equal graph pointers (same Graph), compare lowest node
      size_type this_min = std::min(node_id1_, node_id2_);
      size_type e_min = std::min(e.node_id1_, e.node_id2_);
      if (this_min < e_min) { return true; }
      else if (this_min > e_min) { return false; }
      // If lowest node is equal, compare highest node
      size_type this_max = std::max(node_id1_, node_id2_);
      size_type e_max = std::max(e.node_id1_, e.node_id2_);
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
    size_type node_id1_;
    size_type node_id2_;

    /* Private constructor for Graph to create valid edges */
    Edge(Graph* graph_ptr, size_type node_id1, size_type node_id2):
        graph_ptr_(graph_ptr), node_id1_(node_id1), node_id2_(node_id2) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Make sure that the i-th edge exists and get the associated set
    assert(i < num_edges());
    my_set set_nodes = edges_[i];
    // Make sure edge is well-defined (2 nodes) and get the node numbers
    assert(set_nodes.size() == 2);
    size_type node_id1 = *(set_nodes.begin());
    size_type node_id2 = *(++set_nodes.begin());
    return Edge(const_cast<Graph*>(this), node_id1, node_id2);
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
    if (pairs_.find(a.n_id_) == pairs_.end()) { return false; }
    if (pairs_.find(b.n_id_) == pairs_.end()) { return false; }
    // Check if such edges are between a and b
    my_set edges_a = pairs_.find(a.n_id_)->second;
    my_set edges_b = pairs_.find(b.n_id_)->second;
    if (edges_a.find(b.n_id_) == edges_a.end()) { return false; }
    if (edges_b.find(a.n_id_) == edges_b.end()) { return false; }
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
    // TODO: HW2 changes
    // HW0: YOUR CODE HERE
    size_type old_edge_count = num_edges();
    if (not has_edge(a,b)) {
      // If edge does not exist, we have to add it to the Graph
      // First we add it to 'edges_'
      my_set new_edge ({a.n_id_, b.n_id_});
      edges_.push_back(new_edge);
      // Now we add it to 'pairs_' in duplicate (this helps eval degree)
      pairs_[a.n_id_].insert(b.n_id_);
      pairs_[b.n_id_].insert(a.n_id_);
      // TODO: for now edge value must be specified by user externally
      //size_type min_id = std::min(a.n_id_, b.n_id_);
      //size_type max_id = std::max(a.n_id_, b.n_id_);
      //edge_values_[min_id][max_id] = Edge(this, a.n_id_, b.n_id_).length();
      // Make sure everything went well
      assert(num_edges() == old_edge_count + 1);
    }
    return Edge(this, a.n_id_, b.n_id_);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // Allocated memory for points, so we free the memory now
    for (unsigned i = 0; i < node_Points_.size(); i++) {
        delete node_Points_[i];
    }
    // Simply clear all data to remove all nodes and edges
    node_values_.clear();
    node_Points_.clear();
    edges_.clear();
    pairs_.clear();
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
    Node operator*() const {
      return graph_ptr_->node(index_);
    }
    NodeIterator& operator++() {
      index_++;
      return *this;
    }
    bool operator==(const NodeIterator& iter) const {
      bool eq_graph = (graph_ptr_ == iter.graph_ptr_);
      bool eq_index = (index_ == iter.index_);
      return (eq_graph and eq_index);
    }

   private:
    const Graph* graph_ptr_;
    size_type index_;
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /** Construct a valid NodeIterator. */
    NodeIterator(const Graph* graph_ptr, size_type index):
      graph_ptr_(graph_ptr), index_(index) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator(this, size());
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
    Edge operator*() const {
      if (edge_number_ < neighbours_.size()) {
        return Edge(const_cast<Graph*>(graph_ptr_),
                    n_id_, neighbours_[edge_number_]);
      } else { return Edge(); }
    // if ( iter_ == empty_.end() ) { return Edge(); }
    // return Edge(const_cast<Graph*>(graph_ptr_), n_id_, *iter_);
    }
    IncidentIterator& operator++() {
      edge_number_++;
      return *this;
    }
    bool operator==(const IncidentIterator& i_iter) const {
        bool eq_graph = (graph_ptr_ == i_iter.graph_ptr_);
        bool eq_node = (n_id_ == i_iter.n_id_);
        bool eq_edge = (edge_number_ == i_iter.edge_number_);
        return (eq_graph and eq_node and eq_edge);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_ptr_;
    size_type n_id_;
    size_type edge_number_;
    std::vector<size_type> neighbours_;

    /** Construct a valid IncidentIterator. */
    IncidentIterator(Graph* graph_ptr, size_type n_id,
      size_type edge_number): graph_ptr_(graph_ptr),
      n_id_(n_id), edge_number_(edge_number) {
        // If node has edges, initialize a vector of its neighbours
        if (graph_ptr_->pairs_.find(n_id) != graph_ptr_->pairs_.end()) {
          my_set set_neighbours = graph_ptr_->pairs_[n_id];
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
        bool eq_graph = (graph_ptr_ == iter.graph_ptr_);
        bool eq_edge = (edge_number_ == iter.edge_number_);
        return (eq_graph and eq_edge);
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    const Graph* graph_ptr_;
    size_type edge_number_;
    /** Construct an invalid EdgeIterator. */
    EdgeIterator(const Graph* graph_ptr, size_type edge_number):
      graph_ptr_(graph_ptr), edge_number_(edge_number) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.size());
  }

  // TODO: HW2 node and edge removals
   /** Given a node, remove it from the graph and all the edges with it
    * @post if old has_node(n) == TRUE,
            then new num_nodes() == old num_nodes() - 1;
    * @post has_node(n) == FALSE;
    * @complexity: O(degree(n)) (O(1) for node but degree times O(1) for edge)
    * @invalidates: Node(n.n_id_, this);
                    old node_iterator n_it s.t (*n_it) == n
                    old incident_iterator on node n
    * @return: 0 if not deleted (node inexistent) or 1 if deleted successfully
  */
  size_type remove_node(const Node& n) {
    // Return zero if node is not anymore in graph (or never was)
    if (!has_node(n)) { return 0; }
    // Get index of node to be deleted and of last node
    size_type index = n.index();
    size_type last_index = num_nodes() - 1;
    Node last_node = node(last_index);
    // Copy info from last node to current index and pop back
    node_values_[index] = node_values_[last_index];
    node_Points_[index] = node_Points_[last_index];
    node_values_.pop_back();
    node_Points_.pop_back();
    // Correct the mappings of ID to index and index to ID
    id_to_index_[n.n_id_] = -1; // Invalidated ID
    id_to_index_[last_node.n_id_] = index;
    index_to_id_[index] = last_node.n_id_;
    index_to_id_.erase(last_index);
    // Delete all edges
    if (pairs_.find(n.n_id_) != pairs_.end()) {
        my_set neighs = pairs_[n.n_id_];
        for (auto iter = neighs.begin(); iter != neighs.end(); ++ iter) {
            // Remove the edges with each neighbour of the node
            remove_edge(Node(n.n_id_, this), Node(*iter, this));
        }
    }
    return 1;
  };

  /** Given a node, remove it from the graph and all the edges with it
    * @post if old has_node(n) == TRUE,
            then new num_nodes() == old num_nodes() - 1;
    * @post has_node(n) == FALSE;
    * @complexity: O(degree(n)) (O(1) for node but degree times O(1) for edge)
    * @invalidates: Node(n.n_id_, this);
                    old node_iterator n_it s.t (*n_it) == n
                    old incident_iterator on node n
    * @return: 0 if not deleted (node inexistent) or 1 if deleted successfully
  */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it != node_end()) { remove_node(*n_it); }
    return node_begin();
  };

  /** Given 2 nodes, remove the edge between them if any exist.
    * @post if old has_edge(n1,n2) == TRUE, then
            new num_edges() == old num_edges() - 1;
    * @post has_edge(n1,n2) == FALSE;
    * @complexity: O(1) (deletion in sets and maps)
    * @invalidates: Edge(this, n1.n_id_, n2.n_id_);
                    Edge(this, n2.n_id_, n1.n_id_);
                    old edge_iterator e_it s.t (*e_it) == (either edge above)
    * @return: 0 if not deleted (edge inexistent) or 1 if deleted successfully
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // If edge is not present, just return false
    if (!has_edge(n1, n2)) { return 0; }
    // Continue code to remove edge
    my_set edge_set ({n1.n_id_, n2.n_id_});
    // Delete edge for each of the 3 data structures that hold edge info
    auto iter = std::find(edges_.begin(), edges_.end(), edge_set);
    if (iter != edges_.end()) {
      edges_.erase(iter);
    }
    // Delete twice in pairs_ since edge is stored in two keys: both IDs
    if (pairs_.find(n1.n_id_) != pairs_.end()) {
      auto iter = pairs_[n1.n_id_].find(n2.n_id_);
      if (iter != pairs_[n1.n_id_].end()) {
        pairs_[n1.n_id_].erase(iter);
      }
    }
    if (pairs_.find(n2.n_id_) != pairs_.end()) {
      auto iter = pairs_[n2.n_id_].find(n1.n_id_);
      if (iter != pairs_[n2.n_id_].end()) {
        pairs_[n2.n_id_].erase(iter);
      }
    }
    // Delete in edge_values_ require us to know correct key order
    size_type min_id = std::min(n1.n_id_, n2.n_id_);
    size_type max_id = std::max(n1.n_id_, n2.n_id_);
    if (edge_values_.find(min_id) != edge_values_.end()) {
      auto iter = edge_values_[min_id].find(max_id);
      if (iter != edge_values_[min_id].end()) {
        edge_values_[min_id].erase(max_id);
      }
    }
    assert(!has_edge(n1, n2));
    return 1;
  };

  /** Given 2 nodes, remove the edge between them if any exist.
    * @post if old has_edge(n1,n2) == TRUE, then
            new num_edges() == old num_edges() - 1;
    * @post has_edge(n1,n2) == FALSE;
    * @complexity: O(1) (deletion in sets and maps)
    * @invalidates: Edge(this, n1.n_id_, n2.n_id_);
                    Edge(this, n2.n_id_, n1.n_id_);
                    old edge_iterator e_it s.t (*e_it) == (either edge above)
    * @return: 0 if not deleted (edge inexistent) or 1 if deleted successfully
  */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(), e.node2());
  };

  /** Given 1 edge_iterator, remove the edge it points to if it exists.
    * @post if old has_edge(n1,n2) == TRUE, then
            new num_edges() == old num_edges() - 1;
    * @post has_edge(n1,n2) == FALSE;
    * @complexity: O(1) (deletion in sets and maps)
    * @invalidates: Edge(this, n1.n_id_, n2.n_id_);
                    Edge(this, n2.n_id_, n1.n_id_);
                    old edge_iterator e_it s.t (*e_it) == (either edge above)
    * @return edge_begin() since pointer arithmetic is unsafe
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it != edge_end()) { remove_edge(*e_it); }
    return edge_begin();
  };



 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

};

#endif // CME212_GRAPH_HPP
