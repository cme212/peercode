#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <exception>
#include <unordered_map>
#include <cassert>
#include <utility>

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

  // Predeclare the internal node struct
  struct internal_node;
  // Predeclare the internal edge struct
  struct internal_edge;

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
  Graph() : next_node_index_(0), next_edge_index_(0) {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor. */
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

    // The index value of all invalid nodes
    static const size_type INVALID_ID;

 
     /** Construct an invalid node.
     *
     * Valid nodes are obtained from the Graph class, but it
     * is occasionally useful to declare an @i invalid node, and assign a
     * valid node to it later. For example:
     *
     * @code
     * Graph::Node x;
     * if (...should pick the first node...)
     *   x = graph.node(0);
     * else
     *   x = some other node using a complicated calculation
     * do_something(x);
     * @endcode
     *
     * @post is_valid() == false
     * @post index() == Node::INVALID_ID
     * @post position() throws exception 
     */
    Node() : graph_(nullptr), index_(INVALID_ID) {
      // HW0: YOUR CODE HERE
    }


    /** Determine whether this node is valid. */
    bool is_valid() const {
      // Invalid nodes should have nullptr and index_ = INVALID_ID.
      assert((this->graph_ == nullptr) == (this->index_ == INVALID_ID));

      // To be valid, graph_ should not be nullptr and node_ should be there
      // in graph_.
      if (this->graph_ != nullptr) {
        if (this->graph_->nodes_.find(this->index_)
                       != this->graph_->nodes_.end()) {
          return true;
        }
      }
      // Otherwise, the node is invalid.
      return false;
    }

    /** Return this node's position. 
     *
     * throws std::exception if the node is invalid.
     * TODO Maybe add noexpect(fals)
     */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // invalid node -> TODO costom exception/ runtime_error
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the position of the corresponding internal node.
      return this->graph_->nodes_.at(this->index_).position;
    }

    /** Return this node's index, a number in the range [0, graph_size). 
     *
     * Returns INVALID_ID if the node is invalid.
     **/
    size_type index() const {
      // HW0: YOUR CODE HERE
      return this->get_index_();
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
     * Same graph means the same object in the same address in memory.
     *
     * Invalid nodes are assumed equal, and less than all other nodes.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      // Check if they belong to the same graph and if they have the same index

      // If one of them is invalid, return true only if the other is invalid
      // as well.
      if (not(this->is_valid()) or not(n.is_valid())) {
        return not(this->is_valid()) and not(n.is_valid());
      }
      // Check equality of the graphs' pointers and the indices.
      return ((this->graph_ == n.graph_) and (this->index_ == n.index_));
    }


    /** Test whether this node and @a n are not equal.
     *
     * Unequal nodes are those who have == operator returns false.
     */
    bool operator!=(const Node& n) const {
      // Define in terms of ==
      return not(*this == n);
    }


    /** Test whether this node is less than @a n in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any geometric meaning.
     *
     * Invalid nodes are assumed equal, and less than all other nodes.
     * The nodes in a graph with a smaller address are less than others.
     * Within a graph, the nodes with less indices are less.
     *
     * The node ordering relation obeys trichotomy: For any two nodes x
     * and y, exactly one of x == y, x < y, and y < x is true.
     *
     */
    bool operator<(const Node& n) const {
      // HW0: YOUR CODE HERE

      // If the other node is invalid, this node cannot be less than the
      // other one.
      if (not(n.is_valid())) return false;
      
      // Otherwise, and if this node is invalid, it is less than the other one.
      else if (not(this->is_valid())) return true;

      // Otherwise, compare based on graphs, then on indices.
      else return std::tie(this->graph_, this->index_) < \
             std::tie(n.graph_, n.index_);
    }

    /** The remaining comparison operators
     * 
     * These are defined in terms of <.
     */
    bool operator>(const Node& n) const { return n < *this; }
    bool operator<=(const Node &n) const {return not(*this > n);}
    bool operator>=(const Node &n) const { return not(*this < n);}

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    
    // Allow the testing fixture to access Node's private member data and 
    // functions. This is not meant to test these private members, it is to
    // provide some mock integration test-cases.
    // Also, these members should not be used except occasionally, and their
    // tests are not meant to stay too long, because the internal implementation
    // is not guaranteed to stay the same.
    friend class GraphPointFixture;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    /** Construct a node (may be valid) 
     *  
     *  @pre @a graph_ == nullptr iff @a index_ == Node::INVALID_ID
     *
     *  @post index_ = @a index_
     *  @post graph_ = @a graph_
     *
     *  This constructor is private to cope with the requirement that valid
     *  nodes should only be construct-able from graph methods.
     *  DONT use it with any public method without this requirenment in mind.
     * */
    Node(const Graph* graph_, size_type index_) : graph_(graph_),
                                                        index_(index_) {}
    // A pointer to the graph containing this node, nullptr if invalid.
    const Graph* graph_;
    // The index of the node in the graph, Node::INVALID_ID if invalid.
    size_type index_;

    /** internal getter for the index
     * 
     * this is used to get the updated index for this node which is either it
     * is original index or Node::INVALID_ID depending on whether the node was
     * invalidated smoehow without modifying index_ or not.
     * This can happen because of the use of the proxy design pattern for the
     * Node class.
     *
     * @post result = index_ if the node is still valid
     *              = Node::INVALID_ID if it is no longer valid (invalidated
     *              from the graph)
     *
     */
    size_type get_index_() const {
      // If it is invalid/ was invalidated, retrun sizetype(-1).
      if (not (this->is_valid())) return INVALID_ID;
      // Otherwise, return the current index_.
      return this->index_;
    }

    /** internal getter for the graph pointer
     *
     * this is used to get the updated graph pointer for this node. (nullptr
     * if it was invalidated somehow).
     *
     * @post result =  graph_ if the node is still valid
     *              = nullptr if it is no longer valid (invalidated
     *              from the graph)
     *
     */
    const Graph* get_graph_() const {
      // If it is invalid/ was invalidated, retrun nullptr.
      if (not (this->is_valid())) return nullptr;
      // Otherwise, return the current graph pointer.
      return this->graph_;
    }

  };


  /** A hash object for Node.
   *
   * Uses the graph pointer as well as the index.
   * This is to allow using the node with map.
   *
   * Idealy, this is closely related with the == method in the Node class.
   *
   */
  struct NodeHasher{
    std::size_t operator()(const Node& n) const
    {
      // Compute individual hash values for first,
      // second and combine them using XOR
      // and bit shifting:
      std::hash<Graph const*> hash_graph;
      return ((hash_graph(n.graph_) >> 1)
               ^ (std::hash<size_type>()(n.get_index_()) << 1));
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // The size of the graph's nodes map
    return this->nodes_.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.position() == @a position
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position) {
    // HW0: YOUR CODE HERE
    // Create a new internal node and store it at nodes map with a new index.
    internal_node new_node(position);
    this->nodes_.emplace(this->next_node_index_, new_node);
    // Then, return a representation of this added internal node.
    // And, increment the next index member.
    return this->node(this->next_node_index_++);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    // First, check to see if the node's index falls in the valid range of
    // indices of this graph.
    if ((this->nodes_).find(n.index_) != (this->nodes_).end()) {
      // Then, check to see if this node is the same as the node of the same
      // index in this graph.
      return n == this->node(n.index());
    }
    // Otherwise, the node is not part of the graph.
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * TODO Handle the case when the node was removed?
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Create a node representation of the ith internal node.
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

    // The index of all invalid edges
    static const size_type INVALID_ID;
 
    /** Construct an invalid Edge. */
    Edge() : graph_(nullptr), index_(INVALID_ID), is_ordered_(true) {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge 
     *
     * if the edge is invalid, throws std::exception
     */
    Node node1() const {
      // HW0: YOUR CODE HERE
      // invalid edge -> exception.
      // TODO costom exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the first node of the Edge.
      // If is_order, then the first node of the edge is the first node in
      // the internal representation.
      return this->is_ordered_ ? 
        this->graph_->node(this->graph_->edges_.at(this->index_).node1_id):
        this->graph_->node(this->graph_->edges_.at(this->index_).node2_id);
    }

    /** Return the other node of this Edge 
     *
     * if the edge is invalid, throws std::exception
     */
    Node node2() const {
      // HW0: YOUR CODE HERE
      // invalid edge -> exception.
      // TODO costom exception
      if (not(this->is_valid())) { 
        std::exception e;
        throw e;
      } 
      // Return the second node of the Edge.
      // If is_order, then the second node of the edge is the first node in
      // the internal representation.
      return this->is_ordered_ ? 
        this->graph_->node(this->graph_->edges_.at(this->index_).node2_id):
        this->graph_->node(this->graph_->edges_.at(this->index_).node1_id);
    }

    /** Return true if this edge is valid, false otherwise. */
    bool is_valid() const {
      // Invalid edges should have nullptr and index_ = INVALID_ID.
      assert((this->graph_ == nullptr) == (this->index_ == INVALID_ID));

      // TO be valid, graph_ should not be nullptr and edge_ should be there
      // in graph_.
      if (this->graph_ != nullptr) {
        auto it = this->graph_->edges_.find(this->index_);
        if (it != this->graph_->edges_.end()) {
          // Also, the nodes should be still in the graph.
          size_type node1_id = (it->second).node1_id;
          size_type node2_id = (it->second).node2_id;
          if (this->graph_->has_node(this->graph_->node(node1_id)) and
              this->graph_->has_node(this->graph_->node(node2_id))) {
            return true;
          }
        }
      }
      // Otherwise, the edge is invalid.
      return false;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     *
     * Invalid nodes are equal. They are not equal to any valid edge.
     */
    bool operator==(const Edge& e) const {
      // If one of them is invalid, return true only if the other is invalid
      // as well.
      if (not(this->is_valid()) or not(e.is_valid())) {
        return not(this->is_valid()) and not(e.is_valid());
      }

      // Check using nodes equality
        if ((this->node1() == e.node1() and this->node2() == e.node2())
            or (this->node2() == e.node1() and this->node1() == e.node2())) {
          return true;
        }
        // No match, return false.
        return false; 
      }

    /** Test whether this edge and @a n are not equal.
     *
     * Unequal edges are edges with "not ==" return true.
     */
    bool operator!=(const Edge& e) const {
      // Define in terms of ==
      return not(*this == e);
    }



    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      // HW0: YOUR CODE HERE

      // If the other edge is invalid, this edge cannot be less than the
      // other one.
      if (not(e.is_valid())) return false;

      // Otherwise, and if this edge is invalid, it is less than the other one.
      else if (not(this->is_valid())) return true;

      // Otherwise, compare based on graphs, then on indices.
      else return std::tie(this->graph_, this->index_) < \
             std::tie(e.graph_, e.index_);
    }

    /** The remaining comparison operators
     * 
     * These are defined in terms of <.
     */
    bool operator>(const Edge& e) const { return e < *this; }
    bool operator<=(const Edge &e) const {return not(*this > e);}
    bool operator>=(const Edge &e) const { return not(*this < e);}

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;

    // Allow the testing fixture to access Edge's private member data and
    // functions. This is not meant to test these private members, it is to
    // provide some mock integration test-cases.
    // Also, these members should not be used except occasionally, and their
    // tests are not meant to stay too long, because the internal implementation
    // is not guaranteed to stay the same.
    friend class GraphPointFixture;

    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    /** Construc an edge.
     *
     * @post index_ = @a index_
     * @post graph_ = @a graph_
     * @post is_ordered_ = @a is_ordered_
     *
     */
    Edge(const Graph* graph_, size_type index_, bool is_ordered_ = true) \
            : graph_(graph_), index_(index_), is_ordered_(is_ordered_) {}

    // A pointer to the graph containing this edge, nullptr if invalid.
    const Graph* graph_;
    // The index of the edge in the graph, INVALID_ID if invalid.
    size_type index_;
    
    /** A boolean specifies whether node1 and node2 are in the same ordered
     * used in the internal implementaion of the edges container in the graph.
     * so, if is_ordered: node1() = graph_->edge(index_).node1 
     *                and node2() = graph_->edge(index_).node2
     *     else         : node1() = graph_->edge(index_).node2 
     *                and node2() = graph_->edge(index_).node1
     */
    bool is_ordered_;

    /** internal getter for the index
     * 
     * this is used to get the updated index for this edge which is either it
     * is original index or Edge::INVALID_ID depending on whether the edge was
     * invalidated smoehow without modifying index_ or not.
     * This can happen because of the use of the proxy design pattern for the
     * Edge class.
     *
     * @post result = index_ if the edge is still valid
     *              = Edge::INVALID_ID if it is no longer valid (invalidated
     *              from the graph)
     *
     */
    size_type get_index_() const {
      // If it is invalid/ was invalidated, retrun sizetype(-1).
      if (not (this->is_valid())) return INVALID_ID;
      // Otherwise, return the current index_.
      return this->index_;
    }

    /** internal getter for the graph pointer
     *
     * this is used to get the updated graph pointer for this edge. (nullptr
     * if it was invalidated somehow).
     *
     * @post result =  graph_ if the edge is still valid
     *              = nullptr if it is no longer valid (invalidated
     *              from the graph)
     *
     */
    const Graph* get_graph_() const {
      // If it is invalid/ was invalidated, retrun nullptr.
      if (not (this->is_valid())) return nullptr;
      // Otherwise, return the current graph pointer.
      return this->graph_;
    }


  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    // The size of the graph's edges map
    return this->edges_.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * The order of the nodes in the returned edge matches that of the first
   * time the edge was added.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Create an edge representation of the ith internal node.
    return this->edge_(i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    // get_edge_index returns Edge::INVALID_ID if the edge is not in the graph.
    return this->get_edge_index_(a,b) != Edge::INVALID_ID;
  }

  /** Add an edge to the graph, or return the current edge if it already exists.
   * @pre @a a and @a b are distinct valid nodes of this graph
   * @return an Edge object e with e.node1() == @a a and e.node2() == @a b
   * @post has_edge(@a a, @a b) == true
   * @post If old has_edge(@a a, @a b), new num_edges() == old num_edges().
   *       Else,                        new num_edges() == old num_edges() + 1.
   *
   * Can invalidate edge indices -- in other words, old edge(@a i) might not
   * equal new edge(@a i). Must not invalidate outstanding Edge objects.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge add_edge(const Node& a, const Node& b) {
    // HW0: YOUR CODE HERE
    size_type index = this->get_edge_index_(a,b);
    bool is_ordered;
    if (index == Edge::INVALID_ID) {
      internal_edge new_edge(a.index(), b.index());
      index = (this->next_edge_index_)++;
      this->edges_.emplace(index, new_edge);
      this->add_edge_index_(a,b,index);
      is_ordered = true;
    } else {
      is_ordered = this->node(this->edges_.at(index).node1_id) == a;
    }
    return this->edge_(index, is_ordered);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodes_.clear();
    this->edges_.clear();
    this->clear_edges_ids_();
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

  // Allow the testing fixture to access Graph's private member data and 
  // functions. This is not meant to test these private members, it is to
  // provide some mock integration test-cases.
  // Also, these members should not be used except occasionally, and their
  // tests are not meant to stay too long, because the internal implementation
  // is not guaranteed to stay the same.
  friend class GraphPointFixture;

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /** Internal node struct. */
  struct internal_node {
  
    /** Construct an internal node. */
    internal_node(const Point position) : position(position) {};
  
    const Point position; // The position of the node.
  };

  /** An unordered map that holds the nodes of the graph.
   * Always: (nodes_[index]) == this->nodes(index)
   */
  std::unordered_map<size_type, internal_node> nodes_;

  // The index to be used next time a node is added to this graph.
  size_type next_node_index_; 

  /** Internal edge struct. */
  struct internal_edge {

    /** Construct an internal edge. */
    internal_edge(size_type node1_id, size_type node2_id) : 
                  node1_id(node1_id), node2_id(node2_id) { }

    size_type node1_id; // Pointer to the first node of the edge.
    size_type node2_id; // Pointer to the second node of the edge.
  };

  /** An unordered map that holds the edge objects of the graph.
   * Always edges_[index] == this->edge(index) 
   * and edges_[index].node1 == this->edge(index).node1()
   * */
  std::unordered_map<size_type, internal_edge> edges_;
  
  /** An unordered map that holds the indices of the edges in the graph.
   *
   * This is used to:
   * - provide O(1) access to the index given a pair of nodes
   * - provide O(1) check if a pair of nodes is an edge in this graph
   * - provide direct access to the indices of the edges that contain a
   *   specific node.
   *
   * It stores each edge index twice, once with (n1,n2) and once with (n2,n1)
   * This is necessary to support the 3rd use-case above. It uses nested maps
   * so for every node a, edges_ids_[a] is a map that maps every other node b
   * connected to a - to the index of the edge (a,b). Note that the same is
   * stored for edges_id_[b].
   * So always edges_ids_[a][b] == edges_ids_[b][a] == index of (a,b).
   *
   * This with the map edges_ form a special data structure that maps every
   * index to a pair of nodes with a specific order in O(1), and maps every
   * pair of nodes that represents an edge to the index regardless of the
   * order of the pair.
   *
   * Another abroach was to use CSR or CSC form of sparse matrices.
   *
   * IMPORTANT: In general, this mamber variable is not to be used directly
   * within the Graph class, instead the private methods add_edge_index_, 
   * erase_edge_index, get_edge_index_, and clear_edges_ids_.
   */
  std::unordered_map<Node,
    std::unordered_map<Node, size_type, NodeHasher>,
      NodeHasher> edges_ids_;

  // The index to be used next time a edge is added to this graph.
  size_type next_edge_index_;


  /** Return an edge with index @a i. Whose nodes order depends on @is_ordered
   * and the internal_edge stored at i order
   * @pre 0 <= @a i < num_edges()
   *
   * @post IF is_ordered : result_edge.node1() == this->edges_[i].node1_id
   *                       result_edge.node2() == this->edges_[i].node2_id
   *       ELSE          : result_edge.node1() == this->edges_[i].node2_id
   *                       result_edge.node2() == this->edges_[i].node1_id
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge_(size_type i, bool is_ordered = true) const {
    // HW0: YOUR CODE HERE
    // Create an edge representation of the ith internal node.
    return Edge(this, i, is_ordered);
  }


  /** Return the index of the edge (n1,n2) or Edge::INVALID_ID if not 
   *  found.
   *
   * if one node is invalid or not part of the graph, or if the edge is not
   * part of the graph, this returns Edge::INA:ID_ID. Otherwise it returns
   * the index of the edge (n1,n2). Returns the same index for both (n1,n2) and
   * (n2,n1).
   *
   * @post if has_edge(n1,n2) : edges_[result].node1() = either n1 or n2 
   *                            edges[result].node2() = the other one.
   *       else               : result = Edge::INVALID_ID
   */
  size_type get_edge_index_(Node n1, Node n2) const {
    // The fact that both n1,n2 and n2,n1 are stored is utilized.
    // Search for n1
    auto it = this->edges_ids_.find(n1);
    if (it == this->edges_ids_.end()) return Edge::INVALID_ID;
    // Search for n2 in the edges of n1
    auto nested_edges_ids = it->second;
    auto nested_it = (nested_edges_ids.find(n2));
    if (nested_it == nested_edges_ids.end()) return Edge::INVALID_ID;
    // Return the index
    return nested_it->second;
  }

  /** Add the index of the edge (n1,n2)/(n2,n1) to the graph.
   *  @pre has_node(n1) == true
   *  @pre has_node(n2) == true
   *  @post get_edge_index_(n1, n2) = index
   *  @post if has_edge(n1,n2) : num_edges_ = num_edges_+ 1
   */
  void add_edge_index_(Node n1, Node n2, size_type index) {
    // if (not(this->has_edge(n1,n2))) (this->num_edges_)++;
    // Store the index twic, once for each order.
    this->edges_ids_[n1][n2] = index;
    this->edges_ids_[n2][n1] = index;
  }

  /** Remove the index of the edge (n1,n2) from the graph edge indices.
   *  
   *  @post has_edge(n1,n2) == false
   *  @post if has_edge(n1,n2) : num_edges_ = num_edges_ - 1
   *
   *  Doesn't handle edges_[get_edge_index(n1, n2)] 
   */
  void erase_edge_index_(Node n1, Node n2) {
    // if (not(this->has_edge(n1,n2))) (this->num_edges_)--;
    this->edges_ids_[n1].erase(n2);
    this->edges_ids_[n2].erase(n1);
  }

  /** Clear the indices of the graph edges 
   *
   *  @post has_edge(n1,n2) == false for every n1 and n2
   *  @post num_edges_ == 0
   */
  void clear_edges_ids_() {
    this->edges_ids_.clear();
    // this->num_edges = 0;
  }


};

// When using YouCompleteMe plugin in vim, different compilation settings are
// used, in these compilation settings, these two lines should not be included
// for the program to compile well, another pair of lines are used instead, in a
// different file. The key here is to define the macro YCM when using
// YouCompleteMe compilation and possibly NOT TO DEFINE it when using the
// default Compilation. However, defining it will probably not cause problems in
// general.
// See _Graph.cpp for more about the other path.
#ifndef YCM
constexpr Graph::size_type Graph::Node::INVALID_ID = Graph::size_type(-1);
constexpr Graph::size_type Graph::Edge::INVALID_ID = Graph::size_type(-1);
#endif

#endif // CME212_GRAPH_HPP
