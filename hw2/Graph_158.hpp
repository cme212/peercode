#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cmath>
#include <exception>
#include <type_traits>
#include <unordered_map>
#include <cassert>
#include <utility>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * @tparam T Type of the nodes' value attribute. It is assumed to have default
 *           initializer/constructor defined. default type: int
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V = int, typename E = int>
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
  using graph_type = Graph<V,E>;

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

  /** Type of nodes' value attribute */
  using node_value_type = V;

   /** Type of edges' value attribute */
  using edge_value_type = E; //

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
  class Node : private totally_ordered<Node> {
   public:

    /** predeclare the hasher struct for the node */
    struct Hasher;
    /** Synonym for the node hasher type. */
    using hasher_type = Hasher;

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
    Node() : graph_(nullptr), index_(INVALID_ID)//, user_id(INVALID_ID) 
    {
      // HW0: YOUR CODE HERE
    }


    /** Determine whether this node is valid. */
    bool is_valid() const {
      // Invalid nodes should have nullptr and index_ = INVALID_ID.
      assert((this->graph_ == nullptr) == (this->index_ == INVALID_ID));

      // To be valid, graph_ should not be nullptr and node_ should be there
      // in graph_.
      if (this->graph_ != nullptr) {
        auto it = this->graph_->nodes_.find(this->index_);
        if ( it != this->graph_->nodes_.end()) {
         return true;
        }
      }
      // Otherwise, the node is invalid.
      return false;
    }

    /** Return this node's position. 
     *
     * throws std::exception if the node is invalid.
     * TODO_hw0 Maybe add noexpect(false)
     */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      // invalid node -> TODO_hw0 costom exception/ runtime_error
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the position of the corresponding internal node.
      return this->graph_->nodes_.at(this->index_).position;
    }

 
    /** Return a read/write reference to this node's position. 
     *
     * throws std::exception if the node is invalid.
     * TODO_hw0 Maybe add noexpect(false)
     */
    Point& position() {
      // HW0: YOUR CODE HERE
      // invalid node -> TODO_hw0 costom exception/ runtime_error
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the position of the corresponding internal node.
      return const_cast<Graph*>(this->graph_)->nodes_.at(this->index_).position;
    }

   /** Return this node's index, a number in the range [0, graph_size). 
     *
     * Returns INVALID_ID if the node is invalid.
     */
    size_type index() const {
      // HW0: YOUR CODE HERE
      size_type user_id = this->graph_->nodes_.at(this->get_index_()).user_id;
      return user_id;
    }

    /** Return a read/write reference to the node's value
     *
     * throws std::exception if the node is invalid.
     */
    node_value_type& value() {
      // Invalid node -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }

      // Return a reference to the costom value of the corresponding internal
      // node.
      return const_cast<Graph*>(this->graph_)->nodes_.at(this->index_).value;
    }

    /** Return a read-only reference to the node's value
     *
     * throws std::exception if the node is invalid.
     */
    const node_value_type& value() const {
      // Invalid node -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return a constant reference to the costom value of the corresponding
      // internal node.
      return this->graph_->nodes_.at(this->index_).value;
    }

    /** Return the number of edges in the graph containing the node.
     *
     * throws std::exception if the node is invalid.
     */
    size_type degree () const {
      // Invalid node -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the size of the map that contains the edges containing this
      // node.
      return this->graph_->get_node_edges_ids_(*this).size();
    }

    /** Returns a read/write iterator that points to the first edge incident to
     *  the node.
     *
     * throws an exception if the node is invalid.
     *
     * @post if this->is_valid() : 
     *       for (; result != this->edge_end(); ++result)
     *                        (*result).node1() == *this
     * @post if this->is_valid() :
     *       for (count=0; result != this->edge_end(); ++result, ++count);
     *       count == this->degree()
     *
     */
    incident_iterator edge_begin() const {
      // Invalid node => exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }

      return IncidentIterator(this->graph_, *this, true); 
    }

    /** Return a read/write iterator that points one past the last edge incident
     *  to this node.
     *
     * throws an exception if the node is invalid.
     *
     * @post if this->is_valid() : 
     *       for (count=0, it=this->edge_begin(); result != it; ++it, ++count);
     *       count == this->degree()
     */
    incident_iterator edge_end() const { 
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
        return IncidentIterator(this->graph_, *this, false); 
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

     /** A hash object for Node.
     *
     * Uses the graph pointer as well as the index.
     * This is to allow using the node with maps.
     *
     * Idealy, this is closely related with the == method in the Node class.
     *
     */
    struct Hasher{
      std::size_t operator()(const Node& n) const
      {
        // Compute individual hash values for first,
        // second and combine them using XOR
        // and bit shifting:
        std::hash<graph_type const*> hash_graph;
        return ((hash_graph(n.graph_) >> 1)
                 ^ (std::hash<size_type>()(n.get_index_()) << 1));
      }
    };

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
    Node(const graph_type* graph_, size_type index_) : 
                   graph_(graph_), index_(index_) {
      assert (graph_ == this->graph_);
      assert (graph_ != nullptr);
      assert (this->graph_ != nullptr);
    }
    // A pointer to the graph containing this node, nullptr if invalid.
    const graph_type* graph_;
    // The index of the node in the graph, Node::INVALID_ID if invalid.
    size_type index_;
    
    //TODO (uid_inside the Node to invalidate?)
    // The user index of the node in the graph when this node object is created
    // Should be constant, but can be copied via copy constructors.
    // This is added to invalidate every node object that represent the same
    // node but whose user_id doesn't match the user_id of the graph's internal
    // representation of the node. This happens because the latter changes
    // during the program execution.
    // For example:
    // Assume graph is a graph_type of 8 nodes:
    // node_type n0 = graph.node(4);
    // node_type n1 = graph.node(7);
    // graph.remove(n0);
    // Here, dependong on the implementation, n1 internal representation may 
    // have user_id = 4 for exmple (if we used swap and pop). So:
    // n1.index() will return 4 which is not-intuitive. With this member
    // we can track that n1 has user_id == 7 when it was created. so when the
    // internal representation changes the user_id, and we try n1.index() we can
    // detect the mismatch and then invalidate n1 and throw an exception.
    // size_type user_id;

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
    const graph_type* get_graph_() const {
      // If it is invalid/ was invalidated, retrun nullptr.
      if (not (this->is_valid())) return nullptr;
      // Otherwise, return the current graph pointer.
      return this->graph_;
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
   * @param[in] value The new node's value
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   * @post result_node.position() == @a position
   * @post result_node.value() == @a value
   * @post edges_ids_ contains result_node as a key.
   *
   * Complexity: O(1) amortized operations.
   */
  node_type add_node(const Point& position, 
                const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    // Create a new internal node and store it at nodes map with a new index.
    size_type user_id = num_nodes();
    internal_node new_node(position, user_id, value);
    this->nodes_.emplace(this->next_node_index_, new_node);
    
    assert(user_id == nodes_uids_to_ids_.size());
    
    this->nodes_uids_to_ids_.push_back(this->next_node_index_);
    // Then, initialize the node incident edges map and return a representation
    // of this added internal node. And, increment the next index member.
    size_type index = this->next_node_index_++;
    node_type n = this->node_(index);
    initialize_node_edges_ids_(n);
    return n;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const node_type& n) const {
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
   *
   * Complexity: O(1).
   */
  node_type node(size_type i) const {
    // HW0: YOUR CODE HERE
    // Create a node representation of the ith internal node.
    size_type index = nodes_uids_to_ids_.at(i);
    (void) 1;
    return node_type(this, index);
  }

   /** Return the node with index @a i.
   * @pre nodes_.find(@a i) != nodes_.end()
   * @post result_node.index() == nodes_.at(i).user_id
   *
   *
   * Complexity: O(1).
   */
  node_type node_(size_type i) const {
    // HW0: YOUR CODE HERE
    // Create a node representation of the ith internal node.
    return node_type(this, i);
  }
 
  /** Remove a node from the graph
   *
   *  Remove the node @a n from the graph if it exists. Returns 1 if
   *  the removal succeeds, and 0 otherwise.
   *
   *  For all previously created nodes: We have 3 cases:
   *    (1) The node objects corresponding to @a n: These are invalidated i.e.
   *        They will return false when calling is_valid() on them, and they
   *        will throw an exception when trying to access their attributes.
   *    (2) The nodes whose public index (user index) was changed during the
   *        removal process: Here the problem is that these nodes will result in
   *        different indices when calling index() depending on whether it was
   *        called before or after the removal of @a n.
   *
   *        Note: why the public index changes? because upon removal, the graph
   *        size is lowered by 1, so the range of public indices should be 
   *        adjusted to be in [0,num_nodes()), which means, almost with every
   *        node removal, at leaset a public index will be modified.
   *    (3) Other nodes: Nothing happen, and all public methods will have the
   *        same effect before and after the removal.
   *
   * For all previously created iterators: Nothing changes except what was
   *        mentioned about the nodes, and the fact that any iterator pointing
   *        to @a n will be invalid and the behavior will be undefined.
   *
   *  @post if has_node(@a n) result == 1 else result == 0
   *  @post old num_nodes() == new num_nodes() + result
   *  @post has_node(@a n) == false
   *  @post if old m == @a n then not new m.is_valid()
   *  @post if ai == old node(i), aj == old node(j): 0 <= i,j < old num_nodes(),
   *           ai != @a n, and aj != @a n
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  @post if a == old node(old num_nodes()-1) then not a.is_valid()
   *  @post if old (*it) == old node(i) 0 <= i < old num_nodes()-1,
   *           old (*it) != @a n, (it is node_itrator) and 
   *           new node(i).is_valid() then 
   *           new (*it) == old node(i)
   *  @post if old e.node1() == @a n or old e.node2() == @a n 
   *        then 
   *            (1) new not has_node(e)
   *            (2) new not e.is_valid()
   *  @post old num_edges() == new num_edges() + result*(old n.degree())
   *
   *  Note: not node.is_valid() is equivalent to:
   *             (1) new m.index() throws exception
   *             (2) new m.position() throws exception
   *             (3) new m.value() throws exception
   *             (4) new m == Node()
   *
   *  Complexity: O((@a n).degree() + 1)
   */
  size_type remove_node(const node_type& n) {
    // Only remove nodes in the graph
    if (this->has_node(n)) {

      // First remove the incident edges
      for(auto it = n.edge_begin(); it != n.edge_end();) {
        it = this->remove_edge(it);
      }
      // Get the node's user index and internal index
      size_type user_id = n.index();
      size_type index = this->nodes_uids_to_ids_.at(user_id);

      // Get the user index and internal index of the node with the largest
      // user index
      size_type last_user_id = this->nodes_uids_to_ids_.size()-1;
      size_type last_index = this->nodes_uids_to_ids_[last_user_id];

      // Swap and then pop from all of the data structures
      this->nodes_uids_to_ids_[user_id] = last_index;
      this->nodes_.at(last_index).user_id = user_id;
      this->nodes_uids_to_ids_.pop_back();
      this->nodes_.erase(index);

      // The node was removed => return 1 
      return 1;
    } else {
      // The node is not part of the graph
      return 0;
    }
  } 

  /** Remove a node from the graph
   *
   *  Remove the node *(@a n_it) from the graph if it exists. Returns an
   *  iterator to the next node if the removal succeeds, and the same iterator
   *  otherwise.
   *
   *  For all previously created nodes: We have 3 cases:
   *    (1) The node objects corresponding to *(@a n_it): These are invalidated
   *        i.e.They will return false when calling is_valid() on them, and they
   *        will throw an exception when trying to access their attributes.
   *    (2) The nodes whose public index (user index) was changed during the
   *        removal process: Here the problem is that these nodes will result in
   *        different indices when calling index() depending on whether it was
   *        called before or after the removal of *(@a n_it).
   *
   *        Note: why the public index changes? because upon removal, the graph
   *        size is lowered by 1, so the range of public indices should be 
   *        adjusted to be in [0,num_nodes()), which means, almost with every
   *        node removal, at leaset a public index will be modified.
   *    (3) Other nodes: Nothing happen, and all public methods will have the
   *        same effect before and after the removal.
   *
   * For all previously created iterators: Nothing changes except what was
   *        mentioned about the nodes, and the fact that any iterator pointing
   *        to *(@a n_it) will be invalid and the behavior will be undefined.
   *
   *  @post if has_node(*(@a n_it)) old num_nodes() == new num_nodes() + 1
   *        else old num_nodes() == new num_nodes()
   *  @post has_node(*(@a n_it)) == false
   *  @post if old m == *(@a n_it) then not new m.is_valid()
   *  @post if ai == old node(i), aj == old node(j): 0 <= i,j < old num_nodes(),
   *           ai != *(@a n_it), and aj != *(@a n_it)
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  @post if a == old node(old num_nodes()-1) then not a.is_valid()
   *  @post if old (*it) == old node(i) 0 <= i < old num_nodes()-1,
   *           old (*it) != *(@a n_it), (it is node_itrator) and 
   *           new node(i).is_valid() then 
   *           new (*it) == old node(i)
   *  @post if old e.node1() == *(@a n_it) or old e.node2() == *(@a n_it) 
   *        then 
   *            (1) new not has_node(e)
   *            (2) new not e.is_valid()
   *  @post if has_node(*(@a n_it)) 
   *        then old num_edges() == new num_edges() + (old n.degree())
   *
   *  @post if has_node(*(@a n_it)), then if we define set1 and set2 as:
   *        let it = result
   *        set1 = empty set : while (it != node_end()) set1.add(*it_);++it;
   *        let it = old @a n_it
   *        set2 = empty set : while (it != node_end()) set2.add(*it);++it;
   *        then: set1 == set2
   *        
   *
   *  Note: not node.is_valid() is equivalent to:
   *             (1) new m.index() throws exception
   *             (2) new m.position() throws exception
   *             (3) new m.value() throws exception
   *             (4) new m == Node()
   *  Complexity: O((*(@a n_it)).degree() + 1)
   */
  node_iterator remove_node(node_iterator n_it) {
    // Copy the iterator
    auto it = n_it;
    // Increment the copy
    ++it;
    if (this->remove_node(*n_it)) {
      return it;
    }
    return n_it;
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
    node_type node1() const {
      // HW0: YOUR CODE HERE
      // invalid edge -> exception.
      // TODO_hw0 costom exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the first node of the Edge.
      // If is_order, then the first node of the edge is the first node in
      // the internal representation.

      return this->is_ordered_ ? 
        this->graph_->node_(this->graph_->edges_.at(this->index_).node1_id):
        this->graph_->node_(this->graph_->edges_.at(this->index_).node2_id);
    }

    /** Return the other node of this Edge 
     *
     * if the edge is invalid, throws std::exception
     */
    node_type node2() const {
      // HW0: YOUR CODE HERE
      // invalid edge -> exception.
      // TODO_hw0 costom exception
      if (not(this->is_valid())) { 
        std::exception e;
        throw e;
      } 
      // Return the second node of the Edge.
      // If is_order, then the second node of the edge is the first node in
      // the internal representation.
      return this->is_ordered_ ? 
        this->graph_->node_(this->graph_->edges_.at(this->index_).node2_id):
        this->graph_->node_(this->graph_->edges_.at(this->index_).node1_id);
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
          if (this->graph_->has_node(this->graph_->node_(node1_id)) and
              this->graph_->has_node(this->graph_->node_(node2_id))) {
            return true;
          }
        }
      }
      // Otherwise, the edge is invalid.
      return false;
    }

    /** Return a read/write reference to the edge's value
     *
     * throws std::exception if the edge is invalid.
     */
    edge_value_type& value() {
      // Invalid edge -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }

      // Return a reference to the costom value of the corresponding internal
      // edge.
      return const_cast<Graph*>(this->graph_)->edges_.at(this->index_).value;
    }

    /** Return a read-only reference to the edge's value
     *
     * throws std::exception if the edge is invalid.
     */
    const edge_value_type& value() const {
      // Invalid edge -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return a constant reference to the costom value of the corresponding
      // internal edge.
      return this->graph_->edges_.at(this->index_).value;
    }


    /** Return the Eucleadian length of the edge
     *
     * if the edge is invalid, throws std::exception
     */
    double length() const {
      // Invalid edge -> exception
      if (not(this->is_valid())) {
        std::exception e;
        throw e;
      }
      // Return the length of the edge.
      return norm(node1().position() - node2().position());
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
    Edge(const graph_type* graph_, size_type index_, bool is_ordered_ = true) \
            : graph_(graph_), index_(index_), is_ordered_(is_ordered_) {}

    // A pointer to the graph containing this edge, nullptr if invalid.
    const graph_type* graph_;
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
    const graph_type* get_graph_() const {
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
  edge_type edge(size_type i) const {
    // HW0: YOUR CODE HERE
    // Create an edge representation of the ith internal edge.
    size_type index = this->edges_uids_to_ids_.at(i);
    return this->edge_(index);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const node_type& a, const node_type& b) const {
    // HW0: YOUR CODE HERE
    // get_edge_index returns Edge::INVALID_ID if the edge is not in the graph.
    return this->get_edge_index_(a,b) != edge_type::INVALID_ID;
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
  edge_type add_edge(const node_type& a, const node_type& b,
                     const edge_value_type value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    size_type index = this->get_edge_index_(a,b);
    bool is_ordered;
    if (index == edge_type::INVALID_ID) {
      
      size_type user_id = num_edges();
      internal_edge new_edge(a.index(), b.index(), user_id, value);
      index = (this->next_edge_index_)++;

      assert(edges_uids_to_ids_.size() == user_id);
      
      this->edges_.emplace(index, new_edge);
      this->edges_uids_to_ids_.push_back(index);
      this->add_edge_index_(a,b,index);
      is_ordered = true;
    } else {
      is_ordered = this->node_(this->edges_.at(index).node1_id) == a;
      this->edges_.at(index).value = value;
    }
    return this->edge_(index, is_ordered);
  }

  /** Remove an edge from the graph
   *
   *  Remove the edge (@a n1, @a n2) from the graph if it exists. Returns 1 if
   *  the removal succeeds, and 0 otherwise.
   *
   *  For all previously created edges and iterators, The edge objects
   *  corresponding to (@a n1, @a n2) are invalidated i.e. They will return
   *  false when calling is_valid() on them. For all other edge objects, Nothing
   *  changes. (They have the same ordering, same output for all the public
   *  methods). The iterators that point to (@a n1, @a n2) will be invalid, and
   *  will cause an undefined behavior. For all other iterators, They will have
   *  the same result as if the edge was not part of the graph in the first
   *  place.
   *
   *  @post old num_edges() == new num_edges() + result
   *  @post has_edge(@a n1,@a n2) == false
   *  @post has_edge(@a n2,@a n1) == false
   *  @post if e == old add_edge(@a n1,@a n2) then not e.is_valid()
   *  @post if ai == old edge(i), aj == old edge(j), 
   *           e == old add_edge(@a n1,@a n2), ai != e, and aj != e
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  Note not edge.is_valid() is equivalent to:
   *             (1) e.node1() throws exception
   *             (2) e.node2() throws exception
   *             (3) e.value() throws exception
   *             (4) e.length() throws exception
   *             (5) e == Edge()
   *  Complexity: O( 1 )
   */
  size_type remove_edge(const node_type& n1, const node_type& n2) {
    // Only remove nodes in the graph
    if (this->has_edge(n1,n2)) {
      
      // Get the edge's user index and internal index
      size_type index = this->get_edge_index_(n1,n2);
      size_type user_id = this->edges_.at(index).user_id;
      // Get the user index and internal index of the node with the largest
      // user index
      size_type last_user_id = this->edges_uids_to_ids_.size()-1;
      size_type last_index = this->edges_uids_to_ids_[last_user_id];
      // Swap
      this->edges_uids_to_ids_[user_id] = last_index;
      this->edges_.at(last_index).user_id = user_id;
      // Pop
      this->edges_uids_to_ids_.pop_back();
      this->edges_.erase(index);
      this->erase_edge_index_(n1,n2);
      // Removed the edge => Return 1
      return 1;
    } else {
      // The edge is not part of the graph
      return 0;
    }
  } 

  /** Remove an edge from the graph
   *
   *  Remove the edge (@a e) from the graph if it exists. Returns 1 if
   *  the removal succeeds, and 0 otherwise.
   *
   *  For all previously created edges and iterators, The edge objects
   *  corresponding to (@a e) are invalidated i.e. They will return
   *  false when calling is_valid() on them. For all other edge objects, Nothing
   *  changes. (They have the same ordering, same output for all the public
   *  methods). The iterators that point to (@a e) will be invalid, and
   *  will cause an undefined behavior. For all other iterators, They will have
   *  the same result as if the edge was not part of the graph in the first
   *  place.
   *
   *  @post old num_edges() == new num_edges() + result
   *  @post has_edge((@a e).node1(), (@a e).node2() == false
   *  @post has_edge((@a e).node2(), (@a e).node1() == false
   *  @post if h == old add_edge((@a e).node1(), (@a e).node2() 
   *        then not h.is_valid()
   *
   *  @post if ai == old edge(i), aj == old edge(j), 
   *           h == old add_edge((@a e).node1(), (@a e).node2(), ai != h
   *           and aj != h
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  Note not edge.is_valid() is equivalent to:
   *             (1) e.node1() throws exception
   *             (2) e.node2() throws exception
   *             (3) e.value() throws exception
   *             (4) e.length() throws exception
   *             (5) e == Edge()
  *
   *  Complexity: O( 1 )  
   */
  size_type remove_edge(const edge_type& e) {
    return this->remove_edge(e.node1(), e.node2());
  }

  /** Remove an edge from the graph
   *
   *  Remove the edge *(@a e_it) from the graph if it exists. Returns an
   *  iterator for the next edge if the removal succeeds, and the same iterator
   *  otherwise.
   *
   *  For all previously created edges and iterators, The edge objects
   *  corresponding to *(@a e_it) are invalidated i.e. They will return
   *  false when calling is_valid() on them. For all other edge objects, Nothing
   *  changes. (They have the same ordering, same output for all the public
   *  methods). The iterators that point to *(@a e_it) will be invalid, and
   *  will cause an undefined behavior. For all other iterators, They will have
   *  the same result as if the edge was not part of the graph in the first
   *  place.
   *
   *  @post old num_edges() == new num_edges() + result
   *  @post has_edge(*(@a e_it).node1(), *(@a e_it).node2() == false
   *  @post has_edge(*(@a e_it).node2(), *(@a e_it).node1() == false
   *  @post if e == old add_edge(*(@a e_it).node1(), *(@a e_it).node2() 
   *        then not e.is_valid()
   *
   *  @post if ai == old edge(i), aj == old edge(j), 
   *           h == old add_edge(@a e).node1(), (@a e).node2(), ai != h
   *           and aj != h
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  @post if has_node(*(@a e_it)), then if we define set1 and set2 as:
   *        let it = result
   *        set1 = empty set : while (it != node_end()) set1.add(*it_);++it;
   *        let it = old @a e_it
   *        set2 = empty set : while (it != node_end()) set2.add(*it);++it;
   *        then: set1 == set2
   *
   *  Note not edge.is_valid() is equivalent to:
   *             (1) e.node1() throws exception
   *             (2) e.node2() throws exception
   *             (3) e.value() throws exception
   *             (4) e.length() throws exception
   *             (5) e == Edge()
  *
   *  Complexity: O( 1 )
   */
  edge_iterator remove_edge(edge_iterator e_it) {
    // Copy the iterator
    auto it = e_it;
    // Increment the copy
    ++it;
    if (this->remove_edge(*e_it)) {
      return it;
    }
    return e_it;
  }

   /** Remove an edge from the graph
   *
   *  Remove the edge *(@a e_it) from the graph if it exists. Returns an
   *  iterator for the next edge if the removal succeeds, and the same iterator
   *  otherwise.
   *
   *  For all previously created edges and iterators, The edge objects
   *  corresponding to *(@a e_it) are invalidated i.e. They will return
   *  false when calling is_valid() on them. For all other edge objects, Nothing
   *  changes. (They have the same ordering, same output for all the public
   *  methods). The iterators that point to *(@a e_it) will be invalid, and
   *  will cause an undefined behavior. For all other iterators, They will have
   *  the same result as if the edge was not part of the graph in the first
   *  place.
   *
   *  @post old num_edges() == new num_edges() + result
   *  @post has_edge(*(@a e_it).node1(), *(@a e_it).node2() == false
   *  @post has_edge(*(@a e_it).node2(), *(@a e_it).node1() == false
   *  @post if e == old add_edge(*(@a e_it).node1(), *(@a e_it).node2() 
   *        then not e.is_valid()
   *
   *  @post if ai == old edge(i), aj == old edge(j), 
   *           h == old add_edge(@a e).node1(), (@a e).node2(), ai != h
   *           and aj != h
   *        then
   *             (1) if old ai < old aj => new ai < new aj
   *             (2) if old ai == old aj => new ai == new aj
   *
   *  @post if has_node(*(@a e_it)), then if we define set1 and set2 as:
   *        let it = result
   *        set1 = empty set : while (it != node_end()) set1.add(*it_);++it;
   *        let it = old @a e_it
   *        set2 = empty set : while (it != node_end()) set2.add(*it);++it;
   *        then: set1 == set2
   *
   *  Note not edge.is_valid() is equivalent to:
   *             (1) e.node1() throws exception
   *             (2) e.node2() throws exception
   *             (3) e.value() throws exception
   *             (4) e.length() throws exception
   *             (5) e == Edge()
  *
   *  Complexity: O( 1 )
   */
  incident_iterator remove_edge(incident_iterator e_it) {
    // Copy the iterator
    auto it = e_it;
    // Increment the copy
    ++it;
    if (this->remove_edge(*e_it)) {
      return it;
    }
    return e_it;
  }

 /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->nodes_.clear();
    this->nodes_uids_to_ids_.clear();
    this->edges_.clear();
    this->edges_uids_to_ids_.clear();
    this->clear_edges_ids_();
  }

  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = node_type;                // Element type
    using pointer           = node_type*;               // Pointers to elements
    using reference         = node_type&;               // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    // Internal iterator type. An iterator of a map that maps an index to an
    // internal node.
    typedef typename std::unordered_map<size_type, internal_node> \
                        ::const_iterator nodes_map_iterator;

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr) {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const

    /** Return true if this iterator is valid, false otherwise. */
    bool is_valid() const {
      return this->graph_ != nullptr;
    }

    // Forward iterator requirements
    /** Dereference the iterator and return the corresponding node.
     *
     * If the iterator is invalid this throws an exception.
     *
     * @pre this != this->graph_.edge_end()
     *
     */
    node_type operator*() const {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      // Construct a node representation in graph_ with the index of the
      // iterator. 
      return this->graph_->node_(this->it_->first);
    }


    /** increment the iterator to point to the next Node in the graph.
     *
     *  If the node is invalid, it throws an exception.
     *
     *  @pre this != this->graph_.edge_end()
     *
     */
    NodeIterator& operator++() {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      ++(this->it_);
      return *this;
    }

    /** Compare two iterators 
     *
     *  Two iterators are equal iff one of the following is true:
     *    (1) Both iterators are invalid.
     *    (2) They belong to the same graph object (in the same address) and
     *    they point to the same node.
     *
     */
    bool operator==(const NodeIterator& it) const {
      if (not (this->is_valid()) or not(it.is_valid())) {
        return not (this->is_valid()) and not (it.is_valid());
      }
      return this->it_ == it.it_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    /** Constructor 
     * 
     * @param begin_ whether the constructed iterator be a begin or an end
     * iterator.
     */
    NodeIterator(const graph_type* graph_, bool begin = true) : graph_(graph_) {
      this->it_ = begin ? graph_->nodes_.begin(): graph_->nodes_.end();
    }

    /* Pointer to the graph whose nodes are iterated via the iterator. */
    const graph_type* graph_;

    /* Iterator for the internal nodes map in the graph. */
    nodes_map_iterator it_;
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /** Returns a read/write iterator that points to the first node in the 
   *  graph.
   */
  NodeIterator node_begin() const {
    return NodeIterator(this, true);
  }

  /** Returns a read/write iterator that points one past the last node in the 
   *  graph.
   */
  NodeIterator node_end() const { 
    return NodeIterator(this, false);
  } 


  //
  // Incident Iterator
  //

  /** @class graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = edge_type;                // Element type
    using pointer           = edge_type*;               // Pointers to elements
    using reference         = edge_type&;               // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    
    typedef typename 
      std::unordered_map<node_type, size_type, typename node_type::hasher_type>
                        ::const_iterator node_edges_map_iterator;

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    /** Return true if this iterator is valid, false otherwise. */
    bool is_valid() const {
      return this->graph_ != nullptr and 
        this->node_.index() != node_type::INVALID_ID;
    }

    // Forward iterator requirements
    /** Dereference the iterator and return the corresponding edge.
     *
     * @pre this != this->node_.edge_end()
     *
     * @post result.node1() == this->node_
     *
     * If the iterator is invalid this throws an exception.
     */
    edge_type operator*() const {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      size_type edge_id = this->it_->second;
      size_type node_id = this->graph_->edges_.at(edge_id).node2_id;
      node_type n = this->graph_->node_(node_id);
      bool is_ordered = n != this->node_;
      return edge_type(this->graph_, edge_id, is_ordered);
    }

    /** increment the iterator to point to the next edge in the graph.
     *
     *  If the node is invalid, it throws an exception.
     *
     * @pre this != this->node_->edge_end()
     */
    IncidentIterator& operator++() {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      ++(this->it_);
      return *this;
    }

    /** Compare two iterators for equality 
     *
     *  Two iterators are equal iff one of the following is true:
     *    (1) Both iterators are invalid.
     *    (2) They belong to the same graph object (in the same address) and
     *    they point to the same edge.
     *
     */
    bool operator==(const IncidentIterator& it) const {
      // If one of them is invalid
      if (not (this->is_valid()) or not(it.is_valid())) {
        // Both should be invalid to be ==
        return not (this->is_valid()) and not (it.is_valid());
      }
      return (this->graph_ == it.graph_) and (this->node_ == this->node_) and
        (this->it_ == it.it_);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE

    /** Constructor 
     * 
     * @param begin_ whether the constructed iterator be a begin or an end
     * iterator.
     */
    IncidentIterator(const graph_type* graph_, const node_type node_,
                     bool begin_ = true): graph_(graph_), node_(node_) {
      this->it_ = begin_ ? this->graph_->get_node_edges_ids_(node_).begin() :
                           this->graph_->get_node_edges_ids_(node_).end();
    }

    /* Pointer to the graph whose nodes are iterated via the iterator. */
    const graph_type* graph_;
    
    /* The node whose incident edges are to be iterated by the iterator */
    node_type node_;

    /* Iterator for the internal edges map of the node node_ in the graph. */
    node_edges_map_iterator it_;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>  {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = edge_type;                // Element type
    using pointer           = edge_type*;               // Pointers to elements
    using reference         = edge_type&;               // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy
    typedef typename std::unordered_map<size_type, internal_edge> \
                        ::const_iterator edges_map_iterator;

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    /** Return true if this iterator is valid, false otherwise. */
    bool is_valid() const {
      return this->graph_ != nullptr;
    }

    // Forward iterator requirements
    /** Dereference the iterator and return the corresponding edge.
     *
     * @pre this != this->graph_.edge_end()
     *
     * If the iterator is invalid this throws an exception.
     */
    edge_type operator*() const {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      return edge_type(this->graph_, this->it_->first);
    }

    /** increment the iterator to point to the next edge in the graph.
     *
     *  If the node is invalid, it throws an exception.
     *
     * @pre this != this->graph_->edge_end()
     */
    EdgeIterator& operator++() {
      if (not (this->is_valid())) {
        throw std::exception();
      }
      ++(this->it_);
      return *this;
    }

    /** Compare two iterators for equality.
     *
     *  Two iterators are equal iff one of the following is true:
     *    (1) Both iterators are invalid.
     *    (2) They belong to the same graph object (in the same address) and
     *    they point to the same edge.
     *
     */
    bool operator==(const EdgeIterator& it) const {
      if (not (this->is_valid()) or not(it.is_valid())) {
        return not (this->is_valid()) and not (it.is_valid());
      }
      return this->it_ == it.it_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    /** Constructor 
     * 
     * @param begin_ whether the constructed iterator be a begin or an end
     * iterator.
     */
    EdgeIterator(const graph_type* graph_, bool begin = true) : graph_(graph_) {
      this->it_ = begin ? graph_->edges_.begin(): graph_->edges_.end();
    }

    /* Pointer to the graph whose edges are iterated via the iterator. */
    const graph_type* graph_;

    /* Iterator for the internal edges map in the graph. */
    edges_map_iterator it_;

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  /** Returns a read/write iterator that points to the first edge in the 
   *  graph.
   */
  EdgeIterator edge_begin() const {
    return EdgeIterator(this, true);
  }

  /** Returns a read/write iterator that points one past the last edge in the 
   *  graph.
   */
  EdgeIterator edge_end() const { 
    return EdgeIterator(this, false);
  } 



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
    internal_node(Point position, size_type user_id) : position(position),
                                  user_id(user_id), value() {}

    /** Construct an internal node with value */
    internal_node(Point position, size_type user_id, node_value_type value) :
                  position(position), user_id(user_id), value(value) {}

    Point position; // The position of the node. 
    size_type user_id; // public index of the node.
    node_value_type value; // The value of the node.
  };


  /** An unordered map that holds the nodes of the graph.
   * Always: (nodes_[index]) == this->nodes(index)
   */
  std::unordered_map<size_type, internal_node> nodes_;

  // The index to be used next time a node is added to this graph.
  size_type next_node_index_; 

  /** A vector whose ith element is the index of the node with the ith user id.
   * 
   *  Always: nodes_uids_to_ids_.size() == num_nodes() // == nodes_.size()
   *          nodes_[i].uid == nodes_uids_to_ids_[i]
   *          node(i).position() == nodes_[nodes_uids_to_ids_[i]].position()
   */
  std::vector<size_type> nodes_uids_to_ids_;

  /** Internal edge struct. */
  struct internal_edge {

    /** Construct an internal edge. */
    internal_edge(size_type node1_id, size_type node2_id, size_type user_id) :
                  node1_id(node1_id), node2_id(node2_id),
                  user_id(user_id), value() { }

    /** Construct an internal node with value */
    internal_edge(size_type node1_id, size_type node2_id,
                  size_type user_id,  edge_value_type v) : 
                  node1_id(node1_id), node2_id(node2_id),
                  user_id(user_id), value(v) {}

    size_type node1_id; // Pointer to the first node of the edge.
    size_type node2_id; // Pointer to the second node of the edge.
    size_type user_id; // public index of the edge.
    edge_value_type value; // The value of the edge.
  };

  /** An unordered map that holds the edge objects of the graph.
   * Always edges_[index] == this->edge(index) 
   * and edges_[index].node1 == this->edge(index).node1()
   */
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
   * Another approach was to use CSR or CSC form of sparse matrices.
   *
   * IMPORTANT: In general, this mamber variable is not to be used directly
   * within the Graph class, instead the private methods add_edge_index_, 
   * erase_edge_index, get_edge_index_, get_node_edges_ids, and
   * clear_edges_ids_.
   */
  std::unordered_map<node_type,
    std::unordered_map<node_type, size_type, typename node_type::hasher_type>,
      typename node_type::hasher_type> edges_ids_;

  // The index to be used next time a edge is added to this graph.
  size_type next_edge_index_;

  /** A vector whose ith element is the index of the edge with the ith user id.
   * 
   *  Always: edges_uids_to_ids_.size() == num_edges() // == edges_.size()
   *          edges_[i].uid == edges_uids_to_ids_[i]
   *          edge(i).value() == edges_[edges_uids_to_ids_[i]].value()
   */ 
  std::vector<size_type> edges_uids_to_ids_;

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
  edge_type edge_(size_type i, bool is_ordered = true) const {
    // HW0: YOUR CODE HERE
    // Create an edge representation of the ith internal node.
    return edge_type(this, i, is_ordered);
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
  size_type get_edge_index_(node_type n1, node_type n2) const {
    // The fact that both n1,n2 and n2,n1 are stored is utilized.
    // Search for n1
    auto it = this->edges_ids_.find(n1);
    if (it == this->edges_ids_.end()) return edge_type::INVALID_ID;
    // Search for n2 in the edges of n1
    auto nested_edges_ids = it->second;
    auto nested_it = (nested_edges_ids.find(n2));
    if (nested_it == nested_edges_ids.end()) return edge_type::INVALID_ID;
    // Return the index
    return nested_it->second;
  }

  /** Add the index of the edge (n1,n2)/(n2,n1) to the graph.
   *  @pre has_node(n1) == true
   *  @pre has_node(n2) == true
   *  @post get_edge_index_(n1, n2) = index
   *  @post if has_edge(n1,n2) : num_edges_ = num_edges_+ 1
   */
  void add_edge_index_(node_type n1, node_type n2, size_type index) {
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
  void erase_edge_index_(node_type n1, node_type n2) {
    // if (not(this->has_edge(n1,n2))) (this->num_edges_)--;
    this->edges_ids_[n1].erase(n2);
    this->edges_ids_[n2].erase(n1);

  }

  /** Default initializes the map of adjacent edges to a node.
   *
   * if the node is already initialized, this function doesn't do any thing.
   * Otherwise, it creates a default map of {node_type: size_type}. This map
   * should hold all the edges including the node @a n
   *
   * @post if @a n in edges_ids_ keys this function doesn't modify the internal
   *       state of the graph or the node.
   *       Otherwise, edges_ids_.at(n).size() == 0
   *
   */
  void initialize_node_edges_ids_(node_type n) {
    auto it = this->edges_ids_.find(n);
    if (it == this->edges_ids_.end()) {
      (void) this->edges_ids_[n];
    }
  }

  /** Return the edges connected to a node.
   *
   * Returns a read-only map that maps a node to an index such that each pair
   * (m,idx) in this map means that the edge (n,m) has the index idx in the
   * graph.
   *
   * @pre edges_ids_ contains n (hint: if n was added through add_node)
   * 
   * @post result.size() == n.degree()
   */
  const std::unordered_map<node_type, size_type, 
                           typename node_type::hasher_type>&
                           get_node_edges_ids_(node_type n) const {
    return this->edges_ids_.at(n);
  }

  /** Clear the indices of the graph edges 
   *
   *  @post has_edge(n1,n2) == false for every n1 and n2
   *  @post num_edges_ == 0
   */
  void clear_edges_ids_() {
    this->edges_ids_.clear();
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

template <typename V, typename E>
constexpr typename Graph<V,E>::size_type Graph<V,E>::Node::INVALID_ID = \
                                                    Graph<V,E>::size_type(-1);
template <typename V, typename E>
constexpr typename Graph<V,E>::size_type Graph<V,E>::Edge::INVALID_ID = \
                                                    Graph<V,E>::size_type(-1);
#endif

#endif // CME212_GRAPH_HPP
