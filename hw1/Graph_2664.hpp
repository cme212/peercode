#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <set>
#include <cassert>
#include <unordered_set>

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
template <typename V>
class Graph {
 public:
    /** Predeclaration of Node and edge type. */
    class Node;
    class Edge;
    /** Type of indexes and sizes.
        Return type of Graph::Node::index(), Graph::num_nodes(),
        Graph::num_edges(), and argument type of Graph::node(size_type) */
    using size_type = unsigned;
    using index_type = int;
    using node_value_type = V;
    using edge_custom_type = std::array<size_type, 2>;
    using edge_node1_custom_type = std::array<index_type, 2>;
    using edge_set_type = std::set<edge_custom_type>;
    using edge_node1_set_type = std::set<edge_node1_custom_type>;

private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

  // TX: Add sets of nodes and edges
  std::vector <Point> point_collection;
  std::vector<std::array<size_type, 2> > edge_collection;
  edge_set_type edge_collection_set;
  edge_node1_set_type edge_node1_set;
  std::vector<node_value_type> node_value_collection;
  std::vector<Node*>  node_ptr_collection;

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
    this->edge_collection = std::vector<std::array<size_type, 2> >();
    this->node_value_collection = std::vector<node_value_type>(); //
    this->node_ptr_collection = std::vector<Node*>(); // Stores pointer to node
    this->point_collection = std::vector<Point>(); //
    this->edge_collection_set = edge_set_type(); //
    this->edge_node1_set = edge_node1_set_type();


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
  class Node: private totally_ordered<Node> {
      // Predeclaration
  private:
   // Allow Graph to access Node's private member data and functions.
   friend class Graph;
   Graph* graph_ptr;
   size_type node_idx;

    Node(const Graph* ptr_tmp, size_type idx_tmp) {
      // HW0: YOUR CODE HERE

      //TX: Done
      this->graph_ptr = const_cast<Graph*>(ptr_tmp);
      this->node_idx = idx_tmp;

    }

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




    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      if (this->graph_ptr != nullptr){
          return (*(this->graph_ptr)).point_collection[this->node_idx];
      }
      else{
          return Point();
      }
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
    /* value()
    * input: None
    * output: value of node
    */
    node_value_type& value(){
        return ((*(this->graph_ptr)).node_value_collection)[this->node_idx];
    };

    const node_value_type& value() const{
        return ((*(this->graph_ptr)).node_value_collection)[this->node_idx];
    };

    // Function used for debugging
    size_type get_idx(){
        return this->node_idx;
    }

    size_type degree() const {
        // HW1: YOUR CODE HERE
        size_type deg_count = this->graph_ptr->edge_node1_set.count({index_type(this->node_idx), -1});
        return deg_count;
    };

    incident_iterator edge_begin() const {
        // HW1: Your code here
        return IncidentIterator(this->graph_ptr, this->node_idx, 0);
    };

    incident_iterator edge_end() const {
        // HW1: YOUR CODE HERE
        return IncidentIterator(this->graph_ptr, this->node_idx, -1);
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      return this->node_idx == n.node_idx && this->graph_ptr == n.graph_ptr;
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
        assert_GRAPH_PTR_EXIST_TX;
        assert_GRAPH_PTR_EXIST_TX_general(n);
        return this->node_idx < n.node_idx;
    }


  };// End of node class

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE

    // TX: add a very easy function

    return this->point_collection.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return this->size();

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
    // HW1 update: added null pointer for this method
    Node node_tmp(this, this->size());
    node_value_type tmp_value {};
    (this->point_collection).push_back(position); // Add the point
    (this->node_value_collection).push_back(tmp_value); // Add the point
    (this->node_ptr_collection).push_back(&node_tmp); // Add node pointer

    // TODO: is the above line legal? Do I encounter illegal change to the Point instance?
    return node_tmp;
  }

  Node add_node(const Point& position, const node_value_type& tmp_value) {
    // HW0: YOUR CODE HERE
    // HW1 update: added null pointer for this method
    Node node_tmp(this, this->size());
    (this->point_collection).push_back(position); // Add the point
    (this->node_value_collection).push_back(tmp_value); // Add the point
    (this->node_ptr_collection).push_back(&node_tmp); // Add node pointer

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

  void set_value(size_type i, node_value_type v) {
    // HW0: YOUR CODE HERE
    assert (i < this->size());
      (this->node_value_collection).at(i) = v;
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
      Edge(const Graph* ptr_tmp, size_type idx_tmp) {
          // HW0: YOUR CODE HERE

          //TX: Done
          this->graph_ptr = const_cast<Graph*>(ptr_tmp);
          this->node_idx_1 = (*ptr_tmp).edge_collection[idx_tmp][0];
          this->node_idx_2 = (*ptr_tmp).edge_collection[idx_tmp][1];
      }



   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
      this->graph_ptr = nullptr;
    }

    Edge swap(){
        Edge new_edge;
        new_edge.graph_ptr = const_cast<Graph*>(this->graph_ptr);
        new_edge.node_idx_1 = this->node_idx_2;
        new_edge.node_idx_2 = this->node_idx_1;
        return new_edge;
    }


    bool is_valid(){
        return this->graph_ptr != nullptr;
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

    /* HAS_NODE
    Deterine if the Node is in the edge
    */
    bool has_node(const Node& n) const {
      // HW0: YOUR CODE HERE
        bool cond1 = n == this->node1();
        bool cond2 = n == this->node2();

      // TX: Determine if node is a valid index
      return cond1 || cond2;
    }


    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //TX: added my code here
        bool cond1 = (this->node_idx_1 == e.node_idx_1) && (this->node_idx_2 == e.node_idx_2);
        bool cond2 = (this->node_idx_1 == e.node_idx_2) && (this->node_idx_2 == e.node_idx_1);
        bool cond3 = this->graph_ptr == e.graph_ptr;
      return (cond1 || cond2) && cond3;
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

    unsigned int result = this->edge_collection_set.count(edge_custom_type({node_idx_a, node_idx_b}));
    return result != 0;
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
      (this->edge_collection).push_back({node_idx_a, node_idx_b}); // Add the point
      this->edge_collection_set.insert(edge_custom_type({node_idx_a, node_idx_b}));
      this->edge_collection_set.insert(edge_custom_type({node_idx_b, node_idx_a}));
      this->edge_node1_set.insert(edge_node1_custom_type({int(node_idx_a), -1}));
      this->edge_node1_set.insert(edge_node1_custom_type({int(node_idx_b), -1}));

      return Edge(this , this->num_edges());
  }

  std::vector<size_type> index_total(size_type  n_tmp){
      std::vector<size_type> tmp_vec;
      for (unsigned int idx = 0; idx < this->num_edges(); idx++){
          size_type node_1 = this->edge_collection[idx][0];
          size_type node_2 = this->edge_collection[idx][1];

          bool cond1 = (n_tmp == node_1) || (n_tmp == node_2);
          if (cond1) {
              tmp_vec.push_back(idx);
          }
      }
      return tmp_vec;
  }

  std::vector<bool> swap_total(size_type  n_tmp){
      std::vector<bool> tmp_vec;
      for (unsigned int idx = 0; idx < this->num_edges(); idx++){
          size_type node_1 = this->edge_collection[idx][0];
          size_type node_2 = this->edge_collection[idx][1];

          bool cond1 = (n_tmp == node_1) || (n_tmp == node_2);
          if (cond1) {
              if (n_tmp == node_2){
                  tmp_vec.push_back(true);
              }
              else{
                  tmp_vec.push_back(false);
              }

          }
      }
      return tmp_vec;
  }


  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    this->point_collection = std::vector<Point>(); //TODO: check correctness
    this->edge_collection = std::vector<std::array<size_type, 2> >();
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
    NodeIterator() {
        this->graph_ptr = nullptr;
        this->node_idx = -1;
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    // Return the Node pointing to
    Node operator*() const {
        // HW1 #2: YOUR CODE HERE
        if (this->graph_ptr == nullptr || this->node_idx == -1){
            return Node();
        }
        else{
            return (*this->graph_ptr).node(this->node_idx);
        }
    }

    // Advance and turn pointer to null if out of bound
    node_iterator& operator++() {
        // HW1 #2: YOUR CODE HERE
        // TX: Done
        if (node_idx != -1){
            node_idx++;

            if (node_idx >= int((*this->graph_ptr).size() )){
                node_idx = -1; // Becomes null iterator
            }
        }
        return *this;
    }


    bool operator==(const node_iterator& iter_tmp) const {
        // HW1 #2: YOUR CODE HERE
        // TX: Done
        bool cond1 = this->graph_ptr == iter_tmp.graph_ptr;
        bool cond2 = this->node_idx == iter_tmp.node_idx;
        return cond1 && cond2;
    }

    index_type get_idx(){
        return this->node_idx;
    }

    Graph* get_graph(){
          return this->graph_ptr;
    }

    bool is_valid(){
        bool cond1 = this->graph_ptr != nullptr;
        bool cond2;
        if (cond1) {
            cond2 = node_idx < int((*this->graph_ptr).size());
        }
        else{
            cond2 = false;
        }
        return cond1 && cond2;
    }

   private:
    friend class Graph;
    index_type node_idx;
    Graph* graph_ptr;

    NodeIterator(const Graph* ptr_tmp, index_type idx_tmp) {
      // HW1: YOUR CODE HERE

      this->graph_ptr = const_cast<Graph*>(ptr_tmp);
      this->node_idx = idx_tmp;
    }

    // HW1 #2: YOUR CODE HERE
 }; // End of node iterator

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Initialize with 0 position
  node_iterator node_begin() const {
      return NodeIterator(this, 0);
  }

  // Initialize to null position
  node_iterator node_end() const {
      return NodeIterator(this, -1);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator> {
  private:
      index_type edge_idx;
      std::vector<bool> swp_total;
      Graph* graph_ptr;
      size_type node_idx;
      std::vector<size_type> idx_total;

   public:

    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() {
        this->graph_ptr = nullptr;
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    // Obtain edge pointed to with null exception
    Edge operator*() const {
        if (this->graph_ptr == nullptr || this->edge_idx == -1){
            return Edge();
        }
        else{
            if (swp_total[edge_idx]){
                return (this->graph_ptr->edge(idx_total[edge_idx])).swap();
            }
            else{
                return (this->graph_ptr->edge(idx_total[edge_idx]));
            }
        }
    }

    // Advance and turn iterator to null if necessary

    IncidentIterator& operator++() {
        // HW1 #2: YOUR CODE HERE
        // TX: Done
        if (edge_idx != -1){
            edge_idx++;

            if (edge_idx >= int(idx_total.size() )){
                edge_idx = -1; // Becomes null iterator
            }
        }
        return *this;
    }

    bool is_null(){
        return edge_idx == -1;
    }

    size_type get_idx() {
        return (**this).node2().get_idx();
    }


    bool operator==(const IncidentIterator& iter_tmp) const {
        // HW1 #2: YOUR CODE HERE
        // TX: Done
        bool cond1 = this->graph_ptr == iter_tmp.graph_ptr;
        bool cond2 = this->node_idx == iter_tmp.node_idx;
        bool cond3 = this->edge_idx == iter_tmp.edge_idx;

        return cond1 && cond2 && cond3;
    }

   private:
    friend class Graph;


    IncidentIterator(Graph* ptr, size_type n_tmp, index_type e_tmp) {
        this->graph_ptr = ptr;
        this->edge_idx = e_tmp;
        this->node_idx = n_tmp;
        this->idx_total = (this->graph_ptr->index_total(n_tmp));
        this->swp_total = this->graph_ptr->swap_total(n_tmp);
        if (this->idx_total.size() == 0 ){
            this->edge_idx = -1; // Defaults to null iterator
        }
    }
    // HW1 #3: YOUR CODE HERE
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator> {
  private:
      index_type edge_idx;
      Graph* graph_ptr;
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() {
        this->graph_ptr = nullptr;
    }

    // HW1 #5: YOUR CODE HERE
    Edge operator*() const {
        // HW1 #2: YOUR CODE HERE
        if (this->graph_ptr == nullptr || this->edge_idx == -1){
            return Edge();
        }
        else{
            return (*this->graph_ptr).edge(this->edge_idx);
        }
    }

      EdgeIterator& operator++() {
          // HW1 #2: YOUR CODE HERE
          // TX: Done
          if (edge_idx != -1){
              edge_idx++;

              if (edge_idx >= int((*this->graph_ptr).size() )){
                  edge_idx = -1; // Becomes null iterator
              }
          }
          return *this;
      }

      bool operator==(const EdgeIterator& iter_tmp) const {
          // HW1 #2: YOUR CODE HERE
          // TX: Done
          bool cond1 = this->graph_ptr == iter_tmp.graph_ptr;
          bool cond2 = this->edge_idx == iter_tmp.edge_idx;
          return cond1 && cond2;
      }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const Graph* ptr_tmp, index_type idx_tmp) {
      this->graph_ptr = const_cast<Graph*>(ptr_tmp);
      this->edge_idx = idx_tmp;
    }
 }; //End of edge class

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:

  // Intitialize to 0 position
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  }
  // Intitialize to null position

  edge_iterator edge_end() const {
      return EdgeIterator(this, -1);
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

}; // End of Graph class

#endif // CME212_GRAPH_HPP
