#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <set>
#include <iostream>
#include <map>
#include <unordered_map>
#include <utility>
#include <tuple>
#include <unordered_set>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> // Setting up Graph as a template class.
class Graph {
 //Declare internal structs
 struct internal_node;
 struct internal_edge;
 struct node_adjacent;

    public:
        // Value type assignments.
        typedef V node_value_type;
        typedef E edge_value_type;

        std::vector<internal_edge> int_edges_vec; //Vector of internal edge objects

        std::vector<internal_node> int_nodes_vec; //Vector of internal node objects

        std::unordered_map<unsigned, int> nid_index; // Map from node id's to their indexes
        // in the internal structures

        std::unordered_map<unsigned, int> eid_index; // Map from edges to their indexes

 private:
  // HW0: CODE ADDED HERE
  // Use this space for declarations of important internal types you need
  // Interal type for nodes:
     struct internal_node;
     struct internal_edge;
     struct node_adjacent;

     unsigned size_; //Number of nodes

     unsigned num_edges_; //Number of edges

     struct internal_node{
         unsigned _int_id; // Id of node
	     Point _point; //Point of node
         node_value_type _value; // Value of node
         std::vector<node_adjacent> intadjnodes_; //Vector of adjacent nodes
         Point _initial_point; // Used to track constraint of fixed points.
         internal_node(Point intpoint, unsigned intnidx
                , node_value_type value, Point initial_point)
               : _point(intpoint), _int_id(intnidx),
               _value(value), _initial_point(initial_point) {}

         // This returns the adjacent node index for a given index in the
         // adjecency vector.
         unsigned adj_node(int i) const {
             return intadjnodes_[i].adj_node_id;
         }
      };

     // Internal type of edges:
     struct internal_edge{
         unsigned _edge_id; //Index of edge
         unsigned _node_1_id; //Index of node1
         unsigned _node_2_id; //Index of node2
         edge_value_type _value;
         internal_edge(unsigned inteidx, unsigned intnode1idx,
                 unsigned intnode2idx, edge_value_type value)
             : _edge_id(inteidx), _node_1_id(intnode1idx),
             _node_2_id(intnode2idx), _value(value) {}
     };


     //Internal type of adjacent node
     struct node_adjacent{
         unsigned adj_edge_id; //Index of edge
         unsigned adj_node_id; //Index of node attached
         node_adjacent(unsigned adjeidx, unsigned adjnodeidx)
             : adj_edge_id(adjeidx), adj_node_id(adjnodeidx) {}
         node_adjacent() : adj_edge_id(0), adj_node_id(0) {}
     };

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<node_value_type, edge_value_type>;

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

  // Vector for storage of Edges/Nodes
  std::vector<Node> nodes_vec;
  std::vector<Edge> edges_vec;

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph()
      : size_(0), num_edges_(0) {
    // HW0: CODE ADDED HERE
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
    Node()
        : _gptr(nullptr), _id(0){
      // HW0: CODE ADDED HERE
    }
    /** Return this nodes position. */
    const Point& position() const {
      // HW0: CODE ADDED HERE
      // Reach into vector of internal nodes and return point attribute.
        return _gptr->int_nodes_vec[_gptr->nid_index.at(_id)]._point;
    }

    /** Returns this nodes position*/
    Point& position() {
        return _gptr->int_nodes_vec[_gptr->nid_index.at(_id)]._point;
    }


    Graph* graph() const {
        return _gptr;
    }

    /** Return this nodes unique id. */
    unsigned id() const {
      // HW0: CODE ADDED HERE
      // Return index attribute
        return _id;
    }

    /** Return this nodes index in the internal structures, a number in range
     * (0, graph.size) */
    unsigned index() const{
        return _gptr->nid_index.at(_id);
    }


    /** Returns this nodes internal node representation */
    internal_node& fetch() const{
        return _gptr->int_nodes_vec[_gptr->nid_index.at(_id)];
    }

    /** Returns this nodes starting position */
    const Point& initial_position() const {
        return fetch()._initial_point;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    node_value_type& value() {
        return fetch()._value;
    }

    const node_value_type& value() const {
        return fetch()._value;
    }

    size_type degree() const {
        return this->fetch().intadjnodes_.size();
    }

    incident_iterator edge_begin() const {
        // Goes to the intadjnodes and returns a vector iterator to the first
        // element.
        incident_iterator _adj_edge_it(this, 0);
        return _adj_edge_it;
    }

    incident_iterator edge_end() const {
        // Goes to the intadjnodes and returns a vector iterator to the last
        // element.
        incident_iterator _adj_edge_it(this, degree());
        return _adj_edge_it;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: CODE ADDED HERE
      // Given 2 Nodes, return True if they have same Graph and index.
      if(this-> _gptr == n._gptr && this-> _id == n._id){
          return true;
      }
      (void) n;          // Quiet compiler warning
      return false;
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
      // HW0: CODE ADDED HERE
      // Overload operator to compare node indices
        return (this->_id < n._id);
        (void) n;           // Quiet compiler warning
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: CODE ADDED HERE
    Graph* _gptr; //Pointer to graph
    unsigned _id; //Index of node
    // Private constructor to return a valid Node
    Node(const Graph* ngraph, size_type nidx)
        : _gptr(const_cast<Graph*>(ngraph)), _id(nidx) {
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  unsigned size() const {
    // HW0: CODE ADDED HERE
    // Return size_ (number of nodes) attribute
     return nodes_vec.size();
    // return 0;
  }

  /** Synonym for size(). */
  unsigned num_nodes() const {
    return size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(Point position,
          const node_value_type& _value = node_value_type()) {
    // HW0: CODE ADDED HERE

    // Add a new instance of internal_node to the back of the vector
    int_nodes_vec.push_back(internal_node(position, size_,
                _value, position));

    // Add id to nid_index
    nid_index[size_] = nodes_vec.size();

    // Create Node object
    Node _node(this, size_);
    nodes_vec.push_back(_node);
    size_ = size_ + 1; // Increment number of nodes
    (void) position;      // Quiet compiler warning
    return _node;        // Return node
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: CODE ADDED HERE
    // Check if node is present in the unordered_map of nodes.
    if(nid_index.count(n.id()) == 1) {return true;}

    (void) n;            // Quiet compiler warning
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: CODE ADDED HERE
    assert(i < size_);
    (void) i;             // Quiet compiler warning
    return nodes_vec[i];        // Return node at id i
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
    Edge()
        : egraph_(nullptr), eidx_(0), n1idx_(0), n2idx_(0) {
      // HW0: CODE ADDED HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: CODE ADDED HERE
      return Node(egraph_, n1idx_); //Return the first Node
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: CODE ADDED HERE
    return Node(egraph_, n2idx_); // Return the second Node
    }

    edge_value_type& value() {
        return egraph_->int_edges_vec[egraph_->eid_index.at(eidx_)]._value;
    }

    const edge_value_type& value() const {
        return egraph_->int_edges_vec[egraph_->eid_index.at(eidx_)]._value;
    }

    /** Helper methods to sawp the orientation of the nodes, ensures the passed
    * node is returned from node1()
    */
    void assign_first(const Node* n) {
        if ((*n).id() != n1idx_) {
            unsigned tmp = n1idx_;
            n1idx_ = n2idx_;
            n2idx_ = tmp;
        }
    }

    /** Returns the length of the node */
    double length() {
        return norm(node1().position() - node2().position());
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: CODE ADDED HERE
      // If the graph and index of the edge match, they are ==
      if (egraph_ == e.egraph_ && eidx_ == e.eidx_){
          return true;
      }
      return false;
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: CODE ADDED HERE
      // Have an arbitrary output if they do not belong to the same graph.
      if (egraph_ != e.egraph_) {return true;}
      //Overload operator to compare edge indices
      return (this-> eidx_ < e.eidx_);
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: CODE ADDED HERE
    graph_type* egraph_; // Pointer to graph
    unsigned eidx_; // Index of edge
    unsigned n1idx_; // Index of first node
    unsigned n2idx_; //Index of second node
    //Private constructor for returning valid Edge objects:
    Edge(const graph_type* graph, unsigned eidx, unsigned n1idx, unsigned n2idx)
        : egraph_(const_cast<Graph*>(graph)), eidx_(eidx), n1idx_(n1idx), n2idx_(n2idx){}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: CODE ADDED HERE
      return edges_vec.size(); //Return number of edges
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: CODE ADDED HERE
    (void) i;             // Quiet compiler warning
    unsigned node1idx = int_edges_vec[i]._node_1_id;
    unsigned node2idx = int_edges_vec[i]._node_2_id;
    return Edge(this, i, node1idx, node2idx); // Return Edge at i
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: CODE ADDED HERE
    // For each node in the corresponding adjacent nodes vector,
      // check to see if the index matches node b
    for(unsigned i=0; i< a.fetch().intadjnodes_.size();i++){
        if(a.fetch().intadjnodes_[i].adj_node_id == b.id()){
            return true;
    	}
    }
    (void) a; (void) b;   // Quiet compiler warning
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
  Edge add_edge(const Node& a, const Node& b,
          edge_value_type value = edge_value_type()) {
    // HW0: YOUR CODE HERE
    // Check if edge already exists

    // We loop again (instead of calling has_edge) to get the index
    for(unsigned i=0; i< a.fetch().intadjnodes_.size();i++){
        if(a.fetch().intadjnodes_[i].adj_node_id == b.id()){
            return edges_vec[eid_index.at(a.fetch().intadjnodes_[i].adj_edge_id)];
        }
    }

    // Add two new instances of adjacent nodes (one for a, one for b)
    // Add these instances to the corresponding vector
    a.fetch().intadjnodes_.push_back(node_adjacent(num_edges_,b.id()));
    b.fetch().intadjnodes_.push_back(node_adjacent(num_edges_,a.id()));

    // Add one new instance of internal edge and add to vector
    int_edges_vec.push_back(internal_edge(num_edges_,a.id() ,b.id(), value));

    // Update the eid_index map
    eid_index[num_edges_] = edges_vec.size();

    // Create new edge object and add to vector
    Edge _edge(this, num_edges_, a.id(), b.id());
    edges_vec.emplace_back(_edge);

    // Increment number of edges
    num_edges_ = num_edges_ + 1;

    (void) a, (void) b;   // Quiet compiler warning
    return _edge; // Return New Edge
   }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: CODE ADDED HERE
    nid_index.clear();
    eid_index.clear();
    edges_vec.clear();
    nodes_vec.clear();
    int_edges_vec.clear();
    int_nodes_vec.clear();
    size_ = 0;
    num_edges_ = 0;
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

    // Location in the vector of nodes.
    int _index;
    const graph_type* _gptr;

    Node operator*() const {
        Node node = _gptr->nodes_vec[_index];
        return node; // Returns the nodes at location
        //_index of the vector of nodes.
    }

    NodeIterator& operator++() {
        _index = _index + 1; // Increment vector iterator on nodes_vec
        return *this;
    }

    bool operator==(const NodeIterator& iter) const {
        // Check if the location is the same.
        if ((this->_gptr == iter._gptr) and (this->_index == iter._index)) {
            return true;
        }
        return false;
    }

    bool operator!=(const NodeIterator& iter) const {
        // Check if graph and location are different
        if ((this->_gptr == iter._gptr) and (this->_index == iter._index)) {
            return false;
        }
        return true;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    NodeIterator(const graph_type* _gptr, int _start)
        : _gptr(_gptr), _index(_start) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
      node_iterator node_it(this, 0); // Create an iterator object with index 0.
    return node_it;
  }

  node_iterator node_end() const {
      node_iterator node_it(this, nodes_vec.size()); // Create an iterator object
      // with index size_
    return node_it;
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
    IncidentIterator() {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    int _index;
    const Node* _node;

    Edge operator*() const {
        // Access edge id of edge of edge pointed to by it in adjacent node list
        // of required node.
        Edge e = _node->_gptr->edges_vec
            [_node->fetch().intadjnodes_[_index].adj_edge_id];
        e.assign_first((_node));
        return e;
    }

    IncidentIterator& operator++() {
        _index = _index + 1;
        return *this;
    }

    bool operator==(const IncidentIterator& iter) const {
        if (this->_node == iter._node
                and this->_index == iter._index) {return true;}
        return false;
    }

    int index() const {
        return _index;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    IncidentIterator(const Node* node, int start_idx)
        : _node(node), _index(start_idx) {}
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
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const

    int _index;
    const graph_type* _gptr;

    Edge operator*() const {
        return _gptr->edges_vec[_index];
    }

    EdgeIterator& operator++() {
        _index = _index + 1;
        return *this;
    }

    bool operator==(const EdgeIterator& iter) const {
        if (this->_gptr == iter._gptr and this->_index == iter._index) {return true;}
        return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    EdgeIterator(const graph_type* gptr, int start_idx)
        : _gptr(gptr), _index(start_idx) {}
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin() const {
      edge_iterator edge_it(this, 0);
      return edge_it;
  }

  edge_iterator edge_end() const {
      edge_iterator edge_it(this, edges_vec.size());
      return edge_it;
  }

  /** Removes the specified node in the graph, if it exitst */
  size_type remove_node(const Node& n) {

    // Checks the node is in the graph.
    if(has_node(n) == true) {

        // Stores the information from the current node and the node being
        // swapped out.
        int node_id = n._id;
        int node_index = nid_index.at(node_id);
        int last_node_id = nodes_vec.back().id();

        // Iterates through each of the incident edges and removes them by
        // calling remove_edge()
        int incidents = n.degree();
        for (auto i = 0; i < incidents; i++) {
            int edge_id = n.fetch().intadjnodes_.back().adj_edge_id;
            remove_edge(edges_vec[eid_index.at(edge_id)]);
        }

        // Copies the last element into the location of the one I wish to
        // remove.
        nodes_vec[node_index] = nodes_vec.back();
        // Removes the now duplicate final element.
        nodes_vec.pop_back();

        // Does the same steps for the internal representations
        int_nodes_vec[node_index] = int_nodes_vec.back();
        int_nodes_vec.pop_back();

        // Organises the mapping from id_index.
        nid_index.erase(node_id);
        nid_index[last_node_id] = node_index;

        return true;
    } else {
        return false;
    }
  }

  node_iterator remove_node(node_iterator n_it) {
      // Simply calls the node removel method based upon a node.
    remove_node(*n_it);
    return n_it;
  }

  size_type remove_edge(const Node& n1, const Node& n2) {
      // Checks the edge is in the graph.
    if (has_edge(n1, n2) == true) {

        // Will the information when we locate the edge id in the adjacency
        // list of n1.
        int edge_id;
        int edge_index;
        node_adjacent k;

        // Find the edge in the adjacency list of n1, locates and then erases
        // the adjacency object.
        for (auto i = n1.edge_begin(); i!=n1.edge_end(); ++i) {
            k = i._node->fetch().intadjnodes_[i._index];
            if (k.adj_node_id == n2.id()) {
                edge_id = k.adj_edge_id;
                edge_index = eid_index.at(edge_id);
                n1.fetch().intadjnodes_.erase(
                        n1.fetch().intadjnodes_.begin()+i._index);
                break;
            }
        }

        // Remove from adjacency list of n2, locates and then erases the
        // adjacency object.
        for (auto i = n2.edge_begin(); i!=n2.edge_end(); ++i) {
            k = i._node->fetch().intadjnodes_[i.index()];
            if (k.adj_node_id == n1.id()) {
                n2.fetch().intadjnodes_.erase(
                        n2.fetch().intadjnodes_.begin()+i.index());
                break;
            }
        }

        // Copy the final element and then remove the now duplicate final
        // element.
        edges_vec[edge_index] = edges_vec.back();
        edges_vec.pop_back();

        // Do the same for the internal representations.
        int_edges_vec[edge_index] = int_edges_vec.back();
        int_edges_vec.pop_back();

        // Reorganise the mapping eid_index.
        eid_index.erase(edge_id);
        eid_index[edges_vec[edge_index].eidx_] = edge_index;

        return true;

      } else {
          return false;
      }
  }

  size_type remove_edge(const Edge& e) {
      // Call the removeal method based on two nodes.
      return remove_edge(e.node1(), e.node2());
  }

  edge_iterator remove_edge(edge_iterator e_it) {
      // Call the removal method based on an edge.
      remove_edge(*e_it);
      return e_it;
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.


};

#endif // CME212_GRAPH_HPP
