#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <map>
#include <vector>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


// ############ used code from https://raw.githubusercontent.com/cme212/peercode/master/hw0/Graph_2489.hpp


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  struct internal_node;
  struct internal_edge;

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  /** Type of this graph. */
  using graph_type = Graph<V>;

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
  Graph()
    : nodes_(), node_counts_(0), next_nid_(0),
      edges_(), edge_counts_(0), next_eid_(0), edge_lookup(), map_node_to_node() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
    for (size_type i = 0; i < node_counts_; i++) {
      internal_node* node_ptr = nodes_[i];
      delete node_ptr;
    }
    for (size_type i = 0; i < edge_counts_; i++) {
      internal_edge* edge_ptr = edges_[i];
      delete edge_ptr;
    }
    clear();
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
    

    // return the number of incident edges .
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().coord;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return fetch().nid;
    }
    


    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    //returns the node value, from the graph structure
    node_value_type& value() {
        return graph_->nodes_[this->index()]->value;
    }
    //const returns the node value, from the graph structure
    const node_value_type&  value() const {
        return value_;
    }

    /** Return this node's degree, defined as its number of neighbors. */
    size_type degree () const {
        if (graph_->map_node_to_node.find(*this) != graph_->map_node_to_node.end()) {
            // Node is present in the map, return the size of the adjacency list at that node
              return graph_->map_node_to_node.at(*this).size();
          }
        else {
            //Node absent from the map, can't have incident edges, return 0
            return 0;
          }
    }
    //begin - end iterators for incident edges
    incident_iterator edge_begin() const {
        return IncidentIterator(graph_, *this, graph_->map_node_to_node[*this], 0);
    }
    incident_iterator edge_end() const {
          return IncidentIterator(graph_, *this, graph_->map_node_to_node[*this], graph_->map_node_to_node[*this].size());
      }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == n.graph_) && (index() == n.index()));
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
      return ((index() < n.index()) && (graph_ == n.graph_)) ||  (graph_ < n.graph_);
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    Graph* graph_;
    // This node's id
    size_type nid_;
    // This node's value
    node_value_type value_;
    /** Private Constructor that can be accessed by Graph **/
    Node(const Graph* graph, size_type nid, node_value_type value)
      : graph_(const_cast<Graph*>(graph)), nid_(nid), value_(value) {}
    /** Helper method to return the data associated with node that is
        stored in Graph **/
    internal_node& fetch() const {
      if (graph_->has_node(*this)){
        return *(graph_->nodes_[nid_]);
      }else{
        throw("node is invalid");
      }
    }


  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return node_counts_;
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
  Node add_node(const Point& position, const node_value_type& value = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(new internal_node{next_nid_, position, value});
    next_nid_ += 1;
    node_counts_ += 1;
    return Node(this, next_nid_-1, value);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    return ((n.nid_ < node_counts_) && (n.graph_ == this));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    assert(i < size());
    node_value_type val = nodes_[i]->value;
    return Node(this, i, val);
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
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      if (!nodes_swapped_){
        size_type n1id = fetch().n1;
        node_value_type n1val = graph_->node(n1id).value();
        return Node(graph_, n1id, n1val);
        }
      else{
        size_type n1id = fetch().n2;
        node_value_type n1val = graph_->node(n1id).value();
        return Node(graph_, n1id, n1val);
      }
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      if (nodes_swapped_){
        size_type n1id = fetch().n1;
        node_value_type n1val = graph_->node(n1id).value();
        return Node(graph_, n1id, n1val);
        }
      else{
        size_type n1id = fetch().n2;
        node_value_type n1val = graph_->node(n1id).value();
        return Node(graph_, n1id, n1val);
      }
    }
    //swap the associated nodes
    void swap_nodes() {
      nodes_swapped_ = true;
    }
    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      bool cond1 = (node1() == e.node1()  && node2() == e.node2());
      bool cond2 = (node1() == e.node2() && node2() == e.node1());
      return (cond1 || cond2);
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      return ((eid_ < e.eid_) && (graph_ == e.graph_)) || (graph_ < e.graph_);
    }


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    // Pointer back to the Graph container
    Graph* graph_;
    // This edge's id
    size_type eid_;
    /** Private Constructor that can be accessed by Graph **/
    bool nodes_swapped_ = false;
    Edge(const Graph* graph, size_type eid)
      : graph_(const_cast<Graph*>(graph)), eid_(eid) {}
    /** Helper method to return the data associated with edge that is
        stored in Graph **/
    internal_edge& fetch() const {
      if (eid_ >= graph_->edge_counts_) throw("edge id is invalid");
      return *(graph_->edges_[eid_]);
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return edge_counts_;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
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
    if (has_node(a) && has_node(b)) {
      std::vector<size_type> edge_key{a.nid_, b.nid_};
    if (edge_lookup.find(edge_key) != edge_lookup.end()) {
        return true;
      }
    }
    return false;
  }



  Edge swap_nodes(const Edge& e) const {
    size_type nid1 = edges_[e.eid_]->n1;
    size_type nid2 = edges_[e.eid_]->n2;
    edges_[e.eid_] = {e.eid_, nid2, nid1};
    return Edge(this, e.eid_);

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
    std::vector<size_type> edge_key{a.nid_, b.nid_};
    if (has_edge(a, b)) {
      return Edge(this, edge_lookup[edge_key]);

    }

    internal_edge* new_edge = new internal_edge{next_eid_, a.nid_, b.nid_};
    edges_.push_back(new_edge);
    edge_lookup.insert({edge_key, next_eid_});
    std::vector<size_type> edge_key_rev{b.nid_, a.nid_};
    edge_lookup.insert({edge_key_rev, next_eid_});
    next_eid_ += 1;
    edge_counts_ += 1;
    // update adjacency list
    map_node_to_node[a].push_back(b);
    map_node_to_node[b].push_back(a);
    return Edge(this, next_eid_-1);
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_.clear();
    edges_.clear();
    edge_lookup.clear();
    node_counts_ = 0;
    edge_counts_ = 0;
    next_nid_ = 0;
    next_eid_ = 0;
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
    
    // * operator for NodeIterator, returns the node pointed at
    Node operator*() const {
        return graph_->node(nid_);
    }
    // * operator for NodeIterator, point to the next node, ordered by nid_
    NodeIterator& operator++() {
        ++nid_;
        return *this;
    }
    // == operator for NodeIterator, equality qith @a n means having the same graph pointer and the same current nid
    bool operator==(const NodeIterator& otherNodeIterator) const {
        return ( graph_ == otherNodeIterator.graph_) && (nid_ == otherNodeIterator.nid_);
    }
    bool operator!=(const NodeIterator& otherNodeIterator) const {
        return !(*this == otherNodeIterator);
    }
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph * graph_;
    size_type nid_;
    NodeIterator (const Graph* graph, size_type nid) : graph_(const_cast<Graph*>(graph)), nid_(nid) {} 
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const
    node_iterator node_begin() const {
        return NodeIterator(this, 0);

    }
    node_iterator node_end() const {
        return NodeIterator(this, this->num_nodes());
    }
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

    //dereferencing operator for IncidentIterator
    Edge operator*() const {
        Edge corresponding_edge = graph_->add_edge(query_node_, neighbor_nodes_[neighbor_number_]);
        if (corresponding_edge.node1().nid_ != query_node_.nid_){
          corresponding_edge.swap_nodes();
          return corresponding_edge;
        }
        return corresponding_edge;
        
    }
    //prefix ++ operator for IncidentIterator
    IncidentIterator& operator++() {
        ++neighbor_number_;
        return *this;
    }
    //== operator for IncidentIterator, this and @a otherIncidentIt must have the same query node and current index
    bool operator==(const IncidentIterator& otherIncidentIt) const {
        return (query_node_ == otherIncidentIt.query_node_) && (neighbor_number_ == otherIncidentIt.neighbor_number_);
    }
    bool operator!=(const IncidentIterator& otherIncidentIt) const {
          return !(*this == otherIncidentIt);
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph * graph_;
    // The node whose incident Edges we iterate on
    node_type query_node_;
    // The query node's neighbors involved in incident edges with the query node
    std::vector<node_type> neighbor_nodes_;
    // The index of the neighbor we are currently visiting (different from nid_)
    size_type neighbor_number_;
    // Constructor
    IncidentIterator(const Graph* graph, node_type query_node, std::vector<node_type> neighbor_nodes, size_type neighbor_number): 
            graph_(const_cast<Graph*>(graph)), query_node_(query_node), neighbor_nodes_(neighbor_nodes), neighbor_number_(neighbor_number) {}
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

    //++ prefix operator for EdgeIterator
    EdgeIterator& operator++() {
        ++eid_;
        return *this;
    }
    //* prefix operator for EdgeIterator
    Edge operator*() const {
        return graph_->edge(eid_);

    }

    //Equality operator for EdgeIterator, this and @a otherEdgeIt must have the same graph and same current edge id
    bool operator==(const EdgeIterator& otherEdgeIt) const {
        return (graph_ == otherEdgeIt.graph_) && (eid_ == otherEdgeIt.eid_);
    }
    bool operator!=(const EdgeIterator& otherEdgeIt) const {
        return !(*this == otherEdgeIt);
    }
   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph * graph_;
    size_type eid_;
    EdgeIterator (const Graph* graph, size_type eid) : graph_(const_cast<Graph*>(graph)), eid_(eid) {} 
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  // returns the beginning of edge iterator
  edge_iterator edge_begin() const {
      return EdgeIterator(this, 0);
  }
  // returns the end of edge iterator
  edge_iterator edge_end() const {
      return EdgeIterator(this, this->num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    size_type nid;
    Point coord;
    node_value_type value;
  };

  struct internal_edge{
    size_type eid;
    size_type n1;
    size_type n2;
  };

  std::vector<internal_node*> nodes_;
  size_type node_counts_;
  size_type next_nid_;

  std::vector<internal_edge*> edges_;
  size_type edge_counts_;
  size_type next_eid_;

  /* for the following map, the key is a vector of start node and end node;
     the value is the edge id. It stores each current edge and its
     equivalent edge (size of this map is 2*num_edges)*/
  std::map<std::vector<size_type>, size_type> edge_lookup;

  /* for the following map, the key is a start node and the value is a vector of end nodes;
     it is bidirectional, each edge will be an element of two key's vector values*/
  std::map<Node, std::vector<Node>> map_node_to_node;


};

#endif // CME212_GRAPH_HPP
