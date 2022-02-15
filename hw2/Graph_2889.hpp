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


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E>
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
  Graph()
    : nodes_(), node_counts_(0), next_nid_(0), i2u_n(),
      edges_(), edge_counts_(0), next_eid_(0), i2u_e(),
      edge_lookup(), node_incidents(), eid2val() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
    for (size_type i = 0; i < nodes_.size(); i++) {
      internal_node* node_ptr = nodes_[i];
      delete node_ptr;
    }
    for (size_type i = 0; i < edges_.size(); i++) {
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
    Node() {
      // HW0: YOUR CODE HERE
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE
      return fetch().coord;
    }

    // HW2: Modifialble Node Position
    Point& position() {
      return fetch().coord;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return (*(graph_->nodes_[nid_])).idx_n;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    //public function to write node's value
    node_value_type& value() {
      return fetch().val;
    };
    //public function to read node's value
    const node_value_type& value() const {
      return fetch().val;
    };

    size_type degree() const {
      if (graph_->node_incidents.find(nid_) == graph_->node_incidents.end()) {
        return 0;
      } else {
        return graph_->node_incidents[nid_].size();
      }
    }

    incident_iterator edge_begin() const {
      if (this->degree() == 0) {
        return IncidentIterator();
      } else {
        return IncidentIterator(graph_, nid_, 0);
      }
    }
    incident_iterator edge_end() const {
      return IncidentIterator();
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return ((graph_ == n.graph_) && (index() == n.index())
              && (nid_ == n.nid_) );
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
      return (index() < n.index());
    }

    size_type get_nid() {
      return nid_;
    }

    bool pub_valid() {
      return valid();
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // Allow IncidentIterator to use Node's private member data and functions.
    //friend class IncidentIterator;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // Pointer back to the Graph container
    Graph* graph_;
    // This node's id
    size_type nid_;
    /** Private Constructor that can be accessed by Graph **/
    Node(const Graph* graph, size_type nid)
      : graph_(const_cast<Graph*>(graph)), nid_(nid) {
      }
    /** Helper method to return the data associated with node that is
        stored in Graph **/
    internal_node& fetch() const {
      if (graph_->has_node(*this)){
        return *(graph_->nodes_[nid_]);
      }else{
        throw("node is invalid");
      }
    }
    bool valid() const {
      return nid_ >= 0 && nid_<graph_->nodes_.size()
             && graph_->nodes_[nid_]->idx_n < graph_->i2u_n.size()
             && graph_->i2u_n[graph_->nodes_[nid_]->idx_n] == nid_;
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
  Node add_node(const Point& position,
                const node_value_type& default_val = node_value_type()) {
    // HW0: YOUR CODE HERE
    nodes_.push_back(new internal_node{next_nid_, position, default_val, (size_type)i2u_n.size()});
    i2u_n.push_back(next_nid_);
    next_nid_ += 1;
    node_counts_ += 1;
    return Node(this, next_nid_-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    size_type idx = n.graph_->nodes_[n.nid_]->idx_n;
    return ((idx < i2u_n.size()) && (n.graph_ == this));
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */

  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    return Node(this, i2u_n[i]);
  }


  /** Delete node n from this graph object
   *
   * @param[in] n node object
   * @return 0 if n doesn't exist in graph, 1 if n is successfully removed
   *
   *@pre n is a valid node
   *@pre old node_counts_ > 0
   *@post new node_counts_ = old node_counts - 1
   *@post new i2u_n.size() = old i2u_n.size() -1
   *@post all incident edges of node n are removed (see function remove_edge)
   *@post new edge_lookup.size() <= old edge_lookup.size()
   *@post new graph stores exactly the same internal information
   *      except for those relevant to the node deleted by this function
   *
   * The time complexity is O(max(d^2, d*log(num_edges))). where d is the average
   * degree of nodes in the graph. assumes d is very small for sparse graph
  */

  size_type remove_node(const Node& n) {
    if (!has_node(n)) return 0;
    assert(n.valid());
    /* remove all incident edges. O(max(d^2, d*log(num_edges))). assumes d is
       very small for sparse graph */
    //create a list of edge to remove O(d)
    std::vector<Edge> edge2rm{};
    for (size_type i = 0; i < node_incidents[n.nid_].size(); ++i) {
      Node neighbor = Node(this, node_incidents[n.nid_][i]);
      edge2rm.push_back(Edge(n, neighbor));
    }
    for (size_type i = 0; i < edge2rm.size(); ++i) {
      Edge inc_edge = edge2rm[i];
      remove_edge(inc_edge); //O(max(d, log(num_edges)))
    }
    /* modify idx to unique nid vector, swap and pop. O(1) */
    size_type rmidx = n.index();
    size_type swapidx = i2u_n.size()-1;
    size_type nid_end = i2u_n.back();
    (*nodes_[nid_end]).idx_n = rmidx; //change the idx for edge affected by swap
    (*nodes_[n.nid_]).idx_n = swapidx; //idx_n for nid_ of removed node is i2u_n.size() after the pop
    std::iter_swap(i2u_n.begin()+rmidx, i2u_n.end()-1); //O(1)
    i2u_n.pop_back(); //O(1)

    node_counts_ -= 1;
    return 1;
  }

  /** Delete the node n_it refers to from this graph object
   *
   * @param[in] n_it node iterator
   * @return this node iterator after swapping the internal information
   *
   *@pre n_it is a valid node iterator
   *@pre old node_counts_ > 0
   *@post new node_counts_ = old node_counts - 1
   *@post new i2u_n.size() = old i2u_n.size() -1
   *@post all incident edges of node n are removed (see function remove_edge)
   *@post new edge_lookup.size() <= old edge_lookup.size()
   *@post new graph stores exactly the same internal information
   *      except for those relevant to the node deleted by this function
   *
   * The time complexity is O(max(d^2, d*log(num_edges))). where d is the average
   * degree of nodes in the graph. assumes d is very small for sparse graph
  */
  node_iterator remove_node(node_iterator n_it) {
    if (n_it == node_end()) return n_it;
    Node n = *n_it;
    assert(n.valid());
    /* remove all incident edges. O(max(d^2, d*log(num_edges))). assumes d is
       very small for sparse graph */
   //create a list of edge to remove O(d)
   std::vector<Edge> edge2rm{};
   for (size_type i = 0; i < node_incidents[n.nid_].size(); ++i) {
     Node neighbor = Node(this, node_incidents[n.nid_][i]);
     edge2rm.push_back(Edge(n, neighbor));
   }
   for (size_type i = 0; i < edge2rm.size(); ++i) {
     Edge inc_edge = edge2rm[i];
     remove_edge(inc_edge); //O(max(d, log(num_edges)))
   }
    /* modify idx to unique nid vector, swap and pop. O(1) */
    size_type rmidx = n.index();
    size_type swapidx = i2u_n.size()-1;
    size_type nid_end = i2u_n.back();
    (*nodes_[nid_end]).idx_n = rmidx; //change the idx for edge affected by swap
    (*nodes_[n.nid_]).idx_n = swapidx; //idx_n for nid_ of removed node is i2u_n.size() after the pop
    std::iter_swap(i2u_n.begin()+rmidx, i2u_n.end()-1); //O(1)
    i2u_n.pop_back(); //O(1)

    node_counts_ -= 1;
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

   private:
   // Allow Graph to access Edge's private member data and functions.
   friend class Graph;
   // Allow IncidentIterator to use Edge's private member data and functions.
   //friend class IncidentIterator;
   // HW0: YOUR CODE HERE
   // Use this space to declare private data members and methods for Edge
   // that will not be visible to users, but may be useful within Graph.
   // i.e. Graph needs a way to construct valid Edge objects

   // two nodes that defines this edge
   Node n1_, n2_;

   /** Private Constructor **/
   Edge(const Node& n1, const Node& n2)
     : n1_(n1), n2_(n2) {}

   public:
    /** Construct an invalid Edge. */
    Edge() : n1_(), n2_() {
      // HW0: YOUR CODE HERE
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return n2_;
    }

    // HW2 Edge length
    double length() const {
      return norm(n1_.position() - n2_.position());
    }

    // HW2
    edge_value_type& value() {
      std::vector<size_type> edge_key{n1_.nid_, n2_.nid_};
      size_type eid = n1_.graph_->edge_lookup[edge_key];
      return n1_.graph_->eid2val[eid];
    }

    const edge_value_type& value() const {
      std::vector<size_type> edge_key{n1_.nid_, n2_.nid_};
      size_type eid = n1_.graph_->edge_lookup[edge_key];
      return n1_.graph_->eid2val[eid];
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
      if (*this == e) {
        return false;
      } else {
          if (n1_.graph_ == e.node1().graph_) {
            return ( (n1_.index() + n2_.index() ) <
                     (e.node1().index() + e.node2().index()) );
          } else {
            return (n1_.graph_ < e.node1().graph_);
          }
      }
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
    size_type n1_id = edges_[i2u_e[i]]->n1;
    size_type n2_id = edges_[i2u_e[i]]->n2;
    Node node1(this, n1_id);
    Node node2(this, n2_id);
    return Edge(node1, node2);
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
      return Edge(a, b);
    }
    internal_edge* new_edge = new internal_edge{next_eid_, a.nid_, b.nid_, (size_type)i2u_e.size()};
    edges_.push_back(new_edge);
    edge_lookup.insert({edge_key, next_eid_});
    std::vector<size_type> edge_key_rev{b.nid_, a.nid_};
    edge_lookup.insert({edge_key_rev, next_eid_});
    eid2val[next_eid_] = edge_value_type();
    i2u_e.push_back(next_eid_);
    next_eid_ += 1;
    edge_counts_ += 1;
    // store "b" as a connected node for "a"
    if (node_incidents.find(a.nid_) == node_incidents.end()) {
      //initiailize key node "a" if not already in the map
      std::vector<size_type> incident_nodes;
      node_incidents.insert({a.nid_, incident_nodes});
    }
    node_incidents[a.nid_].push_back(b.nid_);
    // store "a" as a connected node for "b"
    if (node_incidents.find(b.nid_) == node_incidents.end()) {
      //initiailize key node "b" if not already in the map
      std::vector<size_type> incident_nodes;
      node_incidents.insert({b.nid_, incident_nodes});
    }
    node_incidents[b.nid_].push_back(a.nid_);
    return Edge(a, b);
  }

  /** Delete the edge connecting n1 and n2 from this graph object
   *
   * @param[in] n1 and n2 of node type
   * @return 0 if edge doesn't exist, 1 if edge is successfully removed
   *
   *@pre has_edge(n1, n2)
   *@pre old edge_counts_ > 0
   *@post new edge_counts_ = old edge_counts - 1
   *@post new i2u_e.size() = old i2u_e.size() -1
   *@post Let n1id and n2id be the unique node ids for the two nodes of the
   *      edge to be deleted
   *      new node_incidents[n1id].size() < old node_incidents[n1id].size()
   *      new node_incidents[n2id].size() < old node_incidents[n2id].size()
   *@post new edge_lookup.size() < old edge_lookup.size()
   *@post new graph stores exactly the same internal information
   *      except for those relevant to the edge deleted by this function
   *
   * The time complexity is max(O(log(num_edge)), O(d)) where d is the average
   * degree of nodes in the graph
  */
  size_type remove_edge(const Node& n1, const Node& n2) {
    if (!has_edge(n1, n2)) return 0;
    size_type n1id = n1.nid_;
    size_type n2id = n2.nid_;
    /* modify lookup map for the edge existence O(log(num_edges)) */
    std::vector<size_type> edge_key{n1id, n2id};
    std::vector<size_type> edge_key_rev{n2id, n1id};
    size_type eid = edge_lookup[edge_key];
    edge_lookup.erase(edge_key);  //O(log(num_edges))
    edge_lookup.erase(edge_key_rev);//O(log(num_edges))

    /* modify each node's incident nodes
     * O(d) where d is the degree of each node, in general a lot smaller than
     * num_nodes for a sparse graph */
    std::vector<size_type> n1_inc = node_incidents[n1id];
    for (auto it = n1_inc.begin(); it != n1_inc.end(); ){
      if (*it == n2id) {
        it = n1_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n1id] = n1_inc;
    std::vector<size_type> n2_inc = node_incidents[n2id];
    for (auto it = n2_inc.begin(); it != n2_inc.end(); ){
      if (*it == n1id) {
        it = n2_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n2id] = n2_inc;

    size_type idx = (*edges_[eid]).idx_e; //find idx of edge O(1)

    /* modify idx to unique eid vector, swap and pop. O(1) */
    size_type eid_end = i2u_e.back();
    (*edges_[eid_end]).idx_e = idx; //change the idx for edge affected by swap
    std::iter_swap(i2u_e.begin() + idx, i2u_e.end()-1); //O(1)
    i2u_e.pop_back(); //O(1)

    edge_counts_ -= 1;
    return 1;
  }

  /** Delete the edge e_ from this graph object
   *
   * @param[in] e_ edge object
   * @return 0 if edge doesn't exist, 1 if edge is successfully removed
   *
   *@pre e_ is a valid edge
   *@pre old edge_counts_ > 0
   *@post new edge_counts_ = old edge_counts - 1
   *@post new i2u_e.size() = old i2u_e.size() -1
   *@post Let n1id and n2id be the unique node ids for the two nodes of the
   *      edge to be deleted
   *      new node_incidents[n1id].size() < old node_incidents[n1id].size()
   *      new node_incidents[n2id].size() < old node_incidents[n2id].size()
   *@post new edge_lookup.size() < old edge_lookup.size()
   *@post new graph stores exactly the same internal information
   *      except for those relevant to the edge deleted by this function
   *
   * The time complexity is max(O(log(num_edge)), O(d)) where d is the average
   * degree of nodes in the graph
  */
  size_type remove_edge(const Edge& e_) {
    if (!has_edge(e_.node1(), e_.node2() )) return 0;
    size_type n1id = e_.node1().nid_;
    size_type n2id = e_.node2().nid_;
    /* modify lookup map for the edge existence O(log(num_edges)) */
    std::vector<size_type> edge_key{n1id, n2id};
    std::vector<size_type> edge_key_rev{n2id, n1id};
    size_type eid = edge_lookup[edge_key];
    edge_lookup.erase(edge_key);  //O(log(num_edges))
    edge_lookup.erase(edge_key_rev);//O(log(num_edges))

    /* modify each node's incident nodes
     * O(d) where d is the degree of each node, in general a lot smaller than
     * num_nodes for a sparse graph */
    std::vector<size_type> n1_inc = node_incidents[n1id];
    for (auto it = n1_inc.begin(); it != n1_inc.end(); ){
      if (*it == n2id) {
        it = n1_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n1id] = n1_inc;
    std::vector<size_type> n2_inc = node_incidents[n2id];
    for (auto it = n2_inc.begin(); it != n2_inc.end(); ){
      if (*it == n1id) {
        it = n2_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n2id] = n2_inc;

    size_type idx = (*edges_[eid]).idx_e; //find idx of edge O(1)

    /* modify idx to unique eid vector, swap and pop. O(1) */
    size_type eid_end = i2u_e.back();
    (*edges_[eid_end]).idx_e = idx; //change the idx for edge affected by swap
    std::iter_swap(i2u_e.begin() + idx, i2u_e.end()-1); //O(1)
    i2u_e.pop_back(); //O(1)

    edge_counts_ -= 1;
    return 1;
  }



  /** Delete the edge which e_it reference to from this graph object
   *
   * @param[in] e_it edge_iterator
   * @return the same edge iterator after swapping internal information
   *
   *@pre e_it != edge_end() and e_it is a valid edge iterator
   *@pre old edge_counts_ > 0
   *@post new edge_counts_ = old edge_counts - 1
   *@post new i2u_e.size() = old i2u_e.size() -1
   *@post Let n1id and n2id be the unique node ids for the two nodes of the
   *      edge to be deleted
   *      new node_incidents[n1id].size() < old node_incidents[n1id].size()
   *      new node_incidents[n2id].size() < old node_incidents[n2id].size()
   *@post new edge_lookup.size() < old edge_lookup.size()
   *@post new graph stores exactly the same internal information
   *      except for those relevant to the edge deleted by this function
   *
   * The time complexity is max(O(log(num_edge)), O(d)) where d is the average
   * degree of nodes in the graph
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    if (e_it == edge_end()) return e_it;
    edge_type edge_ = *e_it;
    size_type n1id = edge_.node1().nid_;
    size_type n2id = edge_.node2().nid_;
    /* modify lookup map for the edge existence O(log(num_edges)) */
    std::vector<size_type> edge_key{n1id, n2id};
    std::vector<size_type> edge_key_rev{n2id, n1id};
    edge_lookup.erase(edge_key);  //O(log(num_edges))
    edge_lookup.erase(edge_key_rev);//O(log(num_edges))

    /* modify each node's incident nodes
     * O(d) where d is the degree of each node, in general a lot smaller than
     * num_nodes for a sparse graph */
    std::vector<size_type> n1_inc = node_incidents[n1id];
    for (auto it = n1_inc.begin(); it != n1_inc.end(); ){
      if (*it == n2id) {
        it = n1_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n1id] = n1_inc;
    std::vector<size_type> n2_inc = node_incidents[n2id];
    for (auto it = n2_inc.begin(); it != n2_inc.end(); ){
      if (*it == n1id) {
        it = n2_inc.erase(it);
      } else {
        ++it;
      }
    }
    node_incidents[n2id] = n2_inc;

    /* modify idx to unique eid vector, swap and pop. O(1) */
    size_type eid_end = i2u_e.back();
    (*edges_[eid_end]).idx_e = e_it.idx_; //change the idx for edge affected by swap
    std::iter_swap(i2u_e.begin()+e_it.idx_, i2u_e.end()-1); //O(1)
    i2u_e.pop_back(); //O(1)

    edge_counts_ -= 1;
    return e_it;
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
    node_incidents.clear();
    eid2val.clear();
    edge_lookup.clear();
    i2u_n.clear();
    i2u_e.clear();
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
  class NodeIterator : private totally_ordered<NodeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Node;                     // Element type
    using pointer           = Node*;                    // Pointers to elements
    using reference         = Node&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : graph_(nullptr), idx_(0){
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    Node operator*() const {
        return graph_->node(idx_);
    }

    NodeIterator& operator++() {
        if (idx_ < graph_->num_nodes() - 1) {
          ++idx_;
        } else {
          graph_ = nullptr;
          idx_ = 0;
        }
        return *this;
    }

    bool operator==(const NodeIterator& node_iter) const {
      return ((graph_ == node_iter.graph_) &&
              (idx_ == node_iter.idx_));
    }

    size_type get_idx() const {
      return idx_;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph* graph_;
    size_type idx_;
    //Private constructor that can be accessed by the Graph class
    NodeIterator(const Graph* graph , size_type idx)
     : graph_(const_cast<Graph*>(graph)), idx_(idx) {}
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const {
    return NodeIterator(this, 0);
  }
  node_iterator node_end() const {
    return NodeIterator();
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
    IncidentIterator()
     : graph_(nullptr), spawn_node_id(0), iter_idx(0) {
    }

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:

    Edge operator*() const {
        size_type other_nid = graph_->node_incidents[spawn_node_id][iter_idx];
        Node node1(graph_, spawn_node_id);
        Node node2(graph_, other_nid);
        return Edge(node1, node2);
    }

    IncidentIterator& operator++() {
      size_type idx = graph_->nodes_[spawn_node_id]->idx_n;
      size_type num_incident_edges = (graph_->node(idx)).degree();
      if (iter_idx < num_incident_edges - 1) {
        ++iter_idx;
      } else {
        graph_ = nullptr;
        iter_idx = 0;
        spawn_node_id = 0;
      }
      return *this;
    }

    bool operator==(const IncidentIterator& inc_iter) const {
      return ((graph_ == inc_iter.graph_) &&
             (spawn_node_id == inc_iter.spawn_node_id) &&
             (iter_idx == inc_iter.iter_idx));
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* graph_;
    size_type spawn_node_id;
    size_type iter_idx;
    IncidentIterator(const Graph* graph, size_type _nid, size_type _iter_idx)
     : graph_(const_cast<Graph*>(graph)), spawn_node_id(_nid),
       iter_idx(_iter_idx) {}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator : private totally_ordered<EdgeIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : graph_(nullptr), idx_(0) {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return graph_->edge(idx_);
    }
    EdgeIterator& operator++() {
      if (idx_ < graph_->num_edges() - 1) {
        ++idx_;
      } else {
        graph_ = nullptr;
        idx_ = 0;
      }
      return *this;
    }
    bool operator==(const EdgeIterator& edge_iter) const {
      return ((graph_ == edge_iter.graph_) &&
             (idx_ == edge_iter.idx_));
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* graph_;
    size_type idx_;
    EdgeIterator(const Graph* graph , size_type idx)
     : graph_(const_cast<Graph*>(graph)), idx_(idx) {}

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }
  edge_iterator edge_end() const {
    return EdgeIterator();
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  struct internal_node{
    size_type nid; //unique node id
    Point coord;
    node_value_type val;
    size_type idx_n; //user interface node id
  };

  struct internal_edge{
    size_type eid; //unique edge id
    size_type n1;
    size_type n2;
    size_type idx_e; //user interface edge id
  };

  std::vector<internal_node*> nodes_; //indexed by nid
  size_type node_counts_;
  size_type next_nid_;
  std::vector<size_type> i2u_n; //indexed by user idx for node

  std::vector<internal_edge*> edges_; //indexed by eid
  size_type edge_counts_;
  size_type next_eid_;
  std::vector<size_type> i2u_e; //indexed by user idx for edge

  /* for the following map, the key is a vector of start node and end node;
     the value is the edge id. It stores each current edge and its
     equivalent edge (size of this map is 2*num_edges)*/
  std::map<std::vector<size_type>, size_type> edge_lookup;
  /* node_incidents maps from node id to a vector of node ids that the
     key node is connected to by edges */
  std::map<size_type, std::vector<size_type>> node_incidents;
  /* eid2val maps from edge id to edge value */
  std::map<size_type, edge_value_type> eid2val;


};

#endif // CME212_GRAPH_HPP
