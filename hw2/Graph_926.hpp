#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type. Peer code adopted form Graph_1414.hpp"
 */

#include <algorithm>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V, typename E> //Suddenly not passing HW0
class Graph {
  //typedef V node_value_type; 
  //typedef E edge_value_type;

 private:
  /** Predeclaration of internal_node and internal_edge types. */
  struct internal_node;
  struct internal_edge;

 public:
  /** Type of this graph. */
  using graph_type = Graph;

  /** Modifying graph to become class template for edges and nodes*/
  using node_value_type = V;
  using edge_value_type = E;

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
  Graph() : nodes_(), uid_to_idx(), uid_tracker(0), edges_(), 
            eid_tracker(0), eid_to_idx(), map_edges() {}

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
  class Node : private totally_ordered<Node> { // Proxy Node actually
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

    /** Define fetch() similarly to in proxy_example. Fetch helps us get the
     * internal node given a node.
     */
    internal_node& fetch() const {
      return parent_graph_->nodes_.at(parent_graph_->uid_to_idx.at(uid_));
    }

    /** Return this node's position. */
    const Point& position() const {
      return fetch().pos_; 
    }

    /** Return this node's position in a way that is modifiable */
    Point& position() {
      return fetch().pos_;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return parent_graph_->uid_to_idx[fetch().uid_];
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    size_type degree() const { //return the number of incident edges
      return fetch().adj_nodes.size(); 
    }
    
    incident_iterator edge_begin() const {
      //std::cout << "The type of adj_nodes input is: " << typeof(fetch().adj_nodes.begin()) << std::endl;
      return IncidentIterator(parent_graph_, *this, fetch().adj_nodes.begin());
    }
    incident_iterator edge_end() const {
      return IncidentIterator(parent_graph_, *this, fetch().adj_nodes.end());
    }

    //** Return the node's value */
    node_value_type& value() {
      return fetch().value_;
    }

    const node_value_type& value() const {
      return fetch().value_;
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      if (fetch().uid_ == n.uid_ &&
          this->parent_graph_ == n.parent_graph_) {
        return true;
      }

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
      size_type myuid = fetch().uid_;
      bool check_nodes = myuid < n.uid_; 
      bool check_graph = this->parent_graph_ < n.parent_graph_;

      if (check_nodes || check_graph) {
        return true;
      } else {
        return false;
      }
    }

   private:
    friend class Graph;
    size_type uid_;
    graph_type* parent_graph_;

    //This constructor has an empty body!
    Node() : uid_{}, parent_graph_{} {} //everyone that calls this must be friends with the node class
    //uid_ initialized to zero; parent_graph initialized to nullptr
    //If one of these attributes was a standard library container, they have default constructors!!
    //  The default constructors do sensible things

    //this constructor is okay! 
    Node(const graph_type* graph, size_type uid) 
      : uid_(uid), parent_graph_(const_cast<graph_type*>(graph)) {}
  };

  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    return uid_to_idx.size();
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
    Node add_node(const Point& position, const node_value_type& val = node_value_type()) {
    size_type assign_uid = this->uid_tracker;
    uid_tracker += 1;
    size_type idx = nodes_.size();
    uid_to_idx[assign_uid] = idx;

    //This could also be: 
    // internal_node new_internal_node(position, assign_uid, val); //Constructor
    internal_node new_internal_node; //This another way to initialize this struct
    new_internal_node.pos_ = position;
    new_internal_node.uid_ = assign_uid;
    new_internal_node.value_ = val;
    //NodeData new_internal_node;
    new_internal_node.adj_nodes = {}; 
    this->nodes_.emplace_back(new_internal_node); //Adding this new node to the nodes list

    return Node(this, assign_uid);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    size_type n_uid = n.uid_;
    if (uid_to_idx.find(n_uid) == uid_to_idx.end()) {
      return false;
    }
    return true;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * INPUT: Node IDX
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // assert(i < num_nodes());
    return Node(this, (nodes_.at(i)).uid_);
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
    Edge() {}

    /** Help function to get internal_edge given an Edge object. Similar to
     *  fetch() in Node class, and fetch() in proxy_example.cpp.
     */
    internal_edge& edge_fetch() const {
      //std::cout << "edge_fetch() is trying get " << eid_ << std::endl;
      return edge_p_graph_->edges_.at(edge_p_graph_->eid_to_idx.at(eid_));
    }

    /** Return a node of this Edge */
    Node node1() const {
      //std::cout<< "Node 1 calls edge fetch" << std::endl;
      return edge_fetch().n1_;
    }

    /** Return the other node of this Edge */
    Node node2() const {
      //std::cout<< "Node 2 calls edge fetch" << std::endl;
      return edge_fetch().n2_;
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     * Added graph comparison functionality
     */
    bool operator==(const Edge& e) const {
      size_type myeid = edge_fetch().eid_;
      if (this->edge_p_graph_ != e.edge_p_graph_) {
        return false; 
      } else if (myeid == e.eid_) {
        return true;
      } else {
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     * HW2: Added graph comparison functionality
     */
    bool operator<(const Edge& e) const {
      size_type myeid = edge_fetch().eid_;
      bool node_compare = ((node1() < e.node1()) && (node2() < e.node2())); 
      bool graph_compare = this->edge_p_graph_ < e.edge_p_graph_;

      if (this->edge_p_graph_ != e.edge_p_graph_) {
        return graph_compare; //If edges equal
      } else if (myeid == e.eid_) {
        return false;
      } else {
        return node_compare;
      }
    }

    /** Method to return the edge's value */
    edge_value_type& value() {
      return edge_fetch().value_;
    }

    /** Method to return the edge's value */
    const edge_value_type& value() const {
      return edge_fetch().value_;
    }

    /** Calculate the distance between the two nodes that make up an edge */
    double length() const {
      double dist = norm(node1().position() + (-node2().position()));
      return dist; //Return the length
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    graph_type* edge_p_graph_;
    size_type eid_;

    // Initialized for Edge.
    Edge(const graph_type* graph, size_type eid)
      : edge_p_graph_(const_cast<graph_type*>(graph)), eid_(eid) {}
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    return eid_to_idx.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(i < num_edges());
    size_type myeid = edges_[i].eid_;
    return Edge(this, myeid);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // Check if nodes a and b exist
    assert(has_node(a));
    assert(has_node(b));

    std::pair<size_type, size_type> edge_ab = {a.uid_, b.uid_};
    std::pair<size_type, size_type> edge_ba = {b.uid_, a.uid_};

    // if (map_edges.find(edge_ab) != map_edges.end()) {
    //   size_type myidx = eid_to_idx.at(map_edges[edge_ab]);
    //   if ((edges_[myidx].n1_ == a && edges_[myidx].n2_) 
    //     || (edges_[myidx].n1_ == b && edges_[myidx].n2_ == a)) {
    //     return true;
    //   } 
    //   // else if (map_edges.find(edge_ba) != map_edges.end()) {
    //   //   size_type get_eid = map_edges[edge_ba]; 
    //   //   size_type myidx = eid_to_idx.at(get_eid);
    //   //   if ((edges_[myidx].n1_ == a && edges_[myidx].n2_) 
    //   //     || (edges_[myidx].n1_ == b && edges_[myidx].n2_ == a)) {
    //   //     return true;
    //   // } 
    //   else {
    //     return false;
    //   }
    // }
    // but this is more modular and robust!
    if (map_edges.find(edge_ab) == map_edges.end()
     && map_edges.find(edge_ba) == map_edges.end()) {
      return false;
    }
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
    assert(a.uid_ != b.uid_); // a and b are distinct nodes
    assert(has_node(a)); // node a is valid
    assert(has_node(b)); // node b is valid

    //when you add this edge, go into both nodes and add the other adjacent node
    std::pair<size_type, size_type> edge_ab = {a.uid_, b.uid_}; 
    std::pair<size_type, size_type> edge_ba = {b.uid_, a.uid_};

    // std::cout << "I am getting into add_edge" << std::endl;
    //std::cout << "The edge that doesn't exist is: " << a.uid_<< and b.uid_ << std::endl;
    // Check if edge already exists
    if (map_edges.find(edge_ab) != map_edges.end()) { //If this is in the map
      // std::cout << "Does the test get into loop A" << std::endl;
      // size_type get_eid = map_edges.at(edge_ab);
      // std::cout << "The eid of this edge is: " << get_eid << std::endl;
      size_type myidx = eid_to_idx.at(map_edges[edge_ab]);
      return Edge(this, myidx);

    } else if (map_edges.find(edge_ba) != map_edges.end()) {
      // Checking if edge_ba exists but edge_ab doesn't...
      // std::cout << "Does the test get into loop b" << std::endl;
      size_type get_eid = map_edges[edge_ba]; //get edge eid
      size_type ba_idx = eid_to_idx.at(get_eid); //edge the edge index
      map_edges.erase(edge_ba); //Switching the key
      map_edges[edge_ab] = get_eid;
      edges_[ba_idx].n1_ = a;  //Adding a & b to (a,b)
      edges_[ba_idx].n2_ = b;  
      return Edge(this, ba_idx);
    }

    internal_edge new_internal_edge;
    size_type assign_eid = this->eid_tracker;
    eid_tracker += 1;
    size_type idx = edges_.size();
    eid_to_idx[assign_eid] = idx;
    map_edges[edge_ab] = assign_eid; // undirected graph so a->b and b->a are
    //map_edges[edge_ba] = assign_eid; // the same, so same eid_, too! //THIS was causing my erro

    (a.fetch().adj_nodes).push_back(b); //Add b to the adjacent nodes of a
    (b.fetch().adj_nodes).push_back(a); //Add a to the adjacent nodes of b

    new_internal_edge.n1_ = a;
    new_internal_edge.n2_ = b;
    new_internal_edge.eid_ = assign_eid;
    this->edges_.push_back(new_internal_edge);

    return Edge(this, assign_eid);
  }

  /** Remove edge from this graph
   * @post num_edges -= 1
   * @pre Two different nodes
   * 
   * Invalidates the edge
   */
  size_type remove_edge(const Node& a, const Node& b) {
    if (has_edge(a, b) == false) {
      return false; //Edge already removed
    } else {
      size_type a_idx = uid_to_idx.at(a.uid_);
      size_type b_idx = uid_to_idx.at(b.uid_);

      std::vector<Node> &a_neighbors = nodes_[a_idx].adj_nodes;
      std::vector<Node> &b_neighbors = nodes_[b_idx].adj_nodes;

      for(unsigned int i = 0; i < a_neighbors.size(); ++i) {
        if (a_neighbors[i].uid_ == b.uid_) {
          a_neighbors[i] = a_neighbors.back(); 
          a_neighbors.pop_back(); //remove b reference to myuid
        }
      } 

      for(unsigned int i = 0; i < b_neighbors.size(); ++i) {
        if (b_neighbors[i].uid_ == a.uid_) {
            b_neighbors[i] = b_neighbors.back(); 
            b_neighbors.pop_back(); 
        }
      } 

      //Find the nodes in the map
      //  If not (a,b), then switch to (b,a);
      std::pair<size_type, size_type> edge_ab = {a.uid_, b.uid_};
      if (map_edges.find(edge_ab) == map_edges.end()) { 
        edge_ab = {b.uid_, a.uid_};
      }      

      size_type myeid = map_edges.at(edge_ab);
      size_type idx = eid_to_idx.at(myeid);

      // swap and pop internal_edges
      internal_edge last_edge = edges_.back();
      edges_[idx] = edges_.back();
      edges_.pop_back();

      // // swap and pop eid_to_idx map
      eid_to_idx[last_edge.eid_] = idx;
      eid_to_idx.erase(myeid);
      
      //Remove the edges stored in map_edges
      map_edges.erase(edge_ab);

      return true;
    }
  }

  /** Remove the edge from the graph given edge */
  size_type remove_edge(const Edge& e) {
    remove_edge(e.node1(), e.node2());
  }

  /** Remove edge using an edge_iterator 
   * @pre Two different notes
   * @post num_edges -= 1
  */
  edge_iterator remove_edge(edge_iterator e_it) {
    //dereference the iterator and call remove edge
    Edge e = *e_it;
    if(remove_edge(e)){
      return (*this).edge_begin();
    } else{
      return ++e_it; 
    }
  }

  /** Remove a node from this graph *
   * @post num_nodes -= 
   * 
   * @output:
   *  - True: if we remove a node
   *  - False: if we don't remove a node
  */
  size_type remove_node(const Node& n) {
    if (has_node(n) == false) {//Don't remove a node that doesn't exist
      return false;
    } else { //Remove the node
      size_type myuid = n.uid_;
      size_type idx = uid_to_idx.at(myuid);

      std::vector<Node> neighs_copy = nodes_[idx].adj_nodes;

      for (unsigned int i = 0; i < neighs_copy.size(); ++i) {
        size_type b_uid = neighs_copy[i].uid_;
        size_type b_idx = uid_to_idx.at(b_uid);
        const Node& b = node(b_idx);
        remove_edge(n, b);
      }

      internal_node last_node = nodes_.back();
      nodes_[idx] = nodes_.back();
      nodes_.pop_back();

      //Update uid_to_idx
      uid_to_idx.at(last_node.uid_) = idx; 
      uid_to_idx.erase(myuid);
      return true;
    }
  }


  node_iterator remove_node(node_iterator n_it) {
    Node n = *n_it; 
    if(remove(n)) {
      return(*this).node_begin();
    } else {
      return ++n_it;
    }
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    //set trackers to zero
    uid_tracker = 0;
    eid_tracker = 0;

    //clear edges and nodes
    edges_.clear();
    nodes_.clear();
    uid_to_idx.clear();
    eid_to_idx.clear();
    map_edges.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() : parent_graph_{}, node_iter_{} {}

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Node operator*() const{ //Dereference operator
      return node_type(parent_graph_, (*node_iter_).uid_); //Return this node
    }

    NodeIterator& operator++() {
      ++node_iter_;
      return *this; // return the lvalue reference
    }
    
    // Defines equality between two nodes. 
    // Two nodes are equal if the iuds are the same and the node_idx is the same
    bool operator==(const NodeIterator& node_iter_b) const{
      return this->node_iter_ == node_iter_b.node_iter_;
    }

    //bool operator!=(const NodeIterator& node_iter_b) const{
    //  return this->node_iter_ != node_iter_b.node_iter_;
   // }


   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE

    // NodeIterator Data
    Graph* parent_graph_; //parent_graph = pointer to a graph
    typename std::vector<internal_node>::const_iterator node_iter_;

    NodeIterator(const Graph* parent_graph, //Const graph
      typename std::vector<internal_node>::const_iterator node_iter) : 
      parent_graph_(const_cast<graph_type*> (parent_graph)), node_iter_(node_iter) //Initialize list
      {} //The implementation body - could do stuff to setup the node iterator
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  node_iterator node_begin() const { //Return node iterator at the beginning
    return NodeIterator(this, nodes_.begin());
  }

  node_iterator node_end() const { //Return node iterator at the end
    return NodeIterator(this, nodes_.end());

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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator() : incident_parent_graph_{}, root_node_{}, incident_iter_{} {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Dereferences the edge to return it */
    Edge operator*() const{ //return this edge
      Node a = Node(incident_parent_graph_, (*incident_iter_).uid_);
      
      return incident_parent_graph_->add_edge(root_node_, a); 
    }

    IncidentIterator& operator++() {
      ++incident_iter_;
      return *this; // return the lvalue reference
    }

    bool operator==(const IncidentIterator& incident_iter_b) const {
      return this->incident_iter_ == incident_iter_b.incident_iter_;
    }


   private:
    friend class Graph;
    //friend class Edge;
    // HW1 #3: YOUR CODE HERE
    //The incident iterator includes all the adjacent nodes
    Graph* incident_parent_graph_; 
    Node root_node_;
    typename std::vector<Node>::const_iterator incident_iter_;

   IncidentIterator(const Graph* parent_graph, Node root_node,  
   typename std::vector<Node>::const_iterator iter) :
    incident_parent_graph_(const_cast<graph_type*> (parent_graph)), 
    root_node_(root_node), incident_iter_(iter) {}
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator() : parent_graph_{}, edge_iter_{} {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    Edge operator*() const {
      return edge_type(parent_graph_, (*edge_iter_).eid_); //return this edge
    }

    EdgeIterator& operator++() {
      ++edge_iter_; 
      return *this; //return the lvalue reference
    }

    bool operator==(const EdgeIterator& edge_iter_b) const {
      return this->edge_iter_ == edge_iter_b.edge_iter_;
    }
    
    //bool operator!=(const EdgeIterator& edge_iter_b) const {
    //  return this->edge_iter_ != edge_iter_b.edge_iter_;
    //}

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE

    //EdgeIterator Data
    Graph* parent_graph_; //pointer to the graph
    typename std::vector<internal_edge>::const_iterator edge_iter_;

    EdgeIterator(const Graph* parent_graph, //Const graph pointer
      typename std::vector<internal_edge>::const_iterator edge_iter) : 
      parent_graph_(const_cast<graph_type*> (parent_graph)), edge_iter_(edge_iter) {} //The implementation body - could do stuff to setup the node iterator

  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  edge_iterator edge_begin() const {
    return EdgeIterator(this, edges_.begin());
  }

  edge_iterator edge_end() const {
    return EdgeIterator(this, edges_.end());
  }

 private:
  // Nodes data
  std::vector<internal_node> nodes_;
  std::unordered_map<size_type, size_type> uid_to_idx;
  size_type uid_tracker;

  // Each internal_node includes a uid (id), position, and a value
  // Also, the internal node keeps track of all the adjacent nodes
  struct internal_node {
    internal_node() {}
    //internal_node(size_type uid, const Point& pos, node_value_type& value) : uid_(uid), pos_(pos) {}
    size_type uid_;
    Point pos_;
    node_value_type value_; //Adding the value attribute
    std::vector<Node> adj_nodes; //vector of adjacent nodes
    //const node_value_type& value_ const; 
  };

  // Edges data
  std::vector<internal_edge> edges_;
  size_type eid_tracker;
  std::unordered_map<size_type, size_type> eid_to_idx;
  std::map<std::pair<size_type, size_type>, size_type> map_edges; // {(n1, n2): eid, etc}

  // Each internal_edge includes two nodes and an eid (id)
  struct internal_edge {
    internal_edge() {} //
    //internal_edge(size_type eid, Node n1, Node n2)
    //  : eid_(eid), n1_(n1), n2_(n2) {}
    size_type eid_;
    Node n1_;
    Node n2_;
    edge_value_type value_; //Runs the value struct
  };
};

#endif // CME212_GRAPH_HPP