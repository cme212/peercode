#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <cassert>
#include <map>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V=int,typename E=double>
class Graph {
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition.
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)

 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //

  /** Type of this graph. */
  using graph_type = Graph<V,E>;

  /** Type of node value. */
  using node_value_type=V;
  /** Predeclaration of Node type. */
  class Node;
  /** Synonym for Node (following STL conventions). */
  using node_type = Node;

  using edge_value_type=E;
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

  /** Type of node id to edge index unordered map, which is used as
      a nested map to match id of node2 of the edge to index of the
      edge, this map will be matched with id of node1 as the key in
      the outer undered map, i.e., node1's incident edge map.  */
  using u_map=std::unordered_map<size_type,size_type>;
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
  }

  /** Default destructor */
  ~Graph() {
    clear(); // free up all dynamically allocated memory
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
    /** Construct an invalid node
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
      assert(G_!=nullptr);// invalid node does not have a position
      // retrieve position from internal node attribute pointer
      return G_->node_vec[id()]->node_pos;
    }

    /** Return position by reference, modifiable. */
    Point& position(){
      assert(G_!=nullptr);// invalid node does not have a position
      return G_->node_vec[id()]->node_pos;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      return node_idx;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Retrieve the value associated with the node from its attribute,
        modifiable. */
    node_value_type& value(){
      assert(G_!=nullptr);
      return G_->node_vec[id()]->node_val;
    }

    /** Retrieve value of the node
        Return as const to avoid modification of the value. */
    const node_value_type& value() const{
      assert(G_!=nullptr);
      return G_->node_vec[id()]->node_val;
    }

    // Return the size of incident edge vectors stored in node attribute
    size_type degree() const {
      if(G_==nullptr)
        return 0; // Invalid Node have zero degree
      return G_->node_2_edge[id()].size();
    }
    
    // Start the incident iterator from beginning of the node's incident edge map 
    incident_iterator edge_begin() const {
      return IncidentIterator(G_, node_idx, G_->node_2_edge[id()].begin());
    }
    
    // Incident iterator end with the end of the node's incident edge map
    incident_iterator edge_end() const {
      return IncidentIterator(G_, node_idx, G_->node_2_edge[id()].end());
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      /** Check whether they share the same Graph pointer
         to the graph where they belongs to and whether
         the two nodes have the same index. */
      if((G_==nullptr)&&(n.G_==nullptr))
        return true; // two invalid nodes will be considered equal
      if((G_==n.G_) && (node_idx==n.node_idx)) {
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
      // HW0: YOUR CODE HERE
      // order based on node index
      if (*this==n)
        return false;
      if (node_idx<n.node_idx) {
        return true;
      }else if(node_idx==n.node_idx) {
        return G_<n.G_; // same indexed node from different graphs compare graph pointer value
      }
      return false;
    }

    // Retrive the id associated with the node from node_idx_uid
    size_type id() const {
      assert(G_!=nullptr); // invalid node should not have a valid uid
      return G_->n_idx_uid[node_idx];
    }


   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects

    // node idx used to get node id from n_idx_uid
    size_type node_idx=size_type();
    Graph* G_=nullptr; // Pointer to the graph where the node belongs to

    /** Private node constructor to construct a valid Node */
    Node(const Graph* G, size_type idx) {
      G_=const_cast<Graph*>(G);
      node_idx=idx;
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    // # of edges in the graph should be # of ids of active nodes
    return n_idx_uid.size();
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
  Node add_node(const Point& position, const node_value_type& val= node_value_type()) {
    // HW0: YOUR CODE HERE
    // Create new internal node attribute pointer and insert
    Node_attr* n=new Node_attr(position,val);
    node_vec.push_back(n);
    // unique id of new node should be # nodes ever added to the graph
    n_idx_uid.push_back(node_vec.size()-1);
    // create empty incident map storing indices of edges incident to the node
    u_map incid_map;
    node_2_edge.insert(std::make_pair(node_vec.size()-1,incid_map));
    return Node(this,size()-1);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE 
    /* the node should exist in this graph if its index is less than total #
       active nodes in the graph. */
    if (n.G_==this&&n.node_idx<size())
      return true;
    return false;
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    // HW0: YOUR CODE HERE
    // if query an existing node, return a copy of Node with such input index
    if (i<size()) {
      return Node(this,i);
    }
    return Node(); // return invalid node if not exist
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
      assert(G_!=nullptr); // need to query a valid edge
      return *(n1);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      assert(G_!=nullptr); // need to query a valid edge
      return *(n2);
    }

    /** Return the value of an edge. */
    edge_value_type& value(){
      assert(G_!=nullptr);// invalid Edge should not have value
      return G_->edge_vec[id()]->edge_val;
    }

    /** Return the value of an edge as constant ref. */
    const edge_value_type& value() const{
      assert(G_!=nullptr);// invalid Edge should not have value
      return G_->edge_val[id()]->edge_val;
    }

    // Return edge's uid based on its index
    size_type id() const {
      // invalid Edge should not have a valid id
      assert(G_!=nullptr);
      return G_->e_idx_uid[edge_idx];
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      // Check whether two edges are in the same graph and
      // share the same connected nodes
      if (G_==nullptr && e.G_==nullptr)
        return true; // consider both invalid edges as equivalent
      if (G_!=e.G_) { // equal edges should come from the same graph
        return false;
      }else if (((*n1==*e.n1 && *n2==*e.n2)||(*n1==*e.n2 && *n2==*e.n1)) &&
edge_idx==e.edge_idx) { // check end nodes and edge idx
        return true;
      }else{
        return false;
      }
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE
      if (*this==e)
        return false;
      if (edge_idx<e.edge_idx) { // order by edge index
        return true;
      } else if(edge_idx==e.edge_idx) {
        // compare graph pointer if two edges share the same other attribtues
        return G_<e.G_;
      }
      return false;
    }


    // Return the position difference btw node1 and node2
    Point node_diff() const{
      return (*n1).position()-(*n2).position();
    }

    // Return the length of an edge based on Euclidean dist
    // btw node1 and node2
    double length() const{
      return norm_2(node_diff());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    // pointers to two nodes connected by the edge
    const Node* n1;
    const Node* n2;
    // edge index corresponds to its index in e_idx_uid
    size_type edge_idx=size_type();
    Graph* G_=nullptr; // Graph pointer to the graph where the edge belongs to

    /** Private constructor to construct a valid edge. */
    Edge(const Graph* G, const Node* a, const Node* b, size_type idx) {
      G_=const_cast<Graph*>(G);
      n1=const_cast<Node*>(a);
      n2=const_cast<Node*>(b);
      edge_idx=idx;
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return e_idx_uid.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    if (i<num_edges()) {
      return Edge(this,edge_vec[e_idx_uid[i]]->n1,
edge_vec[e_idx_uid[i]]->n2,i); // return deferenced edge pointer
    }
    return Edge(); // return an invalid Edge if cannot find input indexed edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    /* use helper function which returns a internal edge attribute pointer
       with edge index pair to such edge cononecting the two nodes. */
    std::pair<Edge_attr*,size_type> temp=get_edge(a,b);
    if (temp.first!=nullptr)
      return true;
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
    // make sure input valid nodes
    assert(a.index()<size() && b.index()<size());
    std::pair<Edge_attr*,size_type> temp=get_edge(a,b);
    /* Directly return a newly created Edge object with a as node1, b as
       node2. */
    if (temp.first!=nullptr) {
      return Edge(this,&a,&b,temp.second);
    }

    // new edge's index should be # active edges in the graph
    size_type edge_count=num_edges();

    Node* new_a=new Node(this,a.index());
    Node* new_b=new Node(this,b.index());
    Edge_attr* new_edge=new Edge_attr(new_a,new_b);
    Edge new_e=Edge(this,new_a,new_b,edge_count);

    // Assign value to the internal edge attribute
    new_edge->Assign_val(new_e.length());
    edge_vec.push_back(new_edge);
    /* insert neighboring node id to edge idx pairs into each of a,
       and b's incident edge unordered maps. */
    node_2_edge[a.id()].insert(std::make_pair(b.id(),edge_count));
    node_2_edge[b.id()].insert(std::make_pair(a.id(),edge_count));
    e_idx_uid.push_back(edge_vec.size()-1);// insert id of edge into index to id mapping

    return new_e;
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    // free up dynamically allocated memory first
    for (auto& v: node_vec) {
      if (v!=nullptr) {
        delete v;
      }
    }
    for (auto& v: edge_vec) {
      if (v!=nullptr) {
        // clear up node pointers stored as internal edge attributes
        if(v->n1!=nullptr)
          delete v->n1;
        if(v->n2!=nullptr)
          delete v->n2;
        delete v;
      }
    }
    // clear up containers
    node_vec.clear();
    node_2_edge.clear();
    edge_vec.clear();
    n_idx_uid.clear();
    e_idx_uid.clear();

  }


  /**
   *  @brief Remove input node from the graph and all its incident edges.
   *
   *  @param[in] n The node to be removed.
   *  @return 1 if node is removed successfully or 0 if input node is not
   *     an active node in the graph(has been removed before or does not
   *     exist in the graph).
   *
   *  @pre n is a valid active node existing in the graph
   *  @post case I: has_node(@a n) == false
   *     all attributes of the graph remain the same
   *  @post case II: has_node(@a n) == true
   *       (i). graph new size() = old size() - 1.
   *      (ii). graph new num_edges() = old num_edges - n.degree().
   *     (iii). n becomes an inactive node, i.e., its id is erased from n_idx_uid.
   *      (iv). old last active node's index gets swapped with n's index.
   *       (v). edges incident to old last active node update this node's
   *            information (new node index) inside the internal edge attribute.
   *
   *  The complexity of the remove node algorithm is O(num_nodes()).
   */
  size_type remove_node(const Node& n) {
    // Only attempt edge removals and node index updates if n eixsts and active
    if(has_node(n)) {

      /* Clear up edges incident to n
         At most O(num_nodes()) iterations of remove_edge calls which each takes
         O(1) if this node has edges between every other nodes. */
      for(auto iter=n.edge_begin();iter!=n.edge_end();) {
        auto e=*iter;
        ++iter;
        remove_edge(e);
      }

      // swap last node's index with n's index in n_idx_uid mapping
      size_type last=n_idx_uid.back();
      n_idx_uid[n.index()]=last;

      /* Only try to update edges incident to previous last node for new
         node index information.
         At most O(num_nodes()) iterations to update node pointers, which
         takes O(1) if the last node has all other nodes as neighbors. */
      if(n.index()<(size()-1)) {
        for(auto iter=n.edge_begin(); iter!=node(n.index()).edge_end();++iter) {
          auto e=(*iter);
          if(edge_vec[e.id()]->n1->index()==(n_idx_uid.size()-1)) {
            delete edge_vec[e.id()]->n1;
            edge_vec[e.id()]->n1=new Node(this,n.index());
          }else{
            delete edge_vec[e.id()]->n2;
            edge_vec[e.id()]->n2=new Node(this, n.index());
          }
        }
      }
      n_idx_uid.pop_back();
      return 1;
    } else{
      return 0; // return 0 if n does not exist in the graph
    }
  }


  /**
   *  @brief Remove node pointed to by the input node iterator.
   *     Wraps around the remove_node(const Node& n) function by
   *     passing the dereferenced node iterator into the function.
   *
   *  @param[in] n_iter node iterator for the graph
   *  @return A node iterator pointing to end if removing the end
   *     node iterator or node is not removed successfully; otherwise,
   *     return a node iterator pointing to the beginning.
   *
   *  @pre n_iter a valid iterator spawned from the graph
   *  @post case I: @a n_iter==node_end() or remove_node(*@a n_iter) == 0
   *     all attributes of the graph remain the same
   *  @post case II: Else
   *       (i). graph new size() = old size() - 1.
   *      (ii). graph new num_edges() = old num_edges - n.degree().
   *     (iii). node pointed by the iterator becomes an inactive node,
   *            i.e., its id is erased from n_idx_uid.
   *      (iv). old last active node's index gets swapped with pointed
   *            node's index.
   *       (v). edges incident to old last active node update this node's
   *            information (new node index) inside the internal edge
   *            attribute.
   *
   *  The complexity of the remove node algorithm is O(num_nodes())).
   */
  node_iterator remove_node(node_iterator n_iter) {
    if(n_iter==node_end())
      return n_iter; // do nothing and return end iterator if removing from end
    size_type result=remove_node(*n_iter);
    if(result==0)
      return node_end();
    return node_begin();// return begin iterator o.w.
  }


  /**
   *  @brief Remove edge incident to two input nodes.
   *
   *  @param[in] n1, n2 two nodes of the graph
   *  @return 1 if the edge is removed successfully; 0 otherwise (i.e.,
   *     the edge is no longer active or does not exist in the graph.
   *
   *  @pre n1, n2 are valid active nodes of the graph
   *  @post case I: get_edge(@a n1, @a n2).first != nullptr &&
   *   second < num_edges()
   *     all attributes of the graph remain the same
   *  @post case II: Else
   *       (i). graph new num_edges() = old num_edges() - 1
   *      (ii). graph has_edge(@a n1, @a n2) == false &&
   *            graph has_edge(@a n2, @a n1) == false
   *     (iii). last active edge's index gets swapped with the to be removed
   *            edge's
   *      (iv). last active edge's new index gets updated in incident maps
   *            for input nodes
   *
   *  The complexity for the remove edge algorithm is O(1).
   */
  size_type remove_edge(const Node& n1, const Node& n2) {
    // First retrieve the edge with input nodes as two ends
    std::pair<Edge_attr*,size_type> temp=get_edge(n1,n2);
    // Remove information about this edge if it exists in the graph
    if(temp.first!=nullptr&&temp.second<num_edges()) {
      // erase the edge from input nodes' incident maps
      node_2_edge[n1.id()].erase(n2.id());
      node_2_edge[n2.id()].erase(n1.id());

      size_type last=e_idx_uid.back(); // uid of last edge
      Edge last_edge=edge(num_edges()-1); // num_edges()-1 index of last edge

      // if to be removed edge is not the last active edge in the graph
      if(last_edge.edge_idx!=temp.second){
        e_idx_uid[temp.second]=last; // update edge id
        // update new index for last active edge in incident maps from
        // its incident nodes
        node_2_edge[last_edge.node1().id()][last_edge.node2().id()]=temp.second;
        node_2_edge[last_edge.node2().id()][last_edge.node1().id()]=temp.second;

      }
      e_idx_uid.pop_back();
      return 1;
    }else {
      return 0; // return 0 to indicate erasing a non-existing edge
    }
  }

  /**
   *  @brief Remove the input edge from the graph.
   *
   *  @param[in] e an edge in the graph
   *  @return the result from remove_edge(@a e.node1(), @a e.node2())
   *
   *  @pre e is a valid active edge in the graph
   *  @post case I: remove_edge(@a e.node1(), @a e.node2()) == 0
   *     all attributes of the graph remain the same
   *  @post case II: remove_edge(@a e.node1(), @a e.node2()) == 1
   *       (i). graph new num_edges() = old num_edges() - 1
   *      (ii). graph has_edge(@a e.node1(), @a e.node2()) == false &&
   *            graph has_edge(@a e.node2(), @a e.node1()) == false
   *
   *  The complexity for the remove edge algorithm is O(1).
   */
  size_type remove_edge(const Edge& e) {
    return remove_edge(e.node1(),e.node2());
  }

  /**
   *  @brief Remove the edge pointed by the input edge iterator
   *
   *  @param[in] e_iter an edge iterator from the graph
   *  @return the end edge iterator if @a e_iter pointing to the
   *     end or *@a e_iter is not an valid active edge in the graph
   *
   *  @pre e_iter a valid iterator spawned from the graph
   *  @post case I: @a n_iter==node_end() or remove_edge(*@a e_iter) == 0
   *     all attributes of the graph remain the same
   *  @post case II: Else
   *       (i). graph new num_edges() = old num_edges() - 1
   *      (ii). graph has_edge(@a e.node1(), @a e.node2()) == false &&
   *            graph has_edge(@a e.node2(), @a e.node1()) == false
   *
   *  The complexity for the remove edge algorithm is O(1).   
   */
  edge_iterator remove_edge(edge_iterator e_iter) {
    if(e_iter==edge_end())
      return e_iter;
    size_type result=remove_edge(*e_iter);
    if(result==0)
      return edge_end();
    else
      return edge_begin();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    /** Deference a NodeIterator. */
    Node operator*() const{
      assert(G_!=nullptr); // should not deference an invalid node iterator
      return G_->node(iter_idx);
    }

    /** Forward increment the NodeIterator by incrementing iterator index. */
    NodeIterator& operator++(){
      ++iter_idx;
      return *this;
    }

    /** Equivalent node iterators should belong to the same graph and
        point to the node with the same node index. */
    bool operator==(const NodeIterator& iter) const{
      if(G_==iter.G_ && iter_idx==iter.iter_idx){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    size_type iter_idx=size_type();
    Graph* G_=nullptr;
    // Valid private node iterator constructor
    NodeIterator(const Graph* G, size_type idx){
      G_=const_cast<Graph*>(G);
      iter_idx=idx;
    }
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // Node iterator should begin to point to node indexed 0
  node_iterator node_begin() const{
    return  NodeIterator(this,0);
  }

  /** Update iterator index to be one past the total number
      of nodes in the graph. */
  node_iterator node_end() const{
    return NodeIterator(this,size());
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
    /** The edge obtained by deferencing the incident iterator should have the
        the node that spawns the iterator as node1. */
    Edge operator*() const {
      /* Retrieve the existing edge incident to the node where the iterator
         points to. */
      // should not deref invalid incident iterator
      assert(G_!=nullptr);
      //  should not deref the end iterator
      size_type e_idx=outedge_iter->second;
      Edge_attr* temp=G_->edge_vec[G_->e_idx_uid[e_idx]];
      /* If node1 of the existing edge is not the node that spawns the
         iterator, return an edge with flipped node1 and node2. */
      if(temp->n1->index()==n_idx){
        return G_->add_edge(*(temp->n1),*(temp->n2));
      }else{
        Edge temp2=G_->add_edge(*(temp->n2),*(temp->n1));
        return temp2;
      }
    }

    /* Increment the incident iterator by moving onto the next edge
       incident to the node pointed to by the iterator. */
    IncidentIterator& operator++(){
      ++outedge_iter;
      return *this;
    }
    
    /* Equivalent incident iterators should belong to the same graph,
       point to the same node, and currently examine the same edge incident
       to that node. */
    bool operator==(const IncidentIterator& iter) const {
      if(G_==iter.G_ && n_idx==iter.n_idx &&
outedge_iter==iter.outedge_iter){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph* G_=nullptr;
    size_type n_idx=size_type();
    /** index of edge w.r.t out_edge vector from node attribute that
      the IncidentIterator is currently point to. */
    u_map::iterator outedge_iter;
    

    // Private valid incident iterator constructor
    IncidentIterator(const Graph* G, size_type idx_,
u_map::iterator iter){
      G_=const_cast<Graph*>(G);
      n_idx=idx_;
      outedge_iter=iter;
    }
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
    // Dereference operator that returns the edge pointed to by the iterator
    Edge operator*() const {
      assert(G_!=nullptr); // should not dereference an invalid edge iterator
      return G_->edge(edge_idx);
    }

    /** Forward increment operation on EdgeIterator, move onto the next edge
        in the graph with incremented edge index. */
    EdgeIterator& operator++() {
      edge_idx++;
      return *this;
    }

    /** Equivalent edge iterators should belong to the same graph and point
        to the edge with the same edge index. */
    bool operator==(const EdgeIterator& iter) const {
      if(G_==iter.G_ && edge_idx==iter.edge_idx){
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph* G_=nullptr;
    size_type edge_idx=size_type();

    // Valid private edge iterator constructor
    EdgeIterator(const Graph* G, size_type idx) {
      G_=const_cast<Graph*>(G);
      edge_idx=idx;
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // Edge iterator begins to point to the first edge in the graph
  edge_iterator edge_begin() const {
    return EdgeIterator(this, 0);
  }

  // Edge iterator ends one past over all edges in the graph
  edge_iterator edge_end() const {
    return EdgeIterator(this, num_edges());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  // Internal node structure
  struct Node_attr {
      Point node_pos;
      node_value_type node_val;

      Node_attr(const Point& P, const node_value_type& val) {
        node_pos=Point(P.x, P.y, P.z);
        node_val=val;
      }
  };

  // Internal edge structure
  struct Edge_attr {
    const Node* n1; // pointer to node1 of the edge
    const Node* n2; // pointer to node2 of the edge
    E edge_val=E(); // edge value

    // constructor for internal edge attribute structure
    Edge_attr(const Node* a, const Node* b) {
      n1=const_cast<Node*>(a);
      n2=const_cast<Node*>(b);
    }

    // Assign edge value by the length of an edge based on Euclidean dist
    // btw node1 and node2
    void Assign_val(double l) {
      edge_val=l;
    }

  };

  /* vector of internal Node pointers, vector index is consistent with
     node idx. */
  std::vector<Node_attr*> node_vec;

  /* map with the pair of node uids as key and store edge index
     uid should always stay the same whenever removing nodes/edges. */
  std::unordered_map<size_type,u_map> node_2_edge;

  /* vector of internal Edge pointers, vector index is consistent with
     edge idx. */
  std::vector<Edge_attr*> edge_vec;

  /* vector whose each index stores the node id of the node that lives
     at that index in the node vector. */
  std::vector<size_type> n_idx_uid;

  /* vector whose each index stores the edge id of the edge that lives
     at that index in the edge vector. */
  std::vector<size_type> e_idx_uid;

  /** 
   * @brief Get an internal edge attribute pointer with edge index pair for
   *    edge incident to input nodes.
   *
   * @pre @a a and @a b are valid nodes of the graph
   * @return internal edge attribute pointer with edge index pair if
   *    for some @a i, edge(@a i) connects @a a and @a b or nullptr with
   *    deault size_type vlaue if no such edge exists. 
   */
  std::pair<Edge_attr*,size_type> get_edge(const Node& a, const Node& b) const {
    // make sure input valid nodes in the same graph
    assert(a.G_==this&&b.G_==this);
    // search whether node id keys exist
    auto search=node_2_edge.find(a.id())->second.find(b.id());
    if(search!=node_2_edge.find(a.id())->second.end()) {
      return std::make_pair(edge_vec[e_idx_uid[search->second]],search->second);
    }
    return std::make_pair(nullptr,size_type());
  }

};

#endif // CME212_GRAPH_HPP

