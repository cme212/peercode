#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

/*  Added code */
#include <iostream>

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <string>

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
 private:

  // HW0: YOUR CODE HERE
  // Use this space for declarations of important internal types you need
  // later in the Graph's definition..
  // (As with all the "YOUR CODE HERE" markings, you may not actually NEED
  // code here. Just use the space if you need it.)
  
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  /*********Added code HW1***********/
  using node_value_type = V;
  
  using edge_value_type = E;

  /** Type of this graph. */
  using graph_type = Graph<V, E>;

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
    // HW0: YOUR CODE HERE
    
    /********* added code ***********/
    // listinitializing Graph object
    : Nodes_p(), edge_n1(), edge_n2(),node_m(), size_n(0), size_e(0), \
    val_type(),edge_val_type(), adj_n(), uid_edge(), size_e_r(), uid_node(), size_n_r(), real_uid(){
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
    Node() {
      // HW0: YOUR CODE HERE
    }

     Point& position() {
      /************ Added code HW2************/
      int r_id = (graph_p->uid_node)[node_id];
      return graph_p->Nodes_p[r_id];
    }

    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE

      /************ Added code ************/
      int r_id = (graph_p->uid_node)[node_id];
      return graph_p->Nodes_p[r_id];
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE

      /************ Added code ************/
      return node_id;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value();
    // const node_value_type& value() const;
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;
    
    /***********Added code HW1************/
    // return value
    node_value_type& value(){
        int r_id = (graph_p->uid_node)[node_id];
        return graph_p->val_type[r_id];
    }
    // return value with const
    const node_value_type& value() const{
        int r_id = (graph_p->uid_node)[node_id];
        return graph_p->val_type[r_id];
    }
    // calculate degree
    size_type degree() const{
        int r_id = (graph_p->uid_node)[node_id];
        return (graph_p->adj_n[r_id]).size();
    }
    // incident iterator
    incident_iterator edge_begin() const{          

        return incident_iterator(graph_p, node_id, 0);
    }
    incident_iterator edge_end() const{
        int r_id = (graph_p->uid_node)[node_id];
        size_type adj_e_size = (graph_p->adj_n[r_id]).size();
        return incident_iterator(graph_p, node_id, adj_e_size);
    }
    
    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */

    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE

      /************ Added code ************/
      // check equality of graph ptr and node index
      
      if(n.graph_p == graph_p && n.index() == node_id){
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

      /************ Added code ************/
      if(graph_p != n.graph_p){
          return (graph_p < n.graph_p);
      }
      if (node_id < n.index()){
          return true;
      }        
      return false;
    }
    
   
   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    

    /******** Added code **********/
    graph_type* graph_p;
    size_type node_id;
    
    // private constructor
    Node(const graph_type* g, int uid)
        // list initializing node
        :graph_p(const_cast<graph_type*>(g)), node_id(uid){
    }
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {

    /********Added code*****************/
    return size_n;
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
  Node add_node(const Point& position, const node_value_type& V_T = node_value_type()) {
    // HW0: YOUR CODE HERE

    /***************Added code*******************/
    // push new position to list of node
    Nodes_p.push_back(position);
    //push new value type
    val_type.push_back(V_T);
    // create new node
    Node new_node = Node(this, size_n);
    // push to shadow id
    uid_node.push_back(size_n_r);
    // real to uid
    real_uid.insert({size_n_r, size_n});
    // increment real size
    size_n_r++;
    size_n += 1;
    return new_node;
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE

    /************ Added code ************/
    // check wheter index is in the range
    if(n.index() < size_n){
        return true;
    }
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

    /************ Added code ************/
    Node node_i = Node(this, i);
    return node_i;        // Invalid node
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
    /**********Added Code HW2***********/
    edge_value_type& value(){
        int r_id = (graph_e->uid_edge)[edge_id];
        return graph_e->edge_val_type[r_id];
    }
    // return value with const
    const edge_value_type& value() const{
        int r_id = (graph_e->uid_edge)[edge_id];
        return graph_e->edge_val_type[r_id];
    }

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      
      /************ Added code ************/
      int r_id = (graph_e->uid_edge)[edge_id];
      int tmp = (graph_e->edge_n1)[r_id];
      int uid = (graph_e->real_uid)[tmp];
      int tmp1 = (graph_e->edge_n2)[r_id];
      int uid1 = (graph_e->real_uid)[tmp1];
      if(node_pivot >= 0){
          if((int) uid  == node_pivot){
              return Node( graph_e, uid);
          } else {
              return Node( graph_e, uid1);
          }
      }          
      return Node( graph_e, uid);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE

      /************ Added code ************/
      int r_id = (graph_e->uid_edge)[edge_id];
      int tmp = (graph_e->edge_n2)[r_id];
      int uid = (graph_e->real_uid)[tmp];
      int tmp1 = (graph_e->edge_n1)[r_id];
      int uid1 = (graph_e->real_uid)[tmp1];

      if(node_pivot >= 0){
          if((int) uid  == node_pivot){
              return Node( graph_e, uid1);
          } else {
              return Node( graph_e, uid);
          }
      }
      return Node( graph_e, uid);
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
      
      /************ Added code ************/
      // check two nodes that define the edge
      // take into account undirected nature of graph
      int r_id = (graph_e->uid_edge)[edge_id];
      
      if( ( e.node1() == Node(graph_e, (graph_e->edge_n2)[r_id]) && e.node2() == Node(graph_e, (graph_e->edge_n1)[r_id]) ) || \
      ( e.node1() == Node(graph_e, (graph_e->edge_n1)[r_id]) && e.node2() == Node(graph_e, (graph_e->edge_n2)[r_id]))){
          if(e.graph_e == this->graph_e){
              return true;
          }
      }
      return false;
    }

    /*************Added code****************/


    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //HW0: YOUR CODE HERE

      /*************Added code****************/
      if(graph_e != e.graph_e){
          return (graph_e < e.graph_e);
      }
      if((int) edge_id < (int) e.edge_id){
          return true;
      }
      return false;
    }
    /**********Added code HW2***********/
    double length() const{
        int r_id = (graph_e->uid_edge)[edge_id];
        int real_n1 = (graph_e->edge_n1)[r_id];
        int real_n2 = (graph_e->edge_n2)[r_id];
        Point a = (graph_e->Nodes_p)[real_n1];
        Point b = (graph_e->Nodes_p)[real_n2];
        return norm(a - b);
    }
    
    private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects

    /***********Added Codes*************/
    // pointer to graph object
    graph_type* graph_e;
    // id of edge
    int edge_id;
    // pivot node for incident iterator
    int node_pivot;

    //private constructor
    Edge(const graph_type* g, size_type e_id, int pivot = -1)
        :graph_e(const_cast<graph_type*>(g)), edge_id(e_id), node_pivot(pivot){       
    }
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE

    /***********Added Codes*************/
    return size_e;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE

    /***********Added Codes*************/
    if(i < size_e){
        Edge new_edge = Edge(this, i);
        return new_edge;
    }
    return Edge();        // Invalid Edge
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */

  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE

    /***********Added Codes*************/
    // check precondition

    //std::cout << (*(edge_n1[0])).index() <<" " << (*(edge_n2[0])).index();
    if(Graph::has_node(a) && Graph::has_node(b)){
        std::string key;
        
        // construct key in the edge map

        int r_ia = uid_node[a.index()];
        int r_ib = uid_node[b.index()];
        if(r_ia < r_ib){
            key = std::to_string(r_ia) \
            + std::string("+") + std::to_string(r_ib);
        } else {
            key = std::to_string(r_ib) \
            + std::string("+") + std::to_string(r_ia);
        }
        // check whether the key is in the map
        if(node_m.count(key) > 0){
            return true;
        }
    }
    std::string key;

    // construct key in the edge map

    int r_ia = uid_node[a.index()];
    int r_ib = uid_node[b.index()];
    if(r_ia < r_ib){
        key = std::to_string(r_ia) \
        + std::string("+") + std::to_string(r_ib);
    } else {
        key = std::to_string(r_ib) \
        + std::string("+") + std::to_string(r_ia);
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
  Edge add_edge(const Node& a, const Node& b, const edge_value_type& E_T = edge_value_type()) {
    // HW0: YOUR CODE HERE

    /***************Added Code****************/
    // construct the key by concatenating indices of
    // nodes separated by "+"
    
    std::string key;
    int r_ia = uid_node[a.index()];
    int r_ib = uid_node[b.index()];
    if(r_ia < r_ib){
        key = std::to_string(r_ia) \
        + std::string("+") + std::to_string(r_ib);
    } else {
        key = std::to_string(r_ib) \
        + std::string("+") + std::to_string(r_ia);
    }

    if(Graph::has_node(a) && Graph::has_node(b) && (a == b) == false){
        // edge not exist case
        if(Graph::has_edge(a,b) == false){
            edge_n1.push_back(uid_node[a.index()]);
            edge_n2.push_back(uid_node[b.index()]);
            Edge new_edge = Edge(this, size_e);
            // add to map        
            node_m.insert({key, size_e});
            // add to edges to adjacent edges 
            /************ Added Code HW1*******/ 
            if(adj_n.count(uid_node[a.index()]) > 0){
                adj_n[uid_node[a.index()]].push_back(size_e);
            } else {
                std::vector<unsigned> tmp{size_e_r};
                adj_n.insert({uid_node[a.index()],tmp});
            }
            if(adj_n.count(uid_node[b.index()]) > 0){
                adj_n[uid_node[b.index()]].push_back(size_e);
            } else {
                std::vector<unsigned> tmp{size_e_r};
                adj_n.insert({uid_node[b.index()],tmp});
            }
            /**************Added code HW2******************/
            uid_edge.push_back(size_e_r);
            size_e_r++;
            size_e++;
            edge_val_type.push_back(E_T);
            return new_edge;
        } else {
        //edge exists case
            size_type index = node_m[key];
            return Edge(this, index);
        }
    }
    return Edge();        // Invalid Edge
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    /*************Added code****************/
    size_e = 0;
    size_n = 0;
    Nodes_p.clear();
    edge_n1.clear();
    edge_n2.clear();
    node_m.clear();
    val_type.clear();
    edge_val_type.clear();
    adj_n.clear();
    uid_edge.clear(); 
    size_e_r = 0; 
    uid_node.clear();
    size_n_r = 0; 
    real_uid.clear();
  }

  /*********Added code HW2*********/
  std::string key_cal(const Node& a, const Node& b){
    std::string key;
    int r_ia = uid_node[a.index()];
    int r_ib = uid_node[b.index()];
    if(r_ia < r_ib){
        key = std::to_string(r_ia) \
        + std::string("+") + std::to_string(r_ib);
    } else {
        key = std::to_string(r_ib) \
        + std::string("+") + std::to_string(r_ia);
    }
    return key;
  }

  
  /** Remove the edge by two nodes
   *  @post g.node(i).index() == i for all i in 0 to g.num_nodes()
   *  @post g.edge(i).index() == i for all i in 0 to g.num_edges()
   *  @post g.num_edges() (after) = g.num_edges() (before) - 1 if there 
   *  is a node between a and b
   *
   *  The iterator on node, edge, and incident should also function after
   *  this removal
   *  
   *  Runtime: This function only takes O(1) as we use swap and
   *  pop technique on uid array to mask the real id of edges
   */
  size_type remove_edge(const Node& a, const Node& b){
    if(has_edge(a,b) == false){
        return 0;
    }
    std::string key = key_cal(a, b);
    if(size_e == 1){
        node_m.erase(key);
        uid_edge.pop_back();
        size_e--;
        adj_n.erase(uid_node[a.index()]);
        adj_n.erase(uid_node[b.index()]);
        return 1;
    }
    // get edge id
    int e_id = node_m[key];
    int n1_i = uid_node[a.index()];
    int n2_i = uid_node[b.index()];
    int back_id = uid_edge.back();
    Node a_b = Node(this, real_uid[edge_n1[back_id]]);
    Node b_b = Node(this, real_uid[edge_n2[back_id]]);
    std::string key_b = key_cal(a_b, b_b);
    // swap and pop shadow id
    int tmp = uid_edge[e_id];
    int e_id_b = uid_edge.size() - 1;
    uid_edge[e_id] = uid_edge.back();
    uid_edge.back() = tmp;
    uid_edge.pop_back();
    // decrement size_e
    size_e--;
    // change node_m
    if(node_m.count(key_b) > 0){
        node_m[key_b] = e_id;
    }
    node_m.erase(key);
    // change adj_n
    // remove old
    std::vector<size_type> adj_t1 = adj_n[n1_i];
    std::vector<size_type>::iterator pos1 = std::find(adj_t1.begin(), adj_t1.end(), e_id);
    std::vector<size_type> adj_t2 = adj_n[n2_i];
    std::vector<size_type>::iterator pos2 = std::find(adj_t2.begin(), adj_t2.end(), e_id);
    if (pos1 != adj_t1.end() && pos2 != adj_t2.end()){
        adj_t1.erase(pos1);
        adj_t2.erase(pos2);
    }

    adj_n[n1_i] = adj_t1;
    adj_n[n2_i] = adj_t2;

    // change new loc
    std::vector<size_type> adj_t1n = adj_n[uid_node[a_b.index()]];
    std::vector<size_type>::iterator pos1_n = std::find(adj_t1n.begin(), adj_t1n.end(), e_id_b);
    std::vector<size_type> adj_t2n = adj_n[uid_node[b_b.index()]];
    std::vector<size_type>::iterator pos2_n = std::find(adj_t2n.begin(), adj_t2n.end(), e_id_b);
    if (pos1_n != adj_t1n.end() && pos2_n != adj_t2n.end()){
        *pos1_n = e_id;
        *pos2_n = e_id;
    }
    adj_n[uid_node[a_b.index()]] = adj_t1n;
    adj_n[uid_node[b_b.index()]] = adj_t2n;
    return 1;
  }


  /** Remove the edge by iterators
   *  @post g.node(i).index() == i for all i in 0 to g.num_nodes()
   *  @post g.edge(i).index() == i for all i in 0 to g.num_edges()
   *  @post g.num_edges() (after) = g.num_edges() (before) - 1 
   *  
   *
   *  The iterators on node, edge, and incident as well as other
   *  standard functions of edges and nodes from previous
   *  should also function after
   *  this removal
   *
   *  This function will return edge_iterator to next removable
   *  edge.
   *
   *  Runtime: This function only takes O(1) as we use swap and
   *  pop technique on uid array to mask the real id of edges
   */  
  edge_iterator remove_edge(edge_iterator e_it){
    Edge e = *e_it;
    Node a = e.node1();
    Node b = e.node2();

    if(has_edge(a,b) == false){
        return (*this).edge_begin();
    }
    std::string key = key_cal(a, b);
    if(size_e == 1){
        node_m.erase(key);
        uid_edge.pop_back();
        size_e--;
        adj_n.erase(uid_node[a.index()]);
        adj_n.erase(uid_node[b.index()]);
        return (*this).edge_begin();
    }
    // get edge id
    int e_id = node_m[key];
    int n1_i = uid_node[a.index()];
    int n2_i = uid_node[b.index()];
    int back_id = uid_edge.back();
    Node a_b = Node(this, real_uid[edge_n1[back_id]]);
    Node b_b = Node(this, real_uid[edge_n2[back_id]]);
    std::string key_b = key_cal(a_b, b_b);
    // swap and pop shadow id
    int tmp = uid_edge[e_id];
    int e_id_b = uid_edge.size() - 1;
    uid_edge[e_id] = uid_edge.back();
    uid_edge.back() = tmp;
    uid_edge.pop_back();
    // decrement size_e
    size_e--;

    // change node_m
    if(node_m.count(key_b) > 0){
        node_m[key_b] = e_id;
    }
    node_m.erase(key);
    // change adj_n
    // remove old
    std::vector<size_type> adj_t1 = adj_n[n1_i];
    std::vector<size_type>::iterator pos1 = std::find(adj_t1.begin(), adj_t1.end(), e_id);
    std::vector<size_type> adj_t2 = adj_n[n2_i];
    std::vector<size_type>::iterator pos2 = std::find(adj_t2.begin(), adj_t2.end(), e_id);
    if (pos1 != adj_t1.end() && pos2 != adj_t2.end()){
        adj_t1.erase(pos1);
        adj_t2.erase(pos2);
    }

    adj_n[n1_i] = adj_t1;
    adj_n[n2_i] = adj_t2;

    // change new loc
    std::vector<size_type> adj_t1n = adj_n[uid_node[a_b.index()]];
    std::vector<size_type>::iterator pos1_n = std::find(adj_t1n.begin(), adj_t1n.end(), e_id_b);
    std::vector<size_type> adj_t2n = adj_n[uid_node[b_b.index()]];
    std::vector<size_type>::iterator pos2_n = std::find(adj_t2n.begin(), adj_t2n.end(), e_id_b);
    if (pos1_n != adj_t1n.end() && pos2_n != adj_t2n.end()){
        *pos1_n = e_id;
        *pos2_n = e_id;
    }
    adj_n[uid_node[a_b.index()]] = adj_t1n;
    adj_n[uid_node[b_b.index()]] = adj_t2n;
    return (*this).edge_begin();
  }
  

  /** Remove the edge by edge
   *  @post g.node(i).index() == i for all i in 0 to g.num_nodes()
   *  @post g.edge(i).index() == i for all i in 0 to g.num_edges()
   *  @post g.num_edges() (after) = g.num_edges() (before) - 1 if e is
   *  a valid edge
   *
   *  The iterators on node, edge, and incident as well as other
   *  standard functions of edges and nodes from previous
   *  should also function after
   *  this removal
   *
   *  Runtime: This function only takes O(1) as we use swap and
   *  pop technique on uid array to mask the real id of edges
   */
  size_type remove_edge(const Edge& e){
    Node a = e.node1();
    Node b = e.node2();

    if(has_edge(a,b) == false){
        return 0;
    }
    std::string key = key_cal(a, b);
    if(size_e == 1){
        node_m.erase(key);
        uid_edge.pop_back();
        size_e--;
        adj_n.erase(uid_node[a.index()]);
        adj_n.erase(uid_node[b.index()]);
        return 1;
    }
    // get edge id
    int e_id = node_m[key];
    int n1_i = uid_node[a.index()];
    int n2_i = uid_node[b.index()];
    int back_id = uid_edge.back();

    Node a_b = Node(this, real_uid[edge_n1[back_id]]);
    Node b_b = Node(this, real_uid[edge_n2[back_id]]);
    std::string key_b = key_cal(a_b, b_b);
    // swap and pop shadow id
    int tmp = uid_edge[e_id];
    int e_id_b = uid_edge.size() - 1;
    uid_edge[e_id] = uid_edge.back();
    uid_edge.back() = tmp;
    uid_edge.pop_back();
    // decrement size_e
    size_e--;
    // change node_m
    

    if(node_m.count(key_b) > 0){
        node_m[key_b] = e_id;
    }
    node_m.erase(key);
    // change adj_n
    // remove old
    std::vector<size_type> adj_t1 = adj_n[n1_i];
    std::vector<size_type>::iterator pos1 = std::find(adj_t1.begin(), adj_t1.end(), e_id);
    std::vector<size_type> adj_t2 = adj_n[n2_i];
    std::vector<size_type>::iterator pos2 = std::find(adj_t2.begin(), adj_t2.end(), e_id);
    
    if (pos1 != adj_t1.end() && pos2 != adj_t2.end()){
        adj_t1.erase(pos1);
        adj_t2.erase(pos2);
    }

    adj_n[n1_i] = adj_t1;
    adj_n[n2_i] = adj_t2;

    // change new loc
    std::vector<size_type> adj_t1n = adj_n[uid_node[a_b.index()]];
    std::vector<size_type>::iterator pos1_n = std::find(adj_t1n.begin(), adj_t1n.end(), e_id_b);
    std::vector<size_type> adj_t2n = adj_n[uid_node[b_b.index()]];
    std::vector<size_type>::iterator pos2_n = std::find(adj_t2n.begin(), adj_t2n.end(), e_id_b);
    if (pos1_n != adj_t1n.end() && pos2_n != adj_t2n.end()){
        *pos1_n = e_id;
        *pos2_n = e_id;
    }
    adj_n[uid_node[a_b.index()]] = adj_t1n;
    adj_n[uid_node[b_b.index()]] = adj_t2n;
    return 1;
  }
  
  /** Remove the Node by Node
   *  @post g.node(i).index() == i for all i in 0 to g.num_nodes()
   *  @post g.edge(i).index() == i for all i in 0 to g.num_edges()
   *  @post g.num_nodes() (after) = g.num_nodes() (before) - 1 if a is
   *  a valid node
   *
   *  The iterators on node, edge, and incident as well as other
   *  standard functions of edges and nodes from previous
   *  should also function after
   *  this removal
   *
   *  Runtime: This function only takes O(degree of a) 
   *  we need to iterate through all incident edges to a
   *  and remove them by remove_edge.
   */
  // remove node
  size_type remove_node(const Node& a){
      if(has_node(a) == false){
          return 0;
      }
      
      if(size_n == 1){
          uid_node.pop_back();
          adj_n.erase(uid_node[a.index()]);
          size_n--;
          return 1;
      }
      // remove edges incident to this node
      int cur_size = adj_n[uid_node[a.index()]].size();
      while(cur_size != 0){
          remove_edge(Edge(this, adj_n[uid_node[a.index()]].back()));
          cur_size = adj_n[uid_node[a.index()]].size();
      }
      // change real_uid
      // remove node
      int n_id = a.index();
      int tmp = uid_node[n_id];
      uid_node[n_id] = uid_node.back();
      uid_node.back() = tmp;
      uid_node.pop_back();

      // change real_uid
      real_uid[uid_node[n_id]] = n_id;
      real_uid.erase(tmp); 

      // size
      size_n--;
      // pop adj_n
      adj_n.erase(n_id);
      return 1;
  }

  /** Remove the Node by node iterator
   *  @post g.node(i).index() == i for all i in 0 to g.num_nodes()
   *  @post g.edge(i).index() == i for all i in 0 to g.num_edges()
   *  @post g.num_nodes() (after) = g.num_nodes() (before) - 1 if a is
   *  a valid node
   *
   *  The iterators on node, edge, and incident as well as other
   *  standard functions of edges and nodes from previous
   *  should also function after
   *  this removal
   *
   *  Runtime: This function only takes O(degree of a) 
   *  we need to iterate through all incident edges to a
   *  and remove them by remove_edge. 
   *
   * This function will return iterator to the new removable node
   * g.node_begin()                 
   */
  // remove node
  node_iterator remove_node(node_iterator& n){
      Node a = *n;
      if(has_node(a) == false){
          return (*this).node_begin();
      }

      if(size_n == 1){
          uid_node.pop_back();
          adj_n.erase(uid_node[a.index()]);
          size_n--;
          return (*this).node_begin();
      }
      // remove edges incident to this node
      int cur_size = adj_n[uid_node[a.index()]].size();
      while(cur_size != 0){
          remove_edge(Edge(this, adj_n[uid_node[a.index()]].back()));
          cur_size = adj_n[uid_node[a.index()]].size();
      }
      // change real_uid
      // remove node
      int n_id = a.index();
      int tmp = uid_node[n_id];
      uid_node[n_id] = uid_node.back();
      uid_node.back() = tmp;
      uid_node.pop_back();

      // change real_uid
      real_uid[uid_node[n_id]] = n_id;
      real_uid.erase(tmp);

           // size
      size_n--;
      // pop adj_n
      adj_n.erase(n_id);
      return (*this).node_begin();
  }




  //
  // Node Iterator
  //

  /** @class Graph::NodeIterator
   * @brief Iterator class for nodes. A forward iterator. */
  class NodeIterator : private totally_ordered<NodeIterator>{
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
    
    /********** Added code HW1 ***********/
    // Implement operator for node iterator
    Node operator*() const{
        return Node(graph_p, loc);    
    } 
    bool operator==(const NodeIterator& tmp) const{
        if(tmp.graph_p == graph_p and loc == tmp.loc){
            return true;
        }
        return false;
    }
    NodeIterator& operator++(){
        loc++;
        return *this;
    }
    
   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    /********** Added code HW1 ***********/
    graph_type* graph_p;
    size_type loc;
    NodeIterator(const graph_type* g, size_type l)
        // list initiliazing
        :graph_p(const_cast<graph_type*>(g)), loc(l){
    }
    
  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  /********** Added code HW1 ***********/
  // implement node_begin and node_end
  node_iterator node_begin() const{
      return node_iterator(this, 0);
  }
  node_iterator node_end() const{
      return node_iterator(this, size_n);
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator : private totally_ordered<IncidentIterator>{
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
    /************Added code hw1***********/
    // implement standard operator for incident iterator class
    Edge operator*() const{
        
        int real_id = (graph_p->uid_node)[node_inc];
        std::vector<size_type> tmp = graph_p->adj_n[real_id];
        return Edge(graph_p, tmp[edge_inc], node_inc);
    }
    IncidentIterator& operator++(){
        edge_inc++;
        return *this;
    }
    bool operator==(const IncidentIterator& k) const{
        if(k.edge_inc == edge_inc && k.node_inc == node_inc && k.graph_p == graph_p){
            return true;
        } 
        return false;       
    }
   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    /************Added code hw1***********/
    graph_type* graph_p;
    size_type node_inc;
    // order in adjacent edges (adj_n)
    size_type edge_inc;
    IncidentIterator(const graph_type* g, size_type n, size_type e)
        :graph_p(const_cast<graph_type*>(g)), node_inc(n), edge_inc(e){
    }
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
    EdgeIterator() {
    }

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    // standard operator for edge iterator
    Edge operator*() const{
        return Edge(graph_p, edge_loc);
    }
    EdgeIterator& operator++(){
        edge_loc++;
        return *this;
    }
    bool operator==(const EdgeIterator& k ) const{
        if(k.graph_p == graph_p && edge_loc == k.edge_loc){
            return true;
        }
        return false;
    }

   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    /********** Added code hw1 *******/
    graph_type* graph_p;
    size_type edge_loc;
    EdgeIterator(const graph_type* g, size_type l)
        :graph_p(const_cast<Graph*>(g)), edge_loc(l){
    }
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const
  edge_iterator edge_begin() const{
      return edge_iterator(this, 0);
  }
  edge_iterator edge_end() const{
      return edge_iterator(this, size_e);
  }
 
 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  /************Added code*************/
  // vector of nodes
  std::vector<Point> Nodes_p;
  // vector of node 1 for each edge
  std::vector<size_type> edge_n1;
  //std::vector<const Node*> edge_n1;
  // vector of node 2 for each edge
  std::vector<size_type> edge_n2;
  //std::vector<const Node*> edge_n2; 
  // edge map with two nodes as a key
  // and edge index as a value
  std::map<std::string, int> node_m;
  // size of node
  unsigned size_n;
  // size of edge
  unsigned size_e;
  // vector to keep V type of node
  std::vector<node_value_type> val_type;
  //vector to keep E type of edge
  std::vector<edge_value_type> edge_val_type;
  // map of node to adjacent edges
  std::map<size_type, std::vector<size_type>> adj_n;
  // shadow id for edge
  std::vector<int> uid_edge;
  // real_size
  unsigned size_e_r;
  // shadow for node
  std::vector<int> uid_node;
  // real size node
  unsigned size_n_r;
  // real node to uid
  std::map<size_type, size_type> real_uid;
};

#endif // CME212_GRAPH_HPP
