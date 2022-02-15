#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <cassert>
#include <utility> 
#include <iterator>

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

  /** Decleration of node value type, i.e. what it represents */
  using node_value_type = V;

  /** Predeclaration of Edge type. */
  class Edge;
  /** Synonym for Edge (following STL conventions). */
  using edge_type = Edge;

  /** Decleration of edge value type */
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
  Graph() {}

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
    /** Construct an invalid node. */
    Node(){}

    /** Return this node's position. 
     * Complexity O(1) since average lookup time is O(1) for unordered_map
     */
    const Point& position() const {
      size_type uid = graph->n_u2i[idx];
      return graph->nodes.at(uid).point;
    }
 
    /** Return a reference to position so it's modifiable */
    Point& position() {
      size_type uid = graph->n_u2i[idx];
      return graph->nodes.at(uid).point;
    }

    /** Allows Nodes value to be read and modified */
    node_value_type& value(){
      size_type uid = graph->n_u2i[idx];
      return graph->nodes.at(uid).value;
    }

    /** Allows Node's value to be read */
    const node_value_type& value() const {
      size_type uid = graph->n_u2i[idx];
      return graph->nodes.at(uid).value;
    }


    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      return idx;
    }

    /** Test whether this node and @a n are equal.
     * param[in] n, other node we're checking for equality
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      return (graph == n.graph) && (idx == n.idx);
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
      if(graph == n.graph){
          if(idx < n.idx){
              return true;
          } else {
              return false;
          }
      } else {
          //Allows us to compare graphs
          std::vector<graph_type*> pointers = {graph, n.graph};
          std::sort(pointers.begin(), pointers.end());
          if(pointers[0] == graph){
              return true;
          }
      }
      return false;
    }

    /** Returns the number of incident edges
     * Complexity O(1) since average constant time to lookup unordered_map
     */
    size_type degree() const {
        size_type uid = graph->n_u2i[idx];
        return (graph->adj_list[uid]).size();
    }
  
    /** Start of the incident iterator for this specific node
     * @pre uid is a valid key in adj_list, o/w throws key_error 
     * @return Incident_Iterator pointing to beginning of connections vector
     */
    incident_iterator edge_begin() const {
        size_type uid = graph->n_u2i[idx];
        std::unordered_set<size_type>::iterator itr = graph->adj_list.at(uid).begin();
        return IncidentIterator(itr, *this, graph);
    }
   
    /** End of the incident iterator
     * @pre uid is a valid key in adj_list, o/w throws key error
     * @return Incident_iterator pointing to one past end of connections vector
     */
    incident_iterator edge_end() const {
        size_type uid = graph->n_u2i[idx];
        std::unordered_set<size_type>::iterator itr = graph->adj_list.at(uid).end();
        return IncidentIterator(itr,*this, graph);
    }
    

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    graph_type* graph;
    size_type idx;

    /** private constructor for node class */
    Node(const graph_type* graph_, size_type idx_)
        : graph(const_cast<graph_type*>(graph_)), idx(idx_) {
    }
      

  };  


  /** Return the number of nodes in the graph.
   * Complexity: O(1).
   */
  size_type size() const {
    return n_u2i.size();
  }

  /** Synonym for size(). */
  size_type num_nodes() const {
    return n_u2i.size();
  }

  /** Add a node to the graph, returning the added node.
   * @param[in] position The new node's position
   * @param[in] value The new node's value 
   * @post new num_nodes() == old num_nodes() + 1
   * @post result_node.index() == old num_nodes()
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position,
		const node_value_type& value = node_value_type()) {
    internal_node i_node;
    i_node.uid = (size_type) nodes.size();
    i_node.point = position;
    i_node.value = value;
    i_node.idx = (size_type) n_u2i.size();

    nodes[i_node.uid] = i_node;
    n_u2i.push_back(i_node.uid);

    //We initialize an empty adj_list at this uid to relax IncidentItr standards
    adj_list[i_node.uid] = std::unordered_set<size_type>(); 
    return Node(this, i_node.idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   * param[in] n, node we're checking for equality
   * Just makes sure that index is within graph size and that
   * node is referring to this graph
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if(this == n.graph && n.idx < num_nodes()){
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
    assert(i < num_nodes());
    return Node(this, i);
  }

  /** Remove the node with index n.index()
   * @pre n_u2i.size() != 0
   * param[in,out] n, the node we want to remove from the graph
   * Complexity: O(d) where d is the degree of the node
   * 
   * @brief: Removes node from graph, and all edges incident to it
   *         thus invalidating any incident iterators which either 
   *         had n as its root_node, or was in it's adj_list u_set.
   *         Also invalidates node iterators and edge iterators if n
   *         had at least one edge.   
   */
  size_type remove_node(const Node& n){
      if(!has_node(n)){  //can't remove node that's not in the graph
          return 0; 
      }
 
      size_type succ = 1; //flag for success of remove
      //First get incident_node iterator to remove edges, first by removing edge, but also removing connection's
      if(adj_list.count(n_u2i[n.index()])){
          for(auto itr = n.edge_begin(); itr!= n.edge_end(); ){
              Edge e = *itr;
              succ = remove_edge(e);
              itr = n.edge_begin();
          }
      }
   
      //swap and pop the n_u2i vector, need to keep internal nodes index updated
      size_type swapped_idx = n_u2i.size()-1;
      nodes.at(n_u2i[swapped_idx]).idx = n.index();
      std::iter_swap(n_u2i.begin() + n.index(), n_u2i.end() -1);
      
      n_u2i.pop_back();

      return succ;
  }

  /** Remove the node that the node_iterator points to
   * @pre n_it points to a valid iterator
   * param[in,out] n_it , points to node we want to remove
   * Complexity: O(d) where d is the degree of the node
   * @brief: same as remove_node(const Node& n)
   */
  node_iterator remove_node(node_iterator n_it){ 
     Node n = *n_it;
     node_iterator retval = ++n_it;
     remove_node(n);
     return retval;
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
    Edge() { }

    /** Return a node of this Edge 
     * My implementation of edge only retains uid's of the nodes,
     * so need to use u_map of nodes to recover node index from node uid
     * Complexity: O(1) amortized
     */
    Node node1() const {
      internal_node i_node = graph->nodes.at(a_uid);
      return graph->node(i_node.idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      internal_node i_node = graph->nodes.at(b_uid);
      return graph->node(i_node.idx);    
    }

    size_type index() const {
        return idx;
    }

    /** Return the length of this Edge
     * @post e.length() > 0
     * @return double, length of the edge
     */
    double length() const {
	Point displacement = node1().position() - node2().position();
	return norm(displacement);
    }

    /** Return a reference to edge's value
     * @pre edge_uid is a valid index for edge_keys vec and thus for edges u_map
     * Complexity: O(1) amortized (indexing in vector and u_map is avg constant)
     */
    edge_value_type& value(){
        size_type uid = graph->e_u2i[idx];
    	std::pair<size_type,size_type> key = graph->edge_keys[uid];
        return graph->edges.at(key).value;
    }
  
    /** Return a const reference to edge's value that will not modify graph 
     * Same specifications as non_const version of graph
     */
    const edge_value_type& value() const {
        size_type uid = graph->e_u2i[idx];
        std::pair<size_type,size_type> key = graph->edge_keys[uid];
	return graph->edges.at(key).value;
    }

    /** Test whether this edge and @a e are equal.
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      return (node1() == e.node1()) && (node2() == e.node2());
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */
    bool operator<(const Edge& e) const {
      //Just piggyback off ordering of nodes
      if(node1() == e.node1()){
          return (node2() < e.node2());
      }
      return (node1() < e.node1());
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Edge(size_type a_uid_, size_type b_uid_, size_type idx_, \
         const graph_type* graph_)
        : a_uid(a_uid_), b_uid(b_uid_), idx(idx_), \
               graph(const_cast<graph_type*>(graph_)) {
    }
 
    size_type a_uid;
    size_type b_uid;
    size_type idx;
    graph_type* graph; 
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: O(1)
   */
  size_type num_edges() const {
    return e_u2i.size();
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   * @post {Node1.uid, Node2.uid} is in edge_keys
   * Complexity: O(1) amortized
   */
  Edge edge(size_type i) const {
    assert(i < e_u2i.size());
    std::pair<size_type, size_type> key = edge_keys[e_u2i[i]];
    return Edge(key.first, key.second, i, this);       
  }

  /** Helper function that given a tuple of node_uid's, creates the corresponding
   * key pair that enforces the order of {a,b} with a < b
   * @post returns {a,b} with a < b
   */ 
  std::pair<size_type, size_type> get_edge_key(size_type a_idx, size_type b_idx) const {
      size_type a_uid = n_u2i[a_idx];
      size_type b_uid = n_u2i[b_idx];
      if(a_uid < b_uid){
          return {a_uid, b_uid};
      }
      return {b_uid, a_uid};
  }
  

  /** Return the edge with node index @a a and node index @a b
   * @pre (@a a, @a b) elements of edge_keys
   * Complexity: Roughly O(1)
   * @return Edge with node1 index @a a and node2 index @a b
  */
  Edge edge(size_type a_idx, size_type b_idx) const {
    //correct ordering for edge_keys is {a,b} where a < b
    std::pair<size_type, size_type> key = get_edge_key(a_idx, b_idx);
    internal_edge old_edge = edges.at(key);
    return Edge(n_u2i[a_idx],n_u2i[b_idx], old_edge.idx, this);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: O(1) amortized
   */
  bool has_edge(const Node& a, const Node& b) const {
    std::pair<size_type,size_type> key = get_edge_key(a.index(), b.index());
    return edges.count(key);
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
   * Complexity: O(1);
   */
  Edge add_edge(const Node& a, const Node& b) {
    if(has_edge(a, b)){
        return edge(a.index(), b.index());
    } 
    //No edge exists, so we need to create one    
    //First deal with edge indices and interface
    internal_edge e;
    e.idx = e_u2i.size();
    e.uid = edges.size();
    e_u2i.push_back(e.uid);


    //Now add edge_keys
    size_type a_uid = n_u2i[a.index()];
    size_type b_uid = n_u2i[b.index()];
    std::pair<size_type, size_type> new_key = get_edge_key(a.index(), b.index());

    edge_keys.push_back(new_key);
    edges.insert({new_key,e});

    adj_list[a_uid].insert(b_uid);
    adj_list[b_uid].insert(a_uid);

    return edge(a.index(), b.index());   
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    nodes.clear();
    n_u2i.clear();
    edges.clear();
    edge_keys.clear();
    e_u2i.clear();
    adj_list.clear();
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
    using iterator_category = std::forward_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator() {};

    /** Iterate to next node in nodes vector.
      * @return Reference to node iterator
    */
    NodeIterator& operator++() {
        node_idx++;
	return *this;
    }

    /** Equality for NodeIterator class 
     * @return bool, whether or not NodeIter refers to same graph and vector<Node>::iterator
     */
    bool operator==(const NodeIterator& other_iter) const{
        return (node_idx == other_iter.node_idx) &&
		 (graph == other_iter.graph);
    }
    
    /** Dereference operator
     * @pre this != NodeIterator.end()
     * @return Node to which our operator points
    */
    Node operator*(){
	return graph->node(node_idx);
    }

   private:
    friend class Graph;
    //Need to only iterate over active nodes, so iterate over n_u2i and dereference
    size_type node_idx;
    const graph_type* graph;    

    //Private constructor
    NodeIterator(size_type n_u2i_idx, const graph_type* graph_) : 
			node_idx(n_u2i_idx), graph(graph_) {}

  };

  /** begin() function for NodeIterator
   * @return NodeIterator pointing to beginning of nodes vector
  */
  NodeIterator node_begin() const {
      return NodeIterator(0, this);
  }

  /** end() function for NodeIterator
   * @return NodeIterator pointing to one past end of nodes vector
  */
  NodeIterator node_end() const {
      return NodeIterator(n_u2i.size(), this);
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
    IncidentIterator() {}


    /**  Dereference operator for IncidentIterator
     * @pre this->incident_itr != incident_itr.end()
     * @return an Edge that connects root Node and another node
     * @post Edge.node1() == root_node, Edge.node2() == connected_node
     */ 
    Edge operator*() const{
        size_type uid = *incident_itr;
        internal_node incident_node = graph->nodes.at(uid);
        if(graph->has_edge(root_node, graph->node(incident_node.idx))){
            return graph->edge(root_node.idx, incident_node.idx);
        } 
        /** If edge not found, bug exists in program */
        std::cerr << "Couldn't find edge that should exist between nodes" << std::endl;
        exit(1);
    }

    /** Next operator for IncidentIterator
     * @return IncidentIterator pointing to next incident edge's node
     */
    IncidentIterator& operator++(){
        ++incident_itr; 
        return *this;
    }

    /** Equality operator for IncidentIterator
     * @return bool, whether or not two IncidentIterators have same root,
     *         graph, and point to same connection vector element
     */
    bool operator==(const IncidentIterator& other_iter) const{
	return (incident_itr == other_iter.incident_itr) &&
	       (root_node == other_iter.root_node) && 
	        (graph == other_iter.graph);
    }

   private:
    friend class Graph;

    /**Essentially, we have to iterate over the internal nodes member 
     * variable @a connections.This iterates over the incident nodes ids,
     * allowing us to access them in constant time
    */
    std::unordered_set<size_type>::iterator incident_itr;
    Node root_node;
    const Graph* graph;
  
    //Private constructor, takes in an iterator over a vector of size_type
    IncidentIterator(std::unordered_set<size_type>::iterator new_itr, 
		     Node root, const graph_type* graph_) : 
			incident_itr(new_itr), root_node(root), graph(graph_) {};
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
    EdgeIterator() {
    }

    /** Dereference operator for EdgeIterator
     * @pre this->edge_iter != edge_iter.end()
     * @return Edge with index @a i
    */
    Edge operator*() const {
        return graph->edge(edge_iter);
    }
  
    /** Next operator for EdgeIterator
     * @return Reference to Edge Iterator with edge_iter pointing to next elem
    */
    EdgeIterator& operator++() {
	++edge_iter;
	return *this;
    }
   
    /** Equality operator for EdgeIterator
     * @return bool, whether Edge iterator has same graph, points to same edge
    */
    bool operator==(const EdgeIterator& other_iter) const {
	return (edge_iter == other_iter.edge_iter) && (graph == other_iter.graph);
    }

   private:
    friend class Graph;
    size_type edge_iter;
    const graph_type* graph; 

    /** Private constructor for edge iterator */
    EdgeIterator(size_type iter, const graph_type* graph_) :
			edge_iter(iter), graph(graph_) {} ;    
  };


  /** Returns an iterator pointing to beginning of edge keys (set)*/
  edge_iterator edge_begin() const{
      return EdgeIterator(0, this);
  }
  
  /** Returns an iterator pointing to one past end of edge keys (set) */
  edge_iterator edge_end() const {
      return EdgeIterator(e_u2i.size(), this);
  }


  /** Removes the edge from the graph that the edge_iterator points to
   * @pre edge_iterator points to a valid edge
   * @pre map.erase(itr) returns an iterator pointing to next element in map
   * @pre nodes in the edge are valid keys in adj_list
   * param[in,out] e_it points to edge we want to remove from graph
   * note: edge_iterators will become invalidated after this operation as well as
   *       any incident_iterators that had either node in this edge as connections
   * Complexity: O(1) all indexing (or at) operations are constant, and the pop & swap
   *                 for the e_uid vector is a constant operation.
   * @return edge_iterator
   */
  edge_iterator remove_edge(edge_iterator e_it){
      //Need to get nodes so we can remove from connections list
      Edge e = *e_it;
      ++e_it;
      Node n1 = e.node1(); Node n2 = e.node2();
  
      //Now remove key pair from edges and pop out of e_u2i
      std::pair<size_type,size_type> key = get_edge_key(n1.index(), n2.index());
      adj_list.at(key.first).erase(key.second);
      adj_list.at(key.second).erase(key.first);

      //Need to keep internal edge idx updated
      size_type swapped_idx = e_u2i.size()-1;
      std::pair<size_type,size_type> swapped_keys = edge_keys.at(e_u2i[swapped_idx]); 
      edges.at(swapped_keys).idx = e.index();

      //now apply pop and swap to edge u2i vector
      std::iter_swap(e_u2i.begin() + e.index(), e_u2i.end() -1); 
      e_u2i.pop_back();
  
      edges.erase(key);
      return e_it;
  }

  /**
   * same specifications as remove_node(edge_iterator)
   * @return value indicates if removal was succesful
   */
  size_type remove_edge(const Node& n1, const Node& n2){
      //First check that edge exists between the two:
      if(!has_edge(n1,n2)){
          return 0;
      } 
      Edge e = edge(n1.index(), n2.index());
      return remove_edge(e);
  }

  /** 
   * same specfications as remove_node(edge_iterator)
   * @return value indicates if removal was succesful
   */
  size_type remove_edge(const Edge& e){
      size_type succ = 0;
      Node n1 = e.node1();
      Node n2 = e.node2();
      if(!has_edge(n1, n2)){
          return succ;
      } else {
          succ = 1;
           
          //erase the connection in the adj list
          std::pair<size_type,size_type> key = get_edge_key(n1.index(), n2.index());
          adj_list.at(key.first).erase(key.second);
          adj_list.at(key.second).erase(key.first);

          //Need to keep internal edge idx updated
          size_type swapped_idx = e_u2i.size()-1;  
          std::pair<size_type,size_type> swapped_keys = edge_keys.at(e_u2i[swapped_idx]); 
          edges.at(swapped_keys).idx = e.index();  
          
          //now apply pop and swap
          std::iter_swap(e_u2i.begin() + e.index(), e_u2i.end()-1);    
          e_u2i.pop_back();

          edges.erase(key);
          return succ;
      }
  }

 



 private:

  struct internal_node {
    size_type uid;
    size_type idx;
    Point point; 
    node_value_type value;
  };

  std::unordered_map<size_type, internal_node> nodes;
  std::vector<size_type> n_u2i;  //vector responsible for uid <-> idx transform

  //The pair containing the size_types will hold the indices
  //To use the u_map with pairs as keys, need to create hash_func

  struct pair_hash {
      size_type operator() (const std::pair<size_type,size_type> &pair) const {
          return (pair.second << 16) ^ pair.first;
      }
  };



  struct internal_edge {
    size_type uid; //needed to access internal_edge correctly
    size_type idx;
    edge_value_type value; 
  }; 

  std::unordered_map<size_type, std::unordered_set<size_type>> adj_list;
  std::unordered_map<std::pair<size_type,size_type>, internal_edge, pair_hash> edges;
  std::vector<std::pair<size_type, size_type>> edge_keys;
  std::vector<size_type> e_u2i; //vector responsible for uid <-> idx transform
};

#endif // CME212_GRAPH_HPP

