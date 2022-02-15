#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <tuple>
#include <map>
#include <cassert>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"


/** @class Graph
 * @brief A template for 3D undirected graphs.
 *
 * Users can add and retrieve nodes and edges. Edges are unique (there is at
 * most one edge between any pair of distinct nodes).
 */
template <typename V,typename E>
class Graph{
 
 public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type = V;
  using edge_value_type = E;

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

  private:
  
  //attribute: vector of nodes
  // std::vector<std::tuple<Point*,node_value_type>> nodes_;
  
  //Using maps to reduce runtime
  std::map<size_type,std::tuple<Point*,node_value_type>> nodes_;
  //attribute: edge list of vectors of node indices of each edge
  
  // std::vector<std::vector<size_type>> edge_list_;
  //Also using maps to reduce runtime
  //map[a] will contain all nodes that a is connected to
  std::map<size_type,std::vector<std::tuple<size_type,edge_value_type>>> edge_list_;
  // map input_to_curr_: key is index of the node when it is added; 
  //value is index in current graph nodes_ and edge_list_.
  std::map<size_type,size_type> input_to_curr_;
  //map curr_to_input_: key is index of the node in current nodes_ and edge_list_;
  //value is index of the node when it was added
  std::map<size_type,size_type> curr_to_input_;

  public:
  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  /** Construct an empty graph. */

  Graph():nodes_(std::map<size_type,std::tuple<Point*,node_value_type>>()),
          edge_list_(std::map<size_type,std::vector<std::tuple<size_type,edge_value_type>>>()),
          input_to_curr_(std::map<size_type,size_type>()),
          curr_to_input_(std::map<size_type,size_type>()){}

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
  class Node: private totally_ordered <Node> {
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
      //return dereferenced Point
      // return *(graph->nodes_[std::node_idx]);
      return *(std::get<0>(graph->nodes_[input_node_idx]));
    }

    Point& position(){
      // HW0: YOUR CODE HERE
      //return dereferenced Point
      // return *(graph->nodes_[std::node_idx]);
      return *(std::get<0>(graph->nodes_[input_node_idx]));
    }
    size_type num_nodes_in_graph()const{
      return graph->num_nodes();
    }
    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      // if (graph->input_to_curr_.find(input_node_idx) != graph->input_to_curr_.end()){
      // return graph->input_to_curr_.at(input_node_idx);}
      // else{return graph->input_to_curr_.size() +1; }
      return input_node_idx;
    }
    size_type return_curr_idx()const{
      return graph->input_to_curr_.at(input_index());
    }
    size_type input_index() const {
      if (graph->curr_to_input_.find(input_node_idx) != graph->curr_to_input_.end()){
      return graph->curr_to_input_.at(input_node_idx);}
      else{return graph->curr_to_input_.size() +1; }
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    node_value_type& value(){
      return std::get<1>(graph->nodes_.at(input_node_idx));
    };
    //set value function to assign value in shortest_path
    void set_value(node_value_type& new_value){
      std::get<1>(graph->nodes_[input_node_idx]) = new_value;
    };
    const node_value_type& value() const{
      return (std::get<1>(graph->nodes_[input_node_idx]));
    };
    size_type degree() const{
      //return degree of this node
      return graph->edge_list_.at(input_node_idx).size();
    };
    incident_iterator edge_begin() const{
      return IncidentIterator(graph,input_node_idx,0);
    };
    incident_iterator edge_end() const{
      return IncidentIterator(graph,input_node_idx,graph->edge_list_.at(input_node_idx).size());
    };

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      //If two nodes have same graph and same index then they should be equal
      //n.return_curr_idx()
      if (this->graph == n.graph && this->index()==n.return_curr_idx()){
        return true;
      }else{
        return false;
      }
      //(void) n;          // Quiet compiler warning
      //return false;
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
      //If one node is in a smaller graph than the other or 
      //index of one node is smaller than the other
      if ((this->graph<n.graph)||(this->index()<n.index())){
        return true;
      }else{
        return false;
      }
      //(void) n;           // Quiet compiler warning
      //return false;
    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    //Pointer back tp Graph container
    Graph* graph;
    //index of this node
    size_type node_idx;
    size_type input_node_idx;
    //Constructor to construct Node objects from Graph class methods
    Node(const Graph* new_graph, size_type new_node_idx)
        : graph(const_cast<Graph*>(new_graph)),input_node_idx(new_node_idx){
          node_idx = graph->curr_to_input_.at(input_node_idx);
        }

    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  std::map<size_type,std::tuple<Point*,node_value_type>> get_nodes()const{
    return this->nodes;
  }
  std::map<size_type,std::vector<std::tuple<size_type,edge_value_type>>> get_edge_list()const{
    return this->edge_list_;
  }
  size_type size() const {
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
   *
   * Complexity: O(1) amortized operations.
   */
  Node add_node(const Point& position, const node_value_type& type= node_value_type()) {
    //dynamically allocate position
    Point* pos = new Point;
    *pos = position;
    size_type new_node_idx = this->size();
    std::tuple<Point*,node_value_type> new_node(pos,type);
    nodes_[new_node_idx] = new_node;
    input_to_curr_[new_node_idx] = new_node_idx;
    curr_to_input_[new_node_idx] = new_node_idx;
    // this->nodes_.push_back(new_node);
    //check if num_nodes was modified correctly
    assert(this->num_nodes() == new_node_idx +1);
    std::vector<std::tuple<size_type,edge_value_type>> new_vec;
    this->edge_list_[new_node_idx] = new_vec;
    return Node(this,new_node_idx);
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    if (n.graph==this && n.index()<this->size()){
      return true;
    }else{return false;}
  }

  /** Return the node with index @a i.
   * @pre 0 <= @a i < num_nodes()
   * @post result_node.index() == i
   *
   * Complexity: O(1).
   */
  Node node(size_type i) const {
    assert(0<=i && i<num_nodes());
    return Node(this,i);

  }
  /** Remove given node n.
   * @post node n is removed
   * @post all edges incident from node n are removed
   * @post num_nodes = num_nodes-1
   * @post input_to_curr_, curr_to_input_, nodes_, edge_list_ has size = size-1;
   * @post nodes_, edge_list_, curr_to_input_ has keys indexed from 0 to num_nodes-1;
   * @post input_to_curr_ has values from 0 to num_nodes-1;
   * Complexity: O(1) for remove node, O(num edges()) for remove edges
   * Invalidated: node n,
   *              nodes_[num_nodes()-1];
   *              edge_list_[num_nodes()-1];
   *              curr_to_input_[num_nodes()-1];
   *              input_to_curr_[n.index()];
   *              Edge(graph,n.index(),edge_list_.at(a.index())[i]) for all i,
   *              Edge(graph,i,edge_list_.at(i).back()) for all i.
   */
  size_type remove_node(const Node& n){
    std::vector<std::tuple<size_type,edge_value_type>> a_edges = edge_list_.at(n.index());  
      while (n.degree()!=0){
        auto iter = n.edge_begin();
        //remove all edges instanced to node n
        size_type out = remove_edge(*iter);
        if (out==0){break;}
      }
      size_type initial_input_index = n.input_index();
      //swap the nodes_ values between this node index and the last index
      auto it1 = nodes_.find(n.index());
      auto it2 = nodes_.find(num_nodes()-1);
      auto end = nodes_.end();
      if(it1 != end && it2 != end) {
        //delete the pointer to Point
         delete std::get<0>(it1->second);
          std::swap(it1->second, it2->second);
      }

      //swap the edge_list_ values between this node and the last node
      auto edge_it1 = edge_list_.find(n.index());
      auto edge_it2 = edge_list_.find(num_nodes()-1);
      auto edge_end = edge_list_.end();
      if(edge_it1 != edge_end && edge_it2 != edge_end) {
          std::swap(edge_it1->second, edge_it2->second);
      }
      //need to swap the current index between this node and the last node
      //also swap the input index in curr_to_input_
      auto index_it1 = input_to_curr_.find(initial_input_index);
      size_type last_node = curr_to_input_.at(num_nodes()-1);
      auto index_it2 = input_to_curr_.find(last_node);
      auto index_end = input_to_curr_.end();
      auto index_it11 = curr_to_input_.find(n.index());
      auto index_it12 = curr_to_input_.find(num_nodes()-1);
      auto index_end1 = curr_to_input_.end();
      if((index_it1 != index_end && index_it2 != index_end)&&
      (index_it11 != index_end1 && index_it12 != index_end1)) {
          std::swap(index_it1->second, index_it2->second);
          std::swap(index_it11->second, index_it12->second);
      }
      //erase everything indexed to the last node after swapping
      //erase this index in input_to_curr_
      auto index_it = input_to_curr_.find(initial_input_index);
      input_to_curr_.erase(index_it);
      auto it = nodes_.find(num_nodes()-1);
      auto edge_it = edge_list_.find(num_nodes()-1);
      auto index_iter = curr_to_input_.find(num_nodes()-1);
      nodes_.erase(it);
      edge_list_.erase(edge_it);
      curr_to_input_.erase(index_iter);
      return 1;
  }
  /** Remove node n that the node iterator pointed to.
   * @post node n is removed
   * @post n_it points to the node currently at index of n_it
   * @post all edges incident from node n are removed
   * @post num_nodes = num_nodes-1
   * Complexity: O(num edges()) at most
   * Invalidated: node n,
   *              nodes_[num_nodes()-1];
   *              edge_list_[num_nodes()-1];
   *              curr_to_input_[num_nodes()-1];
   *              input_to_curr_[n.index()];
   *              Edge(graph,n.index(),edge_list_.at(a.index())[i]) for all i,
   *              Edge(graph,i,edge_list_.at(i).back()) for all i,
   *              iterator that pointed to the removed node.
   */
  node_iterator remove_node(node_iterator n_it){
    size_type curr_size = num_nodes();
    if((curr_size==0)||(n_it.return_iter_idx()>=curr_size)){
      return node_end();
    }else if(curr_size==1) {
      //if this is the last node in the graph, return node_end()
      Node n = *n_it;
      remove_node(n);
      return node_end();
    }
    else{
      //if not the last node, return iterator pointing to node indexed at iter_idx
      Node n = *n_it;
      size_type iter_idx = n_it.return_iter_idx();
      remove_node(n);
      return node_at_idx(iter_idx);
    }
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
  class Edge: private totally_ordered <Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }
    edge_value_type& value(){
      return std::get<1>(graph->edge_list_.at(node1_idx)[node2_input_vec_idx]);
    }
    const edge_value_type& value() const{
      return std::get<1>(graph->edge_list_.at(node1_idx)[node2_input_vec_idx]);
    }
    double length() const{
      return norm_2(Node(this->graph,node1_idx).position()-Node(this->graph,node2_idx).position());
    }
    Edge reverse_edge()const{
      for (size_type rev_idx = 0;rev_idx<graph->edge_list_.at(node2_idx).size();++rev_idx){
        if(std::get<0>(graph->edge_list_.at(node2_idx)[rev_idx])==node1_idx){
          return Edge(graph,node2_idx,rev_idx);
        }
      }
      return Edge();
    }
    double spring_constant(){
      return (double)graph->edge_list_.at(node1_idx).size();
    }
    /** Return a node of this Edge */
    size_type graph_size(){
      return graph->num_nodes();
    }
    Node node1() const {
      return Node(this->graph,node1_input_idx);
    }

    /** Return the other node of this Edge */
    Node node2() const {
      return Node(this->graph,node2_idx);      
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      //HW0: YOUR CODE HERE
        //If first node are the same and second node are the same.
        //or they are the same but in different order
        //Even though we don't add the edge if the same edge but different
        //node order exists, maybe we should still make this method able
        //to compare for this situation...
        if(((this->node1() == e.node1()) && (this-> node2() == e.node2()))
        ||((this->node1() == e.node2()) && (this-> node2() == e.node1()))){
        return true;
      }
      else{return false;}
    }

    /** Test whether this edge is less than @a e in a global order.
     *
     * This ordering function is useful for STL containers such as
     * std::map<>. It need not have any interpretive meaning.
     */

    std::vector<Node> node_ordered(const Edge& edge0) const {
        std::vector<Node> ordered_nodes;
        if(edge0.node1()<edge0.node2()){
          ordered_nodes.push_back(edge0.node1());
          ordered_nodes.push_back(edge0.node2());
        }else{ordered_nodes.push_back(edge0.node2());
          ordered_nodes.push_back(edge0.node1());}
        return ordered_nodes;
      }

    bool operator<(const Edge& e) const {
      //Same as operator==, we consider ordering in this method
        std::vector<Node> order1 = node_ordered(*this);
        std::vector<Node> order2 = node_ordered(e);

        if((order1[0]<order2[0])||((order1[0]==order2[0])&&(order1[1]<order2[1]))){
              return true;
          }else{return false;}
      }
      


   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    Graph* graph;
    size_type node1_idx;
    size_type node2_input_idx;
    size_type node2_idx;
    size_type node1_input_idx;
    size_type node2_input_vec_idx;
    Edge(const Graph* new_graph,size_type new_node1_idx,size_type new_node2_vec_idx):
      graph(const_cast<Graph*>(new_graph)),node1_input_idx(new_node1_idx),node2_input_vec_idx(new_node2_vec_idx){
        node1_idx = node1_input_idx;
        
        if(node1_idx < graph->edge_list_.size()){
        node2_input_idx = std::get<0>(graph->edge_list_.at(node1_idx)[node2_input_vec_idx]);
        node2_idx = graph->input_to_curr_.at(node2_input_idx);
        }
        else{node2_idx = 0;}
    }

  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    size_type num_edges = 0;
    for (auto const& x : edge_list_){
        num_edges = num_edges + x.second.size();
    }
    num_edges = num_edges/2;
    return num_edges;
  }

  /** Return the edge with index @a i.
   * @pre 0 <= @a i < num_edges()
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  Edge edge(size_type i) const {
    assert(0<=i && i<this->num_edges());
    return *std::next(edge_begin(), i);
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    assert(a.graph==this);
    assert(b.graph==this);
    size_type a_idx = a.return_curr_idx();
    size_type b_idx = b.input_index();
    //check in both orders
    std::vector<std::tuple<size_type,edge_value_type>> a_edges = edge_list_.at(a_idx);
    for (auto v: a_edges){
      if(std::get<0>(v)==b_idx){
        return true;
      }else{continue;}
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
    //(void) a, (void) b;   // Quiet compiler warning
    //return Edge();        // Invalid Edge
    // assert(a.graph==this && b.graph==this);
    assert(a.index() != b.index());
    std::vector<std::tuple<size_type,edge_value_type>> a_vec = edge_list_.at(a.index());

    for (size_type i = 0;i<a.degree();++i){
      if (std::get<0>(a_vec[i]) == b.index()){
        Edge new_edge(this,a.index(),i);
        return new_edge;
      }
    }
    edge_value_type val = edge_value_type();
    std::tuple<size_type,edge_value_type> a_tuple(a.index(),val);
    std::tuple<size_type,edge_value_type> b_tuple(b.index(),val);
    //add new edge to edge_list_
    this->edge_list_[a.index()].push_back(b_tuple);
    this->edge_list_[b.index()].push_back(a_tuple);
    // check if correctly added
    assert(has_edge(a,b));
    assert(has_edge(b,a));
    //check new num_edges
    // assert(old_num_edges+1 == num_edges());
    size_type b_in_vec = edge_list_.at(a.index()).size()-1;
    return Edge(this,a.index(),b_in_vec);
    

  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    
    // delete pointers
    for (auto p: this->nodes_){
      delete std::get<0>(p.second);
    }
    this->nodes_.clear();
    this->edge_list_.clear();
    assert(num_nodes()==0);
    assert(num_edges()==0);

  }
  /** Remove edge between two given nodes.
   * @pre Edge(a,b) still exits in the graph
   * @post Edge(a,b) is removed
   * @post Edge(b,a) also removed
   * @post num_edge = num_edge-1
   * Complexity: O(num edges())
   * Invalidated: Edge(graph,a.index(),edge_list_.at(a.index()).back())
   *              Edge(graph,b.index(),edge_list_.at(b.index()).back())
   */
  size_type remove_edge(const Node& a, const Node& b){
    if (!has_edge(a,b)|| !has_edge(b,a)){return 0;}
    else{
      //node1_idx, node2_idx is the index in current nodes_ and edge_list_
      //node1_idx_input, node2_idx_input is the index of the node when nodes_
      //and edge_list_ were constructed; 
      //Important to match the index in edge_list_.
      size_type node1_idx = a.index();
      size_type node1_idx_input = a.input_index();
      size_type node2_idx_input = b.input_index();
      size_type node2_idx = b.index();
      size_type vec1_size = edge_list_.at(node1_idx).size();
      for (size_type idx = 0;idx<vec1_size;++idx){
        if(std::get<0>(edge_list_.at(node1_idx)[idx])==node2_idx_input){
          //in edge_list_.at(node1_idx), swap the node to be removed in the value vector 
          //to the last and pop back
          std::swap(edge_list_.at(node1_idx)[idx],edge_list_.at(node1_idx).back());
          edge_list_.at(node1_idx).pop_back();
          break;
        }
      }
      //Also remove the reverse edges
      size_type vec2_size = edge_list_.at(node2_idx).size();
      for (size_type idx = 0;idx<vec2_size;++idx){
        if(std::get<0>(edge_list_.at(node2_idx)[idx])==node1_idx_input){
          std::swap(edge_list_.at(node2_idx)[idx],edge_list_.at(node2_idx).back());
          edge_list_.at(node2_idx).pop_back();
          break;
        }
      }
      return ((!has_edge(a,b))&&(!has_edge(b,a)));
    }
  }
  /** Remove input edge.
   * @pre edge still exits in the graph
   * @post edge is removed
   * @post reverse edge with the same two node but different order also removed
   * @post num_edge = num_edge-1
   * Complexity: O(num edges())
   * Invalidated: Edge(graph,edge.node1().index(),edge_list_.at(edge.node1().index()).back())
   *              Edge(graph,edge.node2().index(),edge_list_.at(edge.node2().index()).back())
   */
  size_type remove_edge(const Edge& edge){
    //for removing edge using EdgeIterator, node1() is current index consistent
    //consistent with input in above remove edge;
    //For node2(): from *edge_iterator returns the index of node2 in edge_list_ value vector,
    //Egde constructor is going to find the index of node2 in edge_list_ value vector,
    //which is the index of node2 when edge_list_ was constructed, and node2() returns the 
    //current index of node2, so also consistent with the above remove_node().
    node_type node1 = edge.node1();
    node_type node2 = edge.node2();
    return remove_edge(node1,node2);
  }
  /** Remove edge that edge_iterator points to,
   * and move the edge_iterator to point to a valid edge.
   * @pre edge still exits in the graph
   * @pre edge_iterator is pointing to valid edge
   * @post edge is removed
   * @post reverse edge with the same two node but different order also removed
   * @post num_edge = num_edge-1
   * @post edge_iterator is pointing to a valid edge
   * Complexity: O(num edges())
   * Invalidated: edge_iterator pointed to removed edge
   */
  edge_iterator remove_edge(edge_iterator e_it){
    Edge edge = *e_it;
    remove_edge(edge);
    edge_iterator new_it= e_it.remove();
    return new_it;
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
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid NodeIterator. */
    NodeIterator():graph(nullptr) {
    }

    Node operator*() const{
      return Node(graph,index);
    }
    NodeIterator& operator++(){
      //When index too large or at the end, set it to size of nodes_
      if (index >= graph_size-1){
        index = graph_size;
        return *this;
      }else{
        index = index+1;
        return *this;
      }
    
    }
    bool operator==(const NodeIterator& iter) const{
      return ((this->graph == iter.graph) && (this->index == iter.index));
    }
    size_type return_iter_idx() const{
      return index;
    }
    
    private:
    friend class Graph;
    const Graph* graph;
    size_type index;
    size_type graph_size;
    NodeIterator(const Graph* g, size_type idx):graph(g),index(idx){graph_size = graph->num_nodes();}
  };

  node_iterator node_begin() const{
    return NodeIterator(this,0);
  }
  node_iterator node_end() const{
    return NodeIterator(this,this->num_nodes());
  }
  node_iterator node_at_idx(size_type idx) const{
    assert(idx<this->num_nodes());
    return NodeIterator(this,idx);
  }


  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator:private totally_ordered<IncidentIterator>{
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid IncidentIterator. */
    IncidentIterator():graph(nullptr) {
    }

    Edge operator*() const{
          // return the edge, constructor pass in vector index of node2, not node index
           return Edge(graph,node1_idx,node2_idx);
    }
    IncidentIterator& operator++(){
      //If given node2_idx at the end of vector or probably too large
      //set it to vector size of edge_list_[node1_idx] to prevent seg fault
      size_type vec1_size = graph->edge_list_.at(node1_idx).size();
      if (node2_idx >= vec1_size-1){
        node2_idx = vec1_size;
        return *this;
      }else{
        node2_idx = node2_idx+1;
        return *this;
      }
    }
    bool operator==(const IncidentIterator& iter) const{
      return ((iter.graph == this->graph)&&(iter.node1_idx == this->node1_idx) 
                  && (iter.node2_idx == this->node2_idx));
    }

   private:
    friend class Graph;
    const Graph* graph;
    size_type node1_idx;
    size_type node2_idx;
    IncidentIterator(const Graph* g,size_type node1, size_type node2):
                      graph(g),node1_idx(node1),node2_idx(node2){}
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator:private totally_ordered<EdgeIterator> {
   public:
    // These type definitions let us use STL's iterator_traits.
    using value_type        = Edge;                     // Element type
    using pointer           = Edge*;                    // Pointers to elements
    using reference         = Edge&;                    // Reference to elements
    using difference_type   = std::ptrdiff_t;           // Signed difference
    using iterator_category = std::input_iterator_tag;  // Weak Category, Proxy

    /** Construct an invalid EdgeIterator. */
    EdgeIterator():graph(nullptr){
    }

    Edge operator*() const{
      //return the vector index of node2
      return Edge(graph,node1_idx,node2_idx);
    }
    EdgeIterator& operator++(){
        //When node2_idx is not at the end of its vector, first increment inside this vector
        //When at the end, move to next node1_idx
        node2_idx++;
        //If the node was node removed but all edges it was instanced to were removed
        //move to next node
        size_type vec1_size = graph->edge_list_.size();
        while((node1_idx<vec1_size)&&(graph->edge_list_.at(node1_idx).size()!=0)){
          vec2_size = graph->edge_list_.at(node1_idx).size();
          while(node2_idx<vec2_size){
            //To only visit each edge once, only return the iterator when node1_idx is
            //smaller than the node index of node2
            
            if (node1_idx < graph->input_to_curr_.at(std::get<0>(graph->edge_list_.at(node1_idx)[node2_idx]))){
              return *this;
            }
            ++node2_idx;
          }
          node1_idx = node1_idx+1;
          node2_idx = 0;
        }
        // //If haven't return after while loop, this means the iterator has reach the end
        node1_idx = vec1_size;
        node2_idx = 0;
        return *this;

        
      }
      

    
    bool operator==(const EdgeIterator& iter) const{
      return ((this-> graph == iter.graph)&&(this->node1_idx== iter.node1_idx)
              &&(this->node2_idx== iter.node2_idx));
    }
    
    EdgeIterator& remove(){
        while(node1_idx<vec1_size){
          vec2_size = graph->edge_list_.at(node1_idx).size();
          while(node2_idx<vec2_size-1){
            //To only visit each edge once, only return the iterator when node1_idx is
            //smaller than the node index of node2
            if (node1_idx < std::get<0>(graph->edge_list_.at(node1_idx)[node2_idx])){
              return *this;
            }
            node2_idx = node2_idx+1;
          }
          node1_idx = node1_idx+1;
          node2_idx = 0;
        }
        //If haven't return after while loop, this means the iterator has reach the end
        node1_idx = vec1_size;
        node2_idx = 0;
        return *this;
      }

   private:
    friend class Graph;
    const Graph* graph;
    size_type node1_idx;
    size_type node2_idx;
    size_type vec1_size;
    size_type vec2_size;
    //EdgeIterator Constructor: when passing in node1_idx too large, set vec2_size to 0 to prevent seg fault
    EdgeIterator(const Graph* g,size_type node1,size_type node2):graph(g),node1_idx(node1),node2_idx(node2){
                          vec1_size= graph->edge_list_.size();
                          if (node1_idx >= vec1_size){
                            vec2_size = 0;
                          }else{
                          vec2_size= graph->edge_list_.at(node1_idx).size();};
                        }
  };

  public:
  edge_iterator edge_begin() const{
    return EdgeIterator(this,0,0);
  }
  edge_iterator edge_end() const{
    return EdgeIterator(this,this->edge_list_.size(),0);
  }

 private:
};

#endif // CME212_GRAPH_HPP
