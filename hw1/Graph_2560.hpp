#ifndef CME212_GRAPH_HPP
#define CME212_GRAPH_HPP

/** @file Graph.hpp
 * @brief An undirected graph type
 */

#include <algorithm>
#include <vector>
#include <cassert>
#include <map>
#include <assert.h>

#include "CME212/Util.hpp"
#include "CME212/Point.hpp"

using namespace std;


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
  
  public:

  //
  // PUBLIC TYPE DEFINITIONS
  //
  using node_value_type  = V;
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

  //
  // CONSTRUCTORS AND DESTRUCTOR
  //

  //do a map
  struct adj_nodes{
    size_type edge_index;
    size_type b_index;

    adj_nodes(): edge_index(), b_index() {};
    adj_nodes(size_type index, size_type index_2): 
        edge_index(index), b_index(index_2){}
  };

  //(HW0 comment): constructing a struct storing info of nodes
  //include its position, index, 
  // and a vector of indexes of other nodes it is connected to (edges)
  struct node_details{
    Point the_point;
    size_type the_index;
    node_value_type the_value;
    std::vector<adj_nodes> pair_indexes; 
          //key:index of other nodes, val:index of other edge
    
    //default constructor with proxy design
    node_details(): the_point(), the_index(-1), the_value(), pair_indexes(){}

    //constructor:
    node_details(const Point &position, 
                 size_type index, node_value_type value):
    the_point(position), the_index(index), the_value(value), pair_indexes(){}

  };

  

  //(HW0 comment): construct a struct storing info of edges
  //include its index, the index of the two nodes it connected,
  //so we can quickly connect edge to the two nodes. 
  struct edge_details{
    size_type edge_index;
    size_type a_index;
    size_type b_index;
    //default constructor
    edge_details(): edge_index(), a_index(), b_index(){}
    //constructor with params.
    edge_details(size_type index, size_type index1, size_type index2): 
      edge_index(index), a_index(index1), b_index(index2){}
  };

  std::vector<size_type> proxy_node_ids;
  std::vector<size_type> proxy_edge_ids;


  /** Construct an empty graph. */
  Graph() {
    // HW0: YOUR CODE HERE
    num_of_edges = 0;
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
  class Node : private totally_ordered<Node>{
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

    Node(const Graph *new_graph, size_type idx):
    the_graph(const_cast<Graph*> (new_graph)), curr_index(idx) {}


    /** Return this node's position. */
    const Point& position() const {
      // HW0: YOUR CODE HERE 
      assert(curr_index < the_graph->num_nodes()); //curr_index within range; 
      return the_graph->nodes_map.at(curr_index).the_point;
    }

    /** Return this node's index, a number in the range [0, graph_size). */
    size_type index() const {
      // HW0: YOUR CODE HERE
      assert(curr_index < the_graph->num_nodes()); //curr_index within range; 
      return the_graph->nodes_map.at(curr_index).the_index;
    }

    // HW1: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // node_value_type& value(); okkkk
    // const node_value_type& value() const;  okkkkk
    // size_type degree() const;
    // incident_iterator edge_begin() const;
    // incident_iterator edge_end() const;

    
    node_value_type& value() {
      //return the value that the user input into node
      return the_graph -> nodes_map.at(curr_index).the_value;
    }

    const node_value_type& value() const {
      //return the constant reference to the input value
      return the_graph -> nodes_map.at(curr_index).the_value;
    }

    size_type degree() {
      //return the number of incident edges;
      return the_graph -> nodes_map.at(curr_index).pair_indexes.size();
    }

    incident_iterator edge_begin() const{
      return incident_iterator(the_graph,
                        curr_index, size_type(0));
    }

    incident_iterator edge_end() const{
      return incident_iterator(the_graph,
              curr_index, 
              the_graph -> nodes_map.at(curr_index).pair_indexes.size());
              //number of items in `pair_indexes` for the given node curr_index
    }

    /** Test whether this node and @a n are equal.
     *
     * Equal nodes have the same graph and the same index.
     */
    bool operator==(const Node& n) const {
      // HW0: YOUR CODE HERE
      (void) n;          // Quiet compiler warning
      if(this->the_graph == n.the_graph 
         && this->curr_index == n.curr_index)
      {
        return true;
      }
      else
      {
        return false;
      }
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
      (void) n;           // Quiet compiler warning
      if(this->curr_index < n.curr_index)
      {
        return true;
      }
      else if (this->curr_index > n.curr_index)
      {
        return false;
      }
      else
      {
        if(this->the_graph < n.the_graph)
        {
          return true;
        }

        return false;

      }

    }

   private:
    // Allow Graph to access Node's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Node
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Node objects
    Graph* the_graph;
    size_type curr_index;

  };

  /** Return the number of nodes in the graph.
   *
   * Complexity: O(1).
   */
  size_type size() const {
    // HW0: YOUR CODE HERE
    return proxy_node_ids.size();
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
                const node_value_type& val= node_value_type()) {
    // HW0: YOUR CODE HERE
    (void) position;      // Quiet compiler warning
    nodes_map.insert({num_nodes(), node_details(position, num_nodes(), val)});
    proxy_node_ids.push_back(num_nodes()); 
    //since num_nodes() depende on proxy_node_ids size, 
    //which has not been updated yet
    //so num_nodes() only contain previous number of nodes, 
    //and thus don't need `-1`
    //but when returning the Node, num_nodes() has already been updated
    return Node(this, num_nodes() - 1);  
  }

  /** Determine if a Node belongs to this Graph
   * @return True if @a n is currently a Node of this Graph
   *
   * Complexity: O(1).
   */
  bool has_node(const Node& n) const {
    // HW0: YOUR CODE HERE
    if(n.the_graph == this)
    {
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
    return Node(this, proxy_node_ids[i]); 
    //use proxy_node_ids bc it is possible that we delete/add an node
    //so the index of proxy_node_id != i
  }

  //
  // EDGES
  //

  /** @class Graph::E dge
   * @brief Class representing the graph's edges.
   *
   * Edges are order-insensitive pairs of nodes. Two Edges with the same nodes
   * are considered equal if they connect the same nodes, in either order.
   */
  class Edge: private totally_ordered<Edge> {
   public:
    /** Construct an invalid Edge. */
    Edge() {
      // HW0: YOUR CODE HERE
    }

    Edge(const Graph* new_graph, size_type edge_idx,
         size_type new_index1, size_type new_index2):
    the_graph(const_cast<Graph*>(new_graph)), the_edge_index(edge_idx),
    the_index_1(new_index1), the_index_2(new_index2){};

    /** Return a node of this Edge */
    Node node1() const {
      // HW0: YOUR CODE HERE
      return Node(the_graph, the_index_1);      
    }

    /** Return the other node of this Edge */
    Node node2() const {
      // HW0: YOUR CODE HERE
      return Node(the_graph, the_index_2);  
    }

    /** Test whether this edge and @a e are equal.
     *
     * Equal edges represent the same undirected edge between two nodes.
     */
    bool operator==(const Edge& e) const {
      (void) e;           // Quiet compiler warning
      //HW0: YOUR CODE HERE
      //check graph:
      if(the_graph == e.the_graph && the_edge_index == e.the_edge_index)
      {
        if(the_index_1 == e.the_index_1 &&
           the_index_2 == e.the_index_2)
        {
          return true;
        }
        else if (the_index_1 == e.the_index_2 &&
                 the_index_2 == e.the_index_1)
        { 
          return true;
        }
        
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
      //HW0: YOUR CODE HERE
      if(the_edge_index < e.the_edge_index)
      {
        return true;
      }
      else
      {
        return false;
      } 
    }

   private:
    // Allow Graph to access Edge's private member data and functions.
    friend class Graph;
    // HW0: YOUR CODE HERE
    // Use this space to declare private data members and methods for Edge
    // that will not be visible to users, but may be useful within Graph.
    // i.e. Graph needs a way to construct valid Edge objects
    Graph* the_graph;
    size_type the_edge_index;
    size_type the_index_1;
    size_type the_index_2;
  };

  /** Return the total number of edges in the graph.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  size_type num_edges() const {
    // HW0: YOUR CODE HERE
    return num_of_edges;
  }

  //  Return the edge with index @a i.
  // @pre 0 <= @a i < num_edges()
  //   Complexity: No more than O(num_nodes() + num_edges()), hopefully less
  
  Edge edge(size_type i) const {
    // HW0: YOUR CODE HERE
    (void) i;             // Quiet compiler warning
    size_type e = proxy_edge_ids[i];
    size_type a = edges_map.at(e).a_index;
    size_type b = edges_map.at(e).b_index;
    
    return Edge(this, e, a, b);  
  }

  /** Test whether two nodes are connected by an edge.
   * @pre @a a and @a b are valid nodes of this graph
   * @return True if for some @a i, edge(@a i) connects @a a and @a b.
   *
   * Complexity: No more than O(num_nodes() + num_edges()), hopefully less
   */
  bool has_edge(const Node& a, const Node& b) const {
    // HW0: YOUR CODE HERE
    (void) a; (void) b;   // Quiet compiler warning
    if(a.the_graph == b.the_graph
       && a.curr_index < nodes_map.size()
       && b.curr_index < nodes_map.size())
    {
      if(nodes_map.at(a.curr_index).pair_indexes.size() <
         nodes_map.at(b.curr_index).pair_indexes.size())
      {  
        for (const auto &id: nodes_map.at(a.curr_index).pair_indexes)
        {
          if (id.b_index == b.curr_index)
          {
            return true;
          }
        }
      }
      else
      {
        for (const auto &id: nodes_map.at(b.curr_index).pair_indexes)
        {
          if (id.b_index == a.curr_index)
          {
            return true;
          }
        }
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
    (void) a, (void) b;   // Quiet compiler warning
    if(!has_edge(a,b))
    {
      if(a.curr_index!= b.curr_index)
      {
        nodes_map.at(a.curr_index).\
            pair_indexes.push_back(adj_nodes(num_of_edges, b.curr_index));
        nodes_map.at(b.curr_index).
            pair_indexes.push_back(adj_nodes(num_of_edges, a.curr_index));
        
        
        proxy_edge_ids.push_back(edges_map.size());
        edges_map.insert({num_of_edges, 
            edge_details(num_of_edges, a.curr_index, b.curr_index)});
        
      }
      num_of_edges ++;
      
    }

    // cout << "funcitons beeeing called" << endl;
    // cout << "num of edges now " << num_of_edges << endl; 
    return Edge(this, num_of_edges-1, a.curr_index, b.curr_index);       
  }

  /** Remove all nodes and edges from this graph.
   * @post num_nodes() == 0 && num_edges() == 0
   *
   * Invalidates all outstanding Node and Edge objects.
   */
  void clear() {
    // HW0: YOUR CODE HERE
    nodes_map.clear(); //clear maps
    edges_map.clear();
    num_of_edges = 0;
    cout << "In Graph::clear() method: num of nodes now " 
         << num_nodes() << endl;
    cout << "        num of edges now" << num_edges() << endl;
    cout << "        num of nodes map size " << nodes_map.size() << endl;
    cout << "        num of edges map size " << edges_map.size() << endl;
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
    NodeIterator() {
    }

    // HW1 #2: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Node operator*() const
    // NodeIterator& operator++()
    // bool operator==(const NodeIterator&) const


    //(HW1 Comments)s
      //need to add this constructor that takes in the graph object itself
      //and the index. So it can be used in the begin() and end()
      //methods below
    NodeIterator(const Graph *new_graph, size_type idx):
    the_graph(const_cast<Graph*> (new_graph)), curr_index(idx) {}
      

    Node operator*() const {
      return Node(the_graph, the_graph->proxy_node_ids[curr_index]);
    }

    NodeIterator& operator++() {
      curr_index = curr_index + 1;
      return *this;
    }

    bool operator==(const NodeIterator& itr) const{
      if(this->curr_index == itr.curr_index && 
         this->the_graph == itr.the_graph)
      {
        return true;
      }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #2: YOUR CODE HERE
    Graph *the_graph;
    size_type curr_index;

  };

  // HW1 #2: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // node_iterator node_begin() const
  // node_iterator node_end() const

  node_iterator node_begin() const {
    //(HW1 Comment): see constructor above:
    return node_iterator(this, size_type(0));
    //also since we dont knwo what exactly the size_type is
    //need to use size_type(0) instead of just 0
  }

  node_iterator node_end() const {
    return node_iterator(this, this->proxy_node_ids.size());
  }

  //
  // Incident Iterator
  //

  /** @class Graph::IncidentIterator
   * @brief Iterator class for edges incident to a node. A forward iterator. */
  class IncidentIterator: private totally_ordered<IncidentIterator> {
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

    IncidentIterator(const Graph *newgraph, size_type aidx, size_type pidx):
      the_graph(const_cast<Graph*> (newgraph)), 
      the_index_1(aidx), 
      the_pair_index(pidx) {}

    // HW1 #3: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // IncidentIterator& operator++()
    // bool operator==(const IncidentIterator&) const

    Edge operator*() const{

      //the_index_1 is a proxyid, need to use proxy_node_ids
      //to get the unique id of the node 
      //the_pair_index is an index in the pair_indexes vectors
      //therefore can use iti directly
      size_type id_1 = the_graph->proxy_node_ids[the_index_1];
      size_type edge_id = the_graph->nodes_map.at(id_1).
                          pair_indexes[the_pair_index].edge_index;
      size_type id_2 = the_graph->nodes_map.at(id_1).
                       pair_indexes[the_pair_index].b_index;
       
      return Edge(the_graph, edge_id, id_1, id_2);
    }

    IncidentIterator& operator++(){
      the_pair_index ++;
      return *this;
    }

    bool operator==(const IncidentIterator& itr) const{
      if(this->the_graph == itr.the_graph &&
         this->the_index_1 == itr.the_index_1 &&
         this->the_pair_index == itr.the_pair_index)
         {
           return true;
         }
      return false;
    }

   private:
    friend class Graph;
    // HW1 #3: YOUR CODE HERE
    Graph *the_graph;
    size_type the_index_1;
    size_type the_pair_index;
  };

  //
  // Edge Iterator
  //

  /** @class Graph::EdgeIterator
   * @brief Iterator class for edges. A forward iterator. */
  class EdgeIterator: private totally_ordered<EdgeIterator>{
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

    EdgeIterator(const Graph *new_graph, size_type idx):
    the_graph(const_cast<Graph*> (new_graph)), the_edge_index(idx) {}

    // HW1 #5: YOUR CODE HERE
    // Supply definitions AND SPECIFICATIONS for:
    // Edge operator*() const
    // EdgeIterator& operator++()
    // bool operator==(const EdgeIterator&) const
    
    Edge operator*() const{
      size_type e = the_graph->proxy_edge_ids[the_edge_index];
      size_type a = the_graph->edges_map.at(e).a_index;
      size_type b = the_graph->edges_map.at(e).b_index;

      return Edge(the_graph, e, a, b);
    }

    EdgeIterator& operator++(){
      the_edge_index ++;
      return *this;
    }

    bool operator==(const EdgeIterator& itr) const{
      if(this->the_graph == itr.the_graph && 
         the_edge_index == itr.the_edge_index)
         {
           return true;
         }
      return false;
    }


   private:
    friend class Graph;
    // HW1 #5: YOUR CODE HERE
    Graph *the_graph;
    size_type the_edge_index;
  };

  // HW1 #5: YOUR CODE HERE
  // Supply definitions AND SPECIFICATIONS for:
  // edge_iterator edge_begin() const
  // edge_iterator edge_end() const

  edge_iterator edge_begin() const{
    return edge_iterator(this, size_type(0));
  }

  edge_iterator edge_end() const{
    return edge_iterator(this, this->proxy_edge_ids.size());
  }

 private:

  // HW0: YOUR CODE HERE
  // Use this space for your Graph class's internals:
  //   helper functions, data members, and so forth.

  std::map <size_type, node_details> nodes_map;
  std::map <size_type, edge_details> edges_map;
  size_type num_of_edges;

};

#endif // CME212_GRAPH_HPP
